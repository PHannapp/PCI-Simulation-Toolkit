import numpy as np
from scipy.optimize import minimize_scalar, brentq
from CONSTANTS import *
import matplotlib.pyplot as plt
import GH_data
import copy


#defining Gibbs function of gas phase
def gm_gas(T, p):
    return (GH_data.GHSERH2(T) + T*R_CONST*np.log(0.0000098692327*p) - 4.29e-06*GH_data.FF1(p) - 6.35e-06*GH_data.FF2(p) - 4.25e-06*GH_data.FF3(p) + 1.5e-06*GH_data.FF4(p) + 1.63e-06*GH_data.FF5(p) + 2.479e-06*p + 198428) / 2

    #return (GH_data.GHSERH2(T) + T*R_CONST*np.log(0.0000098692327*p))

# Create function for searching IP
def search_inflection_point(structure, kink):
    gm = structure.gm
    # Calculate approximate first derivative
    dy = np.diff(gm)
    # Find minima
    minimum_indices = np.where(np.diff(np.sign(dy)) > 0)[0] + 1
    # Find maxima
    maximum_indices = np.where(np.diff(np.sign(dy)) < 0)[0] + 1
    # Calculate second derivative (difference of differences)
    d2y = np.diff(dy)
    # Find points where the sign of the second derivative changes
    inflection_points_indices = np.where(np.diff(np.sign(d2y)))[0] + 1
    # Check and remove the first index if it's an inflection point / minima / maxima
    # TODO check if these work!!
    if inflection_points_indices.size > 0 and inflection_points_indices[0] == 1:
        inflection_points_indices = inflection_points_indices[1:]
    if minimum_indices.size > 0 and minimum_indices[0] == 1:
        minimum_indices = minimum_indices[1:]
    if maximum_indices.size > 0 and maximum_indices[0] == 1:
        maximum_indices = maximum_indices[1:]
    # Check and remove the last index if it's an inflection point / minima / maxima
    if inflection_points_indices.size > 0 and inflection_points_indices[-1] == len(d2y):
        inflection_points_indices = inflection_points_indices[:-1]
    if minimum_indices.size > 0 and minimum_indices[-1] == len(dy) - 1:
        minimum_indices = minimum_indices[:-1]
    if maximum_indices.size > 0 and maximum_indices[-1] == len(dy) - 1:
        maximum_indices = maximum_indices[:-1]

    for index in range(4,len(dy)-4):
            # Check for significant change in gradient
        #print(index, x_H[index], dy[index-1], dy[index], np.abs(dy[index] - dy[index-1]))
        if np.abs(dy[index] - dy[index-1]) > 10000:
            print('HEY')
            structure.kink = True
            structure.MINs = minimum_indices
            structure.MAXs = maximum_indices
            return []
    # clean list from consecutive orders
    inflection_points_indices_cleaned = []
    # Track if the previous number was part of a consecutive sequence
    prev_consecutive = False
    for i in range(1, len(inflection_points_indices)):
        if inflection_points_indices[i] != inflection_points_indices[i - 1] + 1 and not prev_consecutive:
            # Current number is not consecutive and previous number was not part of a consecutive sequence
            inflection_points_indices_cleaned.append(inflection_points_indices[i - 1])

        # Update the consecutive flag
        prev_consecutive = inflection_points_indices[i] == inflection_points_indices[i - 1] + 1

    # Handle the last number in the list
    if not prev_consecutive and len(inflection_points_indices) > 0:
        inflection_points_indices_cleaned.append(inflection_points_indices[-1])
    #print(f"Normal: {inflection_points_indices} vs cleaned: {inflection_points_indices_cleaned}")
    structure.kink = False
    structure.MINs = minimum_indices
    structure.MAXs = maximum_indices
    return inflection_points_indices_cleaned
# calculates the y value of a fully defined tangent
def tangent_line(x_H, indice, slope, intercept):
    return slope * x_H[indice] + intercept
# calculate the minimium offset from slope, intercept, and min and max boundaries
def calculate_tangent_offset(structure, slope, intercept, boundary_min, boundary_max):
    gm = structure.gm
    x_H = structure.x_H
    # Iterate over the range of indices and find the maximum negative difference
    max_negative_difference = float('-inf')
    max_negative_indice = 0
    for indice in range(boundary_min, boundary_max + 1):
        diff = (tangent_line(x_H, indice, slope, intercept) - gm[indice]) # Calculate the difference between tangent and function values
        if diff > max_negative_difference:
            max_negative_difference = diff
            max_negative_indice = indice
    return max_negative_indice, abs(max_negative_difference)  # Return the maximum offset and its x-coordinate
# make the minimization function with fixed x boundaries
#def make_calculate_offset_given_point(indice, boundary_min, boundary_max):
def calculate_offset_given_point(slope, structure, indice, boundary_min, boundary_max):
    gm = structure.gm
    x_H = structure.x_H
    intercept = gm[indice] - slope * x_H[indice]  # Calculate the intercept for the given slope
    _, max_offset = calculate_tangent_offset(structure, slope, intercept, boundary_min, boundary_max)  # Calculate the maximum offset for the given slope and intercept
    return max_offset  # Return the absolute maximum offset
#return calculate_offset_given_point
# create function for optimizing the slope to tangent the ab5 curve with fixed point x,y
def optimize_slope_given_point(structure, indice, boundary_min, boundary_max, bounds_slope):
    gm = structure.gm
    x_H = structure.x_H
    #calculate_offset_given_point = make_calculate_offset_given_point(indice, boundary_min, boundary_max)
    result = minimize_scalar(calculate_offset_given_point, bounds = bounds_slope, method='bounded', args=(structure, indice, boundary_min, boundary_max))
    optimal_slope = result.x  # The optimal slope that minimizes the maximum offset
    optimal_intercept = gm[indice] - optimal_slope * x_H[indice] # Calculate the intercept for the optimal slope
    result_x = calculate_tangent_offset(structure, optimal_slope, optimal_intercept, boundary_min, boundary_max) 
    return result_x[0], optimal_slope, optimal_intercept
# create function for searching common tangents below 
def search_common_tangents(structure, IPs, soe, min=0, MAX=1, CTs=None):
    gm = structure.gm
    x_H = structure.x_H
    #print(min, MAX, CTs)
    BOUND_MAX = len(x_H) - 1
    if CTs is None: # if first call --> initialize CT
        CTs = []
    if MAX > len(IPs)-1: # if arrived at last IP -> end recursive function
        return CTs, soe
    if min == 0: # what happens here?
        boundary_min = 0
    else: # what happens else?
        boundary_min = IPs[min - 1]
    if MAX == len(IPs)-1: 
        boundary_max = len(gm) - 1 - 1
    else: 
        boundary_max = IPs[MAX + 1]
    x1, x2 = IPs[min], IPs[MAX]
    #make line through IPs[min] and IPs[MAX] --> get slope
    # set IP[MAX] and optimize slope in boundaries IP[min] and IP[min-1]
    while True:
        initial_slope = max((gm[x2] - gm[x1]) / (x_H[x2] - x_H[x1]), SLOPE_MIN)
        if initial_slope > SLOPE_MAX - ACC:
            CTs = []
            print("WARNING 472: initial_slope > SLOPE_MAX - ACC")
            return CTs, False
        new_x1, optimal_slope, optimal_intercept = optimize_slope_given_point(structure, x2, boundary_min, x1, (initial_slope, SLOPE_MAX))
        initial_slope = max((gm[x2] - gm[new_x1]) / (x_H[x2] - x_H[new_x1]), SLOPE_MIN)
        new_x2, optimal_slope, optimal_intercept = optimize_slope_given_point(structure, new_x1, x2, boundary_max, (SLOPE_MIN, initial_slope))
        if abs(new_x1 - x1) < ACC and abs(new_x2 - x2) < ACC:
            # TODO CHECK
            x1, x2 = new_x1, new_x2
            break
        x1, x2 = new_x1, new_x2
    if x1 < 1 or abs(x2 - BOUND_MAX) < 1:
        return CTs, True
    # TODO Make the next two for loops better, with np.all(...)
    #look left of x1 if it cuts 
    for j in range(SLICES):
        x1_part = int(x1 / (SLICES - 1) * j)
        if gm[x1_part] - tangent_line(x_H, x1_part, optimal_slope, optimal_intercept) < -1:#or max_offset > 
            # if yes --> search_common_tangents(structure, IPs, min+2, MAX+2, firsttry=False, lasttry=False, CTs=None)
            CTs, soe = search_common_tangents(structure, IPs, soe, MAX+1, MAX+2, CTs)
            return CTs, soe
    #look right of x2 if it cuts (x3)
    for j in range(SLICES):
        x2_part = int((BOUND_MAX - x2) / (SLICES - 1) * (j) + x2)
        if gm[x2_part] - tangent_line(x_H, x2_part, optimal_slope, optimal_intercept) < -1:#or max_offset >
            # if yes --> take point left of x3: 
            for m in range(len(IPs)):
                if IPs[m] > x2_part:
                    new_x2 = m - 1
                    # if its indice difference to x2 is uneven --> take next one
                    if ((new_x2 - MAX ) & 1): # 
                        new_x2 = new_x2 + 1
                    CTs, soe = search_common_tangents(structure, IPs, soe, min, new_x2, CTs)
                    return CTs, soe
            # it cuts right of the last 
            return CTs, soe
    # if everywhere no --> append to CTs AND search_common_tangents(structure, IPs, min+2, MAX+2, firsttry=False, lasttry=False, CTs=None) AND return CTs
    CTs.append([x1,x2])
    CTs, soe = search_common_tangents(structure, IPs, soe, MAX+1, MAX+2, CTs)
    return CTs, soe
# finds p of gm_gas of corresponding Gibbs energy
def find_p(T, y):
    lower_bound = 1e-9  # Adjust this as needed.
    upper_bound = 1e19  # Adjust this as needed.
    try:
        return brentq(lambda p: gm_gas(T, p) - y, lower_bound, upper_bound, maxiter=10000)
    except ValueError:
        if y < 0:
            #print("WARNING 109, plateau pressure < 1e-9")
            return lower_bound
        else:
            #print("WARNING 808, plateau pressure > 1e18")
            return upper_bound
# finds plateau information
def plateau(structure, T):
    gm = structure.gm
    x_H = structure.x_H
    IPs = search_inflection_point(structure, kink=False)
    CTs, soe = search_common_tangents(structure, IPs, soe=False)
    #calculate plateau pressure: one point plus tangent slope --> point at x = 1 --> gm_gas
    plateau_slopes = []
    plateau_intercepts = []
    plateau_pressures = []
    for i in range(len(CTs)):
        plateau_slopes.append((gm[CTs[i][1]] - gm[CTs[i][0]]) / (x_H[CTs[i][1]] - x_H[CTs[i][0]]))
        plateau_intercepts.append(gm[CTs[i][0]] - plateau_slopes[i] * x_H[CTs[i][0]])
        g_gas = plateau_slopes[i] * G_X + plateau_intercepts[i]
        p = find_p(T, g_gas)
        plateau_pressures.append(p)
    structure.plateau_pressures = plateau_pressures
    structure.plateau_intercepts = plateau_intercepts
    structure.plateau_slopes = plateau_slopes
    structure.CTs = CTs
    structure.soe = soe

# make the minimization function with fixed x boundaries
#def make_calculate_offset_given_point_gas(x, y, boundary_min, boundary_max):
def calculate_offset_given_point_gas(slope, structure, x, y, boundary_min, boundary_max):
    intercept = y - slope * x  # Calculate the intercept for the given slope
    x1, max_offset = calculate_tangent_offset(structure, slope, intercept, boundary_min, boundary_max)  # Calculate the maximum offset for the given slope and intercept
    return max_offset  # Return the absolute maximum offset
#    return calculate_offset_given_point_gas
# create function for optimizing the slope to tangent the ab5 curve with fixed point x,y
def optimize_slope_given_point_gas(structure, x, y, boundary_min, boundary_max, bounds_slope):
    #calculate_offset_given_point_gas = make_calculate_offset_given_point_gas(x, y, boundary_min, boundary_max)
    result = minimize_scalar(calculate_offset_given_point_gas, bounds = bounds_slope, method='bounded', args=(structure, x, y, boundary_min, boundary_max))
    optimal_slope = result.x  # The optimal slope that minimizes the maximum offset
    optimal_intercept = y - optimal_slope * x # Calculate the intercept for the optimal slope
    result_x = calculate_tangent_offset(structure, optimal_slope, optimal_intercept, boundary_min, boundary_max) 
    #plt.figure()
    #plt.plot(structure.x_H,structure.gm)
    #plt.plot(x,y,'o')
    #plt.plot(np.linspace(0,1,100), [(optimal_slope * i + optimal_intercept) for i in np.linspace(0,1,100)])
    #plt.plot(structure.x_H[result_x[0]], optimal_slope * structure.x_H[result_x[0]] + optimal_intercept, 'o')
    #plt.grid()
    #plt.xticks(np.linspace(0,1,21))
    #plt.show()
    return result_x[0], optimal_slope, optimal_intercept
# optimizer for calculate_x_from_p
def optimize_x_given_p_gas(structure, p, MINIMUM_X, MAXIMUM_X, MINIMUM_SLOPE, MAXIMUM_SLOPE, T):
    x, slope, intercept = optimize_slope_given_point_gas(structure, G_X, gm_gas(T,p), MINIMUM_X, MAXIMUM_X, (MINIMUM_SLOPE, MAXIMUM_SLOPE))
    return x, slope, intercept
# calculates the x value to a corresponding p value from the pop file
def calculate_x_from_p(structure, p, T):
    gm = structure.gm
    x_H = structure.x_H
    CTs = structure.CTs
    plateau_pressures = structure.plateau_pressures
    #find where p is out of plateau pressures
    #calculate x from p by SLOPE_MIN or SLOPE_MAX
    # also X_Min and X_max from CTs
    for i in range(len(plateau_pressures)):
        if plateau_pressures[i] > p: # if p is below the ith plateau-pressure
            if i != 0: # if it is not the first
                MINIMUM_X, MAXIMUM_X = CTs[i-1][1], CTs[i][0]
                #MINIMUM_PRESSURE, MAXIMUM_PRESSURE = plateau_pressures[i-1], plateau_pressures[i]
                MINIMUM_SLOPE = (gm[CTs[i-1][1]] - gm[CTs[i-1][0]]) / (x_H[CTs[i-1][1]] - x_H[CTs[i-1][0]]) if x_H[CTs[i-1][1]] != x_H[CTs[i-1][0]] else 0
                MAXIMUM_SLOPE = (gm[CTs[i][1]] - gm[CTs[i][0]]) / (x_H[CTs[i][1]] - x_H[CTs[i][0]]) if x_H[CTs[i][1]] != x_H[CTs[i][0]] else 0

            else: #otherwise the lower bounds are zero
                MINIMUM_X, MAXIMUM_X = 0, CTs[i][0]
                #MINIMUM_PRESSURE, MAXIMUM_PRESSURE = 0, plateau_pressures[i]
                MINIMUM_SLOPE, MAXIMUM_SLOPE = SLOPE_MIN, max((gm[CTs[i][1]] - gm[CTs[i][0]]) / (x_H[CTs[i][1]] - x_H[CTs[i][0]]) if x_H[CTs[i][1]] != x_H[CTs[i][0]] else 0, SLOPE_MIN)
            break
        elif i == len(plateau_pressures) - 1: # if no plateau-pressure is above:
            MINIMUM_X, MAXIMUM_X = CTs[i][1], len(x_H)-1
            MINIMUM_PRESSURE, MAXIMUM_PRESSURE = plateau_pressures[i], P_MAX
            MINIMUM_SLOPE, MAXIMUM_SLOPE = (
            (gm[CTs[i][1]] - gm[CTs[i][0]]) / (x_H[CTs[i][1]] - x_H[CTs[i][0]]) if x_H[CTs[i][1]] != x_H[CTs[i][0]] else 0,
            SLOPE_MAX
            )
    return optimize_x_given_p_gas(structure, p, MINIMUM_X, MAXIMUM_X, MINIMUM_SLOPE, MAXIMUM_SLOPE, T)

def search_from_left_to_right(gm, MINs, start, end):
    x = np.arange(len(gm))
    connectionline_slope = (gm[end] - gm[start]) / (x[end] - x[start])
    connectionline_intercept = gm[start] - connectionline_slope * x[start]
    connectionline = connectionline_slope * x + connectionline_intercept

    for min_value in reversed(MINs):
        if start < min_value < end and gm[min_value] < connectionline[min_value]:
            return search_from_left_to_right(gm, MINs, start, min_value)

    return start, end

def search_common_tangents2(structure):
    gm = structure.gm
    x_H = structure.x_H
    MINs = structure.MINs
    MAXs = structure.MAXs
    CT2s = []
    # if no minima or maxima
    if len(MINs) == 0 and len(MAXs) == 0:
        # to calculate concavity or convexity i calculate the line between the endpoints and see if the middle point above or below the connection line
        x = np.arange(len(gm))
        connectionline_slope = (gm[-1] - gm[0]) / (x[-1] - x[0])
        connectionline_intercept = gm[0] - connectionline_slope * x[0]
        connectionline = connectionline_slope * x + connectionline_intercept
        # Take the midpoint of the array
        midpoint = len(gm) // 2
        # if concav -> first and last point is CT2s
        if gm[midpoint] > connectionline[midpoint]:
            # if not all points are below curve -> probably cubic function
            if not np.all(connectionline < gm):
                # if the second is higher then the first entry -> set left point, else the right point
                if gm[2] > gm[0]:
                    MIN = [0, 1]
                    MAX = [0, len(gm) - 1]
                    x1, x2 = 0, len(gm) - 1
                else:
                    MIN = [0, len(gm) - 1]
                    MAX = [len(gm) - 2, len(gm) - 1]
                    x1, x2 = 0, len(gm) - 1
                while True:
                    initial_slope = (gm[x2] - gm[x1]) / (x_H[x2] - x_H[x1])
                    if initial_slope > SLOPE_MAX - ACC:
                        CTs = []
                        print("WARNING 472: initial_slope > SLOPE_MAX - ACC")
                        return CTs, False
                    # set right point and open left point: indice = x2, boundaries are MIN and MAX[0]
                    new_x1, optimal_slope, optimal_intercept = optimize_slope_given_point(structure, x2, MIN[0], MAX[0], (initial_slope, SLOPE_MAX))
                    initial_slope = (gm[x2] - gm[new_x1]) / (x_H[x2] - x_H[new_x1])
                    new_x2, optimal_slope, optimal_intercept = optimize_slope_given_point(structure, new_x1, MIN[1], MAX[1], (initial_slope, SLOPE_MAX))
                    if abs(new_x1 - x1) < ACC and abs(new_x2 - x2) < ACC:
                        # TODO CHECK
                        x1, x2 = new_x1, new_x2
                        break
                    x1, x2 = new_x1, new_x2
                CT2s.append([x1, x2])
            else:
                CT2s.append([0,len(gm) - 1])
                return CT2s, 
        # if convex -> no CT2s
        elif gm[midpoint] < connectionline[midpoint]:
            # the lower point (last or first) is the orientation point
            orientation_point = (len(gm) - 1) if gm[-1] < gm[0] else 0
            return [], orientation_point
        else:
            print("HELP")
            CT2s.append([0,len(gm) - 1])
            return CT2s, 
    # if one maximum -> concave -> first and last point is CT2s
    elif len(MINs) == 0 and len(MAXs) == 1:
        x = np.arange(len(gm))
        connectionline_slope = (gm[-1] - gm[0]) / (x[-1] - x[0])
        connectionline_intercept = gm[0] - connectionline_slope * x[0]
        connectionline = connectionline_slope * x + connectionline_intercept
        if not np.all(connectionline < gm):
            if gm[2] > gm[0]:
                MIN = [0, 1]
                MAX = [0, len(gm) - 1]
                x1, x2 = 0, len(gm) - 1
            else:
                MIN = [0, len(gm) - 1]
                MAX = [len(gm) - 2, len(gm) - 1]
                x1, x2 = 0, len(gm) - 1
            while True:
                initial_slope = (gm[x2] - gm[x1]) / (x_H[x2] - x_H[x1])
                if initial_slope > SLOPE_MAX - ACC:
                    CTs = []
                    print("WARNING 472: initial_slope > SLOPE_MAX - ACC")
                    return CTs, False
                # set right point and open left point: indice = x2, boundaries are MIN and MAX[0]
                new_x1, optimal_slope, optimal_intercept = optimize_slope_given_point(structure, x2, MIN[0], MAX[0], (initial_slope, SLOPE_MAX))
                initial_slope = (gm[x2] - gm[new_x1]) / (x_H[x2] - x_H[new_x1])
                new_x2, optimal_slope, optimal_intercept = optimize_slope_given_point(structure, new_x1, MIN[1], MAX[1], (initial_slope, SLOPE_MAX))
                if abs(new_x1 - x1) < ACC and abs(new_x2 - x2) < ACC:
                    # TODO CHECK
                    x1, x2 = new_x1, new_x2
                    break
                x1, x2 = new_x1, new_x2
            CT2s.append([x1, x2])
        else:
            CT2s.append([0,len(gm) - 1])
        return CT2s, 
    # if one minimum -> convex -> minimum as orientation point
    elif len(MINs) == 1 and len(MAXs) == 0:
        orientation_point = MINs[0]
        return [], orientation_point
    # if one minimum and one maximum: 
    elif len(MINs) == 1 and len(MAXs) == 1:
        # if minima left of maxima -> CT2 between minma and gm[-1]
        if MINs[0] < MAXs[0]:
            # optimize slope between gm[-1] and minima: 
            # if minima smaller then gm[-1] --> CT2s[0][0] has to be right of x_H[MINs[0]]
            if gm[MINs[0]] < gm[-1]:
                boundary_min = MINs[0]
                boundary_max = MAXs[0]
            # if minima bigger then gm[-1] --> CT2s[0][0] has to be left of x_H[MINs[0]]
            else:                
                boundary_min = 0
                boundary_max = MINs[0]
            initial_slope = max((gm[-1] - gm[MINs[0]]) / (x_H[-1] - x_H[MINs[0]]), SLOPE_MIN)
            CT2_0, _, _ = optimize_slope_given_point(structure, len(gm) - 1, boundary_min, boundary_max, (SLOPE_MIN, initial_slope))
            CT2s.append([CT2_0, len(gm) - 1])
            return CT2s, 
        # if minima right of maxima -> CT2 between minima and gm[0]
        else:
            # optimize slope between gm[0] and minima: 
            # if minima smaller then gm[0] --> CT2s[0][1] has to be left of x_H[MINs[0]]
            if gm[MINs[0]] < gm[0]:
                boundary_min = MAXs[0]
                boundary_max = MINs[0]
            # if minima bigger then gm[0] --> CT2s[0][1] has to be right of x_H[MINs[0]]
            else:                
                boundary_min = MINs[0]
                boundary_max = len(gm) - 1
            initial_slope = max((gm[0] - gm[MINs[0]]) / (x_H[0] - x_H[MINs[0]]), SLOPE_MIN)
            CT2_1, _, _ = optimize_slope_given_point(structure, 0, boundary_min, boundary_max, (SLOPE_MIN, initial_slope))
            CT2s.append([0, CT2_1])
            return CT2s, 
    
    ## now every gms that have more then two maxima / minima
    # if more MAX then MIN --> has to start at 0 and end at len(gm) - 1. In between have to see
    elif len(MAXs) > len(MINs):        
        # algorithm: 
        # 1. start=0, end=len(gm)-1.
        # 2. draw connection between start and end
        # 3. see (from right to left) if any MIN[?] is lower
        # 4. if yes: end = MIN[?] -> loop to 2.
        # 5. if no: CTs.append([start, end]); start = end, end =len(gm)-1; loop to 2.  
        start = 0  # Starting at the first element
        end = len(gm) - 1  # Ending at the last element

        while start < len(gm) - 1:
            start, end = search_from_left_to_right(gm, MINs, start, end)
            CT2s.append([start, end])
            start = end
            end = len(gm) - 1
    # if more MIN than MAX --> must NOT start at 0 and len(gm) - 1
    elif len(MINs) > len(MAXs):
        start = MINs[0]  # Starting at the first element
        end = MINs[-1] # Ending at the last element

        while start < MINs[-1]:
            start, end = search_from_left_to_right(gm, MINs, start, end)
            CT2s.append([start, end])
            start = end
            end = MINs[-1]
    # if len(MIN) == len(MAX) and more than 1
    else:
        # if the smallest is MIN --> len(gm) - 1 is CT2s[?][1]
        if min(MINs) < min(MAXs):
            start = MINs[0]  # Starting at the first element
            end = len(gm) - 1 # Ending at the last element

            while start < len(gm) - 1:
                start, end = search_from_left_to_right(gm, MINs, start, end)
                CT2s.append([start, end])
                start = end
                end = len(gm) - 1
        else:
            start = 0  # Starting at the first element
            end = MINs[-1] # Ending at the last element
    
            while start < MINs[-1]:
                start, end = search_from_left_to_right(gm, MINs, start, end)
                CT2s.append([start, end])
                start = end
                end = MINs[-1]
    return CT2s, 
                
def plateau_from_minima(structure, T):
    gm = structure.gm
    x_H = structure.x_H
    result = search_common_tangents2(structure)
    if len(result) == 2:
        CT2s, orientation_point = result
    else:
        (CT2s,) = result
    #calculate plateau pressure: one point plus tangent slope --> point at x = 1 --> gm_gas
    plateau2_slopes = []
    plateau2_intercepts = []
    plateau2_pressures = []
    for i in range(len(CT2s)):
        plateau2_slopes.append((gm[CT2s[i][1]] - gm[CT2s[i][0]]) / (x_H[CT2s[i][1]] - x_H[CT2s[i][0]]))
        plateau2_intercepts.append(gm[CT2s[i][0]] - plateau2_slopes[i] * x_H[CT2s[i][0]])
        g_gas = plateau2_slopes[i] * G_X + plateau2_intercepts[i]
        p = find_p(T, g_gas)
        plateau2_pressures.append(p)
    # if no plateau found -> orientation point with zero slope and intercept and gas phase at orientation point
    if len(CT2s) == 0:
        plateau2_slopes.append(0)
        g_gas = gm[orientation_point]
        plateau2_intercepts.append(g_gas)
        p = find_p(T, g_gas)
        plateau2_pressures.append(p)
        structure.plateau_pressures = plateau2_pressures
        structure.plateau_slopes = plateau2_slopes
        structure.plateau_intercepts = plateau2_intercepts
        structure.orientation_point = orientation_point
        structure.CTs = [[structure.orientation_point, structure.orientation_point]]
        return 
    structure.plateau_pressures = plateau2_pressures
    structure.plateau_slopes = plateau2_slopes
    structure.CTs = CT2s
    structure.plateau_intercepts = plateau2_intercepts

def interpolate(x, x1, y1, x2, y2):
    if x2 == x1:
        return y1
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)

def merge_g_values(x_values, g_values):
    x_all = []
    y_all = []
    list_index = []

    # Combine all x-values and sort them
    combined_x = sorted(set(x for sublist in x_values for x in sublist))

    for x in combined_x:
        min_y = float('inf')
        y = float('inf')
        min_idx = -1
        flag = True
        for i, (x_list, y_list) in enumerate(zip(x_values, g_values)):
            if x in x_list:
                y = y_list[np.where(x_list == x)[0][0]]
                flag = True
            else:
                # Find adjacent points for interpolation
                for j in range(len(x_list) - 1):
                    if x_list[j] < x < x_list[j + 1]:
                        y = interpolate(x, x_list[j], y_list[j], x_list[j + 1], y_list[j + 1])
                        flag = False
                        break
            if y < min_y:
                min_y = y
                min_idx = i
                flag_save = flag

        # Add the selected point to the result lists
        if (not x_all or x_all[-1] != x) and (flag_save):
            x_all.append(x)
            y_all.append(min_y)
            list_index.append(min_idx)

    return np.array(x_all), np.array(y_all), np.array(list_index)

def extract_changes(list_index):
    # Convert to numpy array if not already
    list_index = np.array(list_index)
    # Find where changes occur
    change_points = np.where(list_index[:-1] != list_index[1:])[0]
    # Number of changes
    n_changes = len(change_points)
    # List of change points
    list_changes = change_points.tolist()
    return n_changes, list_changes

def plateau_between_phases(global_structure, structure_list, T):
    gms = []
    x_Hs = []
    #structure_list = copy.deepcopy(global_structure_list)
    n_structure = len(structure_list)    
    for i in range(n_structure):
        gms.append(np.array(structure_list[i].gm))
        x_Hs.append(np.array(structure_list[i].x_H))
    # 1. create on gm_all, where each entry is the minimium value of all
    # list_index shows the index, to which structure the corresponding gm and x_H value belongs
    x_all, gm_all, list_index = merge_g_values(x_Hs, gms)
    global_structure.gm = gm_all
    global_structure.x_H = x_all
    # how many changes (and possible plateaus) and where the borders between phase changes are
    n_changes, list_changes = extract_changes(list_index)

    ## check if one gm is fully above the other gm --> its suspended
    #for structure in range(n_structure):
    #    for partner_structure in range(n_structure):
    #        if structure != partner_structure:
    #            # check if upper boundaries are bigger and lower boundaries are smaller then partner structure --> possible to be completely under partner_structure
    #            if x_Hs[structure][0] < x_Hs[partner_structure][0] and x_Hs[structure][-1] > x_Hs[partner_structure][-1]:
    #                # when x_H of structure is bigger then partner_structure -> gm of structure also has to be smaller
    #                if np.all(gms[structure] < gms[partner_structure]):
    #                    global_structure_list[partner_structure].suspended = True
    #                    structure_list[partner_structure].suspended = True
    ## Remove suspended structures and their corresponding gm, x_H, and CTs entries
    #gms = [gm for i, gm in enumerate(gms) if not structure_list[i].suspended]
    #x_Hs = [x_H for i, x_H in enumerate(x_Hs) if not structure_list[i].suspended]
    #structure_list = [structure for structure in structure_list if not structure.suspended]
    #n_structure = len(structure_list)
    #CTs = [[] for _ in range(n_structure)]
    #plateau_pressures = [[] for _ in range(n_structure)]
    
    # check all plateaus if elidgable --> copy to plateau between
    for structure in range(n_structure):
        if not isinstance(structure_list[structure].orientation_point, list):
            continue
        for plateau in range(len(structure_list[structure].plateau_slopes)):
            plateau_flag = True # this flags indicates if the plateau stays (because its below all gms)
            for partner_structure in range(n_structure):
                if structure != partner_structure:
                    # Calculate common tangent line for the range of the partner_structure
                    common_tangent_line = (
                        structure_list[structure].plateau_slopes[plateau] * x_Hs[partner_structure] +
                        structure_list[structure].plateau_intercepts[plateau]
                    )
                    if not np.all(common_tangent_line < gms[partner_structure]):
                        plateau_flag = False
                        break
            if plateau_flag:
                CT_glob = []
                # have to convert the structure_list[structure].CTs[plateau] values to the new x_all, gm_all values
                XH0 = structure_list[structure].x_H[structure_list[structure].CTs[plateau][0]]
                XH1 = structure_list[structure].x_H[structure_list[structure].CTs[plateau][1]]
                #GM0, GM1 = structure_list[structure].CTs[plateau][0], structure_list[structure].CTs[plateau][1]
                for i in range(len(x_all)):
                    if x_all[i] == XH0: #and gm_all[i] == GM0:
                        CT_glob.append(i)
                    if x_all[i] == XH1:# and gm_all[i] == GM1:
                        CT_glob.append(i)
                global_structure.CTs.append(CT_glob)
                global_structure.phase_id.append([list_index[CT_glob[0]], list_index[CT_glob[1]]])
                global_structure.plateau_pressures.append(structure_list[structure].plateau_pressures[plateau])
                global_structure.n_plateau += 1
                global_structure.plateau_slopes.append(structure_list[structure].plateau_slopes[plateau])
                global_structure.plateau_intercepts.append(structure_list[structure].plateau_intercepts[plateau])
    # TODO check if Minima ousiode of plateau! --> postpone, as shouldnt happen, and if its happening  we should discard it anyway
    # concatinate CTs
    CTs_cat = global_structure.CTs#[inner_list for sublist in global_structure.CTs for inner_list in sublist]

    # start with for loop through all possible phase changes
    BOUND_MAX = len(x_all) - 1
    for n in range(n_changes):
        # select boundaries for common tangent search
        # choose first point of n and last point of n + 1 structure.
        # MIN / MAX[0] is for the left phase, MIN / MAX[1] is for the right phase
        if n == 0:
            MIN = [0, list_changes[n]]
            x1 = 0
        else: 
            MIN = [list_changes[n - 1], list_changes[n]]
            x1 = list_changes[n - 1]
        if n == n_changes - 1:
            MAX = [list_changes[n], len(x_all) - 1]
            x2 = len(x_all) - 1
        else:
            MAX = [list_changes[n], list_changes[n + 1]]
            x2 = list_changes[n + 1]
        # check if bounaries change, because of already found plateaus 
        ENCLOSED = False
        for CT in CTs_cat:
            # check if existing CT encloses MIN[0] and MAX[0] or MIN[1] and MAX[1]
            if (CT[0] <= MIN[0] and CT[1] >= MAX[0]) or (CT[0] <= MIN[1] and CT[1] >= MAX[1]):
                ENCLOSED = True
                break
            # if CT[0] (could also be CT[1]) is between MIN[1] and MAX[1]: plateau is in between the right phase: -> change MAX[1] to CT[0]
            if MIN[1] <= CT[0] <= MAX[1]:
                MAX[1] = CT[0]
                x2 = CT[0]
            # if CT[1] (could also be CT[0]) is between MIN[0] and MAX[0]: plateau is in between the left phase: -> change MIN[0] to CT[1]
            if MIN[0] <= CT[1] <= MAX[0]:
                MIN[0] = CT[1]
                x1 = CT[1]
        if ENCLOSED:
            continue
        # search the commmon tangent
        while True:
            initial_slope = (gm_all[x2] - gm_all[x1]) / (x_all[x2] - x_all[x1])
            if initial_slope > SLOPE_MAX - ACC:
                CTs = []
                print("WARNING 472: initial_slope > SLOPE_MAX - ACC")
                return CTs, False
            # set right point and open left point: indice = x2, boundaries are MIN and MAX[0]
            new_x1, optimal_slope, optimal_intercept = optimize_slope_given_point(global_structure, x2, MIN[0], MAX[0], (initial_slope, SLOPE_MAX))
            initial_slope = (gm_all[x2] - gm_all[new_x1]) / (x_all[x2] - x_all[new_x1])
            new_x2, optimal_slope, optimal_intercept = optimize_slope_given_point(global_structure, new_x1, MIN[1], MAX[1], (initial_slope, SLOPE_MAX))
            if abs(new_x1 - x1) < ACC and abs(new_x2 - x2) < ACC:
                # TODO CHECK
                x1, x2 = new_x1, new_x2
                break
            x1, x2 = new_x1, new_x2
        if x1 < 1 or abs(x2 - BOUND_MAX) < 1:
            global_structure.soe = True
        
        # check if below all points of gm_all
        common_tangent_line = (optimal_slope * x_all + optimal_intercept)
        if np.all(common_tangent_line < gm_all):
            # TODO do here another run with new boundaries (... above phase is on left phase -> adjust MIN[0], otherwise adjust MAX[1]
            continue
    # copy results to global list
        global_structure.CTs.append([x1,x2])
        g_gas = optimal_slope * G_X + optimal_intercept
        p = find_p(T, g_gas)
        global_structure.plateau_pressures.append(p)
        global_structure.phase_id.append([list_index[x1], list_index[x2]])
        global_structure.n_plateau += 1
        global_structure.plateau_slopes.append(optimal_slope)
        global_structure.plateau_intercepts.append(optimal_intercept)
    # Step 1: Pair the elements together
    paired = list(zip(
        global_structure.CTs, 
        global_structure.plateau_pressures, 
        global_structure.phase_id,
        global_structure.plateau_slopes, 
        global_structure.plateau_intercepts))
    
    # Step 2: Sort the pairs based on the first element
    paired.sort()
    # Step 3: Separate the pairs back into four arrays
    if paired and all(len(p) == len(paired[0]) for p in paired):  # Check if paired is non-empty and elements have equal lengths
        a_sorted, b_sorted, c_sorted, d_sorted, e_sorted = zip(*paired)
        # Convert them back to lists (optional, if needed as lists)
        global_structure.CTs = list(a_sorted)
        global_structure.plateau_pressures = list(b_sorted)
        global_structure.phase_id = list(c_sorted)
        global_structure.plateau_slopes = list(d_sorted)
        global_structure.plateau_intercepts = list(e_sorted)  
    #else:
    #    print("Error: Paired data is empty or contains elements with unequal lengths.")

  
    
    # if no plateaus -> again orientation point of apparent phase 
    if len(global_structure.CTs) == 0:
        global_structure.CTs.append([structure_list[list_index[0]].orientation_point, structure_list[list_index[0]].orientation_point])
        global_structure.plateau_pressures = (structure_list[list_index[0]].plateau_pressures)
        global_structure.phase_id.append([list_index[0], list_index[0]])
        global_structure.plateau_slopes = (structure_list[list_index[0]].plateau_slopes)
        global_structure.plateau_intercepts = (structure_list[list_index[0]].plateau_intercepts)

    
    
    
    
        