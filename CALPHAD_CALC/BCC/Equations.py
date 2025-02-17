import numpy as np
import types
import inspect
from scipy.optimize import minimize, minimize_scalar, differential_evolution
import matplotlib.pyplot as plt
import time
from itertools import product
import importlib


R_CONST = 8.3144621 # Gas constant

class Structure:
    def __init__(self, solids, multiplicity, site_fractions, name, database_module):
        self.solids = solids
        self.multiplicity = multiplicity
        self.site_fractions_solids = site_fractions.copy() 
        self.n_hydrogen_sites = len(multiplicity) - len(site_fractions)
        self.site_fractions_H = []
        self.site_fractions = site_fractions.copy() 
        self.multiplicity_H = multiplicity[-self.n_hydrogen_sites:]
        self.multiplicity_solids = multiplicity[:len(site_fractions)]
        self.hydrogen_start_value = 0 # if Hydrogen as non HV site -> translation of gm curve
        self.gm = None
        self.x_H = None
        self.name = name
        self.database_module = database_module
        self.plateau_pressures = []
        self.CTs = []
        self.plateau_slopes = []
        self.plateau_intercepts = []
        self.orientation_point = []
        self.MINs = []
        self.MAXs = []
        self.kink = False
        self.soe = False

        for i, site in enumerate(self.site_fractions):
            if abs(sum(site) - 1) > 1e-5:
               raise ValueError("Sum in site ", i ,"is not == 1")
        for i, site in enumerate(self.solids):
            if site == ['HY']:
                self.hydrogen_start_value += self.multiplicity[i]

# if multiple phases -> global structure
class GLOBAL_Structure:
    def __init__(self, solids, multiplicity, site_fractions, name, database_module):
        self.solids = solids
        self.multiplicity = multiplicity
        self.site_fractions_H = []
        self.gm = None
        self.x_H = None
        self.name = name
        self.database_module = database_module
        self.plateau_pressures = []
        self.phase_id = []
        self.n_plateau = 0
        self.CTs = []
        self.plateau_slopes = []
        self.plateau_intercepts = []
        self.soe = False

class Energies:
    def __init__(self, structure, functions, type, HV_Letters, solid_letters, vars, vars_count):
        self.type = type
        self.functions = functions
        self.indices_solids = [] # [[0,1,1,0], [0,0,1,0], [1,0,1,0], [1,1,1,0]]
        self.indices_HV = [] # [[0,1,1,0], [0,0,1,0], [1,0,1,0], [1,1,1,0]]
        self.order = []
        self.HV_Letters = HV_Letters
        self.solid_letters = solid_letters # [[LA, AL], [LA,AL]
        self.variables = vars
        self.vars_count = vars_count

        if self.type == 'reference':
            for i, two_letter in enumerate(solid_letters):
                self.indices_solids.append([])
                self.indices_HV.append([])
                for n, two_let in enumerate(two_letter):
                    self.indices_solids[i].append(structure.solids[n].index(two_let))
                for n, one_letter in enumerate(HV_Letters[i]):
                    self.indices_HV[i].append(0) if one_letter == 'H' else self.indices_HV[i].append(1)

        if self.type == 'interaction':
            for i, function in enumerate(functions):
                self.order.append(int(function[-1]))
                self.indices_solids.append([])
                self.indices_HV.append([])
                for n, solid_letter in enumerate(solid_letters[i]):
                    if len(solid_letter) == 2:
                        self.indices_solids[i].append([structure.solids[n].index(solid_letter)])
                    else:
                        self.indices_solids[i].append([structure.solids[n].index(m) for m in sorted([solid_letter[k:k+2] for k in range(0, (len(solid_letter)), 2)])])
                for n, one_letter in enumerate(HV_Letters[i]):
                    if len(one_letter) == 1:
                        self.indices_HV[i].append([0]) if one_letter == 'H' else self.indices_HV[i].append([1])
                    else:
                        self.indices_HV[i].append([0,1])
        # reference_energies.indices[i]: # every site in endmember indice = [[[0],[1],[0,1],[1]], ...] in [0,1] have to check that alphabetic!!!!

def get_energies(structure, database_module):
    reference_energies, interaction_energies, two_letters_r, one_letters_r, two_letters_i, one_letters_i, ref_vars, int_vars, vars_count  = [], [], [], [], [], [], [], [], 0
    # Loop through all attributes in the module to find functions
    for name in dir(database_module):
        attribute = getattr(database_module, name)
        if isinstance(attribute, types.FunctionType):
            APPEND = True
        # Add the function names to the reference list
            if name[0] == 'G' and name[0:5] != 'GHSER':
                # Split the string by underscores
                split_string = name[1:].split('_')
                # Separate the two-letter substrings
                solid_letters = [s for s in split_string if len(s) == 2]
                HV_letters = [s for s in split_string if len(s) == 1]
                if len(solid_letters) == len(structure.site_fractions_solids) and len(HV_letters) == structure.n_hydrogen_sites:
                    for i, solid in enumerate(solid_letters):
                        if solid in structure.solids[i]:
                            continue
                        else:
                            APPEND = False
                            break
                    if APPEND:
                        params = list(inspect.signature(attribute).parameters) # check if other variables than T
                        variables = []
                        if  len(params) > 1:
                            for param in params[1:]:
                                vars_count += 1
                                variables.append(int(param[1:]))
                        ref_vars.append(variables) 
                        reference_energies.append(name) 
                        two_letters_r.append(solid_letters)
                        one_letters_r.append(HV_letters)
            elif name[0] == 'L':
                # Split the string by underscores
                split_string = name[1:].split('_')
                # Separate the two-letter substrings
                solid_letters = [s for s in split_string if (len(s) == 4) or ((len(s) == 2 and s != 'HV'))] 
                HV_letters = [s for s in split_string if (len(s) == 2 and s == 'HV') or (len(s) == 1 and s.isalpha())] # HV mixing sublattice
                # Multiply the corresponding site_fraction of pure metal sublattices
                if len(solid_letters) == len(structure.site_fractions_solids) and len(HV_letters) == structure.n_hydrogen_sites:
                    for n, solid in enumerate(solid_letters):
                        split_solid = [solid[i:i+2] for i in range(0, (len(solid)), 2)]
                        for two_letter in split_solid:
                            if two_letter in structure.solids[n]:
                                continue
                            else:
                                APPEND = False
                                break
                    if APPEND:
                        params = list(inspect.signature(attribute).parameters) # check if other variables than T
                        variables = []
                        if  len(params) > 1:
                            for param in params[1:]:
                                vars_count += 1
                                variables.append(int(param[1:]))
                        int_vars.append(variables) 
                        interaction_energies.append(name)
                        two_letters_i.append(solid_letters)
                        one_letters_i.append(HV_letters)

    reference_energies = Energies(structure, reference_energies, 'reference', one_letters_r, two_letters_r, ref_vars, vars_count) # initialize reference energy class
    interaction_energies = Energies(structure, interaction_energies, 'interaction', one_letters_i, two_letters_i, int_vars, vars_count) # initialize interaction energy class
    # just a little check of all endmembers being defined
    number_of_endmembers = 1
    for site in structure.site_fractions:
        number_of_endmembers *= len(site)
    if len(reference_energies.functions) < number_of_endmembers:
        print("CAUTION: Not all endmembers defined!")
    # reference gibbs energy
    return reference_energies, interaction_energies

def make_gibbs(T, structure, database_module, vars = None):
    reference_energies, interaction_energies = get_energies(structure, database_module) # initialize the energies for input of calculate_gibbs
    if vars is not None and reference_energies.vars_count != len(vars):
        print(f'WARNING: Amount of optimizing parameters ({len(vars)}) doesnt matches with database ({reference_energies.vars_count}) !')
    n_solid, multiplicity_H = len(structure.multiplicity_solids), structure.multiplicity_H
    Gref_array = []
    for i in range(len(reference_energies.functions)): # every endmember added 
        tmp = 1 
        for n, indice in enumerate(reference_energies.indices_solids[i]): # every site in endmember indice = [[0,1,1,0], [0,0,1,0], [1,0,1,0], [1,1,1,0]]
            tmp *= structure.site_fractions[n][indice]
        args = [T]
        if vars is not None and reference_energies.variables[i] != []: # check if its an optimization run --> has other input params than T
            args += [vars[v] for v in reference_energies.variables[i]]
        Gref_array.append(tmp * getattr(database_module, reference_energies.functions[i])(*args))
    # ideal mixing gibbs energy
    Gid_solid = 0
    for i, site in enumerate(structure.site_fractions): # only solids
        y_mix = 0
        for site_fraction in site:
            y_mix += site_fraction * np.log(max(1e-20, site_fraction))
        Gid_solid += R_CONST * T * (y_mix * structure.multiplicity[i])
    # excess energy
    Gex_array = []
    for i in range(len(interaction_energies.functions)): # every interaction function 
        tmp = 1
        for n, indices in enumerate(interaction_energies.indices_solids[i]): # every site in endmember indice = [[[0],[1],[1,1],[1]], ...]
            if len(indices) == 2:
                tmp *=  (structure.site_fractions[n][indices[0]] - structure.site_fractions[n][indices[1]]) ** interaction_energies.order[i]
            for indice in indices:
                tmp *= structure.site_fractions[n][indice]
        args = [T]
        if vars is not None and interaction_energies.variables[i] != []: # check if its an optimization run --> has other input params than T
            args += [vars[v] for v in interaction_energies.variables[i]]
        Gex_array.append(tmp * getattr(database_module, interaction_energies.functions[i])(*args))
    # n_atoms
    n_solids = sum(structure.multiplicity_solids) # amount of solid atoms
    def function(site_fractions_HV): #site_fractions_HV = [[0.2,0.8],[0,1]]
        Gref = 0
        for i, G in enumerate(Gref_array):
            tmp = 1 
            for n, indice in enumerate(reference_energies.indices_HV[i]): 
                tmp *= site_fractions_HV[n] if indice == 0 else (1 - site_fractions_HV[n])
            Gref += tmp * G
        # ideal mixing gibbs energy
        Gid = Gid_solid
        for i, site_fraction in enumerate(site_fractions_HV):
            y_mix = site_fraction * np.log(max(1e-20, site_fraction)) + (1 - site_fraction) * np.log(max(1e-20, 1 - site_fraction))
            Gid += R_CONST * T * (y_mix * multiplicity_H[i])
        # excess energy
        Gex = 0
        for i, G in enumerate(Gex_array): # every interaction function 
            tmp = 1
            for n, indices in enumerate(interaction_energies.indices_HV[i]): # every site in endmember indice = [[[0],[1],[1,1],[1]], ...]
                if len(indices) == 2:
                    tmp *=  (2 * site_fractions_HV[n] - 1) ** interaction_energies.order[i]
                for indice in indices:
                    tmp *= site_fractions_HV[n] if indice == 0 else (1 - site_fractions_HV[n])
            Gex += tmp * G
        n_atoms = n_solids
        for i, multiplicity in enumerate(structure.multiplicity_H):
            n_atoms += site_fractions_HV[i] * multiplicity
        return (Gref + Gid + Gex) / (n_atoms)
    
    Gref_derarray = [] # DESC: first derivative of reference energy: partial derivative for each hydrogen member 
    # FORMAT: {n_endmember x n_HV_sites} [[G_EM1_PD1, GEM2_PD1 ...],[G_EM1_PD2, ...],...]
    for der in range(len(multiplicity_H)): # -> corresponding to each partial derivative of Gref
        Gref_derarray.append([])
        for i, G in enumerate(Gref_array):
            if reference_energies.indices_HV[i][der] == 1: # VA in the sublattice to derive partially
                Gref_derarray[der].append(-G) # partial derivate of (1 - y_H) * G is -1 * G
            else: # H in the sublattice to derive partially
                Gref_derarray[der].append(G) # partial derivate of (y_H) * G is G
    L_derivative = [ # choose the correct analytical expression depending on interaction order
        lambda x: (1 - 2 * x), # f´(mx(1-x)(x-(1-x))**0) = m*(1-2x)
        lambda x: (- 6 * x**2 + 6 * x - 1), # f´(mx(1-x)(x-(1-x))**1) = m*(-6x^2+6x-1)
        lambda x: (- 16 * x**3 + 24 * x**2 - 10 * x + 1) # f´(mx(1-x)(x-(1-x))**2) = m*(-16x^3+24x^2-10x+1)
    ]
    def derivative(site_fractions_HV):
        grads = []
        for der, Gref_parderarray in enumerate(Gref_derarray): # for each partial derivative: der = PD_der (= 0, 1, 2 ...), Gref_parderarray = for each endmember the precalculated G value
            grads.append(0)
            for i, Gibbs in enumerate(Gref_parderarray): # i is the ith endmember for the current PD der
                tmp = 1 
                for n, indice in enumerate(reference_energies.indices_HV[i]): 
                    if n != der:
                        tmp *= site_fractions_HV[n] if indice == 0 else (1 - site_fractions_HV[n])
                grads[der] += tmp * Gibbs
            # ideal mixing gibbs energy: as Gid = y1*lny1 + ... + .. every summation terms cancels out beside if sublattice = der
            # f´(ln(x)*x + (1-x)*ln(1-x)) = ln(x) + 1 - ln(1-x) - 1 = ln(x) - ln(1-x)
            y_mix = np.log(max(1e-20, site_fractions_HV[der])) - np.log(max(1e-20, 1 - site_fractions_HV[der]))
            grads[der] += R_CONST * T * (y_mix * multiplicity_H[der])
            # excess energy: dont just multiply the base Gex_array with the y or y*(1-y) or (y-1-y)**v but their respective derivates
            # f´(mx) = m; f´(m(1-x)) = -m; f´(mx(1-x)(x-(1-x))**0) = m*(1-2x); f´(mx(1-x)(x-(1-x))**1) = m*(-6x^2+6x-1); f´(mx(1-x)(x-(1-x))**2) = m*(-16x^3+24x^2-10x+1); 
            for i, Gibbs in enumerate(Gex_array): # every interaction function 
                tmp = 1
                for n, indices in enumerate(interaction_energies.indices_HV[i]): # every site in interaction indice = [[[0],[1],[1,1],[1]], ...]
                    if n == der: # if its the site to derivate
                        if len(indices) == 2: # if it is the interaction site
                            tmp *= L_derivative[interaction_energies.order[i]](site_fractions_HV[n]) # choose the correct analytical expression depending on interaction order
                        else: # if its a non-interaction site
                            if indices == 1: # f´(m(1-x)) = -m
                                tmp *= -1 
                    else: # if its the non derivate site
                        if len(indices) == 2:
                            tmp *=  (2 * site_fractions_HV[n] - 1) ** interaction_energies.order[i]
                        for indice in indices:
                            tmp *= site_fractions_HV[n] if indice == 0 else (1 - site_fractions_HV[n])

                grads[der] += tmp * Gibbs
        n_atoms = n_solids
        for i, multiplicity in enumerate(structure.multiplicity_H):
            n_atoms += site_fractions_HV[i] * multiplicity
        # transfer grads to molar grads
        for i in range(len(grads)):
            grads[i] = grads[i] / n_atoms / structure.multiplicity_H[i]
        return grads
    return function, derivative

# Define the constraint function
def constraint(y, m1, m2, y_total):
    y1, y2 = y
    return (m1 * y1 + m2 * y2) / (m1 + m2) - y_total

#def make_objective(derivative, m1, m2, y_total):
#    def objective_function(y1):
#        y2 = ((m1 + m2) * y_total - m1 * y1) / m2
#
#        grad = derivative([y1, y2])
#        return (grad[0] - grad[1])**2
#    return objective_function

def make_objective(calculate, m, y_total):
    def objective_function(y1):
        y2 = ((m[0] + m[1]) * y_total - m[0] * y1) / m[1]

        return calculate([y1, y2])
    return objective_function

def make_minimize_objective(calculate, m, y_total):
    def minimize_objective_function(y):
        #sum_y = sum(y)
        # calculate sum(m[i] * y[i])
        sum_m_times_y = sum(m[i] * y[i] for i in range(len(m) - 1))
        # TODO is check needed if initial guess is correct?
        # rest = m[-1] * y_last
        rest = (sum(m) * y_total - sum_m_times_y)
        # if rest < 0 --> y_last < 0, if rest > m[-1] --> y_last > 1
        #if rest < 0:
        #    #print(y)
        #    #print(rest)
        #    print(y_total)
        #    return 1e9 - 1e9 * rest
        #elif rest > m[-1]:
        #    print("D")
        #    return 1e9 + 1e9 * (rest - m[-1])
        y_last = rest / m[-1]
        y = np.append(y, y_last)
        return calculate(y) # most negative
    return minimize_objective_function

def make_global_objective(calculate, m, y_total):
    def global_objective_function(y):

        sum_m_times_y = sum(m[i] * y[i] for i in range(len(m) - 1))

        rest = (sum(m) * y_total - sum_m_times_y)
        if rest < 0 or rest > m[-1]:
            return 1e9
        y_last = rest / m[-1]
        y = np.append(y, y_last)
        return calculate(y) # most negative
    return global_objective_function

def iterate_over_bounds(bounds, num_points=3):
    # Generate linspace arrays for each bound
    linspace_arrays = [np.linspace(bound[0], bound[1], num_points) for bound in bounds]

    # Use itertools.product to create a cartesian product of these arrays
    for point in product(*linspace_arrays):
        yield point


def gibbs_minimizer(T, structure, vars=None, steps = 100):
    database_module = importlib.import_module(structure.database_module)  # Dynamically import the module
    # faster calculation
    n_solids = sum(structure.multiplicity_solids)
    site_fractions_HV = [] # use array instead of class
    m = structure.multiplicity_H # faster referencing
    sum_mH = np.sum(m)
    for i in range(structure.n_hydrogen_sites):
        site_fractions_HV.append(0)
    # create functions
    calculate, derivative = make_gibbs(T, structure, database_module, vars)
    # Parameters
    n_steps = steps
    def x_from_y(y): # calculate moles from site fractions
        return (sum_mH * y + structure.hydrogen_start_value) / (n_solids + sum_mH * y)
    # this function excludes the non HV interaction Hydrogen sites (full of hydrogen)
    def x_from_y_eq_spaced(y):
        return (sum_mH * y) / (n_solids + sum_mH * y)
    # save results
    g_total, x_total = [calculate(site_fractions_HV)], [x_from_y(0)] 
    fractions = [np.zeros(len(m))]
    # Generate equally spaced x values, instead of equally spaced y values
    def y_from_x(x): # calculate site fractions from moles
        return x * n_solids / (sum_mH * (1 - x))
    x_equally_spaced = np.linspace(x_from_y_eq_spaced(0), x_from_y_eq_spaced(1), n_steps) # Generate equally spaced x values
    y_for_equally_spaced_x = y_from_x(x_equally_spaced) # Calculate the corresponding y values that would give us equally spaced x values
    ROUND = False # if initial guess = bound_max
    RE = False
    result = None
    derivate_old = 0
    # calculate internal equilibria
    if len(m) == 2: # if more than one hydrogen sublattice: Have to calculate internal equilibrium
        for y_t in y_for_equally_spaced_x:
            bound_max = min(sum_mH / m[0] * y_t, 1)
            bound_min = -(min(sum_mH / m[1] * y_t, 1) * m[1] - sum_mH * y_t) / m[0]
            if bound_min < 1e-8:
                bound_min = 0
            #if result is not None:
            #    if ROUND == True:
            #        initial_guess = [bound_max]
            #    else:
            #        initial_guess = [result.x[0]]
            #else:
            #    initial_guess = [0]
            if result is None:
                initial_guess = [0]
            #    old_result = 0
            #bounds = (bound_min, bound_max)  # y1 and y2 must be between 0 and 1 or y_t
            bounds = [(bound_min, bound_max)] 
            # Create an initial population array
            initial_population = np.empty((8, 1))
            # Set the first member of the population to the initial guess
            initial_population[0] = initial_guess
            # Fill the rest of the initial population with variations of the initial guess
            for i in range(1, 8):
                for j, (lower, upper) in enumerate(bounds):
                    variation = np.random.uniform(lower, upper)
                    initial_population[i, j] = variation
            ## Ensure the values are within the boundaries
            #initial_population = np.clip(initial_population, [b[0] for b in bounds], [b[1] for b in bounds])
            #if result is not None:
            #    initial_guess = np.clip(initial_guess, bound_min, bound_max)
            objective_function = make_objective(calculate, m, y_t)
            result = differential_evolution(objective_function, bounds, workers = 1, popsize = 8, maxiter=4)
            initial_guess = result.x
            #result = minimize(objective_function, initial_guess, bounds=bounds, method='L-BFGS-B')
            #derivate = result.x - old_result
            #if derivate_old != 0:
            #    st = max(derivate_old,derivate)/min(derivate_old,derivate)
            #    if st > 10 or st < 0.1:
            #        print("''''''''''''''#####################################################")
            #        print(result.x)
            #        result = differential_evolution(objective_function, bounds, workers = 1, popsize = 8, maxiter=4)
            #        print(result.x)
            #        derivate = result.x - old_result
            #initial_guess = derivate + result.x
            #print("old_result", old_result, "result", result.x, "initial_guess", initial_guess)
            #print(result.x - old_result)
            #old_result = result.x
            #derivate_old = derivate
                #print(f"HEYYYY {y_t}: result.nfev = {result.nfev}, result = {result.fun}, result.x = {result.x}")
            
            #print(f"{y_t}: result.nfev = {result.nfev}, result = {result.fun}, result.x = {result.x}")
            
            #objective_function = make_objective(derivative, m[0], m[1], y_t)
            #result = minimize(objective_function, initial_guess, bounds=bounds, method='L-BFGS-B')
            #for i in np.linspace(bound_min, bound_max, 10):
            #    if objective_function(i) < result.fun:
            #        #print(i)
            #        result = minimize(objective_function, [i], bounds=bounds, method='L-BFGS-B')
            #if result.x > bound_max - 1e-6:
            #    result.x = [bound_max]
            #    ROUND = True
            #else:
            #    ROUND = False
            site_fractions_HV[0] = result.x[0]
            site_fractions_HV[1] = ((m[0] + m[1]) * y_t - m[0] * site_fractions_HV[0]) / m[1]
            g_total.append(calculate(site_fractions_HV))
            x_total.append(x_from_y(y_t))
            #delete for faster calculation, only for debugging
            fractions.append(list(site_fractions_HV))
    elif len(m) > 2: # if more than one hydrogen sublattice: Have to calculate internal equilibrium
        de = 0
        for y_t in y_for_equally_spaced_x:
            de += 1
            bounds = []
            for site in range(len(m) - 1): # create bound_min, max for len(m) - 1 variables, last variable gets calculated through y_total
                bound_max = min(sum_mH / m[site] * y_t, 1) # y1 < sum(m) / m1 * y_total (e.g. y1 < 7/1 * y_total)
                bound_min = -(min(sum_mH / (sum_mH - m[site]) * y_t, 1) * (sum_mH - m[site]) - sum_mH * y_t) / m[site] # if y_total small --> min(...) < 1 --> sum_mH * y_t - sum_mH * y_t = 0
                # if y_t big: min = 1 --> y1 > (sum(m)*y_total - 1) / m1
                if bound_min < 1e-7:
                    bound_min = 0
                if bound_min > bound_max:
                    bound_min = bound_max
                bounds.append((bound_min, bound_max))
            if result is not None:
                initial_guess = result.x
            else:
                initial_guess = np.zeros(len(m) - 1)
            #if de % 5 == 0:
            objective_function = make_global_objective(calculate, m, y_t)
            result = differential_evolution(objective_function, bounds, workers = 1, popsize = 8)
            initial_guess = result.x
                #print(f"HEYYYY {y_t}: result.nfev = {result.nfev}, result = {result.fun}, result.x = {result.x}")
            objective_function = make_minimize_objective(calculate, m, y_t)
            result = minimize(objective_function, initial_guess, bounds=bounds, method='L-BFGS-B')
            #print(f"{y_t}: result.nfev = {result.nfev}, result = {result.fun}, result.x = {result.x}")
            site_fractions_HV[:-1] = result.x
            sum_m_times_y = sum(m[i] * result.x[i] for i in range(len(m) - 1))
            site_fractions_HV[-1] = (sum(m) * y_t - sum_m_times_y) / m[-1]
            g_total.append(calculate(site_fractions_HV))
            x_total.append(x_from_y(y_t))
            fractions.append(list(site_fractions_HV))
    else: # here normal calucaltion 
        for y_t in y_for_equally_spaced_x:
            site_fractions_HV = [y_t]
            g_total.append(calculate(site_fractions_HV))
            x_total.append(x_from_y(y_t))
            # delete for faster calculation, only for debugging
            fractions.append(list(site_fractions_HV))
    structure.gm = g_total[1:]
    structure.x_H = x_total[1:]
    #return  g_total, x_total, fractions

