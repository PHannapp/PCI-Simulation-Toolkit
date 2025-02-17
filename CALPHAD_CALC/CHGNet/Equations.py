import database  # Import custom database module for energy calculations
import numpy as np  # Import NumPy for numerical operations
import types  # Import types module to work with different types
import inspect  # Import inspect module to inspect function signatures
from scipy.optimize import minimize, differential_evolution  # Import optimization functions
from itertools import product  # Import product for cartesian product calculations

R_CONST = 8.3144621  # Define the gas constant in J/(mol*K)

class Structure:
    """
    Represents the structure of a compound, encapsulating information about its solids, sites, and site fractions.
    """
    def __init__(self, solids, multiplicity, site_fractions):
        """
        Initializes the Structure instance with given parameters.
        
        Parameters:
            solids (list): A list of solid phases in the structure.
            multiplicity (list): A list indicating the multiplicity of each site.
            site_fractions (list): A list of lists where each sublist represents the fraction of each solid in a site.
        Raises:
            ValueError: If the sum of site fractions for any site does not equal 1.
        """
        self.solids = solids
        self.multiplicity = multiplicity
        self.site_fractions_solids = site_fractions.copy()
        self.n_hydrogen_sites = len(multiplicity) - len(site_fractions)
        self.site_fractions_H = []
        self.site_fractions = site_fractions.copy()
        self.multiplicity_H = multiplicity[-self.n_hydrogen_sites :]
        self.multiplicity_solids = multiplicity[: len(site_fractions)]

        # Validate site fractions sum to 1 for each site
        for i, site in enumerate(self.site_fractions):
            if sum(site) != 1:
                raise ValueError("Sum in site ", i, "is not == 1")


class Energies:
    """
    Represents the energy calculations for a given structure, including reference and interaction energies.
    """
    def __init__(
        self, structure, functions, type, HV_Letters, solid_letters, vars, vars_count
    ):
        """
        Initializes a new instance of the Energies class.

        Parameters:
            structure (Structure): The structure for which energies are being calculated.
            functions (list): List of function names corresponding to energy calculations.
            type (str): The type of energies being calculated ('reference' or 'interaction').
            HV_Letters (list): List of letters representing hydrogen and vacancies.
            solid_letters (list): List of letters representing solid components.
            vars (list): List of variables used in the energy calculations.
            vars_count (int): The number of variables used in the energy calculations.
        """
        self.type = type
        self.functions = functions
        self.indices_solids = []  # Indices for solid components in the structure.
        self.indices_HV = []  # Indices for hydrogen and vacancies in the structure.
        self.order = []  # Order of interaction for interaction energies.
        self.HV_Letters = HV_Letters
        self.solid_letters = solid_letters
        self.variables = vars
        self.vars_count = vars_count
        
        # Initialize reference energies by mapping solid component letters to their indices in the structure.
        if self.type == "reference":
            for i, two_letter in enumerate(solid_letters):
                self.indices_solids.append([])
                self.indices_HV.append([])
                for n, two_let in enumerate(two_letter):
                    self.indices_solids[i].append(structure.solids[n].index(two_let))
                for n, one_letter in enumerate(HV_Letters[i]):
                    self.indices_HV[i].append(
                        0
                    ) if one_letter == "H" else self.indices_HV[i].append(1)

        # Initialize interaction energies by mapping both solid component letters and interaction orders.
        if self.type == "interaction":
            for i, function in enumerate(functions):
                self.order.append(int(function[-1]))
                self.indices_solids.append([])
                self.indices_HV.append([])
                for n, solid_letter in enumerate(solid_letters[i]):
                    if len(solid_letter) == 2:
                        self.indices_solids[i].append(
                            [structure.solids[n].index(solid_letter)]
                        )
                    else:
                        self.indices_solids[i].append(
                            [
                                structure.solids[n].index(m)
                                for m in sorted(
                                    [
                                        solid_letter[k : k + 2]
                                        for k in range(0, (len(solid_letter)), 2)
                                    ]
                                )
                            ]
                        )
                for n, one_letter in enumerate(HV_Letters[i]):
                    if len(one_letter) == 1:
                        self.indices_HV[i].append(
                            [0]
                        ) if one_letter == "H" else self.indices_HV[i].append([1])
                    else:
                        self.indices_HV[i].append([0, 1])

def get_energies(structure):
    """
    Retrieves reference and interaction energies for a given structure.
    
    Parameters:
        structure (Structure): The structure for which energies are being calculated.
    
    Returns:
        tuple: A tuple containing two Energies instances for reference and interaction energies, respectively.
    """
    (
        reference_energies,
        interaction_energies,
        two_letters_r,
        one_letters_r,
        two_letters_i,
        one_letters_i,
        ref_vars,
        int_vars,
        vars_count,
    ) = [], [], [], [], [], [], [], [], 0
    # Loop through all attributes in the module to find functions
    for name in dir(database):
        attribute = getattr(database, name)
        if isinstance(attribute, types.FunctionType):
            APPEND = True
            # Add the function names to the reference list
            if name[0] == "G" and name[0:5] != "GHSER":
                # Split the string by underscores
                split_string = name[1:].split("_")
                # Separate the two-letter substrings
                solid_letters = [s for s in split_string if len(s) == 2]
                HV_letters = [s for s in split_string if len(s) == 1]
                if (
                    len(solid_letters) == len(structure.site_fractions_solids)
                    and len(HV_letters) == structure.n_hydrogen_sites
                ):
                    for i, solid in enumerate(solid_letters):
                        if solid in structure.solids[i]:
                            continue
                        else:
                            APPEND = False
                            break
                    if APPEND:
                        params = list(
                            inspect.signature(attribute).parameters
                        )  # check if other variables than T
                        variables = []
                        if len(params) > 1:
                            for param in params[1:]:
                                vars_count += 1
                                variables.append(int(param[1:]))
                        ref_vars.append(variables)
                        reference_energies.append(name)
                        two_letters_r.append(solid_letters)
                        one_letters_r.append(HV_letters)
            elif name[0] == "L":
                # Split the string by underscores
                split_string = name[1:].split("_")
                # Separate the two-letter substrings
                solid_letters = [
                    s
                    for s in split_string
                    if (len(s) == 4) or (len(s) == 2 and s != "HV")
                ]
                HV_letters = [
                    s
                    for s in split_string
                    if (len(s) == 2 and s == "HV") or (len(s) == 1 and s.isalpha())
                ]  # HV mixing sublattice
                # Multiply the corresponding site_fraction of pure metal sublattices
                if (
                    len(solid_letters) == len(structure.site_fractions_solids)
                    and len(HV_letters) == structure.n_hydrogen_sites
                ):
                    for n, solid in enumerate(solid_letters):
                        split_solid = [
                            solid[i : i + 2] for i in range(0, (len(solid)), 2)
                        ]
                        for two_letter in split_solid:
                            if two_letter in structure.solids[n]:
                                continue
                            else:
                                APPEND = False
                                break
                    if APPEND:
                        params = list(
                            inspect.signature(attribute).parameters
                        )  # check if other variables than T
                        variables = []
                        if len(params) > 1:
                            for param in params[1:]:
                                vars_count += 1
                                variables.append(int(param[1:]))
                        int_vars.append(variables)
                        interaction_energies.append(name)
                        two_letters_i.append(solid_letters)
                        one_letters_i.append(HV_letters)

    reference_energies = Energies(
        structure,
        reference_energies,
        "reference",
        one_letters_r,
        two_letters_r,
        ref_vars,
        vars_count,
    )  # initialize reference energy class
    interaction_energies = Energies(
        structure,
        interaction_energies,
        "interaction",
        one_letters_i,
        two_letters_i,
        int_vars,
        vars_count,
    )  # initialize interaction energy class
    # just a little check of all endmembers being defined
    number_of_endmembers = 1
    for site in structure.site_fractions:
        number_of_endmembers *= len(site)
    if len(reference_energies.functions) < number_of_endmembers:
        print("CAUTION: Not all endmembers defined!")
    # reference gibbs energy
    return reference_energies, interaction_energies


def make_gibbs(T, structure, vars=None):
    reference_energies, interaction_energies = get_energies(
        structure
    )  # initialize the energies for input of calculate_gibbs
    if vars is not None and reference_energies.vars_count != len(vars):
        print(
            f"WARNING: Amount of optimizing parameters ({len(vars)}) doesnt matches with database ({reference_energies.vars_count}) !"
        )
    multiplicity_H = structure.multiplicity_H
    Gref_array = []
    for i in range(len(reference_energies.functions)):  # every endmember added
        tmp = 1
        for n, indice in enumerate(
            reference_energies.indices_solids[i]
        ):  # every site in endmember indice = [[0,1,1,0], [0,0,1,0], [1,0,1,0], [1,1,1,0]]
            tmp *= structure.site_fractions[n][indice]
        args = [T]
        if (
            vars is not None and reference_energies.variables[i] != []
        ):  # check if its an optimization run --> has other input params than T
            args += [vars[v] for v in reference_energies.variables[i]]
        Gref_array.append(
            tmp * getattr(database, reference_energies.functions[i])(*args)
        )
    # ideal mixing gibbs energy
    Gid_solid = 0
    for i, site in enumerate(structure.site_fractions):  # only solids
        y_mix = 0
        for site_fraction in site:
            y_mix += site_fraction * np.log(max(1e-20, site_fraction))
        Gid_solid += R_CONST * T * (y_mix * structure.multiplicity[i])
    # excess energy
    Gex_array = []
    for i in range(len(interaction_energies.functions)):  # every interaction function
        tmp = 1
        for n, indices in enumerate(
            interaction_energies.indices_solids[i]
        ):  # every site in endmember indice = [[[0],[1],[1,1],[1]], ...]
            if len(indices) == 2:
                tmp *= (
                    structure.site_fractions[n][indices[0]]
                    - structure.site_fractions[n][indices[1]]
                ) ** interaction_energies.order[i]
            for indice in indices:
                tmp *= structure.site_fractions[n][indice]
        args = [T]
        if (
            vars is not None and interaction_energies.variables[i] != []
        ):  # check if its an optimization run --> has other input params than T
            args += [vars[v] for v in interaction_energies.variables[i]]
        Gex_array.append(
            tmp * getattr(database, interaction_energies.functions[i])(*args)
        )
    # n_atoms
    n_solids = sum(structure.multiplicity_solids)  # amount of solid atoms

    def function(site_fractions_HV):  # site_fractions_HV = [[0.2,0.8],[0,1]]
        Gref = 0
        for i, G in enumerate(Gref_array):
            tmp = 1
            for n, indice in enumerate(reference_energies.indices_HV[i]):
                tmp *= (
                    site_fractions_HV[n] if indice == 0 else (1 - site_fractions_HV[n])
                )
            Gref += tmp * G
        # ideal mixing gibbs energy
        Gid = Gid_solid
        for i, site_fraction in enumerate(site_fractions_HV):
            y_mix = site_fraction * np.log(max(1e-20, site_fraction)) + (
                1 - site_fraction
            ) * np.log(max(1e-20, 1 - site_fraction))
            Gid += R_CONST * T * (y_mix * multiplicity_H[i])
        # excess energy
        Gex = 0
        for i, G in enumerate(Gex_array):  # every interaction function
            tmp = 1
            for n, indices in enumerate(
                interaction_energies.indices_HV[i]
            ):  # every site in endmember indice = [[[0],[1],[1,1],[1]], ...]
                if len(indices) == 2:
                    tmp *= (2 * site_fractions_HV[n] - 1) ** interaction_energies.order[i]
                for indice in indices:
                    tmp *= (
                        site_fractions_HV[n]
                        if indice == 0
                        else (1 - site_fractions_HV[n])
                    )
            Gex += tmp * G
        n_atoms = n_solids
        for i, multiplicity in enumerate(structure.multiplicity_H):
            n_atoms += site_fractions_HV[i] * multiplicity
        return (Gref + Gid + Gex) / (n_atoms)


    Gref_derarray = []  # DESC: first derivative of reference energy: partial derivative for each hydrogen member
    # FORMAT: {n_endmember x n_HV_sites} [[G_EM1_PD1, GEM2_PD1 ...],[G_EM1_PD2, ...],...]
    for der in range(
        len(multiplicity_H)
    ):  # -> corresponding to each partial derivative of Gref
        Gref_derarray.append([])
        for i, G in enumerate(Gref_array):
            if (
                reference_energies.indices_HV[i][der] == 1
            ):  # VA in the sublattice to derive partially
                Gref_derarray[der].append(
                    -G
                )  # partial derivate of (1 - y_H) * G is -1 * G
            else:  # H in the sublattice to derive partially
                Gref_derarray[der].append(G)  # partial derivate of (y_H) * G is G
    L_derivative = [  # choose the correct analytical expression depending on interaction order
        lambda x: (1 - 2 * x),  # f´(mx(1-x)(x-(1-x))**0) = m*(1-2x)
        lambda x: (-6 * x**2 + 6 * x - 1),  # f´(mx(1-x)(x-(1-x))**1) = m*(-6x^2+6x-1)
        lambda x: (
            -16 * x**3 + 24 * x**2 - 10 * x + 1
        ),  # f´(mx(1-x)(x-(1-x))**2) = m*(-16x^3+24x^2-10x+1)
    ]

    def derivative(site_fractions_HV):
        grads = []
        for der, Gref_parderarray in enumerate(
            Gref_derarray
        ):  # for each partial derivative: der = PD_der (= 0, 1, 2 ...), Gref_parderarray = for each endmember the precalculated G value
            grads.append(0)
            for i, Gibbs in enumerate(
                Gref_parderarray
            ):  # i is the ith endmember for the current PD der
                tmp = 1
                for n, indice in enumerate(reference_energies.indices_HV[i]):
                    if n != der:
                        tmp *= (
                            site_fractions_HV[n]
                            if indice == 0
                            else (1 - site_fractions_HV[n])
                        )
                grads[der] += tmp * Gibbs
            # ideal mixing gibbs energy: as Gid = y1*lny1 + ... + .. every summation terms cancels out beside if sublattice = der
            # f´(ln(x)*x + (1-x)*ln(1-x)) = ln(x) + 1 - ln(1-x) - 1 = ln(x) - ln(1-x)
            y_mix = np.log(max(1e-20, site_fractions_HV[der])) - np.log(
                max(1e-20, 1 - site_fractions_HV[der])
            )
            grads[der] += R_CONST * T * (y_mix * multiplicity_H[der])
            # excess energy: dont just multiply the base Gex_array with the y or y*(1-y) or (y-1-y)**v but their respective derivates
            # f´(mx) = m; f´(m(1-x)) = -m; f´(mx(1-x)(x-(1-x))**0) = m*(1-2x); f´(mx(1-x)(x-(1-x))**1) = m*(-6x^2+6x-1); f´(mx(1-x)(x-(1-x))**2) = m*(-16x^3+24x^2-10x+1);
            for i, Gibbs in enumerate(Gex_array):  # every interaction function
                tmp = 1
                for n, indices in enumerate(
                    interaction_energies.indices_HV[i]
                ):  # every site in interaction indice = [[[0],[1],[1,1],[1]], ...]
                    if n == der:  # if its the site to derivate
                        if len(indices) == 2:  # if it is the interaction site
                            tmp *= L_derivative[
                                interaction_energies.order[i]
                            ](
                                site_fractions_HV[n]
                            )  # choose the correct analytical expression depending on interaction order
                        else:  # if its a non-interaction site
                            if indices == 1:  # f´(m(1-x)) = -m
                                tmp *= -1
                    else:  # if its the non derivate site
                        if len(indices) == 2:
                            tmp *= (
                                2 * site_fractions_HV[n] - 1
                            ) ** interaction_energies.order[i]
                        for indice in indices:
                            tmp *= (
                                site_fractions_HV[n]
                                if indice == 0
                                else (1 - site_fractions_HV[n])
                            )

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


def make_objective(calculate, m, y_total):
    def objective_function(y1):
        y2 = ((m[0] + m[1]) * y_total - m[0] * y1) / m[1]

        return calculate([y1, y2])

    return objective_function


def make_minimize_objective(calculate, m, y_total):
    def minimize_objective_function(y):
        # sum_y = sum(y)
        # calculate sum(m[i] * y[i])
        sum_m_times_y = sum(m[i] * y[i] for i in range(len(m) - 1))
        rest = sum(m) * y_total - sum_m_times_y
        y_last = rest / m[-1]
        y = np.append(y, y_last)
        return calculate(y)  # most negative

    return minimize_objective_function


def make_global_objective(calculate, m, y_total):
    def global_objective_function(y):
        sum_m_times_y = sum(m[i] * y[i] for i in range(len(m) - 1))

        rest = sum(m) * y_total - sum_m_times_y
        if rest < 0 or rest > m[-1]:
            return 1e9
        y_last = rest / m[-1]
        y = np.append(y, y_last)
        return calculate(y)  # most negative

    return global_objective_function


def iterate_over_bounds(bounds, num_points=3):
    # Generate linspace arrays for each bound
    linspace_arrays = [np.linspace(bound[0], bound[1], num_points) for bound in bounds]

    # Use itertools.product to create a cartesian product of these arrays
    for point in product(*linspace_arrays):
        yield point


def gibbs_minimizer(T, structure, vars=None, steps=100):
    # faster calculation
    n_solids = sum(structure.multiplicity_solids)
    site_fractions_HV = []  # use array instead of class
    m = structure.multiplicity_H  # faster referencing
    sum_mH = np.sum(m)
    for i in range(structure.n_hydrogen_sites):
        site_fractions_HV.append(0)
    # create functions
    calculate, derivative = make_gibbs(T, structure, vars)
    # Parameters
    n_steps = steps
    # save results
    g_total, x_total = [], []
    fractions = []

    # Generate equally spaced x values, instead of equally spaced y values
    def x_from_y(y):  # calculate moles from site fractions
        return sum(m) * y / (n_solids + sum_mH * y)

    def y_from_x(x):  # calculate site fractions from moles
        return x * n_solids / (sum_mH * (1 - x))

    x_equally_spaced = np.linspace(
        x_from_y(0), x_from_y(1), n_steps
    )  # Generate equally spaced x values
    y_for_equally_spaced_x = y_from_x(
        x_equally_spaced
    )  # Calculate the corresponding y values that would give us equally spaced x values
    result = None
    # calculate internal equilibria
    if (
        len(m) == 2
    ):  # if more than one hydrogen sublattice: Have to calculate internal equilibrium
        for y_t in y_for_equally_spaced_x:
            bound_max = min(sum_mH / m[0] * y_t, 1)
            bound_min = -(min(sum_mH / m[1] * y_t, 1) * m[1] - sum_mH * y_t) / m[0]
            if bound_min < 1e-8:
                bound_min = 0
            if result is None:
                initial_guess = [0]
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
            min_result = 1e20
            final_result = 0
            for _ in range(5):
                objective_function = make_objective(calculate, m, y_t)
                result = differential_evolution(
                    objective_function, bounds, workers=1, popsize=8, maxiter=400
                )
                if result.fun < min_result:
                    final_result = result.x[0]
                    min_result = result.fun
            initial_guess = result.x
            site_fractions_HV[0] = final_result
            site_fractions_HV[1] = (
                (m[0] + m[1]) * y_t - m[0] * site_fractions_HV[0]
            ) / m[1]
            g_total.append(calculate(site_fractions_HV))
            x_total.append(x_from_y(y_t))
            # delete for faster calculation, only for debugging
            fractions.append(list(site_fractions_HV))
    elif (
        len(m) > 2
    ):  # if more than one hydrogen sublattice: Have to calculate internal equilibrium
        de = 0
        for y_t in y_for_equally_spaced_x:
            de += 1
            bounds = []
            for site in range(
                len(m) - 1
            ):  # create bound_min, max for len(m) - 1 variables, last variable gets calculated through y_total
                bound_max = min(
                    sum_mH / m[site] * y_t, 1
                )  # y1 < sum(m) / m1 * y_total (e.g. y1 < 7/1 * y_total)
                bound_min = (
                    -(
                        min(sum_mH / (sum_mH - m[site]) * y_t, 1) * (sum_mH - m[site])
                        - sum_mH * y_t
                    )
                    / m[site]
                )  # if y_total small --> min(...) < 1 --> sum_mH * y_t - sum_mH * y_t = 0
                if bound_min < 1e-7:
                    bound_min = 0
                bounds.append((bound_min, bound_max))
            if result is not None:
                initial_guess = result.x
            else:
                initial_guess = np.zeros(len(m) - 1)
            objective_function = make_global_objective(calculate, m, y_t)
            result = differential_evolution(
                objective_function, bounds, workers=1, popsize=8
            )
            initial_guess = result.x
            objective_function = make_minimize_objective(calculate, m, y_t)
            result = minimize(
                objective_function, initial_guess, bounds=bounds, method="L-BFGS-B"
            )
            site_fractions_HV[:-1] = result.x
            sum_m_times_y = sum(m[i] * result.x[i] for i in range(len(m) - 1))
            site_fractions_HV[-1] = (sum(m) * y_t - sum_m_times_y) / m[-1]
            g_total.append(calculate(site_fractions_HV))
            x_total.append(x_from_y(y_t))
            fractions.append(list(site_fractions_HV))
    else:  # here normal calucaltion
        for y_t in y_for_equally_spaced_x:
            site_fractions_HV = [y_t]
            g_total.append(calculate(site_fractions_HV))
            x_total.append(x_from_y(y_t))
            fractions.append(list(site_fractions_HV))
    return g_total, x_total, fractions
