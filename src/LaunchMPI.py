import sys
import imp
import pyipopt
import logging
import time
import os
import numpy as np
from mpi4py import MPI

import cantera as ct

import Global
import Parallelism as Par
import Slave
import Parameters as Param
import OptimFunctions as OptFn
import ReactorCstVolume as RCV
import OneDFlame as Flame


def eval_jac_g_wrapper(x, flag, user_data=None):
    noptim = Global.params.noptim
    if flag:
        logging.debug('Jacobian Initialization')
        tab = []
        for itcon in range(ncon):
            for k in range(noptim):
                tab.append(itcon)
        tab2 = []
        for itcon in range(ncon):
            for k in range(noptim):
                tab2.append(k)
        return (np.array(tab), np.array(tab2))
    else:
        return OptFn.gradient_constraints(x)


if __name__ == '__main__':
    # Main process, common to all tasks (master and slaves)
    # Reading input
    up = imp.load_source('user_param', sys.argv[1])

    # Cantera gas phase
    Global.gas = ct.Solution(up.mechanism)

    # Specific storage for flame computations
    Global.nflames = len(up.P_fl)*len(up.T_fl)*len(up.phi_fl)
    Global.flame_inits = [False for ii in range(Global.nflames)]
    Global.sims = [False for ii in range(Global.nflames)]
    Global.gases = [False for ii in range(Global.nflames)]
    Global.directory = up.directory
    Global.mechanism = up.mechanism

    # Define MPI message tags
    tags = Par.enum('READY', 'DONE', 'EXIT', 'START', 'SLEEP', 'WAKEUP')

    # Initializations and preliminaries
    Global.MPI = MPI
    comm = MPI.COMM_WORLD   # get MPI communicator object
    size = comm.size        # total number of processes
    rank = comm.rank        # rank of this process
    status = MPI.Status()   # get MPI status object

    Global.partitionid = -1

    if rank > 0:
        # Slave Process (waiting for task from Master)
        Slave.slave_execution()
    else:
        # Main Process
        if not os.path.exists(up.directory):
            os.makedirs(up.directory)

        logfile = 'out.log'
        verbosity = 20
        if hasattr(up, 'verbosity'):
            if up.verbosity == 'DEBUG':
                verbosity = 10
            elif up.verbosity == 'INFO':
                verbosity = 20
        logging.basicConfig(filename=logfile, level=verbosity)
        
        logging.info('Start')
        # Initializing gas phase and damping
        gas = Global.gas

        reactants = gas.reactant_stoich_coeffs()
        products = gas.product_stoich_coeffs()

        # Species index build
        species_index_exclude_init = []
        for species in up.species_exclude_init:
            species_index_exclude_init.append(gas.species_index(species))

        # species_index_exclude_zero = []
        # for species in up.species_exclude_zero:
        #     species_index_exclude_zero.append(gas.species_index(species))

        species_exclude_all = up.species_exclude_init + up.species_exclude_zero
        species_index_exclude_all = []
        for species in species_exclude_all:
            species_index_exclude_all.append(gas.species_index(species))

        species_damp_init = ()
        for species in gas.species_names:
            if not(species in up.species_exclude_init) and not(species in up.species_major):
                species_damp_init = species_damp_init + (species,)
        species_index_damp_init = []
        for species in species_damp_init:
            species_index_damp_init.append(gas.species_index(species))

        species_damp = ()
        for species in gas.species_names:
            if not(species in species_exclude_all) and not(species in up.species_major):
                species_damp = species_damp + (species,)
        species_index_damp = []
        for species in species_damp:
            species_index_damp.append(gas.species_index(species))

        logging.info('Species for optimization: %s' % str(species_damp))
        noptim = len(species_damp)
        logging.info('Number of variables for optimization: %s' % noptim)

        cases = []
        quantityrefs = []
        tolerances = []

        # Constructing flame cases
        case_id = -1
        flametab = []
        for P in up.P_fl:
            for T in up.T_fl:
                for phi in up.phi_fl:
                    case_id += 1
                    cur = Flame.OneDFlame(
                        T, P, phi, up.fuel, up.n2_o2_ratio, case_id)
                    cur.case_id = case_id
                    cur.firstflag = True
                    cur.isflame = True
                    flametab.append(cur)
                    cases.append(cur)
                    tolerances.append(up.tolerance_fl)
                    cur.partitionid = case_id

        # Group procs in for flames
        npartition = len(flametab)
        procs = range(0, comm.size-1)
        partitioning = []
        counting = -1
        for k in range(len(procs)):
            if Global.nflames > 0:
                # Assign flames only every two cores
                # TODO parameter for this
                if k % 2 == 0:
                    counting += 1
                    partitioning.append([counting % (npartition)])
                else:
                    partitioning.append([-1])
            else:
                partitioning.append([-1])

        logging.debug('Partitioning: %s' % str(partitioning))

        Global.partitioning = partitioning

        # Constructing AI cases
        nai = len(up.P_ai)*len(up.T_ai)*len(up.phi_ai)
        for P in up.P_ai:
            for T in up.T_ai:
                for phi in up.phi_ai:
                    case_id += 1
                    cur = RCV.ReactorCstVolume(
                        Global.gas, T, P, phi, up.fuel, up.n2_o2_ratio)
                    cur.case_id = case_id
                    cur.isflame = False
                    cases.append(cur)
                    tolerances.append(up.tolerance_ai)

        # Computing reference cases
        tasks = []
        taskindex = -1
        for case in cases:
            taskindex += 1
            tasks.append((case, species_index_exclude_init,
                          species_index_damp_init,
                          np.ones(len(species_index_damp_init)), taskindex))
        quantityrefs = np.array(Par.tasklaunch(tasks))
        for case in flametab:
            case.firstflag = False

        tolerances = np.array(quantityrefs) * np.array(tolerances)

        # TODO remove these global variables
        params = Param.Parameters()
        params.cases = cases
        params.species_index_damp = species_index_damp
        params.species_index_exclude_all = species_index_exclude_all
        params.noptim = noptim
        params.quantityrefs = np.array(quantityrefs)
        params.tolerances = tolerances
        params.species_damp = species_damp
        params.threshold = up.threshold
        Global.params = params

        logging.info('Number of cases: %i' % len(cases))
        logging.info('Reference quantities: %s' % quantityrefs)

        # Setting up IPOpt parameters
        nvar = noptim  # Number of variables
        x_L = np.ones(nvar)*0.0  # Vector of lower bounds
        x_U = np.ones(nvar)*1.0  # Vector of upper bounds
        ncon = len(cases)  # Number of constraints
        logging.info('Number of constraints: %i' % ncon)

        # Building lower and upper bound vectors of tolerance
        g_L = -5.0 * np.ones(len(cases))
        g_U = tolerances

        nnzj = ncon * noptim
        nnzh = 64000

        # Initial coefficient vector
        if hasattr(up, 'x0'):
            x0 = np.array(up.x0)
            logging.info('Initializing x0 from input file')
        else:
            x0 = np.ones(noptim)
            logging.info('Initializing x0 from one')

        if not(len(x0) == noptim):
            logging.error('Wrong length for x0')
            sys.exit(0)

        logging.info('Initial x0 vector')
        for k in range(len(species_damp)):
            logging.info('%s %s' %
                         (gas.species_names[species_index_damp[k]], x0[k]))
        logging.info('x0 passed initial: %s' % str(x0))

        type_calc = 'OPT'
        if hasattr(up, 'type_calc'):
            type_calc = up.type_calc

        if type_calc == 'OPT':
            # Optimization case
            logging.debug('IPOpt initialization')
            nlp = pyipopt.create(nvar, x_L, x_U, ncon, g_L, g_U, nnzj, nnzh,
                                 OptFn.eval_f, OptFn.eval_grad_f, OptFn.eval_g, eval_jac_g_wrapper)

            logging.info('Starting optimization')
            x, zl, zu, constraint_multipliers, obj, status = nlp.solve(x0)
            nlp.close()
            logging.info('End optimization')

            logging.info('Final solution: %s' % x)
            for k in range(len(species_damp)):
                logging.info('%s %s' %
                             (gas.species_names[species_index_damp[k]], x[k]))

        elif type_calc == 'VAL':
            # Validation case
            logging.info('Validation calculation')

            x1 = np.array(x0)
            xreduced = []
            zero_species = []
            remaining_species = []

            for k in range(len(species_damp)):
                if (x1[k] < params.threshold):
                    x1[k] = 0.0
                    zero_species.append(
                        gas.species_names[species_index_damp[k]])
                else:
                    x1[k] = 1.0
                    xreduced.append(x1[k])
                    remaining_species.append(
                        gas.species_names[species_index_damp[k]])
            logging.info('Eliminated species: %s' % str(zero_species))
            for k, species in enumerate(species_damp):
                if species in zero_species:
                    logging.info('%s %s' % (species, x0[k]))
            logging.info('Remaining species: %s' % str(remaining_species))
            for k, species in enumerate(species_damp):
                if species in remaining_species:
                    logging.info('%s %s' % (species, x0[k]))

            # Computing reference quantities
            tasks = []
            taskindex = -1
            for case in cases:
                taskindex += 1
                tasks.append((case, species_index_exclude_init, species_index_damp_init, np.ones(
                    len(species_index_damp_init)), taskindex, 'detailed'))

            for case in cases:
                taskindex += 1
                tasks.append((case, species_index_exclude_all,
                              species_index_damp, x1, taskindex, 'reduced_zero'))

            results = np.array(Par.tasklaunch(tasks))
            refs = results[0:len(cases)]
            optimized = results[len(cases):2*len(cases)]

            errortab = np.abs((optimized-refs)/refs)
            logging.info('Reference quantities: %s' % str(refs))
            logging.info('Errors: %s' % str(errortab))
            if Global.nflames > 0:
                logging.info('Maximum error flame speed: %s' %
                             np.max(errortab[0:len(flametab)]))
            if nai > 0:
                logging.info('Maximum error AI delay: %s' %
                             np.max(errortab[len(flametab):]))
            logging.info('Maximum error: %s' % np.max((optimized-refs)/refs))

        elif type_calc == 'SA':
            logging.info('SA')
            # Sensivity analysis case
            if hasattr(up, 'SA_results'):
                logging.info('Read sensitivities from a file')
                file = open(up.SA_results, 'r')
                maxi = [0 for i in range(len(species_index_damp_init))]
                for line in file:
                    data = line.split()
                    idx = 0
                    for jdx in range(len(species_index_damp_init)):
                        if gas.species_names[species_index_damp_init[jdx]] == data[0]:
                            idx=jdx
                            break
                    maxi[idx] = float(data[1])
            else:
                tasks = []
                taskindex = -1
                perturb = 0.01
                nvars = len(species_index_damp_init)
                for k in range(nvars):
                    x0 = np.ones(nvars)
                    x0[k] = 1.0-perturb
                    for case in cases:
                        taskindex += 1
                        tasks.append((case, species_index_exclude_init,
                                      species_index_damp_init, x0, taskindex))
                results = np.array(Par.tasklaunch(tasks))
                results = results.reshape((nvars, len(cases)))
                results = np.abs(
                    (results[:, :] - quantityrefs[:])/quantityrefs[:]) / perturb
                maxi = np.max(results, axis=1)

            sorting = np.argsort(maxi)
            for indsort in sorting:
                logging.info('%s %s' %
                             (gas.species_names[species_index_damp_init[indsort]], maxi[indsort]))

            nvars = len(species_index_damp_init)
            currentkStart = 1
            if hasattr(up, 'currentkStart'):
                currentkStart = int(up.currentkStart)
            currentkStep = 1
            if hasattr(up, 'currentkStep'):
                currentkStep = int(up.currentkStep)

            error_history = []
            for currentk in range(currentkStart, len(sorting)+1, currentkStep):
                toremove = sorting[0:currentk]
                x1 = np.ones(nvars)
                x1[toremove] = 0.0

                tasks = []
                taskindex = -1
                for case in cases:
                    taskindex += 1
                    tasks.append((case, species_index_exclude_init,
                                  species_index_damp_init, x1, taskindex))

                results = np.array(Par.tasklaunch(tasks))
                refs = quantityrefs
                error_history.append(np.max(np.abs((results-refs)/refs)))
                logging.info('Last species removed: %s' %
                             gas.species_names[species_index_damp_init[toremove[-1]]])
                logging.info('Total error with %s removed species: %s' %
                             (currentk, np.max(np.abs((results-refs)/refs))))
                logging.info('Errors: %s' % str(np.abs((results-refs)/refs)))
                logging.info('Error history: %s' % str(error_history))
