import numpy as np
import Global
import logging
import Parallelism as Par


def getQuantityCase(case, species_index_exclude, species_index_damp, x):
    case.storeMultiplier(species_index_exclude, species_index_damp, x)
    return case.getQuantity()


def eval_f(x, user_data=None):
    logging.debug('F evaluation')
    costfunc = 2.0 * np.sum(x*(1-x)) + np.sum(x*x)
    return costfunc


def eval_grad_f(x, user_data=None):
    logging.debug('Grad F evaluation')
    noptim = Global.params.noptim
    grad_f = 2.0 * (np.ones(noptim)-2*x) + 2*x
    return grad_f


def eval_g(x, user_data=None):
    logging.debug('G evaluation')
    # Retrieve global parameters
    species_index_damp = Global.params.species_index_damp
    species_index_exclude = Global.params.species_index_exclude_all
    cases = Global.params.cases
    quantityrefs = Global.params.quantityrefs

    # Building Tasks
    tasks = []
    taskindex = -1
    for case in cases:
        taskindex += 1
        tasks.append((case, species_index_exclude,
                      species_index_damp, x, taskindex))

    # Retrieve array of results
    quantities = np.array(Par.tasklaunch(tasks))
    return np.abs(quantities - quantityrefs)


def gradient_constraints(x):
    logging.debug('Grad G evaluation')
    # Retrieve global variables
    noptim = Global.params.noptim
    cases = Global.params.cases
    species_index_damp = Global.params.species_index_damp
    species_index_exclude = Global.params.species_index_exclude_all
    quantityrefs = Global.params.quantityrefs
    tolerances = Global.params.tolerances
    species_damp = Global.params.species_damp

    # Building tasks
    tasks = []
    taskindex = -1
    for case in cases:
        taskindex += 1
        tasks.append((case, species_index_exclude,
                      species_index_damp, x, taskindex))

    # Perturbed quantities
    step = 1e-2
    for kcase, case in enumerate(cases):
        for kx in range(noptim):
            xperturb = np.array(x)
            xperturb[kx] += step
            taskindex += 1
            tasks.append((case, species_index_exclude,
                          species_index_damp, xperturb, taskindex))

    # Collecting gradient results
    logging.debug('Collecting Grad G results')

    results = np.array(Par.tasklaunch(tasks))
    initquantities = results[0:len(cases)]
    perturbquantities = results[len(cases):].reshape((len(cases), noptim))

    stack = []
    for kcase, case in enumerate(cases):
        grad = (np.abs(perturbquantities[kcase, :]-quantityrefs[kcase])
                - np.abs(initquantities[kcase]-quantityrefs[kcase]))/step
        stack = np.hstack((stack, grad))

    # Debug information
    logging.info('Gradient G done')
    logging.debug('x')
    logging.debug(str(x))
    logging.debug('sum(x)')
    logging.debug(str(np.sum(x)))


    # Logging info
    logging.info('Current x')
    logging.info(str(x).replace('  ', ','))
    logging.info('Current sum(x)')
    logging.info(str(np.sum(x)))
    logging.info('Current distance vector')
    logging.info(str(np.abs(initquantities-quantityrefs)/quantityrefs))

    candidatecheck = True
    errortab = np.abs(initquantities-quantityrefs)
    for kerror, error in enumerate(errortab):
        # logging.info( 'Error %s %s %s'%( kerror, error, tolerances[kerror] ) )
        if error > tolerances[kerror]*1.1:
            candidatecheck = False
    if candidatecheck:
        logging.info('Current solution satisfies constraint')
        logging.info('Current vector')
        logging.info(str(x).replace('  ', ','))
        logging.info('Number of eliminated species')
        count = np.where(x < Global.params.threshold)
        logging.info(str(len(count[0])))
        logging.info('List of eliminated species')
        for kdamp in range(len(species_damp)):
            if x[kdamp] < Global.params.threshold:
                logging.info(
                    Global.gas.species_names[species_index_damp[kdamp]] + ' ' + str(x[kdamp]))
    else:
        logging.info('Current solution does not satisfy constraint')
        logging.info('Current vector')
        logging.info(str(x).replace('  ', ','))

    logging.info('End Gradient PP')

    return stack
