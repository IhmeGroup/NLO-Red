import numpy as np
import cantera as ct
import logging
import os
import Global


class OneDFlame(object):

    def __init__(self, Tin, P, phi, fuel, n2_o2_ratio, id):
        # Initialize parameters
        self.id = id
        self.phi = phi
        self.P = P
        self.Tin = Tin
        self.fuel = fuel
        self.n2_o2_ratio = n2_o2_ratio
        self.initialized = False

    def storeMultiplier(self, species_index_exclude, species_index_damp, multipliers):
        # Store coefficients and species
        self.species_index_exclude = species_index_exclude
        self.species_index_damp = species_index_damp
        self.multipliers = multipliers

    def putMultiplier(self, curgas):
        # Apply multiplier
        reactants = curgas.reactant_stoich_coeffs()
        products = curgas.product_stoich_coeffs()

        # Damping
        multicoeffs = np.ones(curgas.n_reactions)
        for k in range(curgas.n_reactions):
            for kk, index in enumerate(self.species_index_damp):
                if (not(reactants[index, k] == 0)) or (not(products[index, k] == 0)):
                    multicoeffs[k] *= self.multipliers[kk]
        for k in range(curgas.n_reactions):
            curgas.set_multiplier(multicoeffs[k], k)

        # Species exclusion
        for k in range(curgas.n_reactions):
            for index in self.species_index_exclude:
                if (not(reactants[index, k] == 0)) or (not(products[index, k] == 0)):
                    curgas.set_multiplier(0.0, k)

    def initfromFile(self, filename):
        # Solve flame from initial solution file

        # One Cantera gas phase for each problem, conflict otherwise
        Global.gases[self.id] = ct.Solution(Global.mechanism)
        gas = Global.gases[self.id]

        width = 0.03  # m
        loglevel = 1

        ifuel = gas.species_index(self.fuel)
        io2 = gas.species_index('O2')
        in2 = gas.species_index('N2')

        x = np.zeros(gas.n_species)
        x[ifuel] = self.phi
        x[io2] = gas.n_atoms(self.fuel, 'C') + 0.25 * \
            gas.n_atoms(self.fuel, 'H') - 0.5 * gas.n_atoms(self.fuel, 'O')
        x[in2] = x[io2] * self.n2_o2_ratio
        gas.TPX = self.Tin, self.P, x

        # Set up flame object
        Global.sims[self.id] = ct.FreeFlame(gas, width=width)
        f = Global.sims[self.id]
        f.restore(filename)
        # Solve with mixture-averaged transport model
        f.transport_model = 'Mix'
        # Refinement
        f.set_refine_criteria(ratio=2, slope=0.05, curve=0.05)

        Global.flame_inits[self.id] = True

    def initfromScratch(self):
        # Solve flame from scratch
        Global.gases[self.id] = ct.Solution(Global.mechanism)
        gas = Global.gases[self.id]

        width = 0.03  # m
        loglevel = 1

        ifuel = gas.species_index(self.fuel)
        io2 = gas.species_index('O2')
        in2 = gas.species_index('N2')

        x = np.zeros(gas.n_species)
        x[ifuel] = self.phi
        x[io2] = gas.n_atoms(self.fuel, 'C') + 0.25 * \
            gas.n_atoms(self.fuel, 'H') - 0.5 * gas.n_atoms(self.fuel, 'O')
        x[in2] = x[io2] * self.n2_o2_ratio

        gas.TPX = self.Tin, self.P, x

        # Set up flame object
        Global.sims[self.id] = ct.FreeFlame(gas)
        f = Global.sims[self.id]

        f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
        f.show_solution()

        # Solve with mixture-averaged transport model
        f.transport_model = 'Mix'
        f.solve(loglevel=loglevel, auto=True)

        # Refinement
        f.set_refine_criteria(ratio=2, slope=0.05, curve=0.05)
        f.solve(loglevel)
        f.show_solution()
        logging.debug(
            'Flame speed for the case = {0:7f} m/s'.format(f.u[0]))
        f.save(filename='%s/save_%i.xml' % (Global.directory, self.id),
               name='solution',  description='vamos')

        Global.flame_inits[self.id] = True

    def getSl(self):
        # Handling different cases
        if not(self.firstflag):
            # Not the first time computation
            if not(Global.flame_inits[self.id]):
                self.initfromFile('%s/save_%i.xml' %
                                  (Global.directory, self.id))
            f = Global.sims[self.id]
            gas = Global.gases[self.id]
            self.putMultiplier(gas)
            f.solve(refine_grid=False)
        else:
            # First time computation
            filename = '%s/start_%i.xml' % (Global.directory, self.id)
            # If file exists, restart from file
            if os.path.exists(filename):
                self.initfromFile(filename)
                f = Global.sims[self.id]
                gas = Global.gases[self.id]
                self.putMultiplier(gas)
                f.solve(refine_grid=False)
                f.save(filename='%s/save_%i.xml' % (Global.directory,
                                                    self.id), name='solution',  description='vamos')
            # Otherwise solve the flame from scratch
            else:
                self.initfromScratch()
                f = Global.sims[self.id]
                gas = Global.gases[self.id]
                self.putMultiplier(gas)
                f.solve(refine_grid=False)
        return f.u[0]

    def getQuantity(self):
        # Return relevant metric of the case
        return self.getSl()
