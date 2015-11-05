import numpy as np
import ase.units as unit
from ase.visualize import view
from gpaw.mm_potentials import DipoleQuad
from calc_qmscme import calc_qmscme
from gpaw.mpi import MASTER,rank
class ase_qmscme:
    """ GPAW<->SCME additive interfacer version 0.0001 
        A.O. Dohn,  E. O. Jonsson, November 2015 """
    def __init__(self,atoms,qmidx=0,
                 calc_qm=None,calc_mm=None,
                 qm_cell=None,rcut=None,qm_fixed=False,
                 mm_pbc = True, LJ_qm=None, LJ_mm=None):
        if rank == MASTER:
            print '##############################  IN __init__()'
        self.atoms = atoms
        self.qmidx = qmidx
        self.energy = None
        self.forces = None
        self.mm_forces = None # split in 2 ? 
        self.qm_forces = None # Not used yet...
        self.numatoms = len(atoms) 
        self.dipoles = None
        self.qpoles = None
        self.eF = None
        self.deF = None
        self.calc_qm = calc_qm
        self.calc_mm = calc_mm
        self.qm_cell = qm_cell # full atoms.cell isfor FULL qmmm system
        self.qm_fixed = qm_fixed
        self.rcut = rcut
        self.mm_pbc = mm_pbc
        self.comp_char = None
        self.LJ_qm = LJ_qm
        self.LJ_mm = LJ_mm
        self.mp = 3 # will always be 3.. SCME...
        self.initialized = False # check if initial dipoles need to be made 

    def get_qm_subsystem(self):
        if rank == MASTER:
            print '##############################  IN get_qm_subsystem()'
        """ All atoms before qmidx are considered the qm subsystem.
            Need to create a minimal cell around the qm subsystem and
            keep track of any displacement to fit it neatly in to the
            cell. """
        qmidx = self.qmidx
        pos = self.atoms[:qmidx].get_positions()
        rcut = self.rcut

        """ Can either specify a particular cell and keep it fixed, or
            rcut, where the qm system is then placed in a cell defined
            by min to max + rcut on both sides. """

        if self.qm_cell is None:
            # Find xmin-xmax, ymin-ymax and zmin-zmax: need rcut pos from origin
            # plus xmax, ymax and zmax + rcut from other border.
            C = np.zeros((3,3))
            xmin = pos[:,0].min(); xmax = pos[:,0].max()
            C[0,0] += xmax - xmin + 2 * rcut
            ymin = pos[:,1].min(); ymax = pos[:,1].max()
            C[1,1] += ymax - ymin + 2 * rcut
            zmin = pos[:,2].min(); zmax = pos[:,2].max()
            C[2,2] += zmax - zmin + 2 * rcut
            origin = np.array([xmin, ymin, zmin]) - rcut #(O)
        else:
            C = self.qm_cell
            throw = self.atoms[:qmidx]
            pos_old = self.atoms[0].position
            throw.set_cell(C)
            throw.center()
            pos_new = throw[0].position
            origin = pos_old - pos_new

        if self.qm_fixed:
            origin = 0.0
            C = self.qm_cell

        qm_subsystem = self.atoms[:qmidx]
        pos -= origin
        qm_subsystem.set_positions(pos)
        qm_subsystem.set_cell(C)
        qm_subsystem.set_pbc((0,0,0))

        # Hold on to the shift from origin and qm_subsystem:
        self.origin = origin # (O)
        self.qm = qm_subsystem
        #view(qm_subsystem)

    def get_mm_subsystem(self):
        if rank == MASTER:
            print '##############################  IN get_mm_subsystem()'
        """ All atoms after qmidx are considered as mm subsystem.
            Need to go through get_qm_subsystem to get the origin. """
        if (self.qm is None):
           self.get_qm_subsystem()

        pbc = np.zeros(3)

        if self.mm_pbc is True:
            pbc += 1

        # Only need to shift the positions of the mm subsystem by (O)
        # then find minimum distance as for the minimum image conv.
        qmidx = self.qmidx
        pos = self.atoms[qmidx:].get_positions()
        pos -= self.origin
        mp = self.mp

        # Minimum image relative to the center of the qm cell!
        n = np.zeros(np.shape(pos))
        c_mid = self.qm.cell.diagonal() * 0.5

        n[::mp] = np.rint((c_mid - pos[::mp]) / self.atoms.cell.diagonal())

        # Grab all atoms of this particular molecule
        for i in range(1,mp):
            n[i::mp] += n[::mp]

        pos += n * self.atoms.cell.diagonal() * pbc

        mm_subsystem = self.atoms[qmidx:]
        mm_subsystem.set_positions(pos)
        mm_subsystem.set_pbc((1,1,1))
        self.mm = mm_subsystem
        #view(mm_subsystem)

    def calculate_mm(self):
        if rank == MASTER:
            print '##############################  IN calculate_mm()'
        mm = self.mm
        calc_mm = self.calc_mm
        calc_mm.eF = self.eF
        calc_mm.deF = self.deF

        mm.set_calculator(calc_mm)

        self.mm_energy = 0
        self.mm_forces = np.zeros((len(mm),3))

        self.mm_energy += mm.get_potential_energy()
        self.mm_forces += mm.get_forces()

        # self.mm_eF = calc.mm.eF # get mm eF for torque for FZqm-on-dip

        self.dipoles = calc_mm.get_dipoles()

    def calculate_qm(self):
        if rank == MASTER:
            print '##############################  IN calculate_qm()'
        calc_qm = self.calc_qm
        qm = self.qm
        mm = self.mm

        self.qm_energy = 0
        self.qm_forces = np.zeros((len(qm),3))
        self.mm_forces = np.zeros((len(mm),3))

        dipoles = self.dipoles
        qpoles = self.qpoles
        if rank == MASTER:
            print 'DIPOLES IN CALC_QM:' 
            print dipoles
            print 'QPOLES IN CALC_QM:' 
            print qpoles
        calc_qm.set(external=DipoleQuad(mm,dipoles,qpoles,3))
        qm.set_calculator(calc_qm)

        self.qm_energy += qm.get_potential_energy()
        self.qm_forces += qm.get_forces()

        self.comp_char = calc_qm.get_compensation_charges(index=self.qmidx)
        # self.eF = calc_qm.external.eF # now they are updated
        # self.deF = calc_qm.external.deF # now they are updated

        #self.mm_forces = calc_qm.get_point_charge_forces(mm_subsystem = mm)

    def calculate_qmmm(self, atoms):
        if rank == MASTER:
            print '##############################  IN calculate_qmmm()'
        if self.comp_char is None:
            self.calculate_qm()
        if self.dipoles is None:
            self.calculate_mm()

        qmidx = self.qmidx
        mm = self.mm
        qm = self.qm
        LJ_mm = self.LJ_mm
        LJ_qm = self.LJ_qm

        comp_char = self.comp_char
        dipoles = self.dipoles
        qpoles = self.qpoles

        self.qmmm_energy = 0
        self.qmmm_forces = np.zeros((len(mm)+len(qm),3))

        self.qmmm_energy, self.qmmm_forces = calc_qmscme(mm_subsystem = mm,
            qm_subsystem = qm, qmidx = qmidx, comp_char = comp_char,
            dipoles=dipoles, qpoles=qpoles,LJ_mm=LJ_mm,LJ_qm=LJ_qm).get_energy_and_forces(atoms)

    def initialize(self): # induce poor initial dipoles to take first step
        if rank == MASTER:
            print '##############################  IN initialize()'
        # temporary.. goes in DipoleQuad..
        mp = self.mp
        qmidx = self.qmidx
        atoms = self.atoms

        # Charges for initial POOR field
        charges = np.zeros(len(atoms))
        #charges[:qmidx] += 1.0
        #charges[0]   *= -1.87
        charges[:qmidx] += 0.42
        charges[0]   *= -2.
        #                   eV->D     e_c
        #charges *= 1.275 * 4.80**2 / 14.4 # WHY?

        # New chages
        #ncharges = np.zeros(np.shape(charges))
        #ncharges += 0.42
        #ncharges[::3] *= -2.
        #atoms.set_initial_charges(ncharges)

        # Make initial 'POOR' external field, and gradient
        nMM = (len(atoms) - qmidx) / mp
        eF  = np.zeros((3,nMM))
        deF = np.zeros((3,3,nMM))

        pos = atoms[:qmidx].positions

        for i in range(nMM):
            #
            CM = atoms[(i+1)*mp:(i+2)*mp].get_center_of_mass()
            for j, ps in enumerate(pos):
                d  = np.sqrt(((ps - CM)**2).sum())
                qr = charges[j]*(CM-ps)
                #
                eF[:,i] += qr / d**3
                #
                for k in range(3):
                    val = 3 * qr * (ps[k] - CM[k]) / d**5 
                    deF[:,k,i] -= val
                deF[:,:,i] += 1./3 * np.diag(charges[:3] / d**3) # this is probably wrong

        scme = atoms[qmidx:]
        from SCME.CALC_SCME2 import CALC_SCME
        calc_scme = CALC_SCME(scme, eF, deF)
        scme.set_calculator(calc_scme) 
        scme.get_potential_energy() #now calc has diples

        self.dipoles = calc_scme.dipoles
        self.qpoles = calc_scme.qpoles
        self.eF = eF
        self.deF =deF

        self.initialized = True

    def get_energy_and_forces(self, atoms):
        if rank == MASTER:
            print '##############################  IN get_energy_and_forces()'
        self.positions = atoms.get_positions()
        self.get_qm_subsystem()
        self.get_mm_subsystem()

        qmidx = self.qmidx

        # check if this is step 1:
        if not self.initialized:
            self.initialize()

        self.energy = 0
        self.forces = np.zeros((len(self.mm)+len(self.qm),3))
        
        #view(self.qm + self.mm)

        self.calculate_mm()
        self.energy += self.mm_energy
        self.forces[qmidx:,:] += self.mm_forces
        if rank == MASTER:
            print 'MM ENERGY:'
            print self.mm_energy

        self.calculate_qm()
        self.energy += self.qm_energy
        if rank == MASTER:
            print 'QM ENERGY:'
            print self.qm_energy 
        self.forces[:qmidx,:] += self.qm_forces
        self.forces[qmidx:,:] += self.mm_forces

        self.calculate_qmmm(atoms)
        if rank == MASTER:
            print 'TOTAL QMMM ENERGY'
            print self.qmmm_energy
        self.energy += self.qmmm_energy
        self.forces += self.qmmm_forces

    def get_potential_energy(self, atoms):
        if rank == MASTER:
            print '##############################  IN get_potential_energy()'
        if self.calculation_required(atoms):
            self.get_energy_and_forces(atoms)
        return self.energy

    def get_forces(self, atoms):
        if rank == MASTER:
            print '##############################  IN get_forces()'
        if self.calculation_required(atoms):
            self.get_energy_and_forces(atoms)
        return self.forces

    def calculation_required(self,atoms):
        if rank == MASTER:
            print '##############################  IN get_calculation_required()'
        """ FIX FIX FIX """
        return True 

