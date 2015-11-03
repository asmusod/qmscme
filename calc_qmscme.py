import ase.units as unit
import numpy as np
from math import pi
from SCME.fCMtoAll import fCMtoAll

k_c = 332.1 * unit.kcal / unit.mol

class calc_qmscme:
    def __init__(self, mm_subsystem=None, qm_subsystem=None, 
                 density=None, qmidx=None, mp = 3,
                 comp_char = None, dipoles = None, qpoles = None,
                 origin = None, LJ_mm = None, LJ_qm = None):
        """ Class for calculating the QM/MM terms between gpaw and SCME """
        self.energy = None
        self.forces = None
        self.mm = mm_subsystem
        self.qm = qm_subsystem
        self.density = density 
        self.qmidx = qmidx
        self.mp = mp
        self.comp_char = comp_char
        self.origin = origin
        self.dipoles = dipoles
        self.qpoles = qpoles
        self.LJ_mm = LJ_mm
        self.LJ_qm = LJ_qm


    def calculate_LJ(self):
        """ The combination rules for the LJ terms follow Waldman-Hagler
            J. Comp. Chem. 14, 1077 (1993)
        """
        mm = self.mm
        qm = self.qm
        mp = self.mp

        LJ_mm = self.LJ_mm
        LJ_qm = self.LJ_qm

        energy = 0
        forces = np.zeros((len(qm)+len(mm),3))

        if LJ_qm == None:
            LJ_mm = np.zeros((2,len(mm)))
            LJ_qm = np.zeros((2,len(qm)))

        mm_c = self.mm.get_initial_charges()
        qmidx = self.qmidx
        for i in range(0,len(mm),mp): # only do oxygens!
            for j in range(len(qm)):
                D = qm[j].position - mm[i].position
                d = (D**2).sum()
                if LJ_mm[1,i] * LJ_qm[1,j] == 0:
                    epsilon = 0
                else:
                    epsilon = 2 * LJ_mm[1,i]**3 * LJ_qm[1,j]**3 \
                              * np.sqrt(LJ_mm[0,i] * LJ_qm[0,j]) \
                              / (LJ_mm[1,i]**6 + LJ_qm[1,j]**6)

                sigma = ((LJ_mm[1,i]**6 + LJ_qm[1,j]**6) / 2)**(1./6)

                energy  += (4 * epsilon * (sigma**12 / d**6 - sigma**6 / d**3))     

                f  = (4 * epsilon * (12 * sigma**12 / d**6 - \
                           6 * sigma**6 / d**3)) * D / d

                forces[qmidx+i,:] -= f
                forces[j,:] += f        

        return energy, forces


    def calculate_qmnuclei_dipoles(self):
        qm = self.qm
        mm = self.mm
        mp = self.mp
        qmidx = self.qmidx
        dipoles = self.dipoles
        qpoles = self.qpoles
        comp_char = self.comp_char

        nMM = len(mm)/mp

        # Energy
        energy = 0
        for a, atom in enumerate(qm):
            pos = atom.position
            ch  = comp_char[a]
            # Interact with all SCME:
            for i in range(nMM):
                cm = mm[i*mp:(i+1)*mp].get_center_of_mass()
                r  = (cm - pos)
                d  = np.sqrt(((cm - pos)**2).sum())
                mUr = dipoles[i].dot(r)
                Q   = qpoles[:,:,i]

                energy += mUr / d**3 * ch
                for j in range(3):
                    for k in range(3):
                        energy += Q[j,k]*r[j]*r[k] / d**5 * ch
                        #if j==k:
                        #    energy -= 1./3 * Q[j,k]*r[j]**2 / d**5 * ch

        energy /= 4.8
        """ forces """
        forces = np.zeros((len(qm)+len(mm),3))
        # Forces: MM dipoles on QM nuclei - CHECK UNITS
        f_a = np.zeros((len(qm),3))
        for a, atom in enumerate(qm):
            pos = atom.position
            ch = comp_char[a]
        # Interact with all SCME:
            for i in range(nMM):
                cm = mm[i*mp:(i+1)*mp].get_center_of_mass()
                r = (cm - pos)
                d = np.linalg.norm(r)
                r_u = r/d
                mUr = dipoles[i].dot(r)
                f_a[a,:] += ch * (3 * (mUr / d**4)  - dipoles[i] / d**3) * r_u
                
        # Forces: QM nuclei on MM dipoles    
        f_iCM = np.zeros((len(mm)/mp,3))
        for i in range(0,len(mm),mp):
            cm = mm[i:i+3].get_center_of_mass()
            dip = dipoles[i/mp]
            for a, atom in enumerate(qm):
                ch = comp_char[a]
                pos = atom.position
                r = (pos - cm)  
                d = np.linalg.norm(r)
                r_u = r/d
                mUr = dip.dot(r)
                f_iCM[i/mp,:] += ch *(dip/d**3 - 3*mUr/d**4) * r_u
        
        # Torque mu x Efield
        eT = mm.calc.eT
        tau_cm = np.cross(dipoles,eT) 
        f_i = fCMtoAll(f_iCM=f_iCM, tau_cm=tau_cm,atoms=mm).distribute()

        forces[:qmidx,:] += f_a
        forces[qmidx:,:] += f_i
        return energy, forces

    def update(self, atoms):
        qmidx = self.qmidx
        mm = self.mm
        qm = self.qm
        origin = self.origin
        if self.energy is None:
            self.calculate(atoms)
        elif ((mm.positions + origin != atoms[qmidx:].get_positions()).any() or
            (qm.positions + origin != atoms[:qmidx].get_positions()).any()):
            self.calculate(atoms)

    def calculate(self, atoms):
        qmidx = self.qmidx
        self.energy = 0
        self.forces = np.zeros((len(self.mm)+len(self.qm),3))
        e_lj, f_lj = self.calculate_LJ()
        print 'LJ FORCES:'
        print f_lj
        
        e_zdip, f_zdip = self.calculate_qmnuclei_dipoles()
        print 'Z-dip FORCES:'
        print f_zdip

        self.energy += e_lj 
        self.energy += e_zdip

        self.forces += f_lj
        self.forces[qmidx:,:] += e_zdip

    def get_energy_and_forces(self, atoms):
        self.update(atoms)
        return self.energy, self.forces
