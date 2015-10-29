""" SCME potential as ASE calculator """

import numpy as np
import ase.units as unit
import scme2015 as scme
from ase.visualize import view

class CALC_SCME:
    def __init__(self,efield):
        self.energy = None
        self.forces = None
        # one could make a check for proper numatoms here
        self.dipoles = None
        self.eF      = efield
        self.qpoles  = None

    def ase_to_scme(self, atoms):
        # Reindex in SCME order (HHHH..OO..)
        oxygens = atoms[[atom.index for atom in atoms if atom.symbol =='O']]
        hydrogens = atoms[[atom.index for atom in atoms if atom.symbol =='H']]
        scmeatoms = hydrogens + oxygens
        # Change coordinate system to cell-centered (SCME-Style)
        asecell = atoms.get_cell()
        c_mid = asecell.diagonal() * 0.5
        coords_asecell = scmeatoms.get_positions() 
        coords_scmecell = coords_asecell - c_mid
        scmeatoms.set_positions(coords_scmecell)

        return scmeatoms

    def scme_to_ase(self, atoms):
        # Reindex in ASE order (OHHOHH..)
        # Not used, since the original atoms 
        # object is never updated with SCME stuff
        oidx = self.oidx
        molnum = self.molnum 
        mollist =[]
        for i in range(molnum):
            mollist.append(i+oidx)
            mollist.append(i+i)
            mollist.append(i+i+1)
        
        ase_atoms = atoms[[mollist]] 
        return ase_atoms  

    def f_scme_to_ase(self,f):
        # Convert force array back to ase coordinates
        # numatoms is only mm atoms
        oidx = self.oidx
        nummols = self.molnum # sorry :P 
        numatoms = self.numatoms
        mollist =[]
        aseforces = np.zeros([numatoms,3])
        for i in range(nummols):
            mollist.append(i+oidx)
            mollist.append(i+i)
            mollist.append(i+i+1)

        for i in range(numatoms):
            aseforces[i,0] = f[mollist[i],0]    
            aseforces[i,1] = f[mollist[i],1]   
            aseforces[i,2] = f[mollist[i],2]    
        return aseforces   
 
    def calculate(self,atoms):
        self.numatoms = len(atoms)
        self.oidx = self.numatoms * 2/3
        self.molnum = self.numatoms / 3
        self.positions = atoms.get_positions()
        mm_pos = self.positions[:]
        scmepos = mm_pos.copy()
        self.cell = atoms.get_cell()
        cell = self.cell

        numatoms = self.numatoms
        nummols = self.molnum        
        # Atomwise  MIC PBCs needed for SCME
        n = np.zeros(np.shape(scmepos)) 
        c_mid = cell.diagonal() * 0.5
        n = np.rint((c_mid - mm_pos) / cell.diagonal())
        scmepos += n * cell.diagonal()
        scme_atoms = atoms[:]
        scme_atoms.set_positions(scmepos)
        # Convert to scme coordinates
        scme_atoms = self.ase_to_scme(scme_atoms)
        scme_coords = scme_atoms.get_positions()
        scme_coords = scme_coords.transpose()
        eF = np.reshape(self.eF, np.shape(self.eF), order='F') 
        # Call scme, get mm forces, mm energy, mm dipoles out
        if np.shape(self.eF)[1] != nummols:
            raise ValueError('Dipole array does not match number of MM waters')
        ff,epot,dip = scme.main(scme_coords,cell.diagonal(),self.eF)
        f = np.reshape(ff,[numatoms,3])
        # Convert force array back to ase coordinates
        aseforces = self.f_scme_to_ase(f)
        dip = dip.transpose() # F-> P     
        # et = np.reshape(et, np.shape(et), order='F')
        # CHECK THIS PLEASE! WHERES THE REINDEXING?
        # update atoms with the new dipoles
        self.dipoles =dip.reshape(nummols, 3)
        #self.qpoles  = eq
        self.forces = aseforces
        self.energy = epot 

    def get_potential_energy(self, atoms):
        if self.calculation_required(atoms,[]):
            self.update(atoms)
        return self.energy

    def get_dipoles(self):
        if self.dipoles is None:
            self.update(atoms)
        return self.dipoles

    def get_forces(self, atoms):
        # implement if something changed-check
        self.update(atoms)
        print 'FORCES FROM MM:'
        print self.forces 
        #f_mm = self.mm_forces
        """ all this into calc_qmscme! """
        #f_dip_on_qm = self.get_dipole_on_qm_forces(atoms)
        #print 'DIPOLE ON QM NUCLEI FORCES:'
        #print f_dip_on_qm
        # call the other coupling functions. 
        # f_dip_on_qm = calc_dipole_on_qm_forces(atoms) - me 
        # f_qm_on_dip = calc_qm_on_dipoles_forces(atoms) - elvar
        # add up energies and forces
        return self.forces

    def get_dipole_on_qm_forces(self, atoms):
        if self.dipoles is None:
            self.update(atoms)
        dipoles = self.dipoles
        f_dip = 0 # write shit
        return(f_dip)

    def get_stress(self, atoms):
        raise NotImplementedError

    def calculation_required(self, atoms, quantities):
        return True

    def update(self, atoms):
        self.calculate(atoms)
