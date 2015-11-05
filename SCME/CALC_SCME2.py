""" SCME potential as ASE calculator """

import numpy as np
import ase.units as unit
import scme2015 as scme
from ase.visualize import view

class CALC_SCME:
    def __init__(self,atoms,eF=None,deF=None):
        self.energy = None
        self.forces = None
        # one could make a check for proper numatoms here
        self.numatoms = len(atoms)
        self.oidx = len(atoms) * 2/3
        self.molnum = len(atoms) / 3
        self.dipoles = None
        self.eF      = eF
        self.qpoles  = None
        self.deF = deF
        self.eT = None # total field out of SCME, made under infl of QM

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
        self.positions = atoms.get_positions()
        pos = self.positions[:]
        scmepos = self.positions[:]
        truepos = pos[:]
        self.cell = atoms.get_cell()
        cell = self.cell
        
        numatoms = self.numatoms
        nummols = self.molnum        

        # Atomwise  MIC PBCs needed for SCME
        n = np.zeros(np.shape(scmepos)) 
        c_mid = cell.diagonal() * 0.5
        n = np.rint((c_mid - pos) / cell.diagonal())
        scmepos += n * cell.diagonal()
 
        scme_atoms = atoms[:]
        scme_atoms.set_positions(scmepos)
        
        # Convert to scme coordinates
        scme_atoms = self.ase_to_scme(scme_atoms)
        scme_coords = scme_atoms.get_positions()
        scme_coords = scme_coords.transpose()
        
        # Call SCME
        # with shape tests:
        eF = np.reshape(self.eF, np.shape(self.eF), order='F')
        deF = np.reshape(self.deF, np.shape(self.deF), order='F')
        #print np.shape(deF)
        #eQM = np.zeros([3,nummols])
        testin = np.zeros([3,3,nummols])
        testin = testin.reshape(3,3,nummols,order='F')
        eF = eF #/ unit.Debye
        deF = deF #/ unit.Debye 
        ff,epot,eT,dipole,qpole = scme.main(scme_coords,cell.diagonal(),eF,deF)
        # also get eT out: total field.
        #f = np.reshape(ff,[numatoms,3],order='F')
        f = np.reshape(ff,[numatoms,3])
        #testout = testout.reshape(nummols,3,3,order='F')
        # Convert force array back to ase coordinates
        aseforces = self.f_scme_to_ase(f)
        dipole = dipole*unit.Debye # go back to ase units
        dipole = dipole.transpose() # F-> P     
        eT = -1 * eT * unit.Debye
        eT = eT.transpose() # F-> P     

        # update calc and atoms
        self.dipoles = dipole.reshape(nummols, 3)
        self.qpoles  = qpole
        self.atoms = atoms 
        self.forces = aseforces
        self.energy = epot 
        self.eT = eT
    def get_dipoles(self):
        if self.dipoles is None:
            print 'no dipoles initialized, updating'
            self.update(atoms)
        return self.dipoles

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        raise NotImplementedError

    def update(self, atoms):
        if self.energy is None:
            self.calculate(atoms)
        elif (self.positions != atoms.get_positions()).any():
            self.calculate(atoms)
