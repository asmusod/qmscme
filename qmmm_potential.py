# script-version of ase_qmscme with worked in finland...
from ase import Atoms
from ase.data import s22
import numpy as np
atoms = s22.create_s22_system('Water_dimer')
atoms.center(vacuum=50.0)

# PC test objects
from gpaw.calc_qmmm import calc_qmmm

# rotate the second water mol into the same x plane as
# the first one, to make the potential curve easier to make
center = atoms[0].position
h = atoms[3].position[1]-atoms[0].position[1]
l = np.linalg.norm(atoms[0].position - atoms[3].position)
angle = np.arcsin(h/l)
atoms.rotate('-z',a=angle,center=center)

# Add an additional atom at the corner!
new = atoms[:3].copy()
pp = new[0].position
new.set_positions(new.get_positions() - pp)

# Shift water molecule 1 by 5.0 ang in X

from ase.visualize import view
tot = atoms + new

pos = tot.get_positions()

x = np.arange(-1.0,3.0,0.1)

new = np.zeros_like(pos)
new[3:6,0] += 0.0

tot.set_positions(pos + new)
view(tot)

from gpaw import GPAW

# QMMM GPAW-SCME Specific things
from SCME.CALC_SCME import CALC_SCME
from gpaw.mm_potentials import DipoleQuad

from ase.units import kcal, mol, Bohr, Hartree
k_c = 332.1 * kcal / mol

mp = 3

# NO MIC
# Charges for initial POOR field
charges = np.zeros(len(tot))
charges[:mp] += 1.0
charges[0]   *= -1.87
charges *= 1.275 * 4.80**2 / 14.4 # WHY?

# New charges
ncharges = np.zeros(len(tot))
ncharges += 0.42
ncharges[::3] *= -2.
tot.set_initial_charges(ncharges)

# Make initial 'POOR' external field, and gradient
nMM = (len(tot) - mp) / mp
eF  = np.zeros((3,nMM))
deF = np.zeros((3,3,nMM))

pos = tot[:mp].positions

for i in range(nMM):
    #
    CM = tot[(i+1)*mp:(i+2)*mp].get_center_of_mass()

    for j, ps in enumerate(pos):
        d  = np.sqrt(((ps - CM)**2).sum())
        qr = charges[j]*(CM-ps)
        #
        eF[:,i] += qr / d**3
        #
        for k in range(3):
            val = 3 * qr * (ps[k] - CM[k]) / d**5 
            deF[:,k,i] -= val
        deF[:,:,i] += 1./3 * np.diag(charges[:3] / d**3)

print deF[:,:,0]
print deF[:,:,1]
xxx
scme = tot[mp:]

calc_scme = CALC_SCME(scme, eF, deF)
scme.set_calculator(calc_scme)
scme.get_potential_energy()

#print calc_scme.energy, calc_scme.dipoles[0]

# Make QM BOX:
# Create QM cell, which is to be fixed
pos = tot[:mp].get_positions()
rcut = 4.0 # Well beyong basis set cut-off

C = np.zeros((3,3))
xmin = pos[:,0].min(); xmax = pos[:,0].max()
C[0,0] += xmax - xmin + 2 * rcut
ymin = pos[:,1].min(); ymax = pos[:,1].max()
C[1,1] += ymax - ymin + 2 * rcut
zmin = pos[:,2].min(); zmax = pos[:,2].max()
C[2,2] += zmax - zmin + 2 * rcut

old_pos = pos.copy()
small = tot[:3]
small.set_cell(C)
small.center()
new_pos = small.positions
shift = old_pos - new_pos
scme.set_positions(scme.positions - shift[:,0])

view(small + scme)

# Grab energy of lone water
calc = GPAW(h=.18, mode='lcao', xc='PBE', txt='tst.txt')
# small.set_calculator(calc)
# small.get_potential_energy()

# Grab dipoles and qpoles
dipole = calc_scme.dipoles
quad   = calc_scme.qpoles

calc.set(external=DipoleQuad(scme,dipole,quad,3))
small.set_calculator(calc)
energyQ = small.get_potential_energy()

# Calculate the effect of the Vext on the nuclei
comp_char = calc.get_compensation_charges(index=mp)

energy = 0

for a, atom in enumerate(small):
    pos = atom.position
    ch  = comp_char[a]
    # Interact with all SCME:
    for i in range(nMM):
        cm = scme[i*mp:(i+1)*mp].get_center_of_mass()
        r  = (cm - pos)
        d  = np.sqrt(((cm - pos)**2).sum())
        mUr = dipole[i].dot(r) 
        Q   = quad[:,:,i] 

        energy += mUr / d**3 * ch
        for j in range(3):
            for k in range(3):
                energy += Q[j,k]*r[j]*r[k] / d**5 * ch
                #if j==k:
                #    energy -= 1./3 * Q[j,k]*r[j]**2 / d**5 * ch

energy /= 4.8

# TEST calc_qmmm energy
#energyB, forcesB = calc_qmmm(mm_subsystem=scme,qm_subsystem=small,index=3,mp=mp,
#                             comp_charge=comp_char).get_energy_and_forces(tot)

#print energy, energyB
# Get distance between Oxygen SCME, Oxygen QM
O1 = small[0].position
O2 = scme[0].position
d = np.sqrt(((O1-O2)**2).sum())
print energy
# LJ Energy:
LJene = 4 * 6.59e-3 * (3.15**12 / d**12 - 3.15**6 / d**6)

print d, energyQ + 9.55600 + energy + LJene
