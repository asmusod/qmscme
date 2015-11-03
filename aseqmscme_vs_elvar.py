from ase import Atoms
from ase.data import s22
import numpy as np
from SCME.CALC_SCME2 import CALC_SCME
from ase_qmscme import ase_qmscme
from gpaw import GPAW
from ase.visualize import view
dimer = s22.create_s22_system('Water_dimer')
dimer.center(vacuum=5.0)
# rotate the second water mol into the same x plane as
# the first one, to make the potential curve easier to make
center = dimer[0].position
h = dimer[3].position[1]-dimer[0].position[1]
l = np.linalg.norm(dimer[0].position - dimer[3].position)
angle = np.arcsin(h/l)
dimer.rotate('-z',a=angle,center=center)

nsys = dimer[:3].copy()
nsys.set_positions(dimer[:3].get_positions() + (2.,2.,0))
trimer = dimer + nsys # 3 mols now
trimer_elvar = trimer.copy() # test for comparison

qmidx = 3
mp = 3

calc_qm = GPAW(h=0.20,mode='lcao',xc='PBE',basis={None:'dzp'},txt='classguy.txt')
eF = np.zeros((3,2))
calc_mm = CALC_SCME(trimer[qmidx:])

eps = 6.59e-3 
sig = 3.15
LJ_qm = np.zeros((3,2))
LJ_qm[0,0] = eps
LJ_qm[0,1] = sig
LJ_mm = np.zeros((6,2))
LJ_mm[0::mp,0] = eps
LJ_mm[0::mp,1] = sig

trimer.set_calculator(ase_qmscme(trimer, qmidx=3,calc_qm=calc_qm,
                      calc_mm=calc_mm,qm_cell=trimer.get_cell(),
                      LJ_qm=LJ_qm.T, LJ_mm=LJ_mm.T))
Etot_classes = trimer.get_potential_energy()

print '\n################################## #'
print '########### TOTAL ENERGY classes## #'
print Etot_classes
print '################################## #'

""" NOW IS TIME FOR ELVAR OLD STYLE """

from gpaw import GPAW
# QMMM GPAW-SCME Specific things
from SCME.CALC_SCME2 import CALC_SCME
from gpaw.mm_potentials import DipoleQuad

tot = trimer_elvar.copy()

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

scme = tot[mp:]
calc_scme = CALC_SCME(scme, eF, deF)
scme.set_calculator(calc_scme)
mm_energy = scme.get_potential_energy()

# make QM BOX BUT THE SAME AS for trimer! 
pos = tot[:mp].get_positions()
#rcut = 4.0 # Well beyong basis set cut-off
C = trimer.get_cell()
#
#old_pos = pos.copy()
small = tot[:3]
small.set_cell(C)
#small.center()
#new_pos = small.positions
#shift = old_pos - new_pos
#scme.set_positions(scme.positions - shift[:,0])
view(small)

view(small + scme)
view(trimer.copy())
print 'Classes pos'
print trimer.get_positions()
print 'Elvar pos'
print (small+scme).get_positions()

# Grab energy of lone water
calc = GPAW(h=.20, mode='lcao', xc='PBE', basis={None:'dzp'}, txt='elvarguy.txt')
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
#print energy
# LJ Energy:
LJene1 = 4 * 6.59e-3 * (3.15**12 / d**12 - 3.15**6 / d**6)

O1 = small[0].position
O2 = scme[3].position
d = np.sqrt(((O1-O2)**2).sum())
#print energy
# LJ Energy:
LJene2 = 4 * 6.59e-3 * (3.15**12 / d**12 - 3.15**6 / d**6)

LJene = LJene1 + LJene2

Etot_elvar =  energyQ + 9.55600 + energy + LJene

print "TOTAL MM ENERGY:"
print mm_energy

print "TOTAL QM ENERGY: "
print energyQ

print "TOTAL QMMM ENERGY:"
print LJene+energy

print '\n################################## #'
print '########### TOTAL ENERGY elvar #### #'
print mm_energy+energyQ+LJene+energy
print '################################## #'
