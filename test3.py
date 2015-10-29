from ase import Atoms
from ase.data import s22
import numpy as np
from calc_scme2015 import CALC_SCME
dimer = s22.create_s22_system('Water_dimer')
dimer.center(vacuum=5.0)
# rotate the second water mol into the same x plane as
# the first one, to make the potential curve easier to make
center = dimer[0].position
h = dimer[3].position[1]-dimer[0].position[1]
l = np.linalg.norm(dimer[0].position - dimer[3].position)
angle = np.arcsin(h/l)
dimer.rotate('-z',a=angle,center=center)

# checking if i managed to properly change numdimer to num_mm_dimer
#dimer.set_calculator(CALC_SCME(dimer,np.zeros((3,1)), qmidx = 3))
#dimer.get_potential_energy()
#XXX

nsys = dimer[:3].copy()
nsys.set_positions(dimer[:3].get_positions() + (2.,2.,0))
trimer = dimer + nsys # 3 mols now

#from ase.visualize import view
#view(trimer)
#              xyz,nummols        
eF = np.zeros((3,3))
calc = CALC_SCME(trimer, eF, qmidx=0)
trimer.set_calculator(calc)
trimer.get_potential_energy()

eF2 = np.zeros((3,2))
charges = np.array([-2.00,1.0,1.0]) * 4.80**2 / 14.4

# Geta distance vector and distance to center of mass MM
pos = dimer.get_positions()
mass = dimer.numbers
pos_qm = nsys.positions

for i in range(2):
    #
    pos  = dimer[i*3:(i+1)*3].positions
    CM_i = (mass[:3]*pos.T).T
    CM   = CM_i.sum(axis=0) / mass[:3].sum()
    #
    for j, pos in enumerate(pos_qm):
        d = np.sqrt(((pos - CM)**2).sum())
        qr = charges[j]*(CM - pos)
        #
        eF2[:,i] += qr / d**3

# eF2 to correct units: ??
#eF2 *= 4.803206799**2 / 14.39975841

calc2 = CALC_SCME(dimer, 0*eF2) # still two mols
dimer.set_calculator(calc2)
dimer.get_potential_energy()
calc3 = CALC_SCME(dimer, eF2) # two mols
dimer.set_calculator(calc3)
dimer.get_potential_energy()
# Dipole to field
dip = calc.dipoles[2]
print dip

eF3 = np.zeros((3,2))
CM_dip = ((mass[:3]*pos_qm.T).T).sum(axis=0) / mass[:3].sum()

for i in range(2):
    pos  = dimer[i*3:(i+1)*3].positions
    CM_i = (mass[:3]*pos.T).T
    CM   = CM_i.sum(axis=0) / mass[:3].sum()

    rd  = CM_dip - CM
    d   = np.sqrt(((CM - CM_dip)**2).sum())
        
    eF3[:,i] += 1./d**3 * (3*np.dot(dip,rd)*rd/d**2 - dip)

#eF3 /= np.sqrt(4.803206799**2 / 14.39975841)
qmidx = 3
calc4 = CALC_SCME(dimer, eF3[:,qmidx/3:] ,qmidx=3)
dimer.set_calculator(calc4)
dimer.get_potential_energy()

print eF2
print eF3
print 'SCME tot. dip:'
print calc.dipoles
print 'POOR dip'
print calc3.dipoles
print 'dip. back to'
print calc4.dipoles
print 'No dip'
print calc2.dipoles

true_calc_dipoles  = [[1.50660513,1.72364953,0.], 
                      [1.05891058,-1.84194872,0.],
                      [1.27713458, 1.44286774 , 0.]]
true_calc3_dipoles = [[1.48032539,1.71848724,0.],[1.02122728,-1.63225503,0.]]
true_calc4_dipoles = [[1.48064826,1.74366363,0.],[1.05443859,-1.63582227,0.]]
true_calc2_dipoles = [[1.32856459,1.57615217,0.],[1.27993796,-1.72734717,0.]]
calc4_qmidx3_dipoles = [[ 0.68069779,-1.50795639,0.]]

maxdiff = 1e-7

if not (np.absolute(calc.dipoles-np.array(true_calc_dipoles)) < maxdiff).all():
    print 'SOMETHING CHANGED IN calc.dipoles' 
if not (np.absolute(calc3.dipoles-np.array(true_calc3_dipoles)) < maxdiff).all():
    print 'SOMETHING CHANGED IN calc3.dipoles' 
if not (np.absolute(calc4.dipoles-np.array(true_calc4_dipoles)) < maxdiff).all():
    print 'SOMETHING CHANGED IN calc4.dipoles:'
    if (np.absolute(calc4.dipoles-np.array(calc4_qmidx3_dipoles)) < maxdiff).all(): 
        print 'BUT they are the same as the first time I implemented qmidx ...'
    else:
        print 'Original:'
        print true_calc4_dipoles
        print 'qmidx original:'
        print calc4_qmidx_dipoles
        print 'NOW:'
        print calc4.dipoles
if not (np.absolute(calc2.dipoles-np.array(true_calc2_dipoles)) < maxdiff).all():
    print 'SOMETHING CHANGED IN calc2.dipoles'


dimer.get_forces()
