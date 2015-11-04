from ase import Atoms
from ase.data import s22
import numpy as np
from SCME.CALC_SCME2 import CALC_SCME
from ase_qmscme import ase_qmscme
from gpaw import GPAW
from ase.visualize import view
from gpaw.mpi import rank, MASTER
dimer = s22.create_s22_system('Water_dimer')
# rotate the second water mol into the same x plane as
# the first one, to make the potential curve easier to make
center = dimer[0].position
h = dimer[3].position[1]-dimer[0].position[1]
l = np.linalg.norm(dimer[0].position - dimer[3].position)
angle = np.arcsin(h/l)
dimer.rotate('-z',a=angle,center=center)
dimer.center(vacuum=5.0)

path = '/zhome/c7/a/69784/SCME2015_Nov/QMSCME/'

# and then i need to place another scme water mol FAR away so the program
# actually initializes


qmidx = 3
mp = 3
eps = 6.59e-3 
sig = 3.15
LJ_qm = np.zeros((3,2))
LJ_qm[0,0] = eps
LJ_qm[0,1] = sig
LJ_mm = np.zeros((6,2))
LJ_mm[0::mp,0] = eps
LJ_mm[0::mp,1] = sig


qm_cell = dimer.get_cell()
calc_qm = GPAW(h=0.20,mode='lcao',xc='PBE',basis={None:'dzp'},txt='classguy.txt')

monomer = dimer[:qmidx]
monomer.set_calculator(calc_qm)
monomer.center(vacuum=5.0)
Esingle = monomer.get_potential_energy()

print Esingle

if rank == MASTER:
    f = open(path+'add_gpaw-scme_dimerbinding.ene','a')
    f.write('Single H2O Energy : %24s\n' %Esingle)
    print('Single H2O Energy : %24s\n' %Esingle)
    f.write('%24s%24s%24s\n' % ('Dimer Dist','E_dimer','E_binding'))
    print('%24s%24s%24s\n' % ('Dimer Dist','E_dimer','E_binding'))

# start closer in, and move out
for mol1 in range(3):
    dimer[mol1].x += 1.5

dimer.center(vacuum=5.0)

for step in range(0,40,1): 
    for mol1 in range(3):
        dimer[mol1].x -= 0.05
    for mol2 in range(3,6):
        dimer[mol2].x += 0.05
    if step > 0: # if you accidentally break the script before done
        dimer.center(vacuum=5.0)
        qm_cell = dimer.get_cell()
        dist = np.linalg.norm(dimer[0].position - dimer[3].position)
        extra = dimer[0:3].copy()
        extra.translate([100,0,0])
        full_sys = dimer + extra
        calc_mm = CALC_SCME(full_sys[qmidx:])
        full_sys.set_cell([200,200,200])
        calc_qm = GPAW(h=0.20,mode='lcao',xc='PBE',
                       basis={None:'dzp'},txt=path+'dimer_%2.4f.txt'%dist)
        full_sys.set_calculator(ase_qmscme(full_sys, qmidx=3,calc_qm=calc_qm,
                          calc_mm=calc_mm,qm_cell=qm_cell,
                          LJ_qm=LJ_qm.T, LJ_mm=LJ_mm.T, qm_fixed = True))
        #try:
        Epot = full_sys.get_potential_energy()
        print 'Step: %5.2f, dist: %5.2f, Epot: %5.2f'%(step, dist,Epot)
        #except:
        #    Epot = float('Nan')
        if rank == MASTER:
            f.write('%24s%24s%24s\n' %(dist,Epot,Epot-Esingle))
            print('%24s%24s%24s\n' %(dist,Epot,Epot-Esingle))




