from ase import Atoms
from ase.data import s22
import numpy as np
from calc_scme2015 import CALC_SCME
from ase_qmscme import ase_qmscme
from gpaw import GPAW
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

qmidx = 3

calc_qm = GPAW(h=0.20,mode='lcao',basis={None:'dzp'})
eF = np.zeros((3,2))
calc_mm = CALC_SCME(eF)

dimer.set_calculator(ase_qmscme(trimer, qmidx=3,calc_qm=calc_qm,
calc_mm=calc_mm,qm_cell=trimer.get_cell()))
dimer.get_potential_energy()
