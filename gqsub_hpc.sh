#!/bin/sh
IN=$1

FIL=${IN%.*}
SUBMIT=.qsub.${FIL}
THISDIR=$PWD

cat > $SUBMIT << !EOF
#!/bin/sh
# -------------------------------------------------------------------
#PBS -N $FIL
#PBS -S /bin/bash
#PBS -o $THISDIR/\$PBS_JOBNAME.\$PBS_JOBID.log
#PBS -e $THISDIR/\$PBS_JOBNAME.\$PBS_JOBID.err
#PBS -q hpc
#PBS -l nodes=8:ppn=4,mem=22gb,walltime=2:00:00
# -------------------------------------------------------------------
#module unload mpi/gcc
#module load oldmpi/gcc
#module load oldmpi

cd \$PBS_O_WORKDIR
 
# Setup temporary directory
export PATH=/zhome/c7/a/69784/SCME2015_Nov/QMSCME:$PATH
export PYTHONPATH=/zhome/c7/a/69784/SCME2015_Nov/QMSCME:$PYTHONPATH

export PBS_O_TMPDIR=/SCRATCH/\$PBS_O_LOGNAME/\$PBS_JOBID
#mkdir -p \$PBS_O_TMPDIR
#cd \$PBS_O_TMPDIR
 
cd \$PBS_O_WORKDIR
 
# Run GPAW
mpirun -np 32 gpaw-python $IN
#mpirun -np 32 python $IN

# Copy files back
# NO because this overwrites the .py files to 0 kB if hard limit quota
# is hit. It's not needed since the traj and out files are written
# directly to the original dir due to PATH in the py files
#cp \$PBS_O_TMPDIR/* \$PBS_O_WORKDIR/
 
# Write ending statements to stdout
echo Stop time is \`date\`
echo ========= Job finished ===================
!EOF

qsub $SUBMIT

