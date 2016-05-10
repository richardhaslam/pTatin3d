#!/bin/bash -l

#SBATCH --mail-user=dave.mayhem23@gmail.com
#SBATCH --mail-type=ALL

# Job Name (do not use a space!)
#SBATCH --job-name="pTat3d"

#SBATCH --output=p3davx-%j.out
#SBATCH --error=p3davx-%j.err

# Number of MPI processors (-n )
#SBATCH --ntasks=64

# Number of cores per node - max=12 (-N ) 
#SBATCH --ntasks-per-node=24

# #SBATCH --cpu_bind=verbose,cores
#SBATCH --mem_bind=verbose,local

#SBATCH --time=00:10:00


# no space
EXEC=${PWD}/${PETSC_ARCH}/bin/ptatin_driver_linear_ts.app

aprun -n $SLURM_NTASKS -N $SLURM_NTASKS_PER_NODE $EXEC -options_file src/models/viscous_sinker/examples/sinker-mfscaling.opts -a11_op avx

exit


