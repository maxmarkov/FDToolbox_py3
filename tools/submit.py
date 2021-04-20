import os

def generate_submit_manneback(folders, nprocs = 16, path_to_exec="/home/ucl/modl/mmarkov/soft/vasp/vasp.5.4.1_couplageEM/bin/vasp_std"):
    """
    Generate submit file for Manneback cluster
    """

    submit_string = (
    '#!/bin/bash -l\n'
    '#SBATCH --nodes=1\n'
    '#SBATCH --ntasks='+str(nprocs)+'\n'
    '#SBATCH --ntasks-per-core=1\n'
    '#SBATCH --cpus-per-task=1\n'
    '#SBATCH --time=120:00:00\n'
    '#SBATCH --job-name=relax_job\n'
    '#SBATCH --output=relax_job.out\n'
    '#SBATCH --error=relax_job.error\n'
    '#SBATCH --constraint=Haswell\n'
    '#SBATCH --mem-per-cpu=3900\n\n'

    'ulimit -s unlimited\n\n'

    'module purge\n'
    'module load releases/elic-2017b\n'
    'module load intel/2017b\n\n'

    'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n\n'

    'MPIRUN=""srun""\n'
    'MPIOPT=""--mpi=pmi2 -K1 -n $SLURM_NTASKS""\n'
    'VASP='+path_to_exec+'\n\n'

    '${MPIRUN} ${MPIOPT} ${VASP}\n\n'

    'echo ""--""')

    for folder in folders:
        with open(os.path.join(folder, 'submit.job'), 'w') as writer:
           writer.write(submit_string)
