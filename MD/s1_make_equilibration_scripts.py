import glob
import os

# Find all folders that contain at least one *_nvt_heat.gro file
all_gromacs_dirs = {os.path.dirname(f) for f in glob.glob("*/**/*_nvt_heat.gro", recursive=True)}

if not all_gromacs_dirs:
    print("[ERROR]: No matching gromacs folders found. Check file structure.")
    exit(1)

# Define the content of the Slurm script
slurm_script_template = """#!/bin/bash -e
#SBATCH --job-name={job}_MinEq
#SBATCH --account=<enter, if required>
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=128
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=36
#SBATCH --hint=nomultithread
#SBATCH --output={job}_MinEq.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<enter, if required>

module load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid

# THIS RUNS THE MINIMIZATION
srun -n 1 gmx_serial grompp -f ../templates/minimization_template.mdp -c {job}_nvt_heat.gro -r {job}_nvt_heat.gro -n index.ndx -p topol.top -o {job}_min.tpr -maxwarn 1

srun gmx_mpi mdrun -ntomp ${{SLURM_CPUS_PER_TASK}} -v -x -deffnm {job}_min

# THIS RUNS THE EQUILIBRATION
srun -n 1 gmx_serial grompp -f ../templates/equilibration_template.mdp -c {job}_min.gro -r {job}_min.gro -n index.ndx -p topol.top -o {job}_EQ.tpr -maxwarn 1

srun gmx_mpi mdrun -ntomp ${{SLURM_CPUS_PER_TASK}} -v -x -deffnm {job}_EQ
"""

# Create the scripts
for gromacs_dir in all_gromacs_dirs:
    job_name = os.path.basename(gromacs_dir)  # Use the folder name as the job name
    script_name = f"{job_name}_MinEq.sh"
    script_path = os.path.join(gromacs_dir, script_name)

    # Write the script
    with open(script_path, "w") as slurm_script:
        slurm_script.write(slurm_script_template.format(job=job_name))

    print(f"Created Slurm script {script_name} in {script_path}")

