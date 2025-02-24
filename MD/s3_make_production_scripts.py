import glob
import os

# Search for directories that contain the *_EQ.gro file
all_gromacs_dirs = glob.glob("*/")  # Get all top-level folders (protein folders)

# Filter out directories that don't contain *_EQ.gro files
gromacs_dirs_with_eq = [d for d in all_gromacs_dirs if glob.glob(os.path.join(d, '*_EQ.gro'))]

print(gromacs_dirs_with_eq)
if len(gromacs_dirs_with_eq) == 0:
    print("[ERROR]: no folders with *_EQ.gro files found. Check file structure.")
    exit(1)  # Exit the script with an error code

# Define the content of the Slurm script
slurm_script_content = """#!/bin/bash -e
#SBATCH --job-name={}_Prod
#SBATCH --account=<enter, if required>
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=512
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=12
#SBATCH --hint=nomultithread
#SBATCH --output={}_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<enter, if required>
#SBATCH --array=0-2%1

# Extract base job name; name of the protein
base_job_name="{}"

module load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid

if [ $SLURM_ARRAY_TASK_ID -eq 0 ];then

srun -n 1 gmx_serial grompp -f ../templates/production_template.mdp -c ${base_job_name}_EQ.gro -r ${base_job_name}_EQ.gro -n index.ndx -t ${base_job_name}_EQ.cpt -p topol.top -o ${base_job_name}.tpr -maxwarn 2

srun gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -v -x -deffnm ${base_job_name} -noappend

else

srun gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -v -x -deffnm ${base_job_name} -noappend -cpi ${base_job_name}.cpt

fi
"""

# Traverse the directory structure and create Slurm scripts in the protein folders
for gromacs_dir in gromacs_dirs_with_eq:
    # Get the folder's name (protein name)
    job_name = os.path.basename(os.path.dirname(gromacs_dir))  # protein folder name
    script_name = job_name + "_Prod.sh"  # Create script name based on protein folder name
    script_path = os.path.join(gromacs_dir, script_name)

    # Create the Slurm script in the protein folder
    with open(script_path, "w") as slurm_script:
        slurm_script.write(slurm_script_content.replace("{}", job_name))

    print(f"Created Slurm script {script_name} in {script_path}")

