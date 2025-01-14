import glob
import os

# Search for matching directories using glob
all_gromacs_dirs = glob.glob("*/charmm*/gromacs")  # get all gromacs folders with respective paths
print(all_gromacs_dirs)
if len(all_gromacs_dirs) == 0:
    print("[ERROR]: no gromacs folders found. check file structure.")
    exit(1)  # Exit the script with an error code


# Define the content of the Slurm script
slurm_script_content = """#!/bin/bash -e
#SBATCH --job-name={}_MinEq
#SBATCH --account=<enter, if required>
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=128
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=36
#SBATCH --hint=nomultithread
#SBATCH --output={}_MinEq.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<enter your e-mail address>


module load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid

# Extract base job name without "_MinEq" for file naming
base_job_name="{}"


# THIS RUNS THE MINIMISATION
srun -n 1 gmx_serial grompp -f step4.0_minimization.mdp -c step3_input.gro -r step3_input.gro -n index.ndx -p topol.top -o ${base_job_name}_min.tpr -maxwarn 1

srun gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -v -x -deffnm ${base_job_name}_min

# THIS RUNS THE EQUILIBRATION
srun -n 1 gmx_serial grompp -f step4.1_equilibration.mdp -c {}_min.gro -r {}_min.gro -n index.ndx -p topol.top -o step4.1_${base_job_name}_eq.tpr -maxwarn 1

srun gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -v -x -deffnm step4.1_${base_job_name}_eq 
"""

# Traverse the directory structure and create Slurm scripts in the gromacs folders
for gromacs_dir in all_gromacs_dirs:
    job_name = os.path.basename(os.path.dirname(os.path.dirname(gromacs_dir)))  # Extract the first folder name
#   print(job_name)
    script_name = job_name +"_MinEq.sh"  # Use the first folder name as the script name
#    print(script_name)
    script_path = os.path.join(gromacs_dir, script_name)
    
    # Create the Slurm script in the gromacs directory
    with open(script_path, "w") as slurm_script:
        slurm_script.write(slurm_script_content.replace("{}", job_name))
        
    print("Created Slurm script " + script_name +" in "+ script_path)

