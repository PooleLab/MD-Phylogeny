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
#SBATCH --job-name={}_Prod
#SBATCH --account=<enter if required>
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=512
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=12
#SBATCH --hint=nomultithread
#SBATCH --output={}_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<enter e-mail address>
#SBATCH --array=0-6%1

# Extract base job name without "_MinEq" for file naming
base_job_name="{}"

module load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid

if [ $SLURM_ARRAY_TASK_ID -eq 0 ];then

srun -n 1 gmx_serial grompp -f step5_production.mdp -c step4.1_${base_job_name}_eq.gro -r step4.1_${base_job_name}_eq.gro -n index.ndx -t step4.1_${base_job_name}_eq.cpt -p topol.top -o ${base_job_name}.tpr -maxwarn 2


srun gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -v -x -deffnm ${base_job_name} -noappend

else

srun gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -v -x -deffnm ${base_job_name} -noappend -cpi ${base_job_name}.cpt

fi
"""

# Traverse the directory structure and create Slurm scripts in the gromacs folders
for gromacs_dir in all_gromacs_dirs:
    job_name = os.path.basename(os.path.dirname(os.path.dirname(gromacs_dir)))  # Extract the first folder name
#   print(job_name)
    script_name = job_name +"_Prod.sh"  # Use the first folder name as the script name
#    print(script_name)
    script_path = os.path.join(gromacs_dir, script_name)

    # Create the Slurm script in the gromacs directory
    with open(script_path, "w") as slurm_script:
        slurm_script.write(slurm_script_content.replace("{}", job_name))

    print("Created Slurm script " + script_name +" in "+ script_path)
