import glob
import os
import shutil


# Define the content of the Slurm script with automation for force field and water model selection
slurm_script_template = """#!/bin/bash -e
#SBATCH --job-name={base_job_name}_Sim_setup
#SBATCH --account=<enter, if required>
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=128
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --output={base_job_name}_Sim_setup.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<enter, if required>

module load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid

# Extract base job name 
base_job_name="{base_job_name}"

# THIS RUNS THE FILE CONVERSION FROM PDB TO GRO
# Automatically select CHARMM36 (option 8) and TIP3 water model (option 1)
echo -e "8\\n1" | gmx pdb2gmx -f {base_job_name}.pdb -o {base_job_name}.gro -ignh

# ADD SIMULATION BOX
gmx editconf -f {base_job_name}.gro -o {base_job_name}_box.gro -bt cubic -d 1.2 -c

# ENERGY MINIMISATION
gmx grompp -f ../min_sd.mdp -c {base_job_name}_box.gro -p topol.top -o {base_job_name}_minsd.tpr
gmx mdrun -deffnm {base_job_name}_minsd

# ADD WATER TO THE BOX
gmx solvate -cp {base_job_name}_minsd.gro -cs spc216.gro -p topol.top -o {base_job_name}_solv.gro

# ENERGY MINIMISATION
gmx grompp -f ../min_sd.mdp -c {base_job_name}_solv.gro -p topol.top -o {base_job_name}_solv_minsd.tpr -maxwarn 1
gmx mdrun -deffnm {base_job_name}_solv_minsd 

# ADD IONS
gmx grompp -f ../ions.mdp -c {base_job_name}_solv_minsd.gro -p topol.top -o {base_job_name}_ions.tpr -maxwarn 1
echo -e "13\\n" | gmx genion -s {base_job_name}_ions.tpr -p topol.top -o {base_job_name}_ions.gro -neutral

# MAKE INDEX FILE
echo -e "q\\n" | gmx make_ndx -f {base_job_name}_ions.gro -o index.ndx
# type q and enter

# ENERGY MINIMISATION
gmx grompp -f ../min_sd.mdp -c {base_job_name}_ions.gro -p topol.top -o {base_job_name}_ions_minsd.tpr -maxwarn 1
gmx mdrun -deffnm {base_job_name}_ions_minsd

# INITIALISATION AND HEATING
gmx grompp -f ../charmm-inputs/charmm_nvt_heat.mdp -c {base_job_name}_ions_minsd.gro -r {base_job_name}_ions_minsd.gro -p topol.top -o {base_job_name}_nvt_heat.tpr
gmx mdrun -deffnm {base_job_name}_nvt_heat
"""

# Get all .pdb files in the current directory
pdb_files = glob.glob("*.pdb")

# Loop through each .pdb file
for pdb_file in pdb_files:
    # Get the base name without the .pdb extension
    base_name = os.path.splitext(pdb_file)[0]

    # Define the output folder name
    output_folder = base_name

    # Check if the folder exists
    if not os.path.exists(output_folder):
        # Create the folder if it doesn't exist
        os.makedirs(output_folder)
        print(f"Created folder {output_folder} for {pdb_file}.")
    else:
        print(f"Folder {output_folder} already exists, skipping folder creation.")

    # Check if the .pdb file is already in the folder
    target_pdb_path = os.path.join(output_folder, pdb_file)
    if not os.path.exists(target_pdb_path):
        shutil.copy(pdb_file, target_pdb_path)
        print(f"Copied {pdb_file} to {output_folder}/.")
    else:
        print(f"{pdb_file} already exists in {output_folder}/, skipping copy.")

    # Define the path to the Slurm script inside the folder
    slurm_script_path = os.path.join(output_folder, f"{base_name}_sim_setup.sh")

    # Create the Slurm script content by replacing placeholders with the base_name
    slurm_script_content_filled = slurm_script_template.format(base_job_name=base_name)

    # Write the Slurm script to the file
    with open(slurm_script_path, "w") as slurm_script:
        slurm_script.write(slurm_script_content_filled)

    print(f"Created Slurm script for {pdb_file} in {output_folder}/.")
