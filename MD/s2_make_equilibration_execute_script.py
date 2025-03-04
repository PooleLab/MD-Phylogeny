import glob
import os

# Get all *_nvt_heat.gro files within protein folders
all_paths = glob.glob("*/*nvt_heat.gro")  # Adjusted to match the new file pattern
print(all_paths)

if len(all_paths) == 0:
    print("[ERROR]: no files found. Check the file structure.")
    exit(1)  # Exit the script with an error code

b = open("execute_equilibration.sh", "w")
b.write("#!/bin/bash\nmodule load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid\n")  # Header for bash file

for ff in all_paths:
    # Extract the protein folder name (simname) from the path
    name = ff.split("/")[0]

    # Path to the respective gromacs folder, i.e., protein_name/gromacs/
    path = ff.split("nvt_heat.gro")[0]
    target_folder = ff.split("/nvt_heat.gro")[0]

    # Adjusting the directory structure to use the protein folder name only
    b.write("# Section for simulation: %s\n"
            "cd %s/\n"
            "sbatch %s_MinEq.sh\n"
            "cd ..\n\n"
            % (name, name, name))

b.close()

# After running the script, execute execute_equilibration.sh to submit all jobs

