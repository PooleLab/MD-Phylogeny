import glob
import os

all_paths = glob.glob("*/charmm*/gromacs/step3*gro")  # get all step3_input.gro files with respective paths
print(all_paths)
if len(all_paths) == 0:
    print("[ERROR]: no files found. check file structure.")
    exit(1)  # Exit the script with an error code

b = open("execute_production.sh", "w")
b.write("#!/bin/bash\nmodule load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid\n")  # header needed for indices, declare bashfile

for ff in all_paths:
    # copy the data to a new inputfile, where the folder name serves as simname
    name = ff.split("/")[0]

    # path to the respective gromacs folder, i.e. sim1/gromacs/
    path = ff.split("step3")[0]
    target_folder = ff.split("/gromacs")[0]

    b.write("#Section for simulation: %s\n"
            "cd %s\n"
            "sbatch --partition=milan %s_Prod.sh\n"
            "cd ../../..\n\n"
            % (name, path, name))

b.close()
# after this, run execute.sh to submit everything
