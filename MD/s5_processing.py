import os
import glob

# Get all step3*.gro files
all_paths = glob.glob("*/charmm*/gromacs/step3*gro")
if len(all_paths) == 0:
    print("[ERROR]: no files found. check file structure.")

# Open bash script file for writing
with open("execute_processing.sh", "w") as b:
    b.write("#!/bin/bash\nmodule load GROMACS/2020.2_CPU_single\n")

    # Iterate over each found file
    for ff in all_paths:
        # Extract the base name
        name = ff.split("/")[0]

        # Path to the respective gromacs folder
        path = ff.split("step3")[0]

        # Write processing commands to the bash script
        b.write(f"# Processing files for {name}\n")
        b.write(f"cd {path}\n")

        # Command to concatenate .xtc files
        b.write(f"gmx trjcat -f {name}.part*.xtc -o {name}.all.xtc\n")
        
        # Commands to process the files with Here Document input
        b.write(f"gmx trjconv -f {name}.all.xtc -ur compact -n index.ndx -o {name}.prot.xtc -pbc mol -center -s {name}.tpr <<EOF\n0\n0\nEOF\n")
        b.write(f"gmx trjconv -f {name}.prot.xtc -n index.ndx -o {name}.proc.xtc -s {name}.tpr -fit rot+trans <<EOF\n0\n0\nEOF\n")
        b.write(f"gmx trjconv -f {name}.prot.xtc -b 0 -e 0 -n index.ndx -o {name}.proc.gro -s {name}.tpr <<EOF\n0\nEOF\n")

        # Change back to the main folder
        b.write("cd ../../..\n\n")

