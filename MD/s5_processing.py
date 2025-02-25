import os
import glob

# Get all directories at the top level (assuming you are working from the top directory)
all_dirs = glob.glob("*/")

# Open bash script file for writing
with open("execute_processing.sh", "w") as b:
    b.write("#!/bin/bash\nmodule load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid\n")

    # Iterate over each found directory (one per protein)
    for dir in all_dirs:
        # Debugging: Print the directory being checked
        print(f"Checking directory: {dir}")

        # Check if there are any part*.xtc files in the directory
        part_files = glob.glob(os.path.join(dir, "*part*.xtc"))
        
        # Debugging: Check if part files are found
        if not part_files:
            print(f"No part*.xtc files found in {dir}, skipping.")
            continue  # Skip if no part*.xtc files found

        # Extract the base name from the directory
        name = dir.split("/")[0]

        # Write processing commands to the bash script
        b.write(f"# Processing files for {name}\n")
        b.write(f"cd {dir}\n")

        # Step 1: Concatenate the part*.xtc files into one all.xtc
        b.write(f"gmx trjcat -f {name}.part*.xtc -o {name}.all.xtc\n")

        # Step 2: Extract protein from the concatenated trajectory
        b.write(f"gmx trjconv -f {name}.all.xtc -ur compact -n index.ndx -o {name}.prot.xtc -pbc mol -center -s {name}.tpr <<EOF\n0\n0\nEOF\n")

        # Step 3: Fit the protein rotation and translation
        b.write(f"gmx trjconv -f {name}.prot.xtc -n index.ndx -o {name}.proc.xtc -s {name}.tpr -fit rot+trans <<EOF\n0\n0\nEOF\n")

        # Step 4: Save the processed trajectory as .gro file
        b.write(f"gmx trjconv -f {name}.proc.xtc -b 0 -e 0 -n index.ndx -o {name}.proc.gro -s {name}.tpr <<EOF\n0\nEOF\n")

        # Change back to the main folder
        b.write("cd ..\n\n")

# Debugging: Print a message that the script has finished running
print("Script execution complete.")

