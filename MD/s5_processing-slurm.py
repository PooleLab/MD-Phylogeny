import glob
import os

# Prepare the bash script for writing
bashfile = open("processing.sh", 'w')
bashfile.write('''#!/bin/bash
#SBATCH --job-name=processing
#SBATCH --account=<enter if required>
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=512
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=12
#SBATCH --hint=nomultithread
#SBATCH --output=processing.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<enter e-mail address>

module load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid
''')  # This is the header for the bashfile

# Find the simulation files inside each protein folder
files = sorted([ff for ff in glob.glob("*/" + "*") if ff.endswith(".xtc")])  # Assuming .xtc as the file type
print("Found total of", len(files), "simulation files.")

# Extract unique protein folders (just the top-level directories for each protein)
protein_folders = set(ff.split('/')[0] for ff in files)
print("Found total of", len(protein_folders), "unique protein folders.")

parallel = 0  # Indicator for most efficient parallel threading on the node
for ff in files:
    protein_folder = ff.split("/")[0]  # Protein folder (top level directory)
    file_prefix = ff.split(".")[0]  # Protein folder + protein name

    # Get all part files except those with '_min' or '_MinEq_min'
    extra_files = sorted(glob.glob(os.path.join(protein_folder, "*.part*.xtc")))
    extra_files = [f for f in extra_files if '_min' not in f and '_MinEq_min' not in f]

    if extra_files:  # Only proceed if there are valid extra files
        # Construct the command with only part files
        command = "gmx trjcat -f " + " ".join(extra_files) + " -o " + file_prefix + ".all.xtc"
        parallel += 1
        bashfile.write(command)
        if parallel != 6:
            bashfile.write(" &\n")
        else:
            bashfile.write("\nwait\n")
            parallel = 0

bashfile.write("wait\n\n# Start of processing section: 1) Extract the protein, check PBCs\n\n")

parallel = 0
for ff in files:
    if '_min' in ff or '_MinEq_min' in ff:
        continue  # Skip files with _min or _MinEq_min

    protein_folder = ff.split("/")[0]  # Protein folder
    file_prefix = ff.split(".")[0]  # Protein folder + protein name
    parallel += 1
    bashfile.write("echo -e '0\\n0\\n' | gmx trjconv -f %s.all.xtc -ur compact -n %s/index.ndx -o %s.prot.xtc -pbc mol -center -s %s\n" % (file_prefix, protein_folder, file_prefix, ff))
    if parallel != 6:
        bashfile.write(" &\n")
    else:
        bashfile.write("\nwait\n")
        parallel = 0

bashfile.write("wait\n\n# Start of processing section: 2) Fit the protein rot+trans\n\n")

parallel = 0
for ff in files:
    if '_min' in ff or '_MinEq_min' in ff:
        continue  # Skip files with _min or _MinEq_min

    protein_folder = ff.split("/")[0]  # Protein folder
    file_prefix = ff.split(".")[0]  # Protein folder + protein name
    parallel += 1
    bashfile.write("echo -e '0\\n0\\n' | gmx trjconv -f %s.prot.xtc -n %s/index.ndx -o %s.proc.xtc -s %s -fit rot+trans\n" % (file_prefix, protein_folder, file_prefix, ff))
    if parallel != 6:
        bashfile.write(" &\n")
    else:
        bashfile.write("\nwait\n")
        parallel = 0

bashfile.write("wait\n")

bashfile.write("wait\n\n# Generate the .gro files\n\n")

parallel = 0
for ff in files:
    if '_min' in ff or '_MinEq_min' in ff:
        continue  # Skip files with _min or _MinEq_min

    protein_folder = ff.split("/")[0]  # Protein folder
    file_prefix = ff.split(".")[0]  # Protein folder + protein name
    parallel += 1
    bashfile.write("echo -e '0\\n' | gmx trjconv -f %s.prot.xtc -b 0 -e 0 -n %s/index.ndx -o %s.proc.gro -s %s\n" % (file_prefix, protein_folder, file_prefix, ff))
    if parallel != 6:
        bashfile.write(" &\n")
    else:
        bashfile.write("\nwait\n")
        parallel = 0

bashfile.write("wait\n")

bashfile.close()
print("Finished writing processing.sh.")

