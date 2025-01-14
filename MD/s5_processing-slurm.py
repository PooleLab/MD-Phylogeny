import glob

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
''')  # this is the header for the bashfile

# Find the TPR files to detect the runs
files = sorted([ff for ff in glob.glob("*/charmm*/gromacs/*tpr") if "step" not in ff])
print("Found total of", len(files), "simulation files as TPR run inputs.")

# Extract unique protein folders
protein_folders = set(ff.split('/')[0] for ff in files)
print("Found total of", len(protein_folders), "unique protein folders.")

parallel = 0  # indicator for most efficient parallel threading on the node
for ff in files:
    file_prefix = ff.split(".")[0]  # l1_mut_scan/e774a/gromacs/e774a_1

    # Get all part files except those with '_min' or '_MinEq_min'
    extra_files = sorted(glob.glob(ff.replace(".tpr", ".part*.xtc")))
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

bashfile.write("wait\n\n#Start of processing section: 1) extract the protein, check PBCs\n\n")

parallel = 0
for ff in files:
    if '_min' in ff or '_MinEq_min' in ff:
        continue  # Skip files with _min or _MinEq_min

    file_prefix = ff.split(".")[0]  # l1_mut_scan/e774a/gromacs/e774a_1
    path_parts = ff.split("/")  # ["l1_mut_scan","e774a","gromacs","e774a_1..."]
    file_path = "/".join(path_parts[:-1])  # l1_mut_scan/e774a/gromacs
    parallel += 1
    bashfile.write("echo -e '0\\n0\\n'| gmx trjconv -f %s.all.xtc -ur compact -n %s/index.ndx -o %s.prot.xtc -pbc mol -center -s %s" % (file_prefix, file_path, file_prefix, ff))
    if parallel != 6:
        bashfile.write(" &\n")
    else:
        bashfile.write("\nwait\n")
        parallel = 0
bashfile.write("wait\n\n#Start of processing section: 2) fit the protein rot+trans\n\n")

parallel = 0
for ff in files:
    if '_min' in ff or '_MinEq_min' in ff:
        continue  # Skip files with _min or _MinEq_min

    file_prefix = ff.split(".")[0]  # l1_mut_scan/e774a/gromacs/e774a_1
    path_parts = ff.split("/")  # ["l1_mut_scan","e774a","gromacs","e774a_1..."]
    file_path = "/".join(path_parts[:-1])  # l1_mut_scan/e774a/gromacs
    parallel += 1
    bashfile.write("echo -e '0\\n0\\n' | gmx trjconv -f %s.prot.xtc -n %s/index.ndx -o %s.proc.xtc -s %s -fit rot+trans" % (file_prefix, file_path, file_prefix, ff))
    if parallel != 6:
        bashfile.write(" &\n")
    else:
        bashfile.write("\nwait\n")
        parallel = 0

bashfile.write("wait\n")

bashfile.write("wait\n\n#Generate the .gro files\n\n")

parallel = 0
for ff in files:
    if '_min' in ff or '_MinEq_min' in ff:
        continue  # Skip files with _min or _MinEq_min

    file_prefix = ff.split(".")[0]  # l1_mut_scan/e774a/gromacs/e774a_1
    path_parts = ff.split("/")  # ["l1_mut_scan","e774a","gromacs","e774a_1..."]
    file_path = "/".join(path_parts[:-1])  # l1_mut_scan/e774a/gromacs
    parallel += 1
    bashfile.write("echo -e '0\\n' | gmx trjconv -f %s.prot.xtc -b 0 -e 0 -n %s/index.ndx -o %s.proc.gro -s %s" % (file_prefix, file_path, file_prefix, ff))
    if parallel != 6:
        bashfile.write(" &\n")
    else:
        bashfile.write("\nwait\n")
        parallel = 0

bashfile.write("wait\n")

bashfile.close()
print("Finished writing processing.sh.")

