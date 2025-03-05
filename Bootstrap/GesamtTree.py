#!python
import os
import sys
import glob
import subprocess
import argparse
import shutil
import numpy as np
from Bio import Phylo
from concurrent.futures import ThreadPoolExecutor
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix

def run_gesamt(file1, file2):
    """
    Run the 'gesamt' command on two PDB files and extract the Q-score.
    Returns the Q-score as a string, or '0' if not found or in case of errors.
    """
    print(f"Running 'gesamt' for {file1} and {file2}")

    try:
        result = subprocess.run(
            ['gesamt', file1, file2],
            capture_output=True,
            text=True,
            check=True
        )
        for line in result.stdout.splitlines():
            if "Q-score" in line:
                parts = line.split()
                if len(parts) >= 3:
                    return parts[2]
        return '0'
    except subprocess.CalledProcessError as e:
        print(f"Error running 'gesamt' for {file1} and {file2}: {e}", file=sys.stderr)
        sys.exit()
    except Exception as e:
        print(f"Unexpected error running 'gesamt' for {file1} and {file2}: {e}", file=sys.stderr)
        sys.exit()

def run_gesamt_task(pair):
    """
    Wrapper function for parallel execution.
    """
    file1, file2 = pair
    q_score = run_gesamt(file1, file2)
    return file1, file2, q_score

def generate_resDB(num_threads):
    """
    Generate a 'resDB' file by comparing all PDB files using the 'gesamt' command
    in parallel using the specified number of threads.
    """
    pdb_files = glob.glob('*.pdb')
    if not pdb_files:
        print("No PDB files found in the directory.")
        sys.exit(1)

    print("Preparing tasks for gesamt comparisons...")
    # Create pairs of pdb files (maintaining nested loop order)
    tasks = [(file1, file2) for file1 in pdb_files for file2 in pdb_files]
    print(f"Running gesamt in parallel using {num_threads} threads for {len(tasks)} comparisons...")

    # Run tasks in parallel
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(run_gesamt_task, tasks))

    print("Writing results to resDB file...")
    with open('resDB', 'w') as f:
        for file1, file2, q_score in results:
            f.write('\t'.join([file1.split('.')[0], file2.split('.')[0], q_score]) + '\n')
    print("resDB file generated.")


def build_distance_matrix():
    """
    Build a distance matrix using the data in 'resDB'.
    Returns a tuple (unique_ids, matrix).
    """
    print("Building distance matrix from resDB file...")
    fl = 'resDB'
    d0 = np.genfromtxt(fl, delimiter='\t', dtype=str, usecols=0, comments="~")
    d1 = np.genfromtxt(fl, delimiter='\t', dtype=str, usecols=1, comments="~")
    d2 = np.genfromtxt(fl, delimiter='\t', dtype=str, usecols=2, comments="~")
    
    merged = np.array([f"{i}-{j}" for i, j in zip(d0, d1)])
    unique_ids = np.unique(d0)
    n = len(unique_ids)
    matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            if i > j:
                pat1 = unique_ids[i] + '-' + unique_ids[j]
                pat2 = unique_ids[j] + '-' + unique_ids[i]
                idx1 = np.where(merged == pat1)[0]
                idx2 = np.where(merged == pat2)[0]
                v1 = v2 = 0.0
                if len(idx1) == 1:
                    v1 = float(d2[idx1][0])
                if len(idx2) == 1:
                    v2 = float(d2[idx2][0])
                if len(idx1) == 1 and len(idx2) == 0:
                    v2 = v1
                if len(idx2) == 1 and len(idx1) == 0:
                    v1 = v2
                val = (v1 + v2) / 2
                A = round(1 - val, 5)
                matrix[i, j] = A
                matrix[j, i] = A
    print("Distance matrix built.")
    return unique_ids, matrix

def generate_phylo_data(unique_ids, matrix):
    """
    Generate phylogenetic data and write output files including Data.nex, Taxon.txt,
    tree.nex, and TreeFormatted.nex.
    """
    print("Generating phylogenetic data...")
    ntax = len(unique_ids)
    phy = f'#nexus\nBEGIN taxa;\nDIMENSIONS ntax={ntax};\nTAXLABELS\n'
    label_list = []
    
    for idx, uid in enumerate(unique_ids, start=1):
        label = '_'.join(uid.split('_')[0:2])
        phy += f'[{idx}]{label}\n'
        label_list.append(label)
    phy += ';\nEND [taxa];\nBEGIN distances;\n'
    phy += f'DIMENSIONS ntax={ntax};\nFORMAT labels diagonal triangle=both;\nMATRIX\n'
    
    for idx, row in enumerate(matrix, start=1):
        label = '_'.join(unique_ids[idx-1].split('_')[0:2])
        row_str = '\t'.join(map(str, row))
        phy += f'[{idx}]{label}\t{row_str}\n'
    
    # Prepare distance matrix for NJ tree construction
    tril_matrix = np.tril(matrix)
    dist_list = [list(tril_matrix[x, :x+1]) for x in range(ntax)]
    mm_ = _DistanceMatrix(label_list, dist_list)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(mm_)
    
    # Write tree in Newick format
    Phylo.write(tree, 'tree.nex', 'newick')
    os.system("sed -i 's,:-[0-9\\.]+,:0.0,g' tree.nex")
    
    tree_str = ""
    with open('tree.nex') as f:
        for line in f:
            tree_str += line.strip()
    
    # Write output files
    with open('Data.nex', 'w') as f:
        f.write(phy)
    with open('Taxon.txt', 'w') as f:
        for label in label_list:
            f.write(label + '\n')
    with open('TreeFormatted.nex', 'w') as f:
        f.write("#NEXUS\nBEGIN TREES;\n  Tree tree1 = " + tree_str + "\nEND;\n")
    os.system("sed -i 's,:-[0-9\\.]+,:0.0,g' TreeFormatted.nex")
    print("Phylogenetic data generated and output files created.")

def main():
    parser = argparse.ArgumentParser(
        description="Process PDB files to build a Neighbor-Joining tree using gesamt."
    )
    parser.add_argument("data", help="Path to the directory containing PDB files.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    args = parser.parse_args()
    
    # Check if 'gesamt' is executable
    if not shutil.which("gesamt"):
        print("Error: 'gesamt' command not found. Please ensure it is installed and in your $PATH.", file=sys.stderr)
        sys.exit(1)
    
    original_path = os.getcwd()
    target_dir = os.path.abspath(args.data)
    
    # Change to target directory
    try:
        os.chdir(target_dir)
        print(f"Changed directory to {target_dir}")
    except Exception as e:
        print(f"Error changing directory to {target_dir}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Run processing steps
    generate_resDB(args.threads)
    unique_ids, matrix = build_distance_matrix()
    generate_phylo_data(unique_ids, matrix)
    
    # Return to original directory
    os.chdir(original_path)
    print("Processing completed successfully.")

    os.remove("resDB")

if __name__ == '__main__':
    main()
