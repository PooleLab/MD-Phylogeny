#!python
import os
import sys
import glob
import subprocess
import argparse
import shutil
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import re

def run_gesamt(file1, file2):
    """
    Run the 'gesamt' command on two PDB files and extract the Q-score.
    Returns the Q-score as a string, or '0' if not found or in case of errors.
    """
    # print(f"Running 'gesamt' for {file1} and {file2}") # Commented out for cleaner progress bar output

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
    tasks = [(file1, file2) for file1 in pdb_files for file2 in pdb_files if file1 <= file2]
    num_tasks = len(tasks)
    print(f"Running gesamt in parallel using {num_threads} threads for {num_tasks} comparisons...")

    results = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(run_gesamt_task, task) for task in tasks]
        
        for i, future in enumerate(as_completed(futures)):
            results.append(future.result())
            progress = (i + 1) / num_tasks
            bar_length = 40
            block = int(round(bar_length * progress))
            text = f"\rProgress: [{'#' * block}{'-' * (bar_length - block)}] {progress*100:.1f}% ({i+1}/{num_tasks})"
            sys.stdout.write(text)
            sys.stdout.flush()

    sys.stdout.write('\n')

    print("Writing results to resDB file...")
    with open('resDB', 'w') as f:
        # Results may not be in order, but that's fine for the next step
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

    q_scores = {}
    all_ids = set()

    try:
        # Read resDB into a dictionary for fast lookups
        with open(fl, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("~"):
                    continue
                try:
                    id1, id2, q_score_str = line.split('\t')
                    q_scores[f"{id1}-{id2}"] = float(q_score_str)
                    q_scores[f"{id2}-{id1}"] = float(q_score_str)
                    all_ids.add(id1)
                    all_ids.add(id2)
                except (ValueError, IndexError):
                    print(f"Skipping malformed line: {line}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: '{fl}' file not found.", file=sys.stderr)
        sys.exit(1)

    if not all_ids:
        print("No data found in resDB file.")
        return np.array([]), np.zeros((0, 0))

    label_list = [re.sub(r"[:,()\n\t]","_",x) for x in sorted(list(all_ids))]
    unique_ids = np.array(range(len(label_list)))
    n = len(unique_ids)
    matrix = np.zeros((n, n))
    
    with open("label.list", 'w') as f:
        for i in range(n):
            f.write(f"{label_list[i]}" + "\n")
    print(f"(intermediate) label list: {os.getcwd()}/label.list")

    # Populate the distance matrix
    num_pairs = n * (n - 1) // 2
    if num_pairs > 0:
        print("Populating distance matrix...")
        processed_pairs = 0
        for i in range(n):
            for j in range(i + 1, n):
                uid_i = label_list[unique_ids[i]]
                uid_j = label_list[unique_ids[j]]

                v1 = q_scores.get(f"{uid_i}-{uid_j}", 0.0)
                v2 = q_scores.get(f"{uid_j}-{uid_i}", 0.0)

                # Average the q-scores, handling cases where one might be missing
                if v1 and not v2:
                    v2 = v1
                elif v2 and not v1:
                    v1 = v2

                val = (v1 + v2) / 2.0
                A = round(1 - val, 5)
                matrix[i, j] = A
                matrix[j, i] = A
                
                processed_pairs += 1
                progress = processed_pairs / num_pairs
                bar_length = 40
                block = int(round(bar_length * progress))
                text = f"\rProgress: [{'#' * block}{'-' * (bar_length - block)}] {progress*100:.1f}% ({processed_pairs}/{num_pairs})"
                sys.stdout.write(text)
                sys.stdout.flush()
        sys.stdout.write('\n')

    print(f"(intermediate) distance matrix built: {os.getcwd()}/resDB")

    # Write distance matrix in PHYLIP format for RapidNJ
    phylip_file = 'distances.phy'
    with open(phylip_file, 'w') as f:
        f.write(f"    {n}\n")
        for i in range(n):
            # PHYLIP format requires truncating labels. 10 chars is standard.
            label = unique_ids[i]
            f.write(f"{label:<10} " + " ".join(map(str, matrix[i])) + "\n")
    print(f"(intermediate) distance matrix for RapidNJ created: {os.getcwd()}/{phylip_file}")
    phylip_file = 'distances.phy'

def generate_phylo_data(num_threads):
    """
    Generate phylogenetic data using RapidNJ for fast, parallelized 
    Neighbor-Joining tree construction.
    """
    print("Generating phylogenetic data using RapidNJ...")
    phylip_file = 'distances.phy'
    
    # Run RapidNJ
    output_tree_file = 'tree.nwk'
    print(f"Running RapidNJ with {num_threads} threads...")
    try:
        command = ['rapidnj', phylip_file, '-i', 'pd', '-o', 't', '-x', output_tree_file, '-c', str(num_threads)]
        subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running RapidNJ: {e}", file=sys.stderr)
        print(f"RapidNJ stdout:\n{e.stdout}", file=sys.stderr)
        print(f"RapidNJ stderr:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)

    # Post-process and save final tree files
    os.system(f"sed -i 's,:-[0-9\\.]+,:0.0,g' {output_tree_file}")
    
    label_list=[]
    with open("label.list", 'r') as f:
        label_list = [ line.rstrip() for line in f ]
    n = len(label_list)

    with open(output_tree_file) as f:
        tree_str = f.read().strip()

    def replace_labels_in_newick(tree_str, label_list):
        def repl(match):
            idx = int(match.group(1))
            return(label_list[idx] if 0 <= idx < len(label_list) else match.group(0))
        return re.sub(r'(?<=\(|,)\'(\d+)\'(?=[:\),])', repl, tree_str)

    tree_str_labeled = replace_labels_in_newick(tree_str, label_list)

    # Save updated Newick file
    with open(output_tree_file, 'w') as f:
        f.write(tree_str_labeled + '\n')

    print(f"(result) phylogenetic data output created (newick format): {os.getcwd()}/{output_tree_file}")
    
    with open('tree.nex', 'w') as f:
        f.write("#NEXUS\nBEGIN TREES;\n")
        for idx, label in enumerate(label_list):
            f.write(f'[{idx}]{label}\n')
        f.write("  Tree tree1 = " + tree_str + "\nEND;\n")
    os.system(f"sed -i 's,:-[0-9\\.]+,:0.0,g' tree.nex")
    print(f"(result) phylogenetic data output created (nexus format): {os.getcwd()}/tree.nex")

    # # Also write the Nexus and Taxon files for compatibility/inspection
    # phy = f'#nexus\nBEGIN taxa;\nDIMENSIONS ntax={n};\nTAXLABELS\n'
    # for idx, label in enumerate(label_list):
    #     phy += f'[{idx}]{label}\n'
    # phy += ';\nEND [taxa];\nBEGIN distances;\n'
    # phy += f'DIMENSIONS ntax={n};\nFORMAT labels diagonal triangle=both;\nMATRIX\n'
    # for idx, row in enumerate(label_list):
    #     label = label_list[idx-1]
    #     row_str = '\t'.join(map(str, row))
    #     phy += f'[{idx}]{label}\t{row_str}\n'
    
    # with open('Data.nex', 'w') as f:
    #     f.write(phy)
    # print(f"(intermediate) input data created (nexus format): {os.getcwd()}/Data.nex")

def main():
    parser = argparse.ArgumentParser(
        description="Process PDB files to build a Neighbor-Joining tree using gesamt."
    )
    parser.add_argument("data", help="Path to the directory containing PDB files. All intermediate and result files are written to to this directory as well. Main resulting files are: tree.nwk (newick) and tree.nex (nexus)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    args = parser.parse_args()
    
    # Check if executables are available
    if not shutil.which("gesamt"):
        print("Error: 'gesamt' command not found. Please ensure it is installed and in your $PATH.", file=sys.stderr)
        sys.exit(1)
    if not shutil.which("rapidnj"):
        print("Error: 'rapidnj' command not found. Please ensure it is installed and in your $PATH.", file=sys.stderr)
        print("This script uses RapidNJ for fast, parallelized tree calculation.", file=sys.stderr)
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
    if not os.path.exists("resDB"):
        generate_resDB(args.threads)
    else:
        print(f"(intermedate) found a intermediate {os.getcwd()}/resDB file, I will use this...")

    if not os.path.exists("distances.phy"):
        build_distance_matrix()
    else:
        print(f"(intermedate) found a intermediate {os.getcwd()}/distances.phy file, I will use this...")

    if not os.path.exists("tree.nwk"):
        generate_phylo_data(args.threads)
    else:
        print(f"(result) found phylogenetic data output (nexus format): {os.getcwd()}/tree.nex")
        print(f"(result) found phylogenetic data output (newick format): {os.getcwd()}/tree.nwk")

    # Return to original directory
    os.chdir(original_path)
    print("Processing completed successfully.")

if __name__ == '__main__':
    main()
