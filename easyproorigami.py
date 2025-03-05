#!/usr/bin/python3
#pk
import os
import random
import string
import re
import sys
import subprocess
import requests
from concurrent.futures import ThreadPoolExecutor, TimeoutError, as_completed
import shutil
import tempfile
import argparse

version = "0.0.25" # last updated on 2025-02-24 09:36:12

package_url = 'https://stivalaa.github.io/AcademicWebsite/software/proorigami-cde-package.tar.gz'
package_tar = 'proorigami.tar.gz'
package_dir = 'proorigami'
master_dir = os.path.expanduser("~")  # default; may be overridden via argparse

def execute_process(pdb_file, ptgraph_opts, force):
    bn = os.path.basename(pdb_file)
    dirn = os.path.dirname(pdb_file)
    prefix = random_string()
    bn_base = re.sub(r"[.][^.]+$", "", bn)

    output_png = os.path.join(dirn, f'{bn_base}.png')
    output_svg = os.path.join(dirn, f'{bn_base}.svg')
    if os.path.exists(output_png):
        if force:
            print(f"Force enabled: Overwriting existing output for {pdb_file}")
            try:
                os.remove(output_png)
            except Exception as e:
                print(f"Error removing {output_png}: {e}")
            if os.path.exists(output_svg):
                try:
                    os.remove(output_svg)
                except Exception as e:
                    print(f"Error removing {output_svg}: {e}")
        else:
            print(f"Skipping {pdb_file} because output file exists (use --force to overwrite)")
            return

    print(f"Processing {pdb_file}")

    seed = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
    # Create a temporary directory for this process and copy the package
    with tempfile.TemporaryDirectory() as tmp_dir:
        proorigami_tmp = os.path.join(tmp_dir, f"proorigami_{seed}")
        shutil.copytree(os.path.join(master_dir, package_dir), proorigami_tmp)
        cwd = os.path.join(proorigami_tmp, 'cde-root', 'home', 'proorigami')

        linkp = os.path.join(cwd, f"{prefix}_{bn}")
        os.symlink(os.path.realpath(pdb_file), linkp)

        # Build the command (without appending ptgraph options as arguments)
        cmd = ["bash", "make_cartoon.sh.cde", os.path.basename(linkp)]
        # Prepare environment with PTGRAPH2_OPTIONS
        env = os.environ.copy()
        env["PTGRAPH2_OPTIONS"] = " ".join(ptgraph_opts) if ptgraph_opts else ""
        print("Running command:", " ".join(cmd), " with options: PTGRAPH2_OPTIONS=",env["PTGRAPH2_OPTIONS"])
        p = subprocess.Popen(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env
        )
        stdout, stderr = p.communicate()
        stdout = stdout.decode()
        stderr = stderr.decode()

        # Check for results and move files back to the original directory
        res = re.search(r"Bitmap saved as: (.*)[.]png", stdout + stderr)
        res_base = res.group(1) if res else ""
        if res_base:
            shutil.move(os.path.join(cwd, f'{res_base}.png'), output_png)
            shutil.move(os.path.join(cwd, f'{res_base}.svg'), output_svg)
            if os.path.exists(os.path.join(cwd, f'{res_base}.pml')):
                shutil.move(os.path.join(cwd, f'{res_base}.pml'), os.path.join(dirn, f'{bn_base}.pml'))
            if os.path.exists(os.path.join(cwd, f'{res_base}.M')):
                shutil.move(os.path.join(cwd, f'{res_base}.M'), os.path.join(dirn, f'{bn_base}.M'))
        else:
            raise Exception(f"{stdout}\n{stderr}\n{pdb_file} failed")

    print(f"Finished processing {pdb_file}")


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def download_package(url, dest):
    print(f"Downloading {url}...")
    response = requests.get(url)
    with open(dest, 'wb') as f:
        f.write(response.content)
    print(f"Downloaded {dest}.")

def extract_package(tar_path, extract_to):
    print(f"Extracting {tar_path}...")
    os.makedirs(extract_to, exist_ok=True)
    subprocess.run(['tar', 'zxf', tar_path, '-C', extract_to, '--strip-components=1'])
    print(f"Extracted to {extract_to}.")

def check_cde_exec(exec_path):
    return os.path.isfile(exec_path) and os.access(exec_path, os.X_OK)

def random_string(length=8):
    return ''.join(random.choices(string.ascii_letters + string.digits, k=length))

def modify_make_cartoon_script(master_dir):
    """
    Modify the make_cartoon.sh.cde script so that it first checks whether
    PTGRAPH2_OPTIONS is defined in the environment. If not, it sets it to the default.
    """
    script_path = os.path.join(master_dir, package_dir, 'cde-root', 'home', 'proorigami', 'make_cartoon.sh')
    if os.path.exists(script_path):
        with open(script_path, 'r') as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            if line.startswith("PTGRAPH2_OPTIONS="):
                new_lines.append('if [ -z "$PTGRAPH2_OPTIONS" ]; then\n')
                new_lines.append('    PTGRAPH2_OPTIONS="-r35 -t dssp -k purple -l crossing:black,red,green,navy,blue -b sequential -j -e auto -f auto -o gradient -p ddomain"\n')
                new_lines.append('fi\n')
            else:
                new_lines.append(line)
        with open(script_path, 'w') as f:
            f.writelines(new_lines)
        print(f"Modified {script_path} to use environment PTGRAPH2_OPTIONS")
    else:
        print(f"Warning: {script_path} not found. Cannot modify PTGRAPH2_OPTIONS")


def main(pdb_files, num_threads, proorigami_dir, ptgraph_opts, force):
    global master_dir
    master_dir = proorigami_dir
    cdeexec_bin_path = os.path.join(master_dir, package_dir, 'cde-exec')

    # Ensure the package is downloaded/extracted.
    package_tar_path = os.path.join(master_dir, package_tar)
    package_extract_path = os.path.join(master_dir, package_dir)
    if not os.path.exists(package_tar_path):
        download_package(package_url, package_tar_path)
    if not os.path.exists(package_extract_path):
        extract_package(package_tar_path, package_extract_path)
        # Modify make_cartoon.sh.cde right after extraction
        modify_make_cartoon_script(master_dir)
    else:
        # If the package is already extracted, still try to update the make_cartoon.sh.cde file.
        modify_make_cartoon_script(master_dir)
        
    # Ensure the cde-exec binary is executable.
    subprocess.run(['chmod', '+x', cdeexec_bin_path])

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = {executor.submit(execute_process, pdb_file, ptgraph_opts, force): pdb_file for pdb_file in pdb_files}
        for future in as_completed(futures):
            pdb_file = futures[future]
            try:
                future.result()
            except Exception as e:
                print(f"Error processing {pdb_file}: {e}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Process PDB files using proorigami. This script automatically downloads the proorigami script, executes the script in the desired directory, and adds multithreading for multiple PDB files. The proorigami CDE is stored in the $HOME directory.',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-j', '--num_threads', type=int, default=1, help='Number of threads to use (default: 1)')
    parser.add_argument('-d', '--proorigami_dir', type=str, default=os.path.expanduser("~"), help='Directory to store proorigami, this is the workdirectory for the script (default: $HOME)')
    parser.add_argument('--version', action='version', version=f'%(prog)s {version}', help='Show program version and exit')
    parser.add_argument('--ptgraph-options', dest='ptgraph_opts', nargs='*',  default=["-r35", "-t", "dssp", "-k", "purple", "-l", "crossing:black,red,green,navy,blue",
                                 "-b", "sequential", "-j", "-e", "auto", "-f", "auto", "-o", "gradient", "-p", "ddomain"],
                        help="Extra options to pass to ptgraph2.py via environment variable PTGRAPH2_OPTIONS "
                             """
  -3 include 3_10 helices in diagram
  -5 include pi helices in diagram
  -a use Dunnart automatic graph layout
  -c use HELIX and SHEET cards from PDB file
  -d use GraphViz dot instead of Dunnart SVG
  -n use GraphViz neato instead of Dunnart SVG
  -p domain_prog : use domain_prog to parse domains
       supported is 'none' or 'ddomain' (default)
       or 'cath:cdf_file_name'
  -h graph hydrogen bonds with GraphViz
  -b SSE labelling scheme: 'none', 'sequential', 'separate' (default)
  -t struct_prog : use struct_prog define secondary structure
       supported is 'stride' or 'dssp' (default)
  -r compute angles internally, not with external TableauCreator
  -m write MATLAB M-files to plot strand axes
  -s write PyMOL .pml command file to show SSE definitions
  -v print verbose debugging messages to stderr
  -i use distance matrix information instead of
     heuristic/aesthetic algorithm for helix placement
  -j only valid when not using -i. Don't align helices
     on strand axes if they would push sheets apart
  -k <color> cluster helices, shading them all <color>
  -e <color_list>|auto shade nearby helix clusters the same color
  -x draw connector arrowheads
  -f <color_list>|auto shade each sheet a different color
  -g <separation> set the strand and minimum object separation
  -l connector color scheme: 'all[:<color>]' (default), 'chain[:<color_list>]', 'domain[:<intra_color>,<inter_color>'], crossing:<color_list>
  -o SSE color scheme: 'none' (default), 'simple:sheet=<sheet_colors>.helixcluster=<helixcluster_colors>.alpha=<helix_alpha_colors>.pi=<helix_pi_colors>.310=<helix_310_colors>.terminus=<terminus_colors>', 'gradient', 'sheet', 'fold'
  -u multidomain cartoon: place all domains in the one
     SVG file instead of one per file
  -w interdomain connectors: when using multidomain
     cartoons, draw connectors between domains
     (only in conjunction with -u)
  -q label start and end of helices and strands with
     first and last PDB residue id in that SSE.
  -y use uniform scaling to try to avoid overlaps.
     Ugly and often does not work anyway, use only as
     last resort
  -z print version information and exit

  DEFAULT: '-r35 -t dssp -k purple -l crossing:black,red,green,navy,blue -b sequential -j -e auto -f auto -o gradient -p ddomain'
                               """)
    parser.add_argument('--force', action='store_true', help="Force overwrite of existing output files")
    parser.add_argument('pdb_files', nargs='+', help='List of PDB files to process')

    args = parser.parse_args()

    main(args.pdb_files, args.num_threads, args.proorigami_dir, args.ptgraph_opts, args.force)
