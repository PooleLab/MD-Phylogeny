#!/usr/bin/python3
#pk
import os, random, string
import re
import sys
import subprocess
import requests
from concurrent.futures import ThreadPoolExecutor, TimeoutError, as_completed
import shutil  # Add this import at the top
import tempfile
import random

version="0.0.16" # last updated on 2024-11-13 13:19:56

package_url = 'https://stivalaa.github.io/AcademicWebsite/software/proorigami-cde-package.tar.gz'
master_dir = os.path.expanduser("~")
package_tar = 'proorigami.tar.gz'
package_dir = 'proorigami'
cdeexec_bin_dir = os.path.join(master_dir,package_dir,'cde-exec')
make_cartoon_bin_dir = os.path.join(master_dir,package_dir,'cde-root','home','proorigami','make_cartoon.sh.cde')
make_cartoon_dir = os.path.join(master_dir,package_dir,'cde-root','home','proorigami')

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

timeout=10

def download_package(url, dest):
    print(f"Downloading {url}...")
    response = requests.get(url)
    with open(dest, 'wb') as f:
        f.write(response.content)
    print(f"Downloaded {dest}.")

def extract_package(tar_path, extract_to):
    print(f"Extracting {tar_path}...")
    os.makedirs(os.path.join(master_dir,package_dir),exist_ok=True)
    subprocess.run(['tar', 'zxf', tar_path, '-C', os.path.join(master_dir,package_dir), '--strip-components=1'])
    print(f"Extracted to {extract_to}.")

def check_cde_exec():
    return os.path.isfile('./cde-exec') and os.access('./cde-exec', os.X_OK)

def random_string(length=8):
    return ''.join(random.choices(string.ascii_letters + string.digits, k=length))

def process_pdb(pdb_file):
    print(f"Processing {pdb_file}")
    try:
        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(execute_process, pdb_file)
            future.result(timeout=timeout)
    except TimeoutError:
        print(f"Processing {pdb_file} timed out. Retrying...")
        try:
            with ThreadPoolExecutor(max_workers=1) as executor:
                future = executor.submit(execute_process, pdb_file)
                future.result(timeout=timeout)
        except TimeoutError:
            print(f"Second attempt timed out for {pdb_file}.")
        except Exception as e:
            print(f"Second attempt failed for {pdb_file}: {e}")
    except Exception as e:
        print(f"First attempt failed for {pdb_file}: {e}")
        print(f"Retrying {pdb_file}...")
        try:
            with ThreadPoolExecutor(max_workers=1) as executor:
                future = executor.submit(execute_process, pdb_file)
                future.result(timeout=timeout)
        except TimeoutError:
            print(f"Second attempt timed out for {pdb_file}.")
        except Exception as e:
            print(f"Second attempt failed for {pdb_file}: {e}")

def execute_process(pdb_file):
    bn = os.path.basename(pdb_file)
    dirn = os.path.dirname(pdb_file)
    prefix = random_string()    
    bn_base = re.sub(r"[.][^.]+$","",bn)

    if os.path.exists(os.path.join(dirn,f'{bn_base}.png')):
        print(f"skipping {pdb_file}")
        return

    print(f"processing {pdb_file}")

    seed = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10)) # just to be extra save

    # Create a temporary directory for this process
    with tempfile.TemporaryDirectory() as tmp_dir:
        shutil.copytree(os.path.join(master_dir,package_dir),os.path.join(tmp_dir,f"proorigami_{seed}") )
        cwd = os.path.join(tmp_dir,f"proorigami_{seed}",'cde-root','home','proorigami')
        
        linkp = os.path.join(tmp_dir,f"proorigami_{seed}",'cde-root','home','proorigami', f"{prefix}_{bn}")
        os.symlink(os.path.realpath(pdb_file), linkp)

        # Run the command in the temporary directory
        p = subprocess.Popen(
            ["bash", "make_cartoon.sh.cde", os.path.basename(linkp)],
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = p.communicate()
        stdout = stdout.decode()
        stderr = stderr.decode()

        # Check for results and move files back to the original directory
        res = re.search(r"Bitmap saved as: (.*)[.]png", stdout + stderr)
        res_base = res.group(1) if res else ""

        if res_base:
            shutil.move(os.path.join(cwd, f'{res_base}.png'), os.path.join(dirn, f'{bn_base}.png'))
            shutil.move(os.path.join(cwd, f'{res_base}.svg'), os.path.join(dirn, f'{bn_base}.svg'))
        else:
            raise Exception(f'{stdout}\n{stderr}\n{pdb_file} fail')

    print(f"Finished processing {pdb_file}")
def main(pdb_files, num_threads):
    if not check_cde_exec():
        if not os.path.exists(os.path.join(master_dir,package_tar)):
            download_package(package_url, os.path.join(master_dir,package_tar))
        if not os.path.exists(os.path.join(master_dir,package_dir)):
            extract_package(os.path.join(master_dir,package_tar), os.path.join(master_dir,package_dir))
        subprocess.run(['chmod', '+x', cdeexec_bin_dir ])

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = {executor.submit(process_pdb, pdb_file): pdb_file for pdb_file in pdb_files}
        for future in as_completed(futures):
            pdb_file = futures[future]
            try:
                future.result()
            except Exception as e:
                print(f"Error processing {pdb_file}: {e}")

        #executor.map(process_pdb, pdb_files)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process PDB files using proorigami.')
    parser.add_argument('num_threads', type=int, help='Number of threads to use')
    parser.add_argument('pdb_files', nargs='+', help='List of PDB files to process')

    args = parser.parse_args()
    main(args.pdb_files, args.num_threads)
