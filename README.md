# Structural-Phylogeny-with-MD-Bootstrap
The Structural Phylogeny with MD-Bootstrap method constructs structural phylogenies and provides statistical confidence through molecular dynamics (MD) simulations.

## What does the method do, and how? 

## System Requirements 
Operating system: Ideally, UNIX or Mac. If Windows use Windows Subsystem for Linux (WSL/WSL2), be aware of possible difficulties. 

MD Engine (best to use on an HPC): GROMACS (we’ve used GROMACS version GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid) and CHARMM are common choices for protein simulations, each compatible with different force fields (e.g., CHARMM36 for GROMACS what we have used).


## Installations 

Python for script execution, with NumPy, dendropy, Bio, natsort 
#### Structural alignment
Gesamt (from https://www.ccp4.ac.uk/html/gesamt.html)
#### Phylogenetic Tree Summarization and Annotation: 
sumtrees (https://jeetsukumaran.github.io/DendroPy/programs/sumtrees.html)
#### Visualisation: 
VMD (or Pymol) for structures and simulations 
FigTree or iTOL for phylogenies


## Method

### 1. Create structural phylogeny

- Folder and file setup: Create a main directory containing necessary scripts and a subdirectory specifically for PDB files.
  - PDB Naming Convention: Each PDB file should be a single-chain structure with the format PDBID_chainID.pdb or PDBIDchainID.pdb (e.g., 1abc_A.pdb, 1abcA.pdb). This is critical for script compatibility.
- Create phylogeny: Run the GesamtTree.py script on the specified folder to generate a structural phylogeny based on the PDB files.
  - ``` python GesamtTree.py <FolderName> ```


### 2. MD simulation set-up

In general, this requires: 
1. PDB File Preparation
2. Setup the size and shape of the simulation box
3. Energy minimization
4. Solvation
5. Ionisation
6. Heating of the system

For detailed information about each step, see ((https://link.springer.com/protocol/10.1007/978-1-4939-9869-2_17#Sec16))

These steps can be done:
- manually by using CHARMM GUI Website (https://www.charmm-gui.org/)
- using scripts

#### Option 1: Using CHARMM GUI Website (https://www.charmm-gui.org/) (allows to visually check outputs and intermediate steps visually; very user-friendly interface and handling)

1. Account Access: CHARMM-GUI requires an active account for input generation.
2. Setup Steps:
  - After logging in, navigate to Input Generator > Solution Builder. Solution Builder generates input files for molecular dynamics simulations, allowing you to either solvate your molecule or create a standalone water box for other uses (Alternatively, select a different builder based on the protein type).
  - Load Protein:
    - Option 1: Upload a specific PDB file.
    - Option 2: Enter the PDB ID directly.
  - Model/Chain Selection Option: Select the desired protein chain, then proceed to configure the settings.
3. Configuration Options: 
  - PDB Manipulation Options: Customize manipulation options for the structure based on experimental needs, such as removing certain residues or modifying the structure, if required. Click Next Step to apply the settings. After the calculation, use the CHARMM-GUI viewer to confirm that the desired modifications to the protein structure have been applied accurately before proceeding to the next configuration.
  - Water Box Configuration and ions: Configure the water box size and type to solvate the protein. Add ions to achieve wanted physiological conditions (Note: we’ve used distance placing method and chosen 10 A edge distance and added NaCl ions). Click Next Step to apply the settings. After the calculation, use the CHARMM-GUI viewer to confirm that the desired modifications to the protein structure have been applied accurately before proceeding to the next configuration.
  - Boundary Conditions: Set the periodic boundary conditions (PBC) automatically or adjust manually according to simulation requirements. Click Next Step to apply the settings. After the calculation, use the CHARMM-GUI viewer to confirm that the desired modifications to the protein structure have been applied accurately before proceeding to the next configuration.
  - Force Field Selection: Choose desired force field. CHARMM36m is popular and recommended for protein systems. Ensure that selected force fields and solvent models are compatible. 
  - Equilibration and Simulation Parameters:
    - Select GROMACS as the output format.
    - Set NVT Ensemble for equilibration and NPT Ensemble for dynamics. Specify the temperature of choice for simulations.
    - Click Next Step to apply the settings. After the calculation, use the CHARMM-GUI viewer to confirm that the desired modifications to the protein structure have been applied accurately before proceeding to the next configuration.
4. Download and Organize Input Files:
  - Download the generated files download.tgz and save them in a directory named after the protein’s PDB ID.
  - Upload to HPC: Transfer the folder to the HPC cluster and extract the files

#### Option 2: Using Scripts
> [!CAUTION]
> Scripts only work if pdb-file only contains protein with no engeniered AA, waters, ions,...


1. PDB File Preparation: Manually clean the PDB file to ensure compatibility with the MD setup:
  - Remove non-standard residues (engineered residues, water, ions, ligands, RNA/DNA,... as it can create errors when running gmx and gmx grompp).
  - Only keep the protein of interest.
2. Make force field available: CHARMM27 comes with GROMACS as a default option. But CHARMM27 can also be used. To make the newer version CHARMM36m available, it needs to be downloaded from http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs and once unzipped kept in the resulting force field folder of gromacs (it depends on where gromacs is installed; standard installation location: /usr/share/gromacs/top). Alternatively, when using a HPC it can be placed into the working directory. 
It is recommended to try the simulation setup steps to see which number the force field CHARMM36m has, it has to be adjusted in the scripts used
3. Run Simulation Setup:
  - HPC Version: Use sim_setup_slurm.py to configure system requirements, Slurm header, and other parameters based on the computational resources.
  - Local Version: If running on a local machine, execute sim_setup.py (the non-Slurm version).





### 3. Simulation Run 
Directory Organization: Create a unique directory for each protein and place the charmm-gui.tgz file inside (as mentioned in 3.2 point 4)
Upload to HPC and extract files: Transfer the folder to the HPC cluster and extract the files. Unzip the CHARMM-GUI setup files in each folder (as mentioned in 3.2 point 4)
The folder structure will be <proteinfoldername>/charmm*/gromacs
Equilibration:
Equilibrate the system to stabilize temperature, pressure, and density.
Run s1_make_equilibration_scripts.py to create Slurm scripts in each GROMACS folder for equilibration. (Note: be aware of changing the directory in the script (line 5) and all SBATCH settings to liking and requirements. Also, change the names of files required for simulation. (line 33 (minimization.mdp, input.gro, index.ndx, topol.top) and 38 (equilibration.mdp) )
Initial Testing: Test with one protein to identify potential issues before full batch submission. Change the project name and email address in s2. The email sent when the equilibration started and ended will tell how long the equilibration took, maybe adjust the time as needed.
Execute s2_make_equilibration_execute_script.py to submit all generated equilibration scripts. (Note: Adjust the path in line 3 to your file locations)
Check the trajectories for these parameters for any abnormalities 
Production Run:
Perform the production simulation to under stable, equilibrated conditions.
(If you change dt to 0.002 ps, and you want your simulation to last for 100 ns, then nsteps is 100,000 ps / 0.002 ps = 50,000,000)
The folder structure will be <proteinfoldername>/charmm*/gromacs
Use s3_make_production_scripts.py to create production scripts for each protein.
Submit all production runs with s4_make_production_execute_script.py.
Adjustments: Set the Slurm array size for efficient job management:
#SBATCH --array=0-6%1  
Adjust '6' based on protein size or simulation length
For a protein from roughly  300 AA 6 has proven to be enough
Initial Testing: Test with one protein to identify potential issues before full batch submission. Change the project name and email address in s3. Check Log Files associated with each job ID to estimate run completion time and adjust time/array accordingly.




### 4. Post-Simulation Processing








### 5. Bootstrapping and Structural Phylogeny Analysis
