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


# Method

## 1. Create structural phylogeny

- Folder and file setup: Create a main directory containing necessary scripts and a subdirectory specifically for PDB files.
> [!IMPORTANT]
> **PDB Naming Convention:** Each PDB file should be a single-chain structure with the format ```PDBID_chainID.pdb``` or ```PDBIDchainID.pdb``` (e.g., 1abc_A.pdb, 1abcA.pdb). This is critical for script compatibility.
- Create phylogeny: Run the GesamtTree.py script on the specified folder to generate a structural phylogeny based on the PDB files
  ```
  python GesamtTree.py <FolderName>
  ```

Then, molecular dynamics is used to generate alternative conformations, which can then be sampled randomly to build alternative trees for bootstrap support in phylogenetic analyses

## 2. MD simulation set-up
MD setup prepares a molecular system for simulation, ensuring it is physically and realistic configured.
In general, this requires: 
1. **PDB File Preparation:** Clean the structure (remove non-standard residues, add missing atoms)
2. **Setup of the simulation box:** Define a box around the molecule (size, shape)
3. **Energy minimization:** Resolve steric clashes and optimize geometry.
4. **Solvation:** Add water molecules to mimic a realistic environment.
5. **Ionisation:** Add ions (e.g., Na+, Cl⁻) for charge neutrality or the natural environment of the molecule.
6. **Heating of the system:** Gradually raise system temperature to simulation conditions.

For detailed information about each step, see ((https://link.springer.com/protocol/10.1007/978-1-4939-9869-2_17#Sec16))

These steps can be done:
- manually by using CHARMM GUI Website (https://www.charmm-gui.org/)
- using command line/scripts

### Option 1: Using CHARMM GUI Website (https://www.charmm-gui.org/) 
(allows to visually check outputs and intermediate steps visually; very user-friendly interface and handling)

> [!IMPORTANT]
> CHARMM-GUI requires an active account for input generation. You need to set this up before you can use the CHARMM GUI Solution Builder.
> 
1. Setup Steps:
  - After logging in, navigate to Input Generator > Solution Builder. Solution Builder generates input files for molecular dynamics simulations, allowing you to either solvate your molecule or create a standalone water box for other uses (Alternatively, select a different builder based on the protein type).
  - Load Protein:
    - Option 1: Upload a specific PDB file.
    - Option 2: Enter the PDB ID directly.
  - Model/Chain Selection Option: Select the desired protein chain, then proceed to configure the settings.
3. Configuration Options: 
  - PDB Manipulation Options: Customize manipulation options for the structure based on experimental needs, such as removing certain residues or modifying the structure, if required. Click Next Step to apply the settings. After the calculation, use the CHARMM-GUI viewer to confirm that the desired modifications to the protein structure have been applied accurately before proceeding to the next configuration.
  - Water Box Configuration and ions: Configure the water box size and type to solvate the protein. Add ions to achieve wanted physiological conditions (Note: we’ve used the distance placing method and chosen 10 A edge distance and added NaCl ions). Click Next Step to apply the settings. After the calculation, use the CHARMM-GUI viewer to confirm that the desired modifications to the protein structure have been applied accurately before proceeding to the next configuration.
  - Boundary Conditions: Set the periodic boundary conditions (PBC) automatically or adjust manually according to simulation requirements. Click Next Step to apply the settings. After the calculation, use the CHARMM-GUI viewer to confirm that the desired modifications to the protein structure have been applied accurately before proceeding to the next configuration.
  - Force Field Selection: Choose the desired force field. CHARMM36m is popular and recommended for protein systems. Ensure that selected force fields and solvent models are compatible. 
  - Equilibration and Simulation Parameters:
    - Select GROMACS as the output format.
    - Set NVT Ensemble for equilibration and NPT Ensemble for dynamics. Specify the temperature of choice for simulations.
    - Click Next Step to apply the settings. After the calculation, use the CHARMM-GUI viewer to confirm that the desired modifications to the protein structure have been applied accurately before proceeding to the next configuration.
4. Download and Organize Input Files:
  - Download the generated files download.tgz and save them in a directory named after the protein’s PDB ID.
  - Upload to HPC: Transfer the folder to the HPC cluster and extract the files

### Option 2: Using scripts
> [!CAUTION]
> Scripts only work if pdb-file only contains protein with no engeniered AA, waters, ions,...


1. PDB File Preparation: Manually clean the PDB file to ensure compatibility with the MD setup:
  - Remove non-standard residues (engineered residues, water, ions, ligands, RNA/DNA,... as it can create errors when running gmx and gmx grompp).
  - Only keep the protein of interest.
2. Make force field available: CHARMM27 comes with GROMACS as a default option. But CHARMM27 can also be used. To make the newer version CHARMM36m available, it needs to be downloaded from http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs and once unzipped kept in the resulting force field folder of gromacs (it depends on where gromacs is installed; standard installation location: /usr/share/gromacs/top). Alternatively, an HPC can be placed in the working directory. 
It is recommended to try the simulation setup steps to see which number the force field CHARMM36m has; it has to be adjusted in the scripts used
3. Run Simulation Setup: Use ```MD/sim_setup_slurm.py``` to configure system requirements, Slurm header, and other parameters based on the computational resources.

### Option 3: Using command line
Run the following commands in the command line to set a protein up. Check the output of every step to make sure everything is set up correctly
```
 module load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid
```
THIS RUNS THE FILE CONVERSION FROM PDB TO GRO
```
gmx pdb2gmx -f <protein>.pdb -o <protein>.gro -ignh
```
Choose the force field ( for example, CHARMM36) and TIP3 water model

ADD SIMULATION BOX
```
gmx editconf -f <protein>.gro -o <protein>_box.gro -bt cubic -d 1.2 -c
```
ENERGY MINIMISATION
```
gmx grompp -f min_sd.mdp -c <protein>_box.gro -p topol.top -o <protein>_minsd.tpr
gmx mdrun -deffnm <protein>_minsd
```
ADD WATER TO THE BOX
```
gmx solvate -cp <protein>_minsd.gro -cs spc216.gro -p topol.top -o <protein>_solv.gro
```
ENERGY MINIMISATION
```
gmx grompp -f min_sd.mdp -c <protein>_solv.gro -p topol.top -o <protein>_solv_minsd.tpr -maxwarn 1
gmx mdrun -deffnm <protein>_solv_minsd 
```
ADD IONS
```
gmx grompp -f ions.mdp -c <protein>_solv_minsd.gro -p topol.top -o <protein>_ions.tpr -maxwarn 1
gmx genion -s <protein>_ions.tpr -p topol.top -o <protein>_ions.gro -neutral
```
Choose Iontype 

MAKE INDEX FILE
```
gmx make_ndx -f <protein>_ions.gro -o index.ndx
```
Type 'q' and enter

ENERGY MINIMISATION
```
gmx grompp -f min_sd.mdp -c <protein>_ions.gro -p topol.top -o <protein>_ions_minsd.tpr -maxwarn 1
gmx mdrun -deffnm <protein>_ions_minsd
```
INITIALISATION AND HEATING
```
gmx grompp -f charmm-inputs/charmm_nvt_heat.mdp -c <protein>_ions_minsd.gro -r <protein>_ions_minsd.gro -p topol.top -o <protein>_nvt_heat.tpr
gmx mdrun -deffnm <protein>_nvt_heat
```



## 3. MD Simulation Run 
1. Directory Organization: Create a unique directory for each protein and place the ```charmm-gui.tgz``` file inside (as mentioned in 3.2 point 4)
2. Upload to HPC and extract files: Transfer the folder to the HPC cluster and extract the files. Unzip the CHARMM-GUI setup files in each folder (as mentioned in 3.2 point 4)
The folder structure will be ```<proteinfoldername>/charmm*/gromacs```
3. Equilibration:
  - Equilibrate the system to stabilize temperature, pressure, and density.
  - Run s1_make_equilibration_scripts.py to create Slurm scripts in each GROMACS folder for equilibration. (Note: be aware of changing the directory in the script (line 5) and all SBATCH settings to liking and requirements. Also, change the names of files required for simulation. (line 33 (minimization.mdp, input.gro, index.ndx, topol.top) and 38 (equilibration.mdp) )
  - Initial Testing: Test with one protein to identify potential issues before full batch submission. Change the project name and email address in s2. The email sent when the equilibration started and ended will tell how long the equilibration took, maybe adjust the time as needed.
  - Execute s2_make_equilibration_execute_script.py to submit all generated equilibration scripts. (Note: Adjust the path in line 3 to your file locations)
  - Check the trajectories for these parameters for any abnormalities 
4. Production Run:
  - Perform the production simulation to under stable, equilibrated conditions.
(If you change dt to 0.002 ps, and you want your simulation to last for 100 ns, then nsteps is 100,000 ps / 0.002 ps = 50,000,000)
  - The folder structure will be <proteinfoldername>/charmm*/gromacs
  - Use s3_make_production_scripts.py to create production scripts for each protein.
  - Submit all production runs with s4_make_production_execute_script.py.
  - Adjustments: Set the Slurm array size for efficient job management:
#SBATCH --array=0-6%1  
Adjust '6' based on protein size or simulation length
For a protein from roughly  300 AA 6 has proven to be enough
  - Initial Testing: Test with one protein to identify potential issues before full batch submission. Change the project name and email address in s3. Check Log Files associated with each job ID to estimate run completion time and adjust time/array accordingly.




## 4. Post-Simulation Processing

1. Trajectory Concatenation: Concatenate all trajectory segments into a single trajectory file for continuity. Use GROMACS to concatenate trajectory files for each protein, combining individual parts into a complete trajectory:
 ```
 gmx trjcat -f <name>.part*.xtc -o <name>.all.xtc
 ``` 
3. Protein Extraction and PBC Correction: Transform the concatenated trajectory to remove artefacts caused by periodic boundaries. Centre the protein in the simulation box, ensuring molecules are visually and structurally continuous. Align the protein across frames by removing global rotation and translation, making it easier to analyze internal motions. Use GROMACS to extract the protein coordinates and apply periodic boundary corrections:
```
gmx trjconv -f <name>.all.xtc -ur compact -n index.ndx -o <name>.prot.xtc -pbc mol -center -s <name>.tpr
```
and fit the protein for final analysis: 
``` 
gmx trjconv -f <name>.prot.xtc -n index.ndx -o <name>.proc.xtc -s <name>.tpr -fit rot+trans
``` 
5. Extract final structure: Extract a single representative structure from the trajectory at a specified time point (ideally the starting point) to view the simulation. Save gro file to view the simulation in VMD:
```
gmx trjconv -f <name>.prot.xtc -b 0 -e 0 -n index.ndx -o <name>.proc.gro -s <name>.tpr
``` 






## 5. Bootstrapping and Structural Phylogeny Analysis

This section describes the process for generating structural phylogenies with bootstrapping. Bootstrapping helps assess the stability of phylogenetic relationships by generating multiple replicates, each based on random samples of simulation frames. The following steps require the completed MD simulations with trajectory (.xtc) and structure (.gro) files.
1. Organize Input Files:
Ensure all .xtc and .gro files for each protein are collected in a folder named Trajectories. Each trajectory file should correspond to a unique protein.
2. Prepare Bootstrapping Replicates with BS_master_sort.py:
  - Setup: Place the BS_master_sort.py script outside the Trajectories folder to organize trials and generate the required number of bootstrap samples. This script will select random frames from each trajectory, creating directories with the sampled structures for each trial.
  - Parameter Adjustments:
    - Number of Trials: Define the total number of bootstrap replicates by setting the N_trials variable in BS_master_sort.py. For example:
    ```N_trials = 200  # Set to the desired number of bootstrap trials```
    ```dir_trial = np.arange(0, 200, 1)  # Creates folders for each trial```
    This configuration will generate 200 bootstrap trials, each with randomly selected frames.
    - Frame Selection: Adjust the number of frames per trial based on the length of your MD simulation: ```frame_sel = np.random.randint(0, 10000, N_trials)  # Modify '10000' to match total frames of your simulation```
To determine the correct frame count, load the .gro and .xtc files into a viewer like VMD. Once fully loaded, VMD will display the total frame count for the trajectory.
  - Generate Bootstrapping Directories: Running BS_master_sort.py will create a series of directories named trial_0, trial_1, etc., each containing a unique set of frames selected randomly from the original trajectory. 
3. Run Structural Analysis on Trials:
  - With the trials generated, use GesamtTree.py to analyze each trial folder within the Trajectories/Trials directory. This script will generate phylogenetic trees for each bootstrap replicate.
  - Execution Command:
This command
```
python3 GesamtTree.py Trajectories/Trial/trial_*
``` runs GesamtTree.py on each trial directory, producing a phylogenetic tree output for every bootstrap sample.
4. Generate a Consensus Phylogenetic Tree and Map Bootstrap Values:
  - Once phylogenetic trees have been generated for all bootstrap trials, you can either:
    - Create a Majority-Rule Consensus Tree from the bootstrap replicates, which combines all individual trees to identify the most consistent relationships.
    - Map Bootstrap Values onto a Reference Tree to reflect the stability of nodes in the reference structure.
  - Options:
    - Majority-Rule Consensus Tree: This option compiles the bootstrap trial trees into a single consensus tree, with nodes representing relationships that appear in a specified percentage of the bootstrap samples. For example, setting a 60% threshold (```-f 0.6```) will include only those relationships that appear in at least 60% of the trees.
```
sumtrees -s consensus -o consensus60 -f 0.6 -p -d0 Trajectories/Trials/trial_*/*.nex
```
This produces a consensus tree (consensus60) that highlights the most frequent relationships across bootstrap trials.
    - Mapping Bootstrap Values onto a Reference Tree: If you have a reference tree from Section 3.1, you can add bootstrap support values directly onto it, rather than generating a separate consensus tree. This approach maps the percentage of bootstrap trees supporting each node directly onto the nodes of the reference tree, indicating the stability of each branch.
```
sumtrees -d0 -p -o OutputTree -t Referencetree/<ReferenceTree> Trajectories/Trials/trial_*/*Formatted.nex
```

  - Here:
    - OutputTree will contain the reference tree structure with bootstrap values annotated on each node.
    - Interpretation: Higher bootstrap values on a node in the reference tree suggest strong support for that relationship across bootstrap replicates, while lower values indicate less stability.

