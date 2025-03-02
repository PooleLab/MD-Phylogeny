  # MD-Phylogeny
MD-Phylogeny constructs structural phylogenies and provides statistical confidence through molecular dynamics (MD) simulations.

## Table of content
- [MD-Phylogeny](#md-phylogeny)
  * [System Requirements](#system-requirements)
  * [Installations](#installations)
      - [Visualisation of structures and simulations](#visualisation-of-structures-and-simulations)
      - [Structural superposition and calculation of Qscore](#structural-superposition-and-calculation-of-qscore)
      - [Visualisation of Phylogenies](#visualisation-of-phylogenies)
      - [Summarization and annotation of phylogenies](#summarization-and-annotation-of-phylogenies)
      - [Additional requirements](#additional-requirements)
- [Method](#method)
  * [1. Dataset selection](#1-dataset-selection)
    + [Identify a Suitable Dataset](#identify-a-suitable-dataset)
    + [Evaluate Dataset Suitability](#evaluate-dataset-suitability)
    + [Additional notes](#additional-notes)
    + [Structural comparison and topology maps](#structural-comparison-and-topology-maps)
  * [2. Create structural phylogeny](#2-create-structural-phylogeny)
  * [3. MD simulation set-up](#3-md-simulation-set-up)
  * [4. Post-Simulation Processing](#4-post-simulation-processing)
  * [5. Bootstrapping and Structural Phylogeny Analysis](#5-bootstrapping-and-structural-phylogeny-analysis)
- [Tutorial](#tutorial)




## System Requirements 
Operating system: Ideally, UNIX or Mac (not tested on Mac, yet). If you are using Windows, you can use the Windows Subsystem for Linux (WSL/WSL2), but be aware of possible difficulties. 
For most preparatory and analytical steps it is ok to use laptop or desktop. But for molecular dynamics (MD) simulations an HPC is highly recommended. 

MD Engine (best to use on an HPC, as it gets too computationally intense for a Laptop or Desktop PC): Be aware that our scripts are developed for GROMACS version 2020 and will require significant adaptation for use with other MD engines

GROMACS (we have used GROMACS version GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid) and CHARMM force fields (we have used CHARMM36 for GROMACS). (https://doi.org/10.5281/zenodo.3562512)



## Installations 

#### Visualisation of structures and simulations
Visual Molecular Dynamics (VMD) (https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) or Pymol (https://www.pymol.org/) to visualise structures and VMD to visualise simulation trajectories. 
#### Structural superposition and calculation of Qscore
Gesamt (from https://www.ccp4.ac.uk/html/gesamt.html)
#### Visualisation of Phylogenies 
FigTree (http://tree.bio.ed.ac.uk/software/figtree/) or iTOL (https://itol.embl.de/)
#### Summarization and annotation of phylogenies 
sumtrees (https://jeetsukumaran.github.io/DendroPy/programs/sumtrees.html)
#### Additional requirements 
Python3 with the following packages: 

NumPy (https:/2/numpy.org), dendropy (https://pypi.org/project/DendroPy), Bio Python (https://biopython.org/docs/1.75/api/Bio.html), natsort (https://pypi.org/project/natsort)



# Method
## 1. Dataset selection
### Identify a Suitable Dataset
The dataset should consist of **structures that are likely to share homology**. Consider using tools and databases such as:

Foldseek (https://search.foldseek.com/search , https://github.com/steineggerlab/foldseek), Structome (https://biosig.lab.uq.edu.au/structome_q/), SCOP (https://www.ebi.ac.uk/pdbe/scop/), CATH (https://www.cathdb.info/), InterPro (https://www.ebi.ac.uk/interpro/)

### Evaluate Dataset Suitability
Ensure that the dataset lies within the "twilight zone" of sequence similarity, where traditional sequence-based phylogenetic methods may be unreliable. Therefore e.g. perform an all-against-all pairwise sequence comparisons. Visualize the sequence similarity network using tools such as: CLANS and MMseqs2

:x: Well-connected network → **Use conventional sequence-based phylogenetics**.

:heavy_check_mark: Poorly connected network → **Suitable for MD-phylogeny**.

### Additional notes
Ensure the dataset meets the following criteria:
- Structures should consist of single-chain files
- Identify structures that contain multiple copies of the same motif or fold. This influences the analysis
- Avoid structures that are too large (e.g., ribosomes, entire nucleosomes) as they are impractical for MD simulations. 
- Ensure length homogeneity across structures

### Structural comparison and topology maps
Once the dataset is prepared, perform pairwise structural comparisons to generate a distance matrix for phylogenetic analysis. Careful verification of structural alignments is essential to avoid misinterpretations in phylogenetic reconstruction. Ensure that:
- Only relevant domains are considered to prevent noise in the distance matrix.
- Homologous structures are correctly aligned and superimposed.

While structures may appear similar in 3D, differences in topology may indicate they are unrelated. Checking topology ensures accurate structural comparisons for subsequent phylogenetic analysis.
We've developed a script for a user-friendly application of Pro-origami ```easyproorigami.py```

 ```
  python3 easyproorigami.py <pdb-structure>
  ```


## 2. Create structural phylogeny

- Folder and file setup: Create a main directory containing necessary scripts and a subdirectory specifically for PDB files.
> [!IMPORTANT]
> **PDB Naming Convention:** Each PDB file should be a single-chain structure with the format ```PDBID_chainID.pdb``` or ```PDBIDchainID.pdb``` (e.g., 1abc_A.pdb, 1abcA.pdb). 
- Create phylogeny: Run the GesamtTree.py script on the specified folder that contains all ```.pdb``` files of the structures to generate a structural phylogeny
  ```
  python GesamtTree.py <FolderName>
  ```

## 3. MD simulation set-up
MD setup prepares a molecular system for simulation, ensuring it is physically and realistically configured.
In general, this requires: 
1. **PDB File Preparation:** Clean the structure (remove non-standard residues, add missing atoms)
2. **Setup of the simulation box:** Define a box around the molecule (size, shape)
3. **Energy minimization:** Resolve steric clashes and optimize geometry.
4. **Solvation:** Add water molecules to mimic a realistic environment.
5. **Ionisation:** Add ions (e.g., Na+, Cl⁻) for charge neutrality or the natural environment of the molecule.
6. **Heating of the system:** Gradually raise system temperature to simulation conditions.
7. **Equilibration:** to stabilize temperature, pressure, and density.

For detailed information about each step, see ((https://link.springer.com/protocol/10.1007/978-1-4939-9869-2_17#Sec16))

These steps can be done:
- manually by using CHARMM GUI Website (https://www.charmm-gui.org/)
- using command line/scripts

### Option 1: Using CHARMM GUI Website 
(https://www.charmm-gui.org/)

(allows to visually check outputs and intermediate steps visually; very user-friendly interface and handling)

> [!Note]
> CHARMM-GUI requires an active account for input generation. You need to set this up before using the CHARMM GUI Solution Builder.

Follow the steps and adjust to your requirements, then download the files for simulation. You need  ```<protein>.gro```, ```index.ndx``` and ```topol.top``` files to continue with the MD simulation run. 
> [!IMPORTANT]
> Rename files according to the scripts in the next step or change file names in the script according to the obtained files.

### Option 2: Using scripts
> [!Note]
> Scripts only work if the PDB files contain only standard protein residues, with no engineered amino acids, waters, ions, etc.

1. Make the force field available: GROMACS has CHARMM27 as a default option. While this force field is old, it can be used. However, we recommend using the newer CHARMM36 or CHARMM36m versions. To make these available to GROMACS, they need to be downloaded from http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs and, once unzipped, the resulting folder should be kept in the resulting force field folder of gromacs (it depends on where gromacs is installed; the standard installation location is /usr/share/gromacs/top). Alternatively, e.g. on an HPC, the force field folder can be placed in the working directory. 
It is recommended to try the simulation setup steps to see which number your desired force field has when listed by pdb2gmx; it has to be adjusted in the scripts.

3. Run Simulation Setup and use:
   - ```MD/s0_sim_setup_slurm.py``` on an HPC with Slurm scheduler. Adjust the Slurm header and other parameters based on the computational resources. Run this script in the same folder as you ```.pdb``` files. This will create a folder for each protein and places a ```<protein>_setup.sh``` script in each protein folder. Run each of the scripts. 
   - ```MD/s0_sim_setup.py``` to run on a Laptop or Desktop PC. Run this script in the same folder as you ```.pdb``` files. This will create a folder for each protein and places a ```<protein>_setup.sh``` script in each protein folder. Run each of the scripts.
  
>[!Important]
>You then need the ```<protein>_nvt_heat.gro```, ```index.ndx``` and ```topol.top``` files to continue with the MD simulation run. (Rename them according to the scripts in the next step or change file names in the script according to your files.)

### Option 3: Using command line
Run the following commands in the command line to set a protein up. Check every step's output to ensure everything is set up correctly.

THIS RUNS THE FILE CONVERSION FROM PDB TO GRO
```
gmx pdb2gmx -f <protein>.pdb -o <protein>.gro -ignh
```
Choose the force field ( for example, CHARMM36) and TIP3 water model.

ADD SIMULATION BOX
```
gmx editconf -f <protein>.gro -o <protein>_box.gro -bt cubic -d 1.2 -c
```
ENERGY MINIMISATION
```
gmx grompp -f min_sd.mdp -c <protein>_box.gro -p topol.top -o <protein>_minsd.tpr
```
```
gmx mdrun -deffnm <protein>_minsd
```
ADD WATER TO THE BOX
```
gmx solvate -cp <protein>_minsd.gro -cs spc216.gro -p topol.top -o <protein>_solv.gro
```
ENERGY MINIMISATION
```
gmx grompp -f min_sd.mdp -c <protein>_solv.gro -p topol.top -o <protein>_solv_minsd.tpr -maxwarn 1
```
```
gmx mdrun -deffnm <protein>_solv_minsd 
```
ADD NaCl 
```
gmx grompp -f ions.mdp -c <protein>_solv_minsd.gro -p topol.top -o <protein>_ions.tpr -maxwarn 1
```
```
gmx genion -s <protein>_ions.tpr -p topol.top -o <protein>_ions.gro -neutral
```
Choose 13 (for SOL) and press enter.

MAKE INDEX FILE
```
gmx make_ndx -f <protein>_ions.gro -o index.ndx
```
Type 'q' and enter

ENERGY MINIMISATION
```
gmx grompp -f min_sd.mdp -c <protein>_ions.gro -p topol.top -o <protein>_ions_minsd.tpr -maxwarn 1
```
```
gmx mdrun -deffnm <protein>_ions_minsd
```
INITIALISATION AND HEATING
```
gmx grompp -f charmm-inputs/charmm_nvt_heat.mdp -c <protein>_ions_minsd.gro -r <protein>_ions_minsd.gro -p topol.top -o <protein>_nvt_heat.tpr
```
```
gmx mdrun -deffnm <protein>_nvt_heat
```




## 3. Equilibration and MD Simulation Run 

>[!IMPORTANT]
> The scripts provided for the following steps are designed to be used on an HPC with Slurm scheduler.

>[!NOTE]
>You need an ```minimization.mdp```, ```input.gro```, ```index.ndx```, ```topol.top```, ```equilibration.mdp``` and ```production.mdp``` file to do a MD simulation. There are templates for ```.mdp``` in ```MD/templates``` folder. Adjust to your liking. The following MD-scripts use these templates. If you want to use other ```.mdp``` files, change the scripts accordingly  (```s1_make_equilibration_scripts.py``` and ```s3_make_production_scripts.py```)

**1. Equilibration:**
  - Equilibrate the system to stabilize temperature, pressure, and density.
  - Run ```s1_make_equilibration_scripts.py``` to create Slurm scripts in each protein folder for equilibration.
  - Initial Testing: Test with one protein to identify potential issues before full batch submission. Change the project name and email address in s2. The email sent when the equilibration started and ended will tell how long the equilibration took, maybe adjust the time as needed.
  - Execute ```s2_make_equilibration_execute_script.py``` to submit all generated equilibration scripts. Then run ```execute_equilibration.sh``` for bulk submission of all equilibrations.

**2. Production Run:**
  - Perform the production simulation under stable, equilibrated conditions.
(If you change dt to 0.002 ps, and you want your simulation to last for 100 ns, then nsteps is 100,000 ps / 0.002 ps = 50,000,000)
  - Use ```s3_make_production_scripts.py``` to create production scripts for each protein.
>[!IMPORTANT]
> Adjust Slurm header in the script!

> [!TIP]
> Initial Testing: Test with one protein to identify potential issues before full batch submission. Change the project name and email address in ```s3```. Check Log Files associated with each job ID to estimate run completion time and adjust time/array accordingly.

  - Create a script to bulk submit all MD-simulations by running ```s4_make_production_execute_script.py```. Then run ```execute_production.sh``` to bulk submit all production runs. 
  - Adjustments: Set the Slurm array size for efficient job management:
      - ```#SBATCH --array=0-6%1```  Adjust '6' based on protein size or simulation length. For a protein from roughly  300 AA 6 has proven to be enough





## 4. Post-Simulation Processing
This can be done either using the slurm script on an HPC ```s5_processing-slurm.py``` or on your local laptop or desktop PC using ```s5_processing.py```. Both of these scripts create ```execute_processing.sh``` which performs the processing on each protein. 
Alternatively, you can run the following commands manually on each protein: 
1. Trajectory Concatenation: Concatenate all trajectory segments into a single trajectory file for continuity. Use GROMACS to concatenate trajectory files for each protein, combining individual parts into a complete trajectory:
 ```
 gmx trjcat -f <name>.part*.xtc -o <name>.all.xtc
 ``` 
2. Protein Extraction and PBC Correction: Transform the concatenated trajectory to remove artefacts caused by periodic boundaries. Centre the protein in the simulation box, ensuring molecules are visually and structurally continuous. Align the protein across frames by removing global rotation and translation, making it easier to analyze internal motions. Use GROMACS to extract the protein coordinates and apply periodic boundary corrections:
```
gmx trjconv -f <name>.all.xtc -ur compact -n index.ndx -o <name>.prot.xtc -pbc mol -center -s <name>.tpr
```
and fit the protein for final analysis: 
``` 
gmx trjconv -f <name>.prot.xtc -n index.ndx -o <name>.proc.xtc -s <name>.tpr -fit rot+trans
``` 
3. Extract final structure: Extract a single representative structure from the trajectory at a specified time point (ideally the starting point) to view the simulation. Save ```.gro``` file to view the simulation in VMD:
```
gmx trjconv -f <name>.prot.xtc -b 0 -e 0 -n index.ndx -o <name>.proc.gro -s <name>.tpr
``` 






## 5. Bootstrapping and Structural Phylogeny Analysis

This section describes the process for generating structural phylogenies with bootstrapping. Bootstrapping helps assess the stability of phylogenetic relationships by generating multiple replicates, each based on random samples of simulation frames. The following steps require the completed MD simulations with trajectory (```.xtc```) and structure (```.gro```) files.
1. Organize Input Files:
Ensure all ```.xtc``` and ```.gro``` files for each protein are collected in a folder named ```Trajectories```. Each trajectory file should correspond to a unique protein.
2. Prepare Bootstrapping Replicates with ```BS_master_sort.py```:
  - Setup: Place the ```BS_master_sort.py``` script outside the ```Trajectories``` folder to organize trials and generate the required number of bootstrap samples. This script will select random frames from each trajectory, creating directories with the sampled structures for each trial.
  - Parameter Adjustments:
    - Number of Trials: Define the total number of bootstrap replicates by setting the N_trials variable in BS_master_sort.py. For example:
    ```N_trials = 200  # Set to the desired number of bootstrap trials```
    ```dir_trial = np.arange(0, 200, 1)  # Creates folders for each trial```
    This configuration will generate 200 bootstrap trials, each with randomly selected frames.
    - Frame Selection: Adjust the number of frames per trial based on the length of your MD simulation: ```frame_sel = np.random.randint(0, 10000, N_trials)  # Modify '10000' to match total frames of your simulation```
>[!Tip]
>To determine the correct frame count, load the .gro and .xtc files into a viewer like VMD. Once fully loaded, VMD will display the total frame count for the trajectory. Alternatively, use gmx check on the .xtc file, or calculate the number of frames from your .mdp file as ```nsteps/nstxout-compressed```).
  - Generate Bootstrapping Directories:
    Run
    ```
    python3 BS_master_sort.py
    ```
    Thus will create a series of directories named trial_0, trial_1, etc., each containing a unique set of frames selected randomly from the original trajectory. 
3. Run Structural Analysis on Trials:
  - With the trials generated, use ```GesamtTree.py``` to analyze each trial folder within the ```Trajectories/Trials``` directory. This script will generate phylogenetic trees for each bootstrap replicate.
  - Execution Command:
This command
```
python3 GesamtTree.py Trajectories/Trial/trial_*
```
runs ```GesamtTree.py``` on each trial directory, producing a phylogenetic tree output for every bootstrap sample.

4. Generate a Consensus Phylogenetic Tree and Map Bootstrap Values:
  - Once phylogenetic trees have been generated for all bootstrap trials, you can either:
    - Create a Majority-Rule Consensus Tree from the bootstrap replicates, which combines all individual trees to identify the most consistent relationships.
    - Map Bootstrap Values onto a Reference Tree to reflect the stability of nodes in the reference structure.
  - Options:
    - **Majority-Rule Consensus Tree:** This option compiles the bootstrap trial trees into a single consensus tree, with nodes representing relationships that appear in a specified percentage of the bootstrap samples. For example, setting a 60% threshold (```-f 0.6```) will include only those relationships that appear in at least 60% of the trees.
```
sumtrees -s consensus -o consensus60 -f 0.6 -p -d0 Trajectories/Trials/trial_*/Formatted.nex
```
This produces a consensus tree (consensus60) highlighting the most frequent relationships across bootstrap trials.

  - **Mapping Bootstrap Values onto a Reference Tree:** If you have a reference tree from Section 3.1, you can add Bootstrap support values directly onto it rather than generating a separate consensus tree. This approach maps the percentage of bootstrap trees supporting each node directly onto the reference tree nodes, indicating each branch's stability. 
```
sumtrees -d0 -p -o OutputTree -t Referencetree/<ReferenceTree> Trajectories/Trials/trial_*/Formatted.nex
```
  - Here:
    - ```OutputTree``` will contain the reference tree structure with bootstrap values annotated on each node.
    - Interpretation: Higher bootstrap values on a node in the reference tree suggest strong support for that relationship across bootstrap replicates, while lower values indicate less stability.


# Tutorial 

Coming soon
