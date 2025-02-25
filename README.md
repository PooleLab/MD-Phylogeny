  # MD-Phylogeny
MD-Phylogeny constructs structural phylogenies and provides statistical confidence through molecular dynamics (MD) simulations.

## Publication
Read about this method: 

## System Requirements 
Operating system: Ideally, UNIX or Mac (not tested on Mac, yet). If you are using Windows, you can use the Windows Subsystem for Linux (WSL/WSL2), but be aware of possible difficulties. 
For most preparatory and analytical steps it is ok to use laptop or desktop. But for molecular dynamics (MD) simulations an HPC is highly recommended. 

MD Engine (best to use on an HPC, as it gets too computationally intense for a Laptop or Desktop PC): Be aware that our scripts are developed for GROMACS version 2020 and will require significant adaptation for use with other MD engines

GROMACS (we have used GROMACS version GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid) and CHARMM force fields (we have used CHARMM36 for GROMACS). (https://doi.org/10.5281/zenodo.3562512)



## Installations 

#### Visualisation of structures and simulations
Visual Molecular Dynamics (VMD) (https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) or Pymol (https://www.pymol.org/) to visualise structures and VMD to visualise simulation trajectories. 
#### Structural superposition and calculation of Qscore:
Gesamt (from https://www.ccp4.ac.uk/html/gesamt.html)
#### Visualisation of Phylogenies: 
FigTree (http://tree.bio.ed.ac.uk/software/figtree/) or iTOL (https://itol.embl.de/)
#### Summarization and annotation of phylogenies: 
sumtrees (https://jeetsukumaran.github.io/DendroPy/programs/sumtrees.html)
#### Additional requirements: 
Python3 with the following packages: 

NumPy (https:/2/numpy.org), dendropy (https://pypi.org/project/DendroPy), Bio Python (https://biopython.org/docs/1.75/api/Bio.html), natsort (https://pypi.org/project/natsort)



# Method
## 1. Dataset selection
### Identify a Suitable Dataset:
The dataset should consist of **structures that are likely to share homology**. Consider using tools and databases such as:

Foldseek (https://search.foldseek.com/search , https://github.com/steineggerlab/foldseek), Structome (https://biosig.lab.uq.edu.au/structome_q/), SCOP (https://www.ebi.ac.uk/pdbe/scop/), CATH (https://www.cathdb.info/), InterPro (https://www.ebi.ac.uk/interpro/)

### Evaluate Dataset Suitability:
Ensure that the dataset lies within the "twilight zone" of sequence similarity, where traditional sequence-based phylogenetic methods may be unreliable. Therefore e.g. perform an all-against-all pairwise sequence comparisons. Visualize the sequence similarity network using tools such as: CLANS and MMseqs2

:x: Well-connected network → **Use conventional sequence-based phylogenetics**.

:heavy_check_mark: Poorly connected network → **Suitable for MD-phylogeny**.

### Additional notes:
Ensure the dataset meets the following criteria:
- Structures should consist of single-chain files
- Identify structures that contain multiple copies of the same motif or fold. This influences the analysis
- Avoid structures that are too large (e.g., ribosomes, entire nucleosomes) as they are impractical for MD simulations. 
- Ensure length homogeneity across structures

### Structural comparison and topology maps:
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
- Create phylogeny: Run the GesamtTree.py script on the specified folder to generate a structural phylogeny based on the PDB files
  ```
  python GesamtTree.py <FolderName>
  ```

Then, molecular dynamics is used to generate alternative conformations, which can then be sampled randomly to build alternative trees for bootstrap support in phylogenetic analyses.

## 3. MD simulation set-up
MD setup prepares a molecular system for simulation, ensuring it is physically and realistically configured.
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
> CHARMM-GUI requires an active account for input generation. You need to set this up before using the CHARMM GUI Solution Builder.
> 

Follow the steps and adjust to your requirements, then download the files for simulation. You need  ```<protein>_nvt_heat.gro```, ```index.ndx``` and ```topol.top``` files to continue with the MD simulation run. (Rename them according to the scripts in the next step or change file names in the script according to your files.)

### Option 2: Using scripts
> [!CAUTION]
> Scripts only work if the PDB files contain only standard protein residues, with no engineered amino acids, waters, ions, etc.

1. Make the force field available: GROMACS has CHARMM27 as a default option. While this force field is old, it can be used. However, we recommend using the newer CHARMM36 or CHARMM36m versions. To make these available to GROMACS, they need to be downloaded from http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs and, once unzipped, the resulting folder should be kept in the resulting force field folder of gromacs (it depends on where gromacs is installed; the standard installation location is /usr/share/gromacs/top). Alternatively, e.g. on an HPC, the force field folder can be placed in the working directory. 
It is recommended to try the simulation setup steps to see which number your desired force field has when listed by pdb2gmx; it has to be adjusted in the scripts.
2. Run Simulation Setup: Use ```MD/sim_setup_slurm.py``` to configure system requirements, Slurm header, and other parameters based on the computational resources.

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
You then need the ```<protein>_nvt_heat.gro```, ```index.ndx``` and ```topol.top``` files to continue with the MD simulation run. (Rename them according to the scripts in the next step or change file names in the script according to your files.)


## 3. MD Simulation Run 
> [!IMPORTANT]
> The scripts used in this step are adjusted to the set-up using CHARMM GUI website. Adjust paths and file names according to how you set them up.
> You need an ```minimization.mdp```, ```input.gro```, ```index.ndx```, ```topol.top```, ```equilibration.mdp``` and ```production.mdp``` file to do a MD simulation. There are templates for ```.mdp``` in ```MD/templates``` folder. Adjust to your liking, and make sure you adjust the names of the files themselves or adjust filenames in the scripts.

1. Directory Organization: Create a unique directory for each protein and place the ```charmm-gui.tgz``` file inside (as mentioned in 3.2 point 4)
2. Upload to HPC and extract files: Transfer the folder to the HPC cluster and extract the files. Unzip the CHARMM-GUI setup files in each folder (as mentioned in 3.2 point 4)
The folder structure will be ```<proteinfoldername>/charmm*/gromacs```
3. Equilibration:
  - Equilibrate the system to stabilize temperature, pressure, and density.
  - Run ```s1_make_equilibration_scripts.py``` to create Slurm scripts in each GROMACS folder for equilibration.
> [!NOTE]
> Be aware of changing the directory in the script (line 5) and all SBATCH settings to liking and requirements. Also, change the names of files required for simulation. (line 33 (```minimization.mdp```, ```input.gro```, ```index.ndx```, ```topol.top```) and 38 (```equilibration.mdp```) )
  - Initial Testing: Test with one protein to identify potential issues before full batch submission. Change the project name and email address in s2. The email sent when the equilibration started and ended will tell how long the equilibration took, maybe adjust the time as needed.
  - Execute ```s2_make_equilibration_execute_script.py``` to submit all generated equilibration scripts.
> [!NOTE]
> Adjust the path in line 3 to your file locations
  - Check the trajectories for these parameters for any abnormalities 
4. Production Run:
  - Perform the production simulation under stable, equilibrated conditions.
(If you change dt to 0.002 ps, and you want your simulation to last for 100 ns, then nsteps is 100,000 ps / 0.002 ps = 50,000,000)
  - The folder structure will be ```<proteinfoldername>/charmm*/gromacs```
  - Use ```s3_make_production_scripts.py``` to create production scripts for each protein.
  - Submit all production runs with ```s4_make_production_execute_script.py```.
  - Adjustments: Set the Slurm array size for efficient job management:
      - ```#SBATCH --array=0-6%1```  Adjust '6' based on protein size or simulation length. For a protein from roughly  300 AA 6 has proven to be enough
> [!TIP]
> Initial Testing: Test with one protein to identify potential issues before full batch submission. Change the project name and email address in ```s3```. Check Log Files associated with each job ID to estimate run completion time and adjust time/array accordingly.




## 4. Post-Simulation Processing

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
3. Extract final structure: Extract a single representative structure from the trajectory at a specified time point (ideally the starting point) to view the simulation. Save gro file to view the simulation in VMD:
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

In the folder ```tutorial``` you find 5 ```.pdb``` structures of histone folds. 
Take the time to have a look at the structures in VMD or Pymol. 

CLANS
Sequence alignment

## Build a phylogeny: 
```
mkdir pdbs
```
Copy all pdbs in this folder
```
cp *.pdb pdbs/.
```
Build Phylogeny:
```
python3 Gesamt.py pdbs/
```


## Molecular dynamics simulation:
 Setting up the MD simulation can be done on you local machine. 

Navigate to the ```tutorial``` folder. Then run 
```
Python3 sim_setup.py
```
This creates a folder for each protein, places the pdb in that folder and a .sh script for the setup for this protein. 
Run every ```_sim_setup.sh``` script in every folder. (This step can also already be run on an HPC. Copy the tutorial folder to an HPC and run the ```s0_sim_setup-slurm.py``` script. This script is designed for a slurm sceduler and the scripts uses GROMACS (module load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid). Anything else would need adjustments of this and all following scripts. Adjust script (account, mem-per-cpu,.. e-mail address) and run)
For example: 
```
cd 1p3m_A_H3/
```
```
bash 1p3m_A_H3_sim_setup.sh
```
Repeat for the other folders. 
Be aware that this may take some time. The last step can take 30-60 min depending on you mashine. 
Check outcome: therefore, open VMD, ``` file```>```new molecule```>```browse``` select ```1p3m_A_H3_nvt_heat.gro``` and again and select ```1p3m_A_H3_nvt_heat.xtc```
For the next step you need the last ```.xtc``` file that has been created plus the ```.gro``` ```*_nvt_heat.gro/.xtc```

Then run Equilibration (We found best to do this on a HPC. So if you have access, then copy all folders and file there and continue with the equilibration scripts.)
```
pyhton3 s1_make_equilibration_scripts.py
```
This creates a script for the Equilibration in every protein folder. Run one script to see if errors occur. If not, you can use the next script to submit all at once.
```
python3 s2_make_equilibration_execute_script.py
```
Then run: 
```
bach execute_equilibration.sh
```
The files created are called ```_EQ```. You can check the equilibration again in VMD by loading the ```*.xtc``` and ```*.gro``` file for each protein. 
Then the proteins are ready for simulation. Create simulation scripts by running: 
```
python3 s3_make_production_scripts.py
```
This creates ```_Prod.sh``` scripts in every protein folder. Run one single production run to check if there are errors. Then you can submit all other producton runs at the same time by running: 
```
python3 s4_make_production_execute_script.py
```
```
bash execute_production.sh
```
These simulations are set to 5 ns (you can check this in the ```templates/production_template.mdp``` file. ```nsteps     = 2500000   ; 2 * 2500000 = 5 ns```)
This is a very short simulation for this tutorial. For an actual analysis you would want to set this number up. We have found 100ns to work best. But we also recommend to try with one protein first and check the outcome.) 
Wait for the simulations to finish. This may take 1-2 days. If you don't want to wait, we have provided all files the folder ```full_simulation_files```
Be aware of the size of this folder (2.6 GB)


