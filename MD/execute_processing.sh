#!/bin/bash
module load GROMACS/2020.5-intel-2020a-cuda-11.0.2-hybrid
# Processing files for 4uuz_B_H4
cd 4uuz_B_H4/
gmx trjcat -f 4uuz_B_H4.part*.xtc -o 4uuz_B_H4.all.xtc
gmx trjconv -f 4uuz_B_H4.all.xtc -ur compact -n index.ndx -o 4uuz_B_H4.prot.xtc -pbc mol -center -s 4uuz_B_H4.tpr <<EOF
0
0
EOF
gmx trjconv -f 4uuz_B_H4.prot.xtc -n index.ndx -o 4uuz_B_H4.proc.xtc -s 4uuz_B_H4.tpr -fit rot+trans <<EOF
0
0
EOF
gmx trjconv -f 4uuz_B_H4.proc.xtc -b 0 -e 0 -n index.ndx -o 4uuz_B_H4.proc.gro -s 4uuz_B_H4.tpr <<EOF
0
EOF
cd ..

# Processing files for 4kud_A_H3
cd 4kud_A_H3/
gmx trjcat -f 4kud_A_H3.part*.xtc -o 4kud_A_H3.all.xtc
gmx trjconv -f 4kud_A_H3.all.xtc -ur compact -n index.ndx -o 4kud_A_H3.prot.xtc -pbc mol -center -s 4kud_A_H3.tpr <<EOF
0
0
EOF
gmx trjconv -f 4kud_A_H3.prot.xtc -n index.ndx -o 4kud_A_H3.proc.xtc -s 4kud_A_H3.tpr -fit rot+trans <<EOF
0
0
EOF
gmx trjconv -f 4kud_A_H3.proc.xtc -b 0 -e 0 -n index.ndx -o 4kud_A_H3.proc.gro -s 4kud_A_H3.tpr <<EOF
0
EOF
cd ..

# Processing files for 4kud_D_H2B
cd 4kud_D_H2B/
gmx trjcat -f 4kud_D_H2B.part*.xtc -o 4kud_D_H2B.all.xtc
gmx trjconv -f 4kud_D_H2B.all.xtc -ur compact -n index.ndx -o 4kud_D_H2B.prot.xtc -pbc mol -center -s 4kud_D_H2B.tpr <<EOF
0
0
EOF
gmx trjconv -f 4kud_D_H2B.prot.xtc -n index.ndx -o 4kud_D_H2B.proc.xtc -s 4kud_D_H2B.tpr -fit rot+trans <<EOF
0
0
EOF
gmx trjconv -f 4kud_D_H2B.proc.xtc -b 0 -e 0 -n index.ndx -o 4kud_D_H2B.proc.gro -s 4kud_D_H2B.tpr <<EOF
0
EOF
cd ..

# Processing files for 1p3m_A_H3
cd 1p3m_A_H3/
gmx trjcat -f 1p3m_A_H3.part*.xtc -o 1p3m_A_H3.all.xtc
gmx trjconv -f 1p3m_A_H3.all.xtc -ur compact -n index.ndx -o 1p3m_A_H3.prot.xtc -pbc mol -center -s 1p3m_A_H3.tpr <<EOF
0
0
EOF
gmx trjconv -f 1p3m_A_H3.prot.xtc -n index.ndx -o 1p3m_A_H3.proc.xtc -s 1p3m_A_H3.tpr -fit rot+trans <<EOF
0
0
EOF
gmx trjconv -f 1p3m_A_H3.proc.xtc -b 0 -e 0 -n index.ndx -o 1p3m_A_H3.proc.gro -s 1p3m_A_H3.tpr <<EOF
0
EOF
cd ..

# Processing files for 6m2m_E_H2A
cd 6m2m_E_H2A/
gmx trjcat -f 6m2m_E_H2A.part*.xtc -o 6m2m_E_H2A.all.xtc
gmx trjconv -f 6m2m_E_H2A.all.xtc -ur compact -n index.ndx -o 6m2m_E_H2A.prot.xtc -pbc mol -center -s 6m2m_E_H2A.tpr <<EOF
0
0
EOF
gmx trjconv -f 6m2m_E_H2A.prot.xtc -n index.ndx -o 6m2m_E_H2A.proc.xtc -s 6m2m_E_H2A.tpr -fit rot+trans <<EOF
0
0
EOF
gmx trjconv -f 6m2m_E_H2A.proc.xtc -b 0 -e 0 -n index.ndx -o 6m2m_E_H2A.proc.gro -s 6m2m_E_H2A.tpr <<EOF
0
EOF
cd ..

