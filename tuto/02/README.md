# Tutorial 02
## Overview
We are going to analyse the deactivation simulation produced by all-atom classical molecular dynamic of a β2-adrenergic receptor (ADRB2). The starting structure of this G protein-coupled receptor is in the active state. During the simulation, the receptor state move to an intermediate and then to the inactive state. The most remarkable structural changes involve the lower part of TM6 that shifts inward the receptor center and the lower part of TM7 that slightly shifts outward the receptor center.

This simulation is described in the following article:
-   Dror _et al._ Activation mechanism of the β2-adrenergic receptor. _PNAS_. **2011**, 61(9), 18684-18689. DOI: [10.1073/pnas.1110499108](https://www.pnas.org/content/108/46/18684.long)

## Aim
The aim of this tutorial is to evaluate the conformational changes of ADRB2 during its deactivation.

The helix ends of the 7 helices will be superimposed and will be projected with ATOLL.

## Inputs
### Entries
Simulation is divided into 3 entries that correspond to the ADRB2 structural state (active, intermediate and inactive). Dror _et al._ proposed the TM6-TM3 distance and the NPxxY region RMSD to distinguish the three states. 

Only Cα are represented in these structures.

- active: 195 frames
- intermediate: 467 frames
- inactive: 139 frames

### Reference structure
The reference structure was downloaded from the OPM database (entry id: 4LDE). This structure describes the ADRB2 receptor in the active state bound to an agonist and a nanobody. The structures from OPM database were already in a correct coordinate frame (lipid bilayer aligned onto the XY plane). Unnecessary molecular objects of the structure were removed (water, bilayer representation, lipids, cofactor, ions, ligand and nanobody).

### Sequence alignment
In the sequence alignment file (sequences.sto) the ADRB2 full sequence from the Uniprot database (Uniprot AC: P07550) is given and is set as reference. Due to mutation in the ADRB2 structures, the amino acid sequence of these structures was aligned onto the native one with MOE 2019.01. The aligned structure sequence was triplicated for each structural state. It's noteworthy that it is not necessary to give a sequence for each entry (see tutorial 01).

### Residues used for structural superimposition
The superimposition is done during the ATOLL procedure.
The selected residues for the superimposition are those inserted in the membrane and are assumed te be rigid. These residues are (residue position of the reference sequence in the sequences.sto file):
- TM1: 16 to 31
- TM2: 51 to 67
- TM3: 86 to 103
- TM4: 132 to 148
- TM5: 177 to 194
- TM6: 255 to 269
- TM7: 289 to 302

### Residues used for projection
For this study, the helices definition is given by Uniprot (P51681). However, we reduced the range of helices depending on the secondary structure of ADRB2 in the reference structure. We computed the secondary structure of CCR5 using the DSSP algorithm. If a residue adopts a alpha-helix structure then this residue is considered as part of the helix. Otherwise, the residue is excluded. The residues used for the projection are (residue position of the reference sequence in the sequences.sto file):
- TM1: 10 to 38
- TM2: 46 to 74
- TM3: 82 to 113
- TM4: 126 to 147
- TM5: 175 to 203
- TM6: 246 to 276
- TM7: 284 to 302
