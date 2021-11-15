# Tutorial 01
## Overview
This example described two simulations produced by all-atom classical molecular dynamic of a human chemokine type 5 receptor (CCR5), a HIV-1 glycoprotein 120 (gp120) and a human cluster of differentiation 4 (CD4) complex. The only difference between the two simulations is the amino acid sequence of gp120. The first complex has the #25 gp120 variant and the second has the #34 gp120 variant. A third entry described a single structure of an unliganded CCR5.
For more details about these systems, please refer to the following articles:

-   Colin _et al._ CCR5 structural plasticity shapes HIV-1 phenotypic properties. _PLOS Pathog_. **2018**, _14_, e1007432. DOI: [10.1371/journal.ppat.1007432](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1007432)
-   Jacquemard _et al._ Modeling of CCR5 Recognition by HIV-1 gp120: How the Viral Protein Exploits the Conformational Plasticity of the Coreceptor. _Viruses_. **2021**, _13_, 1395. DOI: [10.3390/v13071395](https://www.mdpi.com/1999-4915/13/7/1395/htm)

## Aim
The gp120#25 and gp120#34 variants show distinct phenotypic properties and tropism. The aim of this tutorial is to evaluate the conformational response of CCR5 to the binding of different gp120 variants.

The helix ends of the CCR5 7 helices will be projected with ATOLL.

## Inputs
### Entries

- CCR5-**gp120#25**-CD4 structures produced by an all-atom molecular dynamic simulation (5 × 100 ns, 100 frames). Waters and ions has been removed.
Path : inputs/25/align-dry.prmtop and inputs/25/sim.nc
- CCR5-**gp120#34**-CD4 structures produced by an all-atom molecular dynamic simulation (5 × 100 ns, 100 frames). Waters and ions has been removed.
Path : inputs/34/align-dry.prmtop and inputs/34/sim.nc
- **unliganded** CCR5 structure from the last frame of an all-atom molecule dynamic simulation (100 ns). Waters and ions has been removed.
Path : inputs/free.pdb

### Reference structure
The reference structure was downloaded from the PDB database (entry id: 4MBS). This structure describes an engineered CCR5 receptor with rubredoxin in the inactive state bound to an inverse agonist. The structures was manually placed to align lipid bilayer onto the XY plane. Unnecessary molecular objects of the structure were removed (water, lipids, ligand and rubredoxin).

### Sequence alignment
The studied protein is CCR5 in the three entries. In the sequence alignment file (sequences.sto) only the CCR5 full sequence is given and is set as reference. The sequence comes from the Uniprot database (Uniprot AC: P51681).

### Residues used for structural superimposition
The superimposition of structures is already done. No alignment procedure is used with ATOLL in this tutorial.
The selected residues for the superimposition are those inserted in the membrane and are assumed te be rigid. These residues are (Uniprot numbering P51681):
- TM1: 31 to 57
- TM2: 64 to 88
- TM3: 99 to 129
- TM4: 143 to 164
- TM5: 190 to 219
- TM6: 235 to 256
- TM7: 277 to 299

### Residues used for projection
For this study, the helices definition is given by Uniprot (P51681). However, we reduced the range of helices depending on the secondary structure of CCR5 during molecular dynamic simulations. We computed the secondary structure of CCR5 using the DSSP algorithm. If a residue adopts a alpha-helix structure in more than 50% of the frames, then this residue is considered as part of the helix. Otherwise, the residue is excluded. The residues used for the projection are (Uniprot numbering P51681):
- TM1: 26 to 57
- TM2: 64 to 89
- TM3: 98 to 131
- TM4: 142 to 165
- TM5: 187 to 223
- TM6: 228 to 259
- TM7: 269 to 300