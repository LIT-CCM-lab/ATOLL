# Tutorial 01
## Overview
In this example, we are going to analyse two simulations produced by all-atom classical molecular dynamic of the complex formed between the human chemokine type 5 receptor (CCR5), HIV-1 glycoprotein 120 (gp120) and human cluster of differentiation 4 (CD4). The two simulations differ in the amino acid sequence of gp120 (#25 or #34 gp120 variant). The MD frames are compared to the structure of unliganded CCR5.
For more details about these systems, please refer to the following articles:

-   Colin _et al._ CCR5 structural plasticity shapes HIV-1 phenotypic properties. _PLOS Pathog_. **2018**, _14_, e1007432. DOI: [10.1371/journal.ppat.1007432](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1007432)
-   Jacquemard _et al._ Modeling of CCR5 Recognition by HIV-1 gp120: How the Viral Protein Exploits the Conformational Plasticity of the Coreceptor. _Viruses_. **2021**, _13_, 1395. DOI: [10.3390/v13071395](https://www.mdpi.com/1999-4915/13/7/1395/htm)

## Aim
The gp120#25 and gp120#34 variants show distinct phenotypic properties and tropism. The aim of this tutorial is to evaluate the conformational response of CCR5 to the binding of different gp120 variants.

The helix ends of the 7 transmembrane helices of CCR5 will be projected with ATOLL.

## Inputs
### Entries

- CCR5-**gp120#25**-CD4 structures produced by an all-atom molecular dynamic simulation (5 × 100 ns, 100 frames). Waters and ions has been removed.
Path : inputs/25/align-dry.prmtop and inputs/25/sim.nc
- CCR5-**gp120#34**-CD4 structures produced by an all-atom molecular dynamic simulation (5 × 100 ns, 100 frames). Waters and ions has been removed.
Path : inputs/34/align-dry.prmtop and inputs/34/sim.nc
- **unliganded** CCR5 structure from the last frame of an all-atom molecule dynamic simulation (100 ns). Waters and ions has been removed.
Path : inputs/free.pdb

### Reference structure
The reference structure is downloaded from the PDB database (entry id: 4MBS). This structure describes an engineered CCR5 receptor (rubredoxin insertion) in the inactive state bound to an inverse agonist. The structures was manually oriented to align the transmebrane domains perpendicular to the XY plane. Water, lipids, ligand and rubredoxin were removed from the coordinate file.

### Sequence alignment
The sequence alignment file (sequences.sto) shows CCR5 full sequence (Uniprot AC: P51681), tagged as reference.

### Residues used for structural superimposition
Input structures are already superimposed in this tutorial. No alignment is needed.
The domains used for the superimposition are transmembrane domains, defined as follows (Uniprot numbering P51681):
- TM1: 31 to 57
- TM2: 64 to 88
- TM3: 99 to 129
- TM4: 143 to 164
- TM5: 190 to 219
- TM6: 235 to 256
- TM7: 277 to 299

### Residues used for projection
Since helix lenght can change during simulation, we re-defined tranmembrane helices boundaries according to the secondary structure, as evaluated using the DSSP algorithm. An alpha-helix thus include a continuous stretch of residue which adopts an alpha-helix conformation in more than 50% of the frames. The residues used for the projection are (Uniprot numbering P51681):
- TM1: 26 to 57
- TM2: 64 to 89
- TM3: 98 to 131
- TM4: 142 to 165
- TM5: 187 to 223
- TM6: 228 to 259
- TM7: 269 to 300
