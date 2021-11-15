# Tutorial 03
## Overview
This example described 26 structures of ligand-gated P2X receptor from the Protein Data Bank. These structures describe 4 structural states:

- open receptor with bound ATP
- closed unliganded receptor
- closed receptor with bound ligand
- desensitized receptor

## Aim
The aim of this tutorial is to evaluate the conformational change of P2X receptors depending on the structural state.

A P2X receptor is a homotrimeric protein containing 2 transmembrane helices per monomer. In this tutorial, the structures will be superimposed and the 6 helices will be projected with ATOLL.

## Inputs
### Entries
The ATOLL entry structures should contain a single chain. Hence, the three chains of P2X receptor was merged into single one and residues were renumbered consecutively. At the early stage of ATOLL development this merging was done by hand. Now it is possible to do this automatically with the "sanitize" tool of ATOLL software.

A detailed list of entries is providen in the table.docx file.

All input structures are stored in the "structures" directory.

### Reference structure
The reference structure was downloaded from the OPM database (entry id: 6U9V). This structure describes the P2X7 receptor of rat (Uniprot id: P2X7_RAT) in the close state. The structures from OPM database were already in a correct coordinate frame (lipid bilayer aligned onto the XY plane). Unnecessary molecular objects of the structure were removed (water, bilayer representation, lipids, cofactor and ions).

### Sequence alignment
The 26 structures of P2X receptor describe 6 proteins (amP2X_GULFCOASTTICK, P2X3_HUMAN, P2X4_ZEBRAFISH, P2X7_CHICKEN, P2X7_GIANTPANDA and P2X7_RAT). The 6 sequences obtained from Uniprot and those from the 26 structures were aligned using the multiple sequence alignment tool of MOE 2019.01. It's noteworthy that the sequences from Uniprot were triplicated to represent the three chains of the P2X receptors. The reference sequence selected is the P2X7_RAT. Then only the aligned sequences from Uniprot were saved in Stockholm format (sequences.sto)

### Residues used for structural superimposition
The superimposition is done during the ATOLL procedure. The selected residues are located in the extracellular part of the proteins near disulfide bridges.

These residues are (residue position of the reference P2X7_RAT sequence in the sequences.sto file):
- 64 to 68
- 126 to 128
- 160 to 164
- 175 to 178
- 181 to 183
- 201 to 206
- 208 to 210
- 246 to 248
- 269 to 274
- 280 to 283
- 291 to 295
- 309 to 319
- 334 to 344
- 462 to 466
- 524 to 526
- 558 to 562
- 573 to 576
- 579 to 581
- 599 to 604
- 606 to 608
- 644 to 646
- 667 to 672
- 678 to 681
- 689 to 693
- 709 to 717
- 732 to 742
- 860 to 864
- 922 to 924
- 956 to 960
- 971 to 974
- 977 to 979
- 997 to 1002
- 1004 to 1006
- 1042 to 1044
- 1065 to 1070
- 1076 to 1079
- 1087 to 1091
- 1105 to 1115
- 1130 to 1140

### Residues used for projection
The projected residues are those inserted in the membrane and are the extremities of helices. In some structures, Intracellular part of the receptor is not resolved until lower part of helices inserted in the membrane. Hence only the upper and the middle of helices are projected.

These residues are (residue position of the reference P2X7_RAT sequence in the sequences.sto file):

- TM1a: 33 to 47
- TM2a: 352 to 368
- TM1b: 431 to 445
- TM2b: 750 to 766
- TM1c: 829 to 843
- TM2c: 1148 to 1164