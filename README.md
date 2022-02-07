## Protein analysis workflow

Analsing a protein family of the users choice without using BioPython.
The user can:
-Download sequences from a user-defined family of proteins and taxonomic group
-Produce a phylogenetic tree
-Scan the sequences for motifs of interests from the PROSITE database (using EMBOSS patmatmotifs)
-Predict and plot transmembrane sequences (using EMBOSS tmap)
-Calculate the isoelectric point of the sequences (using EMBOSS iep)
-Calculate the hydrophobic moment for the chosen proteins (using EMBOSS hmoment)
-Plot a superimposed hydropathy plot (using EMBOSS pepwindowall)
The code is written to prevent the user when making mistakes, facilitating analysis 
