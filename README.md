# Probabilistic Recombinant Finder (Precfinder)

## Description
Precfinder is a tool to identify SARS-CoV-2 recombinants. Rrecfinder is still on a very early stage and it could suffer drastic changes in the future.  

## Instructions

1. Clone the Precfinder repo 

git clone https://github.com/garcia-nacho/Precfinder

2. Save a multifasta file with the sequences that you want to test inside the folder Precfinder/Inference

cp Sequences.fa Precfinder/Inference 

3. Run the script Inference.sh inside the Precfinder folder

cd Precfinder && ./Inference.sh

Inside the /Precfinder/Inference folder you will find three new files
XXX
XXX
XXX

The file XXX includes the output of nextclade and pangolin for each sequence
The file XXX includes recombinant prediction. The score ranges from 0 to 1, the closer it is to one the most likely for the sequence to be a recombinant.
The file XXX includes a plot that visualizes the different mutations found in the sequence and the conditional probability for the different lineages.


## Advanced Instructions


## Under the hood

