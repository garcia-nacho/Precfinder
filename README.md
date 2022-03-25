# Probabilistic Recombinant Finder (Precfinder)

## Description
Precfinder is a tool to identify SARS-CoV-2 recombinants.    
Note that Rrecfinder is still on a very early stage and it could suffer drastic changes in the future.  

## Instructions

1. Clone the Precfinder repo 

<code>git clone https://github.com/garcia-nacho/Precfinder</code>

2. Copy a multifasta file with the sequences that you want to test inside the folder Precfinder/Inference

<code>cp Sequences.fa Precfinder/Inference </code>

3. Run the script Inference.sh inside the Precfinder folder

<code>cd Precfinder && ./Inference.sh</code>


## Output

After finishing the script will save three new files inside the Precfinder/Inference folder
XXX
XXX
XXX

The file XXX includes the output of nextclade and pangolin for each sequence
The file XXX includes recombinant prediction. The score ranges from 0 to 1, the closer it is to one the most likely for the sequence to be a recombinant.
The file XXX includes a plot that visualizes the different mutations found in the sequence and the conditional probability for the different lineages.


## Advanced Instructions


## Under the hood

