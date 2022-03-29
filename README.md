# Probabilistic Recombinant Finder (Precfinder)

## Description
Precfinder is a tool to identify SARS-CoV-2 recombinants.
Precfinder is distributed as docker image with some wrapping scripts to install it, run it or to train additional models.
Note that Rrecfinder is still on a very early stage, it will be buggy and it could suffer drastic changes in the future.

## Instructions

1. Clone the Precfinder repo 

<code>git clone https://github.com/garcia-nacho/Precfinder</code>

2. Install it  

<code>Precfinder/Install.sh</code>

3. Copy your multifasta file inside the Precfinder/Inference folder and run the Inference script inside Precfinder's folder *ModelName*

<code>./Inference.sh ModelName</code>

## Output

Precfinder will save five new files inside the Precfinder/Inference folder:   

*Recombinant_Prediction.xlsx*   
*Recombinant_PredictionwPango.xlsx*   
*Inference_dataset.csv*   
*Date_RecombinantPlots.pdf*   
*Date_RecombinantSelected.pdf*   

*Recombinant_Prediction* contains the recombinant prediction for each sequence. Including two classes (Recombinant/Non Recombinant) and a prediction score. The closer the score to 1 or to 0, the most sure the model about the classification. **It is NOT a probabilistic score.**. There is a third column to display some warnings (e.g. the absence of observed mutations on the training set, the presence of too many mutations, etc)        

*Recombinant_PredictionwPango.xlsx* same as above but it includes pangolin lineage.  

*Inference_dataset.csv* is a temporary file that contains nextclade output alonside pangolin lineage. Note that the insertions and deletions have been reencoded (i.e. D123C is a deletion of 3 nucleotides on the position 123, I1200F is the insertion of 6 nucleotides on the position 1200) and included under the column substitutions.   

*Date_RecombinantPlots.pdf* contains a plot for each of the sequences present on the multifasta file.

*Date_RecombinantSelected.pdf* contains a plot for each sequence predicted to be a recombinant.  

## Interpreting the plots

The pdfs contains plots similar to this one: 

## Current models
We distribute Precfinder with a set of pretrained models:
*Norway10K*   
*Norway20K*   

We expect to increase the number of avaible models and to 

## Advanced Instructions

## Under the hood

