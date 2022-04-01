# Probabilistic Recombinant Finder (Precfinder)

## Description
Precfinder is a tool to identify SARS-CoV-2 recombinants.   
Precfinder is distributed as docker image with some wrapping scripts to install it, run it and train additional models.   
**Note that Rrecfinder is still on a very early stage, it will be buggy and it could suffer drastic changes in the future.**

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

*Recombinant_Prediction.xlsx*    
It contains the recombinant prediction for each sequence. Including two classes (Recombinant/Non Recombinant) and a prediction score. The closer the score to 1 or to 0, the most sure the model about the classification. **It is NOT a probabilistic score.**. There is a third column to display some warnings (e.g. the absence of observed mutations on the training set, the presence of too many mutations, etc)        

*Recombinant_PredictionwPango.xlsx*    
Same as above but it includes pangolin lineage.  

*Inference_dataset.csv*    
It is a temporary file that contains nextclade output alonside pangolin lineage. Note that the insertions and deletions have been reencoded and included under the column substitutions. (e.g. A deletion of 3 nucleotides on the position 123 (*123-125*) would be reencoded as *D123C*. An insertion of 6 nucleotides at position 1200 (*1200:AAAAAA*) will become *I1200F* under the substitutions column)    

*Date_RecombinantPlots.pdf*   
It contains a plot for each of the sequences present on the multifasta file.

*Date_RecombinantSelected.pdf*   
It contains a plot for each sequence predicted to be a recombinant.  

## Interpreting the plots

The pdfs contains plots similar to this one: 

![EX](./XE.jpg)

## Current models
We distribute Precfinder with a set of pretrained models:
*Norway10K*   
*Norway20K*   

We expect to increase the number of avaible models and to update them regularly. 

## Advanced Instructions

## Under the hood

