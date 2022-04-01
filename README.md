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

<img src="/XERecombinant.png" width="300">   
The position on the genome is represented on the X axis and the conditional probability *P(Lineage|Mutation)* is represented on the Y axis (i.e. probability the virus belongs to a certain lineage given a certain mutation). Each point represents a mutation and the colours represent the most likely lineage.   
As you can see, there are mutations which are *"very exclusive"* (i.e. those with very high conditional probability), while others are more *"promiscuous"* (i.e. those with low predictive power / conditional probability). Based on the mutations plot, we can infer that the first part of the genome belongs to a *B.1* (or *B.1.1* more precisely), while the second part of the genome (after position 12K belongs to a *B.2* lineage).  
This *BA.1-BA.2* recombinant has been already circulating around the world and the new *pango-lineage* *XE* has been proposed for it. 

>*XE is a BA.1/BA.2 recombinant, with the majority of the genome including the S gene belonging
to BA.2. XE shows evidence of community transmission within England, although it is currently
<1% of total sequenced cases. Early growth rates for XE were not significantly different from
BA.2, but using the most recent data up to 16 March 2022, XE has a growth rate 9.8% above
that of BA.2. As this estimate has not remained consistent as new data have been added, it
cannot yet be interpreted as an estimate of growth advantage for the recombinant. Numbers
were too small for the XE recombinant to be analysed by region.*
[SARS-CoV-2 variants of concern and variants under investigation in England (#39/25.03.2022)](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1063424/Tech-Briefing-39-25March2022_FINAL.pdf)

## Current models
We distribute Precfinder with a set of pretrained models:
*Norway10K*   
*Norway20K*   

We expect to increase the number of avaible models and to update them regularly. 

## Advanced instructions

## Under the hood

## Understanding the limitations of the method

## Going further
