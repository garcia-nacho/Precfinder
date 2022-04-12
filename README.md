# Probabilistic Recombinant Finder (Precfinder)

:warning:  **IF YOU INSTALLED PRECFINDER BETWEEN 04/04/2022 AND 05/04/2022 YOU MUST INSTALL IT AGAIN**
         

## Description
Precfinder is a tool to identify SARS-CoV-2 recombinants and/or coinfections.   
Precfinder is distributed as docker image with some wrapping scripts to install it, run it and train additional models.   
**Note that Rrecfinder is still on a very early stage, it will be buggy and it could suffer drastic changes in the future.**

## Motivation   
During 2022 we have seen a high *covid-19* prevalence and the coexistance of large viral lineages (e.g. Delta and Omicron, Omicron BA.1 and Omicron BA.2); therefore, the probability of two differente viruses infecting the same person are non-neglectable as it was before.    
In this new context, we decided to create a tool to identify recombinant viruses.   Quickly, we discovered that finding SARS-CoV-2 recombinants is a very challenging task due to the relative low number of mutations between the different lineages, the small viral genome and the relative heterogeneity between samples classified under the same lineage. Moreover, we found that samples from coinfections in which the concentration of two viruses is similar, can have mutations from both lineages on the consensus confounding ever more the hunt for recombinants.   To minimize the number of false positive recombinants, we decided to include a detector of coinfections alongside the recombinant finder.   
The packed tools to find recombinants and coinfections can be used alone or in combination.  

## Instructions

1. Clone the Precfinder repo 

<code>git clone https://github.com/garcia-nacho/Precfinder</code>

2. Install it  

<code>cd Precfinder && ./Install.sh</code>

**Recombinant finder:**

3.A. Copy your multifasta file inside the Precfinder/Inference folder and run the Inference script inside Precfinder's folder *ModelName*   

<code>./Inference.sh ModelName</code>   

**Coinfection finder:**   

3.B. Copy your bam files inside the Precfinder/Coinfections folder. The coinfection finder currently supports two technologies, Illumina and Nanopore, include the flag Illumina or Nanopore accordingly by running:   

<code>./Coinfections.sh Illumina</code>   

or   
 
<code>./Coinfections.sh Nanopore</code>      

   
## Output recombinant finder

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

## Interpreting the recombinant plots

The pdfs contains plots similar to this one: 

<img src="/XERecombinant.png" width="300">   
The position on the genome is represented on the X axis and the conditional probability *P(Lineage|Mutation)* is represented on the Y axis (i.e. probability the virus belongs to a certain lineage given a certain mutation). Each point represents a mutation and the colours represent the most likely lineage.   
As you can see, there are mutations which are *"very exclusive"* (i.e. those with very high conditional probability), while others are more *"promiscuous"* (i.e. those with low predictive power / conditional probability). Based on the mutations plot, we can infer that the first part of the genome belongs to a *B.1* (or *B.1.1* more precisely), while the second part of the genome (after position 12K belongs to a *B.2* lineage).  
The E=0.369 is the entropy ([as in information theory entropy](https://en.wikipedia.org/wiki/Entropy_(information_theory))) of the mutations with a conditional probability larger than 0.8. As a rule of thumb, the larger the more likely to be a recombinant.   

This *BA.1-BA.2* recombinant has been already circulating around the world and the new *pango-lineage* *XE* has been proposed for it.  
   
>*XE is a BA.1/BA.2 recombinant, with the majority of the genome including the S gene belonging
to BA.2. XE shows evidence of community transmission within England, although it is currently
<1% of total sequenced cases. Early growth rates for XE were not significantly different from
BA.2, but using the most recent data up to 16 March 2022, XE has a growth rate 9.8% above
that of BA.2. As this estimate has not remained consistent as new data have been added, it
cannot yet be interpreted as an estimate of growth advantage for the recombinant. Numbers
were too small for the XE recombinant to be analysed by region.*   
[SARS-CoV-2 variants of concern and variants under investigation in England (#39/25.03.2022)](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1063424/Tech-Briefing-39-25March2022_FINAL.pdf)

## Current models for recombinant identification   
We distribute Precfinder with a set of pretrained models:   
**Norway10K**  (Trained using the last 10K Norwegian sequences deposited in [GISAID](https://www.gisaid.org/))    
**Norway22K**  (Trained using the last 20K Norwegian sequences deposited in [GISAID](https://www.gisaid.org/))   
**NorwayMarch2022**  (Trained using all the Norwegian sequences deposited in [GISAID](https://www.gisaid.org/) during March 2022)   
**Europe7K**  (Trained using the last 7K European sequences deposited in [GISAID](https://www.gisaid.org/))   
**Europe15K**  (Trained using the last 15K European sequences deposited in [GISAID](https://www.gisaid.org/))   
**Europe40K**  (Trained using the last 15K European sequences deposited in [GISAID](https://www.gisaid.org/))   

We expect to increase the number of avaible models and to update them regularly. 

## Updating Precfinder
To update precfinder you just need to run
<code>./Update.sh</code>   

## Training your own recombinant models
Precfinder offers the possibility of training your own custom models.   
To train your own model, you need to save a multifasta file on the training folder and run 
<code>./Training.sh MyOwnModel</code>   
A new model named *MyOwnModel* will be trained, saved on the Models folder and transferred to the precfinder docker image.     


## Output coinfection finder   
Precfinder will save a pdf plot containing the calculated noise across the SARS-CoV-2 genome. Samples in which two (or more) distinct viruses are present will have noise on the positions where the two viruses are different:
Example of a regular sample:   

Example of a coinfection:   

Note that Precfinder can't distinguish between in-vivo and ex-vivo contaminations.



## Advanced instructions       
I tried to write precfinder as modular as possible, so it is possible to create you own ways to use precfinder.   
Precfinder contains six core functions that can be used to create you own custom scripts, they are stored in the *Functions.R* file.   
### Core functions for recombinant finding       

***nc.pg.run(collapse.lineage="AY", clean.bad=TRUE, freq.co=0.01, cov.co=0.95, mode="Training")***  
*nc.pg.run()* runs nextclade and pangolin for a multifasta file.   
*collapse.lineage* allows you to collapse all sublineages under the same root lineage (e.g. collapse.lineage="AY" will collapse all AY. pangolin lineages under the AY.NN class) It is possible to pass several lineages to collapse (e.g collapse.lineage="AY,BA").   
clean.bad=TRUE will delete the low-quality sequences during training.
*freq.co* is the minimun frequency of a lineage (or collapsed lineage) on the training set to be included in the model
*cov.co* is the minimun coverage of a sequence to be used during training.   
mode is the mode in which the model is used. It can be set to *"Training"* or *"Inference"*, if the model is set to *"Training"*, it will expect a fasta file on the training folder and the results will be saved there, same for *"Inference"*.    

***table.generator(df="/Training/Training_dataset.csv", cores=10)***   
It will generate a conditional probability table using the pangolin and nextclade results produced by nc.pg.run().
*df* states the path to a csv file   
*cores* defines the number of CPU threads used during the table generation.   
Note that table.generator expects a csv file with at least two columns, "substitutions" (on the A123C format) and "lineage". That means that as long as you have a csv file containing "substitutions" and "lineages", you can generate a probability table for any other organism and/or to use any other classification method.   

***P.calculator(input.data, mutation.table)***   
This function generates a list containing plots and statistics for each sample. It expects input.data and mutation tables to be data frames.   

***trainset.prep(data, sample.mult=5, max.number=150000)***  
It generates the training set that will be used to train the model.
*data* are the results from nc.pg.run() on a data frame format (i.e. *data<-read.csv("/Training/Training_dataset.csv")*)   
*sample.mult* and *max.number* define the number of artificial recombinant sequences that the function will generate to train the model. *sample.mult* sets how many recombinant sequences will be generated for each non-recombinant found on the dataset and *max.number* defines the maximum number of recombinant sequences that will be generated, if *max.number* is larger than the number of non-recombinant by *sample.mult* max.number will prevail.  

***cnn.training(training.set, modelid)***   
It trains and saves the classification models. The function will also save a confusion matrix using a validation set extracted from the training and a  training progress plot.   
*training.set* refers to the traininset generated by the *trainset.prep()* function.   
*modelid* is the name of the new model.  

***ml.inference(inference.set, model, model.id)***   
It predicts the presence of recombinants using the *model* named *model.id*   

### Contaminant finder:   
The identification of contaminants relies on a tool called *FINex* that scans the entire genome looking for noise. It generates two fasta sequences, one for the virus with higher number of reads (major) and one for the virus with smaller number of reads (minor) and it runs pangolin classification on both.   Stay tuned for FINex's own github repo.    


## Under the hood.

### Conditional probability calculation    
The conditional probability for every mutation *P(Lineage|Mutation)* is calculated using Bayes' rule over the training set, and the most likely lineage for each mutation is assigned alongside the conditional probability. The mutations are stored in the *MutationTable.csv* where the colums represent the mutation, the total count, the total count for each lineage, the frequency of that mutation on the different lineages, the most likely lineage and the conditional probability.

### Sample classification
To classify the different samples as recombinant or non-reconbinant, the training process generates a set of artificial recombinants containg one or two breakpoints from two distinct lineges. The training set containing the artificial and real sequences is then parsed into a 3D array in which the first dimension is the sample dimension, the second dimension is where the different mutations are stored and the third dimensions works a a channel dimension that contains the conditional probability for the different lineages. An extra step of 3X agumentation is preformed to increase the number of samples by shuffling the channel dimension.   
Next, a 1D convolutional neural network is trained using 90% the training data, the remaining 10% is used to valiate the model and stop the training before overfitting.   
A confusion-matrix and a plot ilustrating the training process for each model is stored alongside the .hd5f files. 
All the models that we provide have a balaced accuracy >90% being the most accurate one the Europe40K with a balanced accuracy of 97.7%.   

### Understanding the limitations of Precfinder
The main limitation of Precfinder is that it does not know what it has not seen (e.g. Using a Norwegian model to classify Spanish samples could result on incorrect conditional probabilities on the Spanish samples because of Spanish specific sublineages). It is, therefore, important to keep your models updated so that they capture the current status of the pandemic in the context of the samples that you are analyzing.    
Trying to identify recombiations between sublineages with a small amount of unique mutations (e.g BA.1.1 vs BA.1) could result on incorrect predictions.   
Low quality samples (e.g. samples in which big chunks of the genome are missing) could lead to aberrant patterns that can't be understood by the 1-D CNN model.   
Co-infections/contaminants. Highly mixed samples (those in which the DNA of the different viruses have similar yield) can be understood by Precfinder as recombinants since they might contain a mixture of mutations from the two parental samples. Be sure that you run the contaminant detection script on your samples before concluding that it is a recombinant.
Contaminations. Be aware that Precfinder can't distingue between contaminants and coinfections.  
New recombinant lineages. We expect that the recombinants will be included as new independent lineages (e.g. "XE" lineage is classified as "None" by 4th of April 2022). This might confound the first models after the inclusion of those recombinant lineages in the training set. However, after certain time and the adquisition of new XE specific mutations, Precfinder will be able to find recombinants again (including recombinants between the recombinant lineages).     
