# Nextflow Pipeline for Estimating Antimicrobial Resistance Gene Abundance in Metagenomic Samples
This analytical framework is discussed extensively here:  

    Ackerson,Leland K.,,IV. Analytical Framework for Estimating Antimicrobial Resistance Gene Abundance in Metagenomic Samples  
        of Animal Agriculture Origin, Michigan State University, United States -- Michigan, 2023. ProQuest

## Quantification and Estimation of Antimicrobial Resistance Genes from Shotgun Metagenomic Sequence Data
### Methodology
#### 1.  Quality Control + Pre-Processing
~FastQC is performed on both the raw reads and post-processing clean reads.  
~Results from quality control analysis are deposited in the 'QCmetrics' folder upon completion of the workflow proccesses.  
~Adapter trimming and read quality filtering are performed using BBDuk (BBTools).

#### 2.  Antibotic Resistant Gene Mapping
~Reference Database: The Comprehensive Antibiotic Resistance Database (CARD).  
~Alignement Software: BWA MEM

#### 3.  16s rRNA Gene Mapping
Necessary for metagenome taxonomical quantification, and downstream normalization.  
~Reference Database: GreenGenes Database  
~Alignement Software: BWA MEM

#### 4.  Gene Quantification
~ Gene counts are weighted by number of mappings for a read of interest.  
For example, 3 gene mappings by one read results in the allocation of 1/3 of a count for each gene.



