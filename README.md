# Epigenetic Instability Index (EII)

### Protocol from publicly available data to MHIC/VarCGI:
---------------------------------------------------------

#### From Illumina 450k or EPIC array data: 
-------------------------------------------

[Click here for array manifest details](https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/documentation.html)

1)	Process data to obtain a beta value matrix (rows equal CG sites, columns equal samples). You can use any currently available program within the literature to do this, popular programs include minfi, RnBeads, ChAMP. 
2)	If your beta matrix contains multiple cancer types and their subsequent healthy samples, parse your data by cancer type and by disease state of the sample (normal, tumor). 
3)	Beta value matrix table can be further transformed to remove non-cancerous cell type methylation patterns using various deconvolution approaches. 
4)	Now we can calculate our measures of epigenetic instability (MHIC/VarCGI). MHIC stands for Methylation Heterogeneity in CpG sites. It calculates the proportion of CG sites across the genome in a sample which register a beta value of between 0.2 – 0.8. Since cancer samples are characterized as a gain of methylation across time, normal samples should have a lower proportion of CG sites with a beta value of 0.2 – 0.8 than tumor samples. We can pick any beta value range greater than 0.2 to any value less than 1; such as 0.3 to 0.7. The idea is to get a set of genomic regions/probes/CpG sites that are subject to alteration in methylation in any given tumor.  MHIC is calculated within the sample and produces a sample specific score of epigenetic instability. Within sample methylation analyses are not common and none are done in this manner which is what makes MHIC unique. MHIC is represented in the figure below.

![Figure1](https://github.com/sjthursby/EII/blob/main/images/ReadMeFigure1.png)

Here is an example of how to calculate MHIC for an array in R:

``` mhic <- function(meltedDataNormal, meltedDataTumor, Cgid){
  
  meltedDataNormal %>%
    filter(cgid %in% unique(Cgid$cgid)) %>%
    group_by(L1, variable) %>%
    mutate(total=n()) %>%
    filter(value >=0.2 & value <= 0.8) %>%
    mutate(interm=n(), percentinterm=interm/total) %>%
    summarise(percent=unique(percentinterm)) -> mhicNormal
  
  meltedDataTumor %>% 
    filter(cgid %in% unique(Cgid$cgid)) %>%
    group_by(L1, variable) %>%
    mutate(total=n()) %>%
    filter(value >=0.2 & value <= 0.8) %>%
    mutate(interm=n(), percentinterm=interm/total) %>%
    summarise(percent=unique(percentinterm)) -> mhicTumor
  
  rbind(mhicNormal, mhicTumor) -> mhicAll
  
  return(mhicAll)
}
```

5)	VarCGI stands for variance in methylation at regulatory regions, such as CpG Islands (CGI) and is another within sample metric of epigenetic instability. In the human genome, the majority of CpG Islands are unmethylated, it is generally only in disease states such as cancer that we see gains in methylation within the CGI itself. If all sites in the CGI are unmethylated the variance in that CGI’s methylation will be 0. However, if CG sites within the CGI start to gain methylation, the variance of that CGI will increase. Minor alterations to CGI methylation are acceptable and happen during aging, so to correct for this, VarCGI measures the variance in methylation across each CGI in the human genome that we have data for in an acceptable range. CGI vary in length throughout the human genome so in order to characterize the CGI effectively, a criterion of 4 or more CG sites reporting methylation is imposed for inclusion of that CGI into the VarCGI calculation. Once the variance for each acceptable CGI has been computed per sample, we form this into a sample specific score by taking the root mean variance of those individual variance scores. VarCGI is characterized in the figure below - it is novel because there are no other current methylation metrics that use the variability of the methylation within the sample itself to assess that sample’s relation to cancerous state. Most measures compare absolute methylation values in CG sites across all cancerous samples in comparison to all healthy samples and try to work out CG sites that exhibit different methylation status from normal (typically called differential methylation analyses). Our methodology seeks to utilize the variability in the epigenetic state of the sample to our benefit instead of it being something that you have to correct for, as it is in common practice. We denote the latter approach as metrics for epigenetic perturbation or epigenetic instability.


![Figure2](https://github.com/sjthursby/EII/blob/main/images/ReadMeFigure2.png)

Here is an example of how to calculate VarCGI within arrays in R:

```varCGI <- function(meltedDataNormal, meltedDataTumor, Cgid){
  
  meltedDataNormal %>%
    filter(cgid %in% unique(Cgid$cgid)) %>%
    group_by(L1, variable, Arbitrary_CpGI_Name) %>%
    summarise(varianceCGI=var(value, na.rm=T), numProbes=n()) -> a
  
  a <- na.omit(a)
  
  a <- a[a$numProbes >= 10,]
  
  a %>%
    group_by(L1, variable) %>%
    summarise(rmVar=sqrt(sum(varianceCGI/n()))) -> normalRmVar
  
  meltedDataTumor %>%
    filter(cgid %in% unique(Cgid$cgid)) %>%
    group_by(L1, variable, Arbitrary_CpGI_Name) %>%
    summarise(varianceCGI=var(value, na.rm=T), numProbes=n()) -> a
  
  a <- na.omit(a)
  
  a <- a[a$numProbes >= 10,]
  
  a %>%
    group_by(L1, variable) %>%
    summarise(rmVar=sqrt(sum(varianceCGI/n()))) -> tumorRmVar
  
  rbind(tumorRmVar, normalRmVar) -> result
  
  return(result)
}
```

#### From reduced representation bisulfite sequencing (RRBS) or from whole genome bisulfite sequencing (WGBS) on tissues
-------------------------------------------------------------------------------------------------

MHIC and VarCGI work very similar to how they work in array technologies. The only differences are that in RRBS or WGBS there is a higher genome-wide coverage of CG sites that we can take advantage of and there is no beta matrix that hails from RRBS or WGBS analysis. The typical representation of methylation in this case, is a BED file (chromosome, start site, end site, methylation as a beta value) all BED files are read into R (any recent version is acceptable, we are currently using R 4.3.3) and a long format of a beta matrix is created. 

#### From reduced representation bisulfite sequencing (RRBS) or from whole genome bisulfite sequencing (WGBS) on cell-free DNA (cfDNA)
-------------------------------------------------------------------------------------------------

While it is possible to calculate MHIC and VarCGI on cfDNA, and we observe that these metrics yield differences between different stages of cancers and the normal samples, these metrics can be challenging to use on the small amount of sequencing reads derived from circulating tumor DNA (ctDNA) from the cfDNA. Therefore, we developed an alternative metric of EII from cfDNA  to try and extract as much information as possible from the cfDNA sequences. We call this technique (currently) cfEII or cell free epigenetic instability index. 

1)	To calculate cfEII you need to extract the methylation of individual sequencing reads within CGI of interest. We do this with the RLM package. 

Typically, after processing RRBS or WGBS in Bismark (or any other methylation aligner/caller) you obtain BED files containing the chromosome, start, end and beta value representing the methylation of every CG that was analyzed in the dataset. However, that value is made up of a summary of all the cfDNA reads that aligned to that specific CG site within the genome. There could be for example, 10 reads aligning to that area, 5 reporting methylation in that CG site, and 5 reporting lack of methylation for that CG site (this is in the context of one sample). This would result in that site obtaining a beta value of 0.5 (i.e. 50% methylation), which in this context is an over-generalization of what the methylation is like in that specific genomic region; especially when circulating tumor DNA (ctDNA) only makes up ~0.01% of all cfDNA and we are trying to detect that ~0.01% of ctDNA. Additionally, splitting the beta value into its constituent reads will help us identify what tissue/organ any ctDNA came from since methylation is tissue specific and those cfDNA reads with differing methylation values could have originated from any dying cell in the human body. 

2)	Exclude any sequencing read that has less than 5 CG sites on it. 

CfDNA fragments are very short in length, up to ~160bp. Therefore, in order to properly characterize the methylation of the cfDNA fragment via cfEII, only reads with at least 5CG are accepted for further processing. 

3)	Obtain reads that map to CGIs of interest and calculate the variance across the reads that map to those CGI. These variance measures per sample represent that sample’s cfEII. The figure below represents the process of calculating cfEII.

4)	Plot these values as a Cumulative Distribution Function with comparison to the normal samples.

We show that when compared to the normal/healthy samples, even early-stage cancer samples display higher values for cfEII and reach statistical significance. 

![Figure3](https://github.com/sjthursby/EII/blob/main/images/ReadMeFigure3.png)

Data used in the paper can be found [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA494975) and [here](https://ncbi.nlm.nih.gov/bioproject/534206). The former is the breast cancer dataset and the latter is the lung cancer dataset.

