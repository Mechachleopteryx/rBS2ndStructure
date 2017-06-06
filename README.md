
<!-- README.md is generated from README.Rmd. Please edit that file -->
Instructions of the usage of rBS2ndStructure
--------------------------------------------

This package provides thermal stable RNA 2ndary structures predicted on full length transcripts mm10 and hg19. The RNA secondary structures are used to filter the bisulfite sequening sites that potentially have incomplete bisulfite conversion on thermal stable RNA secondary structures.

-   First, install the package in R

``` r
devtools::install_github("ZhenWei10/rBS2ndStructure")
library(rBS2ndStructure)
```

-   Then, you may check and use the following data in the package.

``` r
library(GenomicRanges)

?rBS2ndStructure::Struc_hg19

?rBS2ndStructure::Struc_mm10
```

-   You could directly use the provided GRanges object to filter your sites.

``` r
rBS_gr = GRanges(seqnames = rBS_df$`#SeqID`, strand = rBS_df$refStrand, ranges = IRanges(start = rBS_df$refPos, width = 1))
rBS_gr_filtered <- rBS_gr[!rBS_gr %over% Struc_hg19]
```

-   Alternatively, you could try to filter the sites overlapped with RNA secondary structure with the provided function, this function also helps you mask all the sites on introns and intergenic sequences.

``` r
rBS_gr = GRanges(seqnames = rBS_df$`#SeqID`, strand = rBS_df$refStrand, ranges = IRanges(start = rBS_df$refPos, width = 1))
rBS_gr$mcols = p.adjust(rBS_df$`p-value_mState`,method = "BH") < .05
rBS_gr_filtered <- rBS_2ndStructure_Filter(rBS_gr,"hg19")
```

-   Finally, you could try to predict RNA 2ndary structures of your own transcriptome with RNAfold of viennaRNA package installed on your system.

``` r
#Extract transcript sequences in fasta files from BSgenome and txdb
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
tx_seq_extraction(MMusculus,txdb,getwd()) 

#Run RNAfold in multiple computational clusters 

RNAfold(1:30,"Small","Fsm") 
RNAfold(1:30,"Middle","Fmd") 
RNAfold(1:30,"Large","Flg") 

# Assemble the predicted transcripts

rfold_assembly_tx(txdb,RfdNames = c("Fsm","Fmd","Flg"))
```

-   Notice that you could run RNAfold on command line by yourself, the function RNAfold only works if you implements `qsub` command on your computer system.
