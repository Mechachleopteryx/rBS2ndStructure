#' @title Extract the full length transcripts sequences
#'
#' @description \code{tx_seq_extraction} is used to extract the fasta files of the full length transcripts.
#' @param BS_genome the \code{BSgenome} object of your organism.
#' @param TXDB the \code{txdb} object of your organism.
#' @param Write_PATH the path to write the transcript sequences.
#' @param small_threashold the transcript length threshold to put into the small group.
#' @param middle_threashold the transcript length threshold to put into the middle group.
#' @param trim_wsize the window size for the trimming.
#' 
#' By default, only the transcripts in the large group are trimmed.
#' @param trim_ssize the step size for the trimming.
#' @param divide_num the number of fasta files divided for each of the 3 groups. 
#' 
#' This argument is useful when conducting the parallel computation of the \code{RNAfold} step.
#' 
#' @return This function will create 3 folders named Small, Middle, and Big. 
#' Each of them will contains a number of fasta files containing the full length transcripts of the sizes of each group. 
#' 
#' The number of the fasta files is determined by the argument \code{divide_num}.
#' 
#' @seealso \code{\link{RNAfold}}, \code{\link{Single_RNAfold}}, and \code{\link{rfold_assembly_tx}}
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' tx_seq_extraction(MMusculus,txdb,getwd())
#' }
tx_seq_extraction <-
function(
  BS_genome,
  TXDB,
  Write_PATH = ".",
  small_threashold = 4500,
  middle_threashold = 8000,
  trim_wsize = 2000,
  trim_ssize = 1000,
  divide_num = 30
){
  #Renaming IDs
  exbytx <- exonsBy(TXDB, by = "tx")
  names(exbytx) = paste0("tx_",names(exbytx))
  
  Get_full_tx_sequences <- function(IDs,exBytx){
    
    exBytx_s <- exBytx[IDs]
    
    Views(BS_genome,unlist(exBytx_s)) %>% DNAStringSet %>% split(.,names(.)) %>% lapply(.,function(x) unlist(x) %>% as.character) %>% unlist %>% DNAStringSet -> TX_sequences
    
    return(TX_sequences)
  }
  
  tx_lengths = sum(width(exbytx))
  IDs_small = which(tx_lengths <= small_threashold)
  IDs_middle = which(tx_lengths > small_threashold & tx_lengths < middle_threashold)
  IDs_large = which(tx_lengths >= middle_threashold)
  
  system("mkdir ./mm10")
  setwd("./mm10")
  #Quite Slow! ~5min
  Sequences_small <- Get_full_tx_sequences(IDs_small,exbytx)
  Sequences_middle <- Get_full_tx_sequences(IDs_middle,exbytx)
  Sequences_large <- Get_full_tx_sequences(IDs_large,exbytx)
  
  Sset_cut <- function(Sset,window_size = trim_wsize,step_size = trim_ssize) {
    Sset_list <- suppressWarnings( lapply(Sset,function(x) Views(x,IRanges(start = seq(1,length(x),by = step_size),width = window_size)) %>% DNAStringSet) %>% DNAStringSetList )
    Sset_cut <- unlist(Sset_list)
    elenrows <- elementNROWS(Sset_list)
    idx <- sapply(elenrows,seq_len) %>% unlist
    names(Sset_cut) = paste(names(Sset_cut),idx,sep = ".")
    return(Sset_cut)
  }
  
  Sequences_large <- Sset_cut(Sequences_large)
  
  dividing_saver <- function(x,group_num = 30,title) {
    idx <- cut(seq_len(length(x)),group_num)
    Sequences_lst <- split(x,idx)
    
    mapply(function(x, i) {
      
      file_name = paste0(title,"_d",i,".fasta")
      writeXStringSet(x,file_name,append=F)
    }, x = Sequences_lst, i = seq_len(length(Sequences_lst))
    )
    
  }
  
setwd(Write_PATH)

dir.create("./Small")
dir.create("./Middle")
dir.create("./Large")

setwd("./Small")
dividing_saver(Sequences_small,divide_num,"Small")
setwd("../Middle")
dividing_saver(Sequences_middle,divide_num,"Middle")
setwd("../Large")
dividing_saver(Sequences_large,divide_num,"Large")
}

