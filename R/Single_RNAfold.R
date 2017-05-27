#' @title Generate the bash command of RNAfold
#'
#' @description RNAfold is a RNA secondary structure prediction software in ViennaRNA package. 
#' This R function is used to generate the single bash command of RNAfold, which is used to predict the MEA RNA secondary structures of the full length transcripts.
#'  
#' @param Group_num The number of the group of the fasta file.
#' This value corresponds to the \code{divide_num} argument of the \code{\link{tx_seq_extraction}}.
#' 
#' @param Size The size of the transcript. 
#' Can be one of the "Small", "Middle", and "Large".
#' @param outPutNames The name of the output files of RNA fold.
#' 
#' @return This function will return the code of RNAfold that can directly run in the terminal.
#' 
#' @seealso \code{\link{tx_seq_extraction}}, \code{\link{RNAfold}}, and \code{\link{rfold_assembly}}
#' 
#' @examples 
#' cat(Single_RNAfold(1,"Small","gsmall"))
#' 
#' \dontrun{
#' tx_seq_extraction(MMusculus,txdb,getwd()) 
#' system(Single_RNAfold(1,"Small","gsmall"))
#'}
#'
Single_RNAfold <-
function(Group_num,Size = "Small",outPutNames="Foo") {
    
    RNAfold = c(
      'RNAfold',
      paste0('--infile=','./',Size,"_d",Group_num,".fasta"), 
      paste0('--outfile=',outPutNames),
      '--maxBPspan=150', 
      '--MEA=.1',
      '--temp=70',
      '--noconv',
      '--noPS'
    )
    
    paste(RNAfold, collapse = " \\\n")
  }
