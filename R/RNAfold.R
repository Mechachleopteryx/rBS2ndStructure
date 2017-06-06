#' @title Parallel computation of RNAfold
#'
#' @description \code{RNAfold} is used to submit the RNAfold command in ViennaRNA package into multiple bach jobs.
#'  This function is only useful when you have the qsub command and the ViennaRNA package installed on your linux system. 
#'  
#'  Other wise, please run RNAfold individually on your terminal based on the fasta files extracted by 
#'  \code{\link{tx_seq_extraction}}. 
#'  You can use \code{\link{Single_RNAfold}} to get the code of your RNAfold command.
#'  
#' @param Group_nums This value corresponds to the \code{divide_num} argument of the \code{\link{tx_seq_extraction}}.
#'
#' @param size The group name for the transcripts. Can be one of the "Small", "Middle", and "Large".
#'
#' @param OutPutNames The output names of the RNAfold.
#' 
#' @seealso \code{\link{tx_seq_extraction}}, \code{\link{Single_RNAfold}}, and \code{\link{rfold_assembly_tx}}
#' 
#' @export
#' 
#' @examples 
#' \dontrun{ 
#' tx_seq_extraction(MMusculus,txdb,getwd()) 
#' RNAfold(1:30,"Small","gsmall") 
#' }
RNAfold <-
function(Group_nums = 1:30,size = "Small",OutPutNames="Foo") {
  
  setwd(size)  
  
    qsub_submit <-
      function(command = "./test.sh", output_name = "test.out", q_name = "test", bigMem = F) {
        qsub_command <- paste0('echo "\n', command, '" | qsub -cwd',ifelse(bigMem," -q bigmem.q -l bigmem",""),' -m abe -o ',output_name,' -N ',q_name)
        record_file = paste0(q_name,"_qsub.sh")
        write(paste("#!/bin/bash",qsub_command,sep = "\n"),record_file)
        system(paste0("chmod a+rwx ",record_file))
        system(qsub_command)
      }
    
    N = length(Group_nums)
    
    Commands = vector("character",N)
    
    for (i in 1:N){
      Commands[i] = Single_RNAfold(Group_nums[i],size,OutPutNames)
    }
    
    qsub_M = data.frame(Commands = Commands,
                        output_name = paste0("Rfd_",size,"_",Group_nums,".out"),
                        q_name = paste0("Rfd_",Group_nums,"_",size))
    
    apply(qsub_M,1,function(x) qsub_submit(x[1],x[2],x[3]))
}
