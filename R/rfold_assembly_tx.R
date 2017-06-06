#' @title Assembly of RNA secondary structures
#'
#' @description \code{rfold_assembly_tx} is used to assemble the output files of the RNAfold into a GRangesList object.
#'  
#' @param TXDB The txdb object used to extract the transcripts.
#' @param Rfold_dir The directory containing 3 sub-directories: "Small", "Middle", and "Large".
#' @param RfdNames A character vector of length 3,
#'  which indicates the output names of the RNAfold used in 3 groups. 
#' @param DirNames A character vector of length 3,
#' which indicates the directory names of the RNAfold used in 3 groups. 
#' @param small_threashold the transcript length threshold to put into the small group.
#' @param middle_threashold the transcript length threshold to put into the middle group.
#' @param trim_wsize the window size used for the trimming of the large group.
#' @param trim_ssize the step size used for the trimming of the large group.
#' 
#' @return A \code{GRangesList} object, each element of it represents the MEA RNA structures in a transcript.
#' 
#' @seealso \code{\link{tx_seq_extraction}}, \code{\link{RNAfold}}, and \code{\link{Single_RNAfold}}
#' @export
#'  
#' @examples 
#' 
#' \dontrun{ 
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' tx_seq_extraction(MMusculus,txdb,getwd()) 
#' RNAfold(1:30,"Small","Fsm") 
#' RNAfold(1:30,"Middle","Fmd") 
#' RNAfold(1:30,"Large","Flg") 
#' rfold_assembly_tx(txdb,RfdNames = c("Fsm","Fmd","Flg"))
#' }
#'
rfold_assembly_tx <-
function(TXDB,
           Rfold_dir="./",
           RfdNames = c("Fsm","Fmd","Flg"),
           DirNames = c("Small","Middle","Large"),
           small_threashold = 4500,
           middle_threashold = 8000,
           trim_wsize = 2000,
           trim_ssize = 1000) {
    
    rfold_assembly <-
      function(IDs_read,title,exbytx) {
        Structured_grl <- GRangesList()
        missing_list <- c()
        for(i in IDs_read) {
          file_name = paste0("./",title,"_",title,"_",i,".fold")
          if (file.exists(file_name)) {
            Structure_MEA <- BStringSet(gsub(" .*","",readLines(file_name)[6]))
            Structured_ir <- gaps(reduce(vmatchPattern(".",Structure_MEA)[[1]]))
            if (length(Structured_ir) > 0) {
              Structured_gr <- mapFromTranscripts(GRanges(seqnames = i,
                                                          strand = strand(exbytx[[i]])[1],
                                                          ranges = Structured_ir), 
                                                  exbytx[i], ignore.strand = FALSE)
            }else{
              Structured_gr <- GRanges()
            }
            Structured_grl[[as.character(i)]] = Structured_gr 
          }else{
            missing_list[as.character(i)] = i
          }
        }
        
        write(missing_list,paste0("../",title,"_missing.txt"))
        return(Structured_grl)
      }
    
    rfold_trimmed_assembly <-
      function(dir = "./",title,window_size = 2000,step_size = 1000,exbytx) {
        
        single_trimmed_assembly <-
          function(ID,NUM_cut,title,exbytx,wsize,ssize) {
            Structured_ir <- IRanges()
            for(i in 1:NUM_cut) {
              file_name = paste0("./",title,"_",title,"_",ID,".",i,".fold")
              Structure_MEA <- BStringSet(gsub(" .*","",readLines(file_name)[6]))
              Structured_ir_i <- gaps(reduce(vmatchPattern(".",Structure_MEA)[[1]]))
              Structured_ir = c(Structured_ir,shift(Structured_ir_i,ssize*(i-1)))
            }
            Structured_ir <- reduce(Structured_ir)
            if (length(Structured_ir) > 0) {
              Structured_gr <- mapFromTranscripts(GRanges(seqnames = ID,
                                                          strand = strand(exbytx[[ID]])[1],
                                                          ranges = Structured_ir), 
                                                  exbytx[ID], ignore.strand = FALSE)
            }else{
              Structured_gr <- GRanges()
            }
            return(Structured_gr)
          }
        
        setwd(dir) 
        idx <- grep(".fold$",list.files(dir),value=T) %>% gsub(paste0(title,"_"),"",.)  %>% gsub(".fold","",.) %>% gsub("tx_","",.)
        tx_idx <- gsub("\\.[0-9]*$","",idx) %>% as.numeric
        splice_idx <- gsub("^[0-9]*\\.","",idx) %>% as.numeric
        NUM_each_tx <- tapply(splice_idx,tx_idx,max)
        TX_id <- paste0("tx_",names(NUM_each_tx))
        Structured_grl <- GRangesList()
        
        for(i in 1:length(TX_id)) {
          Structured_grl[[as.character(TX_id[i])]] <- single_trimmed_assembly(TX_id[i],NUM_each_tx[i],title,window_size,step_size,exbytx)
        }
        
        return(Structured_grl)
      }
    
    exbytx = exonsBy(TXDB,by = "tx")
    names(exbytx) = paste0("tx_",names(exbytx))
    tx_lengths = sum(width(exbytx))
    IDs_small = names(exbytx)[which(tx_lengths <= small_threashold)]
    IDs_middle = names(exbytx)[which(tx_lengths > small_threashold & tx_lengths < middle_threashold)]
    
    setwd(paste0(Rfold_dir,"/",DirNames[1]))
    cat("Assembling small transcripts...\n")
    grl_small <- rfold_assembly(IDs_small,RfdNames[1],exbytx)
    setwd(paste0(Rfold_dir,"/",DirNames[2]))
    cat("Assembling middle transcripts...\n")
    grl_middle <- rfold_assembly(IDs_middle,RfdNames[2],exbytx)
    setwd(paste0(Rfold_dir,"/",DirNames[3]))
    cat("Assembling large transcripts...\n")
    grl_Large <- rfold_trimmed_assembly(dir = "./",title = RfdNames[3],exbytx = exbytx,window_size = trim_wsize, step_size = trim_ssize)
    
    
    cat("Combining 3 groups of transcripts...\n")
    gr_small <- unlist(grl_small)
    gr_middle <- unlist(grl_middle)
    gr_large <-  unlist(grl_large)
    
    grl_All <- split(c(gr_small,gr_middle,gr_large),c(names(gr_small),names(gr_middle),names(gr_large)))
   return(grl_All)
 }
