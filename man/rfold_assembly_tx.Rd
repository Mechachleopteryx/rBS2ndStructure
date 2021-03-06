% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rfold_assembly_tx.R
\name{rfold_assembly_tx}
\alias{rfold_assembly_tx}
\title{Assembly of RNA secondary structures}
\usage{
rfold_assembly_tx(TXDB, Rfold_dir = "./", RfdNames = c("Fsm", "Fmd", "Flg"),
  DirNames = c("Small", "Middle", "Large"), small_threashold = 4500,
  middle_threashold = 8000, trim_wsize = 2000, trim_ssize = 1000)
}
\arguments{
\item{TXDB}{The txdb object used to extract the transcripts.}

\item{Rfold_dir}{The directory containing 3 sub-directories: "Small", "Middle", and "Large".}

\item{RfdNames}{A character vector of length 3,
which indicates the output names of the RNAfold used in 3 groups.}

\item{DirNames}{A character vector of length 3,
which indicates the directory names of the RNAfold used in 3 groups.}

\item{small_threashold}{the transcript length threshold to put into the small group.}

\item{middle_threashold}{the transcript length threshold to put into the middle group.}

\item{trim_wsize}{the window size used for the trimming of the large group.}

\item{trim_ssize}{the step size used for the trimming of the large group.}
}
\value{
A \code{GRangesList} object, each element of it represents the MEA RNA structures in a transcript.
}
\description{
\code{rfold_assembly_tx} is used to assemble the output files of the RNAfold into a GRangesList object.
}
\examples{

\dontrun{ 
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
tx_seq_extraction(MMusculus,txdb,getwd()) 
RNAfold(1:30,"Small","Fsm") 
RNAfold(1:30,"Middle","Fmd") 
RNAfold(1:30,"Large","Flg") 
rfold_assembly_tx(txdb,RfdNames = c("Fsm","Fmd","Flg"))
}

}
\seealso{
\code{\link{tx_seq_extraction}}, \code{\link{RNAfold}}, and \code{\link{Single_RNAfold}}
}
