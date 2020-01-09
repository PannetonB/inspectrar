#' Trouve_ErrEchID
#' 
#' Function to identify non matching items in the first columns of Y and X files.
#' 
#'
#' @param fromConsole a boolean. TRUE for running from console. FALSE for running from InSpectoR
#'
#' @return
#'   \itemize{
#'      \item a character string ("Compatible files") when files are compatible 
#'      \item when files are not compatible, a list having three elements:
#'      \itemize{
#'          \item a data frame with three columns: 1-line where non matchings appear, 
#'                                                 2-Y file value in column 1 for line in item 1 above,
#'                                                 3-X file value in column 1 for line in item 1 above.
#'          \item the x full file name
#'          \item the y full file name
#'      }
#'    }      
#'     
#' @export
#'
#' @examples
#' Trouve_ErrEchID()
#' 
Trouve_ErrEchID <- function(fromConsole=T)
{
  yfile <- choose.files(default = "Y*.txt",
                        filters = Filters["txt"],
                        caption = "Choisir un fichier Y",
                        multi=F)
  lepath <- dirname(yfile)
  xfile <- choose.files(default = file.path(lepath,"*.txt"),
                        filters = Filters["txt"],
                        caption = "Choisir un fichier de spectres",
                        multi=F)
  
  ydat <- read.delim(yfile)
  xdat <- read.delim(xfile, header=FALSE)
  if (any(as.character(xdat$V1[-1])!=as.character(ydat$ECHID))){
    indi=which(as.character(xdat$V1[-1])!=as.character(ydat$ECHID))
    dumdf=as.data.frame(cbind(indi,as.character(ydat$ECHID[indi]),as.character(xdat$V1[indi+1])))
    colnames(dumdf) <- c("Lignes", basename(yfile), basename(xfile))
    if (fromConsole)
      cat(paste0(
        "\n**********************************",
        "\nTable of lines with errors: ",
        "\n**********************************\n"))
    return(list(df=dumdf,xfile=xfile,yfile=yfile))
  }else
  {return("Compatibles files!")}
}