#' get_graph_limits
#' 
#' A helper function for InSpectoR. To get plotting limits from a list of spectral 
#' data files. Not exported. 
#'
#' @param x a list of spectral data files. Usually XData or XData_p
#'
#' @return
#' A list with two elements
#' \itemize{
#' \item xlimits: a vector c(xmin,xmax)
#' \item ylimits: a vector c(ymin,ymax)}
#' 
#' @export
#' 
#' @examples 
#' dfile <- system.file("insect_powder","Y_insect.txt",package="InSpectoR")
#' InSpectoR(dfile,0,1600,1024)
#' #**************************************************
#' #GUI may wait for your answer for creating a backup!
#' #**************************************************
#' leslimites <- get_graph_limits(XData_p)
#' #******************************************
#' #Don't forget that an active GUI is running.
#' #******************************************
#' 
get_graph_limits <- function(x)
# To get plotting limits from a list of spectral
# data files (x) for use with InSpectoR. 
# This returns a list with 2 components:
#   1. xlimits: a vector c(xmin,xmax)
#   2. ylimits: a vector (c(ymin,ymax))
#   
#*********************************************************************************************************** 
# B. Panneton, June 2017
#***********************************************************************************************************   
{
  xmin=min(unlist(lapply(x,function(x) min(x[1,]))))
  xmax=max(unlist(lapply(x,function(x) max(x[1,]))))
  xmax=xmax+(xmax-xmin)*0.15  #adds 15% to provide space for legends
  ymin=min(unlist(lapply(x,function(x) min(x[-1,]))))
  ymax=max(unlist(lapply(x,function(x) max(x[-1,]))))
  return(list(xlimits=c(xmin,xmax),ylimits=c(ymin,ymax)))
}  