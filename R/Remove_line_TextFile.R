#' Remove_line_TextFile
#' 
#' Function to remove lines from a text file.
#'
#' @param oldfile text file to modify by removing lines
#' @param line_2_del a vector of line numbers to remove
#' @param newfile text file. This is \strong{oldfile} after line removal
#'
#' @return NOTHING
#' @export
#'
#' 
Remove_line_TextFile <- function(oldfile,line_2_del,newfile=oldfile)
# ********************************************************************  
# Function to remove lines from a text file. 
# Inputs:
#   oldfile    : file to modify
#   newfile    : new file created with modified content 
#                (default to oldfile for overwrite)
#   line_2_del : a vector of line numbers
# ********************************************************************
# B. Panneton, Agriculture and Agri-Food Canada
# July 2017
# bernard.panneton@agr.gc.ca
# ********************************************************************
{
  dum=readLines(oldfile)
  dum=dum[-line_2_del]
  writeLines(dum,newfile)
}