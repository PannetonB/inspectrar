#' Ref_2_Lab
#' 
#' Function to convert reflectance spectra to CIELab values
#' 
#'
#' @param wS wavelenth vector
#' @param S  matrix of reflectance spectra. One spectrum per line. Number of columns should match length(wS)
#'
#' @return
#'  A matrix of dimension (nrow(S),3) where the columns are L, a and b value, one row per sample
#'     
#' @export
#'
#' 
Ref_2_Lab <- function(wS,S){
# Converts reflectance spectra to CIELab.
# Inputs:
#   wS : wavelenth vector
#    S : matrix of reflectance spectra. One spectrum per line.
#        Number of columns should match length(wS)
       
  hasIt <- require(colorSpec) 
  if (hasIt){
    library(colorSpec)
  }else
  {
    install.packages("colorSpec")
    library(colorSpec)
  }
  library(spacesXYZ)
  
  
  I <- D65.1nm    #D65 illuminant, 1 nm
  wI <- wavelength(I)
  I <- as.numeric(coredata(I))
  
  xyz <- xyz1931.1nm  #XYZ 2° observer, 1 nm
  wxyz <- wavelength(xyz)
  xyz <- as.matrix(coredata(xyz))
  
  if(is.vector(S)) S <- t(as.matrix(S))
  
  doIt <- function(S,wS,I,wI,xyz,wxyz)
  {
    rI <- range(wI)
    rxyz <- range(wxyz)
    rS <- range(wS)
    wmin <- max(c(rI[1],rS[1],rxyz[1]))
    wmax <- min(c(rI[2],rS[2],rxyz[2]))
    
    w <- seq(wmin,wmax,1)
    i1 <- which(wI==wmin); i2 <- which(wI==wmax)
    I <- I[i1:i2]
    i1 <- which(wS==wmin); i2 <- which(wS==wmax)
    St <- S[i1:i2]
    i1 <- which(wxyz==wmin); i2 <- which(wxyz==wmax)
    xyz <- xyz[i1:i2,]
    x <- xyz[,1]
    y <- xyz[,2]
    z <- xyz[,3]
    
    N <- length(x)
    X <- 1/sum(I*x)*sum(x*St*I)
    Y <- 1/sum(I*y)*sum(y*St*I)
    Z <- 1/sum(I*z)*sum(z*St*I)
    
    CieLab <- LabfromXYZ(t(rbind(X,Y,Z)), standardXYZ('D65'))
  }
  
  res <- apply(S,MARGIN = 1,doIt, wS=wS,I=I,wI=wI,xyz=xyz,wxyz=wxyz)
  res=t(res)
  colnames(res) <- c("L","a","b")
  
  return(res)
}
