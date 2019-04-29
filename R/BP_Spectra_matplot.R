#' BP_Spectra_matplot
#'
#' BP_matplot is matplot with parameters tailored for use in InSpectoR. Called
#' by InSpectoR. Remember that data for plotting x-axis (e.g. wavelengths) are
#' in the first row of the matrices
#'
#' @param x list of matrices containing spectra. Each matrix is a different type
#'   of spectra.
#' @param yID name of ASCII data files containing spectral matrices. A different
#'   line type will be use to identify data source in plots.
#' @param col_code vector of color IDs. Usually, this will be a column (a
#'   factor) in the Ys_df of the data set.
#' @param subs vector of row number in Y file of selected samples.
#' @param xlimits vector of xmin and xmax
#' @param ylimits vector of ymin and ymax
#' @param preT TRUE for pretreated data. FALSE otherwise.
#'
#' @return 
#' A plot on the active graphic device
#' flags    : list(xdefault, ydefault) indicating caller that x and y axis limits were set to default (TRUE). 
#' @export
#' 
BP_Spectra_matplot<-function(x,yID,col_code,subs,xlimits,ylimits,preT=FALSE){
#***********************************************************************************************************  
#BP_matplot is matplot with parameters tailored for use in InSpecteuR.
# INPUTS:
#     x        : list of matrices containing spectra. Each matrix is a different type of spectra. Remember
#                that data for plotting x-axis (e.g. wavelengths) are in the first row of the matrices  
#     yID      : name of ASCII data files containing spectral matrices. A different line type will be use 
#                to identify data source in plots.
#     col_code : vector of color IDs. Usually, this will be a column (a factor) in the Ys_df of the data set. 
#     subs     : vector of row number in Y file of selected samples.  
#     xlimits  : vector of xmin and xmax  
#     ylimits  : vector of ymin and ymax
#        preT  : TRUE for pretreated data. FALSE otherwise.
# OUTPUTS:
#     flags    : list(xdefault, ydefault) indicating caller that x and y axis limits were set to default.  
#*********************************************************************************************************** 
# B. Panneton, June 2017
#*********************************************************************************************************** 
  #Define colors and appropriate axis limits
  mescols=c("darkred","blue","green3","salmon","yellow3","black","red3","magenta","gray70","cyan") 
  meslty=c("solid","88","26","2686","E6","44C4")
  plotbg="gray94"
  legbg="gray99"
  #Verify that all spectra can be plotted on same graph. If not, quit with a message
  nbwv=as.numeric(lapply(x,"ncol"))
  # if(!all(nbwv[1]==nbwv)){
  #   gmessage("Incompatible spectral data - Will clear plot.")
  #   plot.new()
  #   flags=list(xdefault=TRUE,ydefault=TRUE)
  #   return(flags)
  # }
 
  
  leslimites=list(xlimits=NULL,ylimits=NULL)
  xs_tmp<-lapply(x,function(x) x[c(1,subs+1),])
  leslimites <- get_graph_limits(xs_tmp)
  xmin=xlimits[1]
  xmax=xlimits[2]
  ymin=ylimits[1]
  ymax=ylimits[2]
  flags=list(xdefault=FALSE,ydefault=FALSE)
  if ((xmin==0 && xmax==0) | (xmin>xmax)){
    xmin=leslimites$xlimits[1]
    xmax=leslimites$xlimits[2]
    flags$xdefault=TRUE
  }
  if ((ymin==0 && ymax==0) | (ymin>ymax)){
    ymin=leslimites$ylimits[1]
    ymax=leslimites$ylimits[2]
    flags$ydefault=TRUE
  }
  #Add some extra space to the right by controlling xmax
  N_sets=length(x)
  dum=as.character(yID)
  datatypeIDs=substr(dum,start=c(1,1),regexpr("_",dum)-1)
  
  #Do a dummy plot to set up plotting area
  par(bg=plotbg)
  titre="Raw data"
  if (preT) titre="Pre-treated data"
  if (length(subs)==1){
    matplot(t(x[[1]][1,]),t(x[[1]][subs+1,]),type="n",
            xlim=c(xmin,xmax),ylim=c(ymin,ymax), bty="n",
            xlab="Wavelength or Wavenumber",ylab="Intensity (A.U.)",
            main=titre)
  }else
  {
    matplot(x[[1]][1,],t(x[[1]][subs+1,]),type="n",
            xlim=c(xmin,xmax),ylim=c(ymin,ymax), bty="n",
            xlab="Wavelength or Wavenumber",ylab="Intensity (A.U.)",
            main=titre)
  }

  #Some decoration
  box(lwd=1,lty=1,col="gray35")
  grid(lwd=3,lty=1,col="white")
  
  
  #Plot all curves
  for (k in 1:length(yID)){
    lt= k %% 6
    if (lt==0) lt=6
    if (length(subs)==1){
      matlines(x[[k]][1,],x[[k]][subs+1,],type="l",
               col=mescols,lwd=2,lty=meslty[lt])
    }else
    {
      matlines(x[[k]][1,],t(x[[k]][subs+1,]),type="l",
               col=mescols,lwd=2,lty=meslty[lt])
    }
  }
  
  #Plot legends
  l1=legend("topright", inset=0.01, legend=datatypeIDs, 
         bty='o', box.col = legbg, bg=legbg,
         lty=meslty[c(1:N_sets)],lwd=2,seg.len=4,
         cex=0.8,pt.cex=1,
         title=expression(bold("SOURCE")))

  legend(l1$rect$left+l1$rect$w,l1$rect$top-(l1$rect$h*1.05), inset=0.02, xjust=1,
         legend=col_code[subs],
         bty='o', box.col = legbg, bg=legbg,
         lty=1, col=mescols, lwd=2, seg.len=4,
         cex=0.8,pt.cex=1,
         title=expression(bold("SAMPLE")))
  
  return(flags)
}