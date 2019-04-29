#' BP_Plot_Scores
#' 
#' Generates a biplot of scores from PCA with color coding of points, point labels and ellipses around data groups.
#' Returns some parameters for easy replotting when zooming for example. If ax1 and ax2 are identical, plots
#' a histogram of scores.
#'
#' @param ax1 horizontal axis value of scores 
#' @param ax2 vertical axis value of scores
#' @param pcs 2 element vector giving the PC number for ax1 and ax2 respectively for axis labelling
#' @param frac 2 element vector giving the \% variance explained for PC on ax1 and ax2 respectively. For axis labelling.
#' @param xlimits 2 element vector giving the minimum and maximum value of the horizontal axis.
#' @param ylimits 2 element vector giving the minimum and maximum value of the vertical axis.
#' @param colorby class for color coding of points. Default is all blue. 
#'                When default, no legend is plotted.
#' @param pt.label label for points. NA is used to skip a labels.
#' @param ell_option option for plotting confidence ellipse around groups of points. "None", "0.50" or "0.95"
#'                   Will be done only if there is more than one group defined by colorby.
#'
#' @return
#'  A list with the following elements:
#'  \itemize{
#'         \item Xsc,Ysc  : x and y values of all scores plotted.
#'              \item pcs  : number of x and y axis principal comp.
#'     \item X_Lim, Y_Lim  : xlim and ylim of current plot.
#'          \item colorby  : class for color coding of points.
#'          \item pt.label : label for points. NA is used to skip a labels.
#'       \item  ell_option : option for plotting confidence ellipse around groups of points. "None", "0.50" or "0.95"
#'                     Will be done only if there is more than one group.
#'                     
#' }
#' @export
#'
#' @examples
#'  BP_Plot_Scores(ax1=c(rnorm(20,10,2),rnorm(20,20,3)), ax2=rnorm(40,5,3),
#'                 pcs=c(1,2), frac= c(81.2,9.8),
#'                 colorby=as.factor(c(rep("Cl1",20),rep("Cl2",20))),
#'                 ell_option="0.5")
#'                 
BP_Plot_Scores <- function(ax1,ax2=NULL,pcs,frac,xlimits=NULL,ylimits=NULL,colorby=rep(2,length(ax1)),pt.label=NULL,ell_option)
# Biplot of PCA scores.
# Inputs
#   ax1     : scores on horizontal axis 
#   ax2     : scores on vertical axis
#   pcs     : 2 element vector giving the PC number for ax1 and ax2 respectively
#  frac     : 2 element vector giving the % variance explained for PC on ax1 and ax2 respectively.
#   colorby : class for color coding of points. Default is all blue. 
#             When default, no legend is plotted.
#   pt.label: label for points. NA is used to skip a labels. 
# ell_option: option for plotting confidence ellipse around groups of points. "None", "0.50" or "0.95"
#             Will be done only if there is more than one group.
# Output
#   a list with the following elements:
#          Xsc,Ysc  : x and y values of all scores plotted.
#              pcs  : number of x and y axis principal comp.
#     X_Lim, Y_Lim  : xlim and ylim of current plot.
#          colorby  : class for color coding of points.
#          pt.label : label for points. NA is used to skip a labels.
#        ell_option : option for plotting confidence ellipse around groups of points. "None", "0.50" or "0.95"
#                     Will be done only if there is more than one group.
#     
#*********************************************************************************************************** 
# B. Panneton, November 2018
#*********************************************************************************************************** 
{
  
  mescols=c("darkred","blue","green3","salmon","yellow3","black","red3","magenta","gray70","cyan") 
  #To define corresponding lighter transparent colors for symbol fill
  mescols_fill=col2rgb(mescols,alpha=TRUE)
  mescols_fill[4,]=145
  mescols_fill=rgb(mescols_fill[1,],mescols_fill[2,],mescols_fill[3,],alpha=mescols_fill[4,],maxColorValue = 255)
  recycling=ceiling(length(unique(colorby))/10)  #Only 10 colors defined, so we recycle if necessary
  mescols_fill=rep(mescols_fill,recycling)
  
  
  plotbg="gray94"  # main plot area background color
  legbg="gray99"   # legend background
  par(bg=plotbg)
  

  if (identical(ax1,ax2)){  #Plot histograms per group. Does not work if too many groups
    #ldahist(acp$x[,ax1],colorby,sep=TRUE,col="darkred")
    dum_df<-data.frame(scores=ax1,class=colorby)
    p <- ggplot2::ggplot(dum_df, ggplot2::aes(x=scores)) + ggplot2::geom_histogram(bins=20,colour="white",fill="blue")
    p <- p + ggplot2::facet_wrap( ~ class,ncol=1)
    print(p)
    
    
    outp <- list(Xsc=ax1, Ysc=ax2, X_Lim=NULL, Y_Lim=NULL)
    
  }else
  {
    if (is.null(xlimits)){
      xlimits=range(ax1)
      xlimits[2]=xlimits[2]+diff(xlimits)*0.25   #room for legends
      ylimits=range(ax2)
    }
    
    if (length(unique(colorby))>1 & (!(ell_option=="None"))){
      #add room for ellipse labels
      dx=diff(xlimits)
      xlimits[1]=xlimits[1]-dx*0.25
      xlimits[2]=xlimits[2]+dx*0.25
      dy=diff(ylimits)
      ylimits[1]=ylimits[1]-dy*0.25
      ylimits[2]=ylimits[2]+dy*0.25
      
    }
    plot(ax1,ax2,type="n",
         xlim=xlimits,ylim=ylimits, bty="n",
         xlab=paste("PC",pcs[1]," - ", round(frac[1],1),"%",sep=""),
         ylab=paste("PC",pcs[2]," - ", round(frac[2],1),"%",sep=""))
    
    #Some decoration
    box(lwd=1,lty=1,col="gray35")
    grid(lwd=3,lty=1,col="white")
    
    #Plot the points
    points(ax1,ax2,type="p",
         xlim=xlimits,ylim=ylimits, bty="n",
         xlab=paste("PC",pcs[1],sep=""),ylab=paste("PC",pcs[2],sep=""),
         pch=21,col="black",bg=mescols_fill[colorby],cex=1.5,lwd=1)
    #mescols[colorby]
    #Add label to points
    text(ax1,ax2,pt.label,pos=4,cex=0.8)
    if (length(unique(colorby))>1 & (!(ell_option=="None"))){
        car::dataEllipse(ax1,ax2,as.factor(colorby),col=mescols_fill,
                    plot.points = FALSE,
                    levels=as.numeric(ell_option))
    }
    
    #Plot the legend
    if (diff(range(as.numeric(colorby)))>0){  #points are not same color -> do legend
      l1=legend("topright", inset=0.01, legend=levels(unique(colorby)), 
                bty='o', box.col = legbg, bg=legbg,
                pch=21,pt.bg=mescols_fill,seg.len=4,
                cex=0.8,pt.cex=1.5,
                title=expression(bold("CLASS")))
    }
    outp <- list(Xsc=ax1, Ysc=ax2, pcs=pcs, frac=frac, X_Lim=xlimits, Y_Lim=ylimits, colorby=colorby, pt.label=pt.label, ell_option=ell_option)
  }
  return(outp)
}
