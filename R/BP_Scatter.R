#' BP_Scatter
#' 
#' Produces a scatter plot and return parameters for easy replotting when zooming for example.
#'
#' @param ax1 vector of values on the horizontal axis
#' @param ax2 vector of values on the vertical axis
#' @param labs 2 element list giving x and y axis labels respectively
#' @param titre main figure title
#' @param xlimits minimum and maximum values for the horizontal axis
#' @param ylimits minimum and maximum values for the vertical axis
#' @param colorby class for color coding of points. Default is all blue. 
#'                     When default, no legend is plotted.
#' @param pt.label label for points. NA is used to skip a labels.
#'
#' @return
#'  A list with the following elements:
#'  \itemize{
#'          \item Xsc,Ysc  : x and y values of all scores plotted.
#'              \item  labs : 2 element list giving x and y axis labels respectively
#'         \item  titre     : main figure title
#'     \item  X_Lim, Y_Lim  : xlim and ylim of current plot.
#'          \item  colorby  : class for color coding of points.
#'          \item  pt.label : label for points. NA is used to skip a labels.
#'          }
#' @export
#'
#' @examples
#' BP_Scatter(ax1=c(rnorm(20,10,2),rnorm(20,20,3)), ax2=rnorm(40,5,3),
#'                 labs=list("X-axis","Y-axis"),
#'                 titre = "A scatter plot",
#'                 xlimits=c(5,30), ylimits=c(2,8),
#'                 colorby=as.factor(c(rep("Cl1",20),rep("Cl2",20))))
#'                 
BP_Scatter <- function(ax1,ax2,labs=list("X","Y"),titre="",xlimits=NULL,ylimits=NULL,colorby=rep(2,length(ax1)),pt.label=NULL)
# Biplot of PCA scores.
# Inputs
#   ax1     : scores on horizontal axis 
#   ax2     : scores on vertical axis
#  labs     : 2 element list giving x and y axis labels respectively
# titre     : main figure title
#   colorby : class for color coding of points. Default is all blue. 
#             When default, no legend is plotted.
#   pt.label: label for points. NA is used to skip a labels.
#  
# Output
#   a list with the following elements:
#          Xsc,Ysc  : x and y values of all scores plotted.
#              labs : 2 element list giving x and y axis labels respectively
#         titre     : main figure title
#     X_Lim, Y_Lim  : xlim and ylim of current plot.
#          colorby  : class for color coding of points.
#          pt.label : label for points. NA is used to skip a labels.
#     
#*********************************************************************************************************** 
# B. Panneton, January 2019
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
  

  
  if (is.null(xlimits)){
    xlimits=range(ax1)
    xlimits[2]=xlimits[2]+diff(xlimits)*0.25   #room for legends
    ylimits=range(ax2)
  }
  
  plot(ax1,ax2,type="n",
       xlim=xlimits,ylim=ylimits, bty="n",
       xlab=labs[[1]],
       ylab=labs[[2]])
  
  #Some decoration
  box(lwd=1,lty=1,col="gray35")
  grid(lwd=3,lty=1,col="white")
  
  #Plot the points
  points(ax1,ax2,type="p",
       xlim=xlimits,ylim=ylimits, bty="n",
       xlab=labs[[1]],ylab=labs[[2]],
       pch=21,col="black",bg=mescols_fill[colorby],cex=1.5,lwd=1)
  #mescols[colorby]
  #Add label to points
  text(ax1,ax2,pt.label,pos=4,cex=0.8)
  
  
  #Plot the legend
  if (diff(range(as.numeric(colorby)))>0){  #points are not same color -> do legend
    l1=legend("topright", inset=0.01, legend=levels(unique(colorby)), 
              bty='o', box.col = legbg, bg=legbg,
              pch=21,pt.bg=mescols_fill,seg.len=4,
              cex=0.8,pt.cex=1.5,
              title=expression(bold("CLASS")))
  }
  title(main = titre)
  
  outp <- list(Xsc=ax1, Ysc=ax2, labs=labs, titre=titre, X_Lim=xlimits,
               Y_Lim=ylimits, colorby=colorby, pt.label=pt.label)
}
