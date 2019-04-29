#' Baseline_Correction
#' 
#' A simple GUI to remove a baseline from spectra in a file.
#' 
#' The file format must conform to the one required by InSpectoR.
#' See \link{InSpectoR} for details. The method to calculate the
#' baseline is from the \link{baseline} package, namely the baseline
#' function with \emph{method='modpolyfit'}. A GUI helps the user
#' to select parameters to perform the fit.
#'
#' @param lefichier a filename to work on.
#'
#' @return
#'    Baseline corrected spectra and baseline spectra are stored in files
#'    with the same name as the input file but with "Corrected" or "Baseline"
#'    keyword added.
#' @export
#'
#' @examples
#' dfile <- system.file("raman","Raman_I.txt",package="inspectrar")
#' Baseline_Correction(dfile)
Baseline_Correction <- function(lefichier=NULL)
{
  #library(gWidgets2)
  #library(gWidgets2RGtk2)
  options("guiToolkit"="RGtk2")
  # ok=require(baseline)
  # if (!ok){
  #   install.packages("baseline",dependencies = TRUE)
  #   library(baseline)
  # }
  # library(graphics)
  
  scriptDir=utils::getSrcDirectory(function(x) {x})
  try(setwd(scriptDir),silent=TRUE) #Try because does not work when calling from shortcut/package
  
  
  MakeSpinner <- function()
    #This creates a spinner slightly below the middle of the main GUI window.
    #It returns the gwindow object containing the spinner. Just use dispose
    #on this return object to remove the spinner.  
  {
    wspin<-gWidgets2::gwindow(parent=main,width = 50,height=50,visible=FALSE)
    wspin$widget$decorated<-FALSE
    wspin$widget$modifyBg(GtkStateType["normal"],"white")
    dum<-wspin$widget$getPosition()
    wspin$widget$move(dum$root.x,dum$root.y+100)
    spin<-RGtk2::gtkSpinner()
    gWidgets2::add(wspin,spin)
    gWidgets2::visible(wspin)<-TRUE
    spin$start()
    return(wspin)
  }
  
  
  
  Plot_Baseline <- function(wn,sp){
    degr <- as.numeric(gWidgets2::svalue(flyt)$'Polynomial degree')
    repet <- as.numeric(gWidgets2::svalue(flyt)$'Repetitions')
    tole <- as.numeric(gWidgets2::svalue(flyt)$'Tolerance')
    spin=MakeSpinner()
    bc <- baseline::baseline(sp, method='modpolyfit', deg=degr,t=wn, rep=repet, tol =tole)
    mean_raw <- colMeans(bc@spectra)
    mean_base <- colMeans(bc@baseline)
    mean_corrected <- colMeans(bc@corrected)
    par(mfrow=c(3,1))
    plot(wn,mean_raw,type="l",lwd=2,cex=0.8)
    plot(wn,mean_base,type="l",lwd=2,cex=0.8)
    lines(wn,mean_base,type="l",lwd=2)
    plot(wn,mean_corrected,type="l",lwd=2,cex=0.8)
    lines(wn,mean_corrected,type="l",lwd=2)
    lines(wn,rep(0,length(wn)),lty=2)
    gWidgets2::dispose(spin)
  }
  
  
  if (is.null(lefichier)) lefichier=gWidgets2::gfile("Select a Raman data file",initial.dir = scriptDir)
  dum1=utils::read.table(lefichier,sep="\t", dec=".",header=FALSE)
  wn=as.numeric(dum1[1,-1])
  sp=as.matrix(dum1[-1,-1])
  ids=as.character(dum1[,1])
  
  main=gWidgets2::gwindow(visible=FALSE,width = 1000,height=700,
               title="Raman baseline remover - Polyfit")
  gg=gWidgets2::ggroup(cont=main, horizontal = FALSE)
  file_label=gWidgets2::glabel(paste("File name: ",lefichier,sep=""),cont=gg)
  
  flyt=gWidgets2::gformlayout(container=gg)
  gWidgets2::gspinbutton(label="Polynomial degree",from=1,to=20,by=1,value=5,container=flyt,
          handler=function(h,...){
            Plot_Baseline(wn,sp)
          })
  gWidgets2::gedit("50",label="Repetitions",container=flyt,
        handler=function(h,...){
          Plot_Baseline(wn,sp)
        })
  gWidgets2::gedit("0.005",label="Tolerance",container=flyt,
        handler=function(h,...){
          Plot_Baseline(wn,sp)
        })
  
  save_but=gWidgets2::gbutton("Save data",cont=gg,handler=function(h,...){
    degr <- as.numeric(gWidgets2::svalue(flyt)$'Polynomial degree')
    repet <- as.numeric(gWidgets2::svalue(flyt)$'Repetitions')
    tole <- as.numeric(gWidgets2::svalue(flyt)$'Tolerance')
    spin=MakeSpinner()
    bc <- baseline::baseline(sp, method='modpolyfit', deg=degr,t=wn, rep=repet, tol =tole)
    basedata <- cbind(ids,rbind(wn,bc@baseline))
    correcteddata <- cbind(ids,rbind(wn,bc@corrected))
    lepath=dirname(lefichier)
    lenom=basename(lefichier)
    endname<-substr(lenom,start=regexpr("_",lenom),stop=nchar(lenom))
    startname <- substr(lenom,1,regexpr("_",lenom)-1)
    basefilename <- file.path(lepath, paste(startname,"BaseLine",endname,sep=""))
    correctedfilename <- file.path(lepath, paste(startname,"Corrected",endname,sep=""))
    utils::write.table(basedata,file=basefilename,sep="\t",row.names=FALSE,col.names=FALSE)
    utils::write.table(correcteddata,file=correctedfilename,sep="\t",row.names=FALSE,col.names=FALSE)
    gWidgets2::dispose(spin)
  })
  
  leplot<-gWidgets2::ggraphics(container = gg, expand=TRUE, fill=TRUE)
  
  gWidgets2::visible(main)<-TRUE
  
  gWidgets2::visible(leplot)<-TRUE
  Plot_Baseline(wn,sp)
  
  
  
  
}
