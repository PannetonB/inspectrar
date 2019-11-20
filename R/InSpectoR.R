#' InSpectoR
#' 
#' InSpectoR is a tool for exploring spectroscopic data sets and to develop
#' models. It was developped as a companion to a custom spectrometric data
#' acquisition system for gathering induced fluorescence with up to 6 
#' excitation wavelengths, Raman spectra excited at ~785 nm and transmittance
#' data or reflectance data. It is general enough to accomodate any spectroscopic
#' data sets provided that the format of the required files conformed to the
#' description below. The program is structured around a top level window with tabs giving
#' access to other interfaces dedicated to specialized tasks such as selecting
#' and viewing raw data, performing data normalization, doing some
#' pretreament on the data, performing principal component analysis,
#' developping PLS models... These are referred to as "task interfaces".
#' 
#' The data sets are made of a number of files defined as follows:
#' \itemize{
#'   \item Y_somename.txt  : Y data file.
#'     \itemize{
#'      \item This is a tab delimited text file with a header line and decimal symbol is ".". 
#'     \item After header, each following lines correspond to samples.
#'     \item Columns contain user defined data such as a sample
#'     identification number, some data or factor associated
#'     with the sample (e.g. color, pH, concentration of
#'                  something...). 
#'     \item The filename MUST start with "Y_". 
#'     It is followed by "somename" which is a user defined
#'     name (anything will do) and followed by the ".txt"
#'     extension. 
#'     \item The GUI offers access to a support function
#'     that can merge Y data in an other file. This is useful when
#'     spectral data and accompanying data are acquired in separate
#'     processes that might result in data stored in different files.
#'     The only requirement for successfull merging is that both
#'     files share a common column.
#'     }
#'  \item *_somename*.txt : Spectral data files. Tab delimiter, header and decimal=".".
#'     \itemize{
#'     \item These are files containing spectroscopic data.
#'     \item The first line contains wavelengths or wavenumbers
#'     associated with the spectra. 
#'     \item The first column is a sample ID (same as first column of
#'                     associated Y file). 
#'     \item One sample per line.
#'     \item The first cell [1,1] is not used (could be empty).
#'     \item Name: The first star to the left is a placeholder for
#'     describing the data. For example, EX350 could be used
#'     to indicate fluorescence spectra obtained under
#'     excitation at 350 nm. These data type are used throughout
#'     the program as labels for options, buttons... Should not
#'     start with a number or special character.
#'     It is up to the user to define
#'     and recognize the names. The * before the extension
#'     could be anything. The key is that "somename" matches
#'     the "somename" of the associated Y_somename.txt file.
#'     This is how Y data are connected to spectral data.
#'     }
#'  }
#'  
#'     It is mandatory to group Ys and spectral data files
#'     in the same directory. There could be as many spectral
#'     data files as the user wants. 
#'     
#'     
#'     
#'
#' @param yfile to load data associated with this Y data file.
#'              Default is NULL: no data loaded on startup.
#' @param parcomp to set up for parallel computation with caret package
#'                (default to TRUE) 
#' @param MainWidth width of GUI in pixels (default to 1200)
#' @param MainHeight height of GUI in pixels (default to 800)
#'
#' @return NOTHING
#' @export
#' 
#' @import pls
#'
#' @examples 
#' dfile <- system.file("foodstuff_powder","Y_foodstuff.txt",package="inspectrar")
#' InSpectoR(dfile,0,1600,1024)
InSpectoR <- function(yfile=NULL,parcomp=TRUE,MainWidth=1200,MainHeight=800) 
#*************************************************************************
# InSpectoR is a tool for exploring spectroscopic data sets and to develop
# models. It was developped as a companion to PolySpecteur, an R function 
# for the acquisition of data from SpectrAAC, a custom spectrometric data
# acquisition system for gathering induced fluorescence with up to 6 
# excitation wavelengths, Raman spectra excited at ~785 nm and transmittance
# data or reflectance data. It is general enough to accomodate any spectroscopic
# data sets provided that the format of the required files conformed to the
# description below.
#
# INPUTS:
#       parcomp :   TRUE to enable parallel computing where available (caret package)
#     MainWidth :   width of the main GUI in pixels. Forced to 700 if less than 701
#    MainHeight :   Height of the main GUI in pixels. Forced to 700 if less than 701.
#
# The data sets are made of a number of files defined as follows:
#   Y_somename.txt  : Y data file.
#                     This is a tab delimited text file with a header line and
#                     decimal symbol is ".". After header,
#                     each following lines correspond to samples.
#                     Columns contain user defined data such as a sample
#                     identification number, some data or factor associated
#                     with the sample (e.g. color, pH, concentration of
#                     something...). The filename MUST start with "Y_". It
#                     is followed by "somename" which is a user defined
#                     name (anything will do) and followed by the ".txt"
#                     extension. The GUI offers access to a support function
#                     that can merge Y data in an other file. This is useful when
#                     spectral data and accompanying data are acquired in separate
#                     processes that might result in data stored in different files.
#                     The only requirement for successfull merging is that both
#                     files share a common column.
#   *_somename*.txt : Spectral data files. Tab delimiter, header and decimal=".".
#                     These are files containing spectroscopic data.
#                     The first line contains wavelengths or wavenumbers
#                     associated with the spectra. The first column
#                     is a sample ID (same as first column of
#                     associated Y file). One sample per line.
#                     The first cell [1,1] is not used (could
#                     be empty).
#                     Name: The first star to the left is a placeholder for
#                     describing the data. For example, EX350 could be used
#                     to indicate fluorescence spectra obtained under
#                     excitation at 350 nm. These data type are used throughout
#                     the program as labels for options, buttons... Should not
#                     start with a number or special character.
#                     It is up to the user to define
#                     and recognize the names. The * before the extension
#                     could be anything. The key is that "somename" matches
#                     the "somename" of the associated Y_somename.txt file.
#                     This is how Y data are connected to spectral data.
#                     It is mandatory to group Ys and spectral data files
#                     in the same directory. There could be as many spectral
#                     data files as the user wants.
#                    
#                     
# The program is structured around a top level window with tabs giving
# access to other interfaces dedicated to specialized tasks such as selecting
# and viewing raw data, performing data normalization, doing some
# pretreament on the data, performing principal component analysis,
# developping PLS models... These are referred to as "task interfaces".
# 
# Program by
# Bernard Panneton
# bernard.panneton@agr.gc.ca
# pannetonb@gmail.com
# March 2019
#*************************************************************************
# TO DO list
# 
#*************************************************************************
{  
#-------------------------                                        
  #For RGtk2 interface
  options("guiToolkit"="RGtk2")
  
  #-------------------------  
  #******************************************************************************
  #Create cluster for parallel processing with caret package
  #******************************************************************************
  #-------------------------
  if (parcomp){
    cat("*******************\n")
    cat("Creating cluster of cpus for parallel processing. Wait.\n")
    cat("*******************\n")
    mycluster<<-parallel::makeCluster(parallel::detectCores()-1)
    doParallel::registerDoParallel(mycluster)
    foreach::registerDoSEQ()
  }
  
#Local variables-------------------------  
  #******************************************************************************
  #Local variables
  #******************************************************************************
#-------------------------
  #Load logo image in package
  logoimg <- system.file("InSpectraR_Logo_600_130V.png",package="inspectrar")
  if (!file.exists(logoimg)) 
    logoimg<-"InSpectraR_Logo_600_130V.png"  #Load logo image
  
#Global variables-------------------------
  #******************************************************************************
  #Global variables
  #******************************************************************************
#-------------------------
  #ISR_env is used to store global variables for internal use by InSpectoR
  calling_enviro <<- parent.env(environment())
  ISR_env <- new.env(parent=emptyenv())
  ISR_env$laversion<-"0.1.1" #Version of InSpectoR
  ISR_env$ladate <- "November 2019"    #Date of current version
  index_data <<- 1          #gnotebook page number for data_tab
  index_applymods <<- 2     #gnotebook page number for apply_tab
  index_prepro <<- 3        #gnotebook page number for prepro_tab
  index_acp <<- 4           #gnotebook page number for pca_tab
  index_plsda <<- 5         #gnotebook page number for plsda_tab
  index_pls <<- 6           #gnotebook page number for pls_tab
  Ys_df<<-NULL              #data frame holding Y data file content
  Ys_df_sub<<-NULL          #data frame holding Y data file content after subsetting
  subset_ind <<- NULL       #logical vector as indices for samples to include in subset
  DataDir<<-NULL            #where to find the data set
  XDatalist<<-list()        #list of selected spectral data files
  All_XData <<- list()      #List containing all valid spectral data associated with Y's file.
  XData <<- list()          #List of matrices containing the spectral data
                            #Order in the list matches the one in XDatalist
  XData_p<<-list()          #List of matrices containing the pre-processed spectral data
                            #Order in the list matches the one in XDatalist
  N_ech<<-NULL              #Number of samples in the data set.
  prepro_params <<-list()   #list of preprocessing parameters
  PreProDone <<- FALSE      #Flag to TRUE when preprocessing was applied. Flip  to FALSE when selecting
                            #new spectra type in the Data tab.
  PlotIsByFactor <<- FALSE  #flag for plot type in the data tab
  raw_sort_column <<-NULL   #To remember column used to sort raw data.
  doplot_lesX <<- TRUE      #flag to indicate lesX "changed" handler to plot or not
  isSelectMode<<-FALSE       #flag indicating if interaction with ggraphics is for 
                            #selecting data. See AddPop_2_raw_ggraphics.R
  isZoomMode<<-TRUE         #flag indication if interaction with ggraphics is for
                            #zooming on data. See AddPop_2_raw_ggraphics.R
  isZoomAll<<-FALSE         #flag indication if interaction with ggraphics is for
                            #zooming on all data. See AddPop_2_pca_ggraphics.R
  lesACPs<<-list()          #PCA results for each data type
  lesNCPs<<-list()          #Estimated number of PCs for each data type.
  Sc_plot_params<<-list()   #parameters of the current score plot for zoom and select actions.
  ISR_env$SelectedScores <-NULL  #vector of the index of points selected from the PCA score plot.
  plsda_txt_output<<-""     #Capture console output from plsda analysis for display
  plsda_set<<-list()        #Data sets used for computing plsda
  plsda_inTrain<<-as.numeric() #Indices of elements in plsda_set used for training
                               #-plsda_inTrain used for testing.
  PLS_scatter_params<<-list()   #parameters of the current scatter plot on PLS 
                                #for zoom and select actions.
  plsdaFit <<- NULL         # list of plsda models
  plsFit <<- NULL           # list of pls models
  
#Generic functions-------------------------  
  #******************************************************************************
  #Generic functions
  #******************************************************************************
#-------------------------
  #Make sure we are in the directory of InSpectoR.R where all accompanying
  #R function files are located so sourcing will work.
  #Trick to ID the directory of InSpectoR.R. Then we can work with relative paths.
  scriptDir<-utils::getSrcDirectory(function(x) {x})
  try(setwd(scriptDir),silent=TRUE) #Try because does not work when calling from shortcut
  
  get_DataTypes <- function(datalist)
  # Helper function to extract data type from XDatalist.
  # To be used with InSpectoR() 
  #*********************************************************************************************************** 
  # B. Panneton, June 2017
  #*********************************************************************************************************** 
  
  {
    dum=as.character(datalist)
    datatypeIDs=substr(dum,start=c(1,1),regexpr("_",dum)-1)
  }
  
  #***********************************************************************
  MakeSpinner <- function()
  #This creates a spinner slightly below the middle of the main GUI window.
  #It returns the gwindow object containing the spinner. Just use dispose
  #on this return object to remove the spinner.  
  {
    wspin<-gWidgets2::gwindow(parent=mymain,width = 50,height=50,visible=FALSE)
    wspin$widget$decorated<-FALSE
    wspin$widget$modifyBg(GtkStateType["normal"],"white")
    dum<-wspin$widget$getPosition()
    wspin$widget$move(dum$root.x,dum$root.y+100)
    spin<-RGtk2::gtkSpinner()
    gWidgets2::add(wspin,spin)
    visible(wspin)<-TRUE
    spin$start()
    return(wspin)
  }
      
  #***********************************************************************
  Clear_graph <- function(legraph)
  #Clears graphics of legraph ggraphics object.  
  {
    dum=class(try(visible(legraph)<-TRUE,silent = TRUE))
    if (dum!="try-error") plot.new()
  }
  
  #***********************************************************************
  get_DataType_Names <- function(laliste)
  #Return a character array containing the strings in the spectral data files.
  #i.e. , the string corresponding to the star on the left in  *_somename*.txt.
  #See comments at the beginning of InSpectoR for details about file names.  
  {
    dum<-as.character(laliste)
    labs<-substr(dum,start=c(1,1),regexpr("_",dum)-1)
    return(labs)
  }
  
  #***********************************************************************
  build_truncation_widget <-function(sp_list,gf_truncation,le_r)
  #Build and fills a glayout widget inside the gf_truncation
  #gframe on the pca_tab for defining min and max value for truncating spectra.
  #INPUTS:  
  #     sp_list       : the list of selected files (XDatalist)  
  #     gf_truncation : the parent gframe
  #     le_r          : list (one element by spectrum type) of default ranges (c(MIN,MAX))  
    
  {
    #Empty lists for min and max values of the range
    spin_max <- list()
    spin_min <- list()
    
    #Array of gspinbutton for mini and maxi
    lyt<-gWidgets2::glayout(cont=gf_truncation)   
    #Cheap trick to control the width : gWidgets2::add extra space
    #Column labels
    lyt[1,2] <- "         MINIMUM         "
    lyt[1,3] <- "         MAXIMUM         "
    
    #labels each row. In InSpectoR, get list from XDatalist
    labs<-get_DataType_Names(sp_list)
    
    N<-length(labs)
    
    #Layout gspinbuttons and row labels
    for (k in 1:N){
      lyt[k+1,1] <- labs[[k]]
      lyt[k+1,2] <- (spin_min[[k]] <- gWidgets2::gspinbutton(le_r[[k]][1],le_r[[k]][2],10,value=le_r[[k]][1],cont=lyt))
      lyt[k+1,3] <- (spin_max[[k]] <- gWidgets2::gspinbutton(le_r[[k]][1],le_r[[k]][2],10,value=le_r[[k]][2],cont=lyt))
      
    }
    
    #Handlers to make sure there is a 50 units difference between min and max
    #for a spectra. When acting on mini, maxi is adjusted to suit the rule. Reverse
    #applies when acting on maxi. Absolute mini and maxi are set by the definition
    #of the gspinbutton.
    for (k in 1:N){
      gWidgets2::addHandlerChanged(spin_min[[k]], action=c(spin_max[[k]],k,N,le_r), handler=function(h,..) { 
        PreProDone<<-FALSE
        mini <- gWidgets2::svalue(h$obj)
        maxi <- gWidgets2::svalue(h$action[[1]])
        lek <- h$action[[2]]
        NN <- h$action[[3]]
        TOP <- h$action[[4]][2]
        if ((mini+50)>TOP){
          gWidgets2::svalue(h$action[[1]]) <- TOP
          gWidgets2::svalue(h$obj) <- TOP-50
        }
        if ((mini+50) > maxi)
          gWidgets2::svalue(h$action[[1]])<-mini+50  #force maxi at least 50 over mini
        
        #Force scaling to "none" if it was "Value" (i.e. i=2) and sets new range for value
        dum <- gWidgets2::svalue(gf1$children[[1]][k+1,2],index=TRUE)
        if (dum==2)
          gWidgets2::svalue(gf1$children[[1]][k+1,2],index=TRUE)<-1
        #Get min and max wavelength (number) from gf_truncation
        dum2<-sapply(gf_truncation$children[[1]][-1,-1],gWidgets2::svalue)
        dum2<-matrix(dum2,nrow=N,ncol=2,byrow=FALSE)
        #Build the list of min and max for builg_byvalue_widget
        le_r2<-lapply(seq_len(nrow(dum2)), function(i) dum2[i,])
        #Modify list with user selection
        le_r2[[lek]]<-c(mini,maxi)
        # #Rebuild the gf_byspectra_scaling_data table
        gWidgets2::delete(gf1,gf1$children[[1]])
        gf1 <- build_byvalue_scaling_widget(XDatalist,gf1,le_r2)
      })
      gWidgets2::addHandlerChanged(spin_max[[k]], action=c(spin_min[[k]],k,N,le_r), handler=function(h,..) { 
        PreProDone<<-FALSE
        maxi <- gWidgets2::svalue(h$obj)
        mini <- gWidgets2::svalue(h$action[[1]])
        lek <- h$action[[2]]
        NN <- h$action[[3]]
        BOTTOM <- h$action[[4]][1]
        if ((maxi-50)<BOTTOM){
          gWidgets2::svalue(h$action[[1]]) <- BOTTOM
          gWidgets2::svalue(h$obj) <- BOTTOM+50
        }
        if ((maxi-50)<mini)
          gWidgets2::svalue(h$action[[1]])<-maxi-50  #force mini at least 50 below maxi
        
        #Force scaling to "none" if it was "Value" (i.e. i=2) and sets new range for value
        dum <- gWidgets2::svalue(gf1$children[[1]][k+1,2],index=TRUE)
        if (dum==2)
           gWidgets2::svalue(gf1$children[[1]][k+1,2],index=TRUE)<-1
        #Get min and max wavelength (number) from gf_truncation
        dum2<-sapply(gf_truncation$children[[1]][-1,-1],gWidgets2::svalue)
        dum2<-matrix(dum2,nrow=N,ncol=2,byrow=FALSE)
        #Build the list of min and max for builg_byvalue_widget
        le_r2<-lapply(seq_len(nrow(dum2)), function(i) dum2[i,])
        #Modify list with user selection
        le_r2[[lek]]<-c(mini,maxi)
        # #Rebuild the gf_byspectra_scaling_data table
        delete(gf1,gf1$children[[1]])
        gf1 <- build_byvalue_scaling_widget(XDatalist,gf1,le_r2)
      })
    }
    return(gf_truncation)
  }
  
  #***********************************************************************
  build_byvalue_scaling_widget <-function(sp_list,gf_byspectra_scaling_data,le_r)
    #Build and fills a glayout widget inside the gf_byspectra_scaling_data
    #gframe on the prepro tab for defining wavelength and bandwidth fo spectral
    #data used to compute scaling factor (i.e. mean value over bandwidth).
    #INPUTS:  
    #     sp_list       : the list of selected files (XDatalist)  
    #     gf_byspectra_scaling_data : the parent gframe
    #     le_r          : list (one element by spectrum type) of cneter and BW (c(Ctre,BW))  
    
  {
    #Empty lists for center and range values for the scaling area.
    spin_centre <- list()
    spin_bw <- list()
    byspectra_type <- list()
    
    #Array of gspinbutton for mini and maxi
    lyt2<-gWidgets2::glayout(cont=gf_byspectra_scaling_data)   
    #Cheap trick to control the width : gWidgets2::add extra space
    #Column labels
    lyt2[1,2] <- "         Type            "
    lyt2[1,3] <- "         Centre value         "
    lyt2[1,4] <- "         Bandwith         "
    
    #labels each row. In InSpectoR, get list from XDatalist
    labs<-get_DataType_Names(sp_list)
    
    N<-length(labs)
    
    #Layout gspinbuttons and row labels
    for (k in 1:N){
      lyt2[k+1,1] <- labs[[k]]
      lyt2[k+1,2] <- (byspectra_type[[k]] <- gWidgets2::gradio(c("None      ","Value          ",
                                               "Mean=1    "),
                                             horizontal = TRUE,
                                             cont=lyt2,
                                             handler=function(h,...){ 
                                               PreProDone <<- FALSE
                                             }
      ))
      lyt2[k+1,3] <- (spin_centre[[k]] <- gWidgets2::gspinbutton(le_r[[k]][1]+1,le_r[[k]][2],1,value=le_r[[k]][1],cont=lyt2))
      lyt2[k+1,4] <- (spin_bw[[k]] <- gWidgets2::gspinbutton(0,100,1,value=0,cont=lyt2))
      
    }
    for (k in 1:N){
      gWidgets2::addHandlerChanged(spin_centre[[k]], action=list(spin_bw[[k]],N,le_r,k), handler=function(h,..){
        PreProDone<<-FALSE
        bw<-gWidgets2::svalue(h$action[[1]])
        NN<-h$action[[2]]
        mon_r<-h$action[[3]]
        lek<-h$action[[4]]
        ctr<-gWidgets2::svalue(h$obj)
        if ((ctr-bw)<mon_r[[lek]][1]){  #start of bw below range, adjust centre value
          gWidgets2::svalue(h$obj)<-mon_r[[lek]][1]+bw
        }
        if ((ctr+bw)>mon_r[[lek]][2]){   #end of bw over range, adjust centre value
          gWidgets2::svalue(h$obj)<-mon_r[[lek]][2]-bw
          
        }
      })
      gWidgets2::addHandlerChanged(spin_bw[[k]], action=list(spin_centre[[k]],le_r,k), handler=function(h,..){
        PreProDone<<-FALSE
        bw<-gWidgets2::svalue(h$obj)
        ctr<-gWidgets2::svalue(h$action[[1]])
        lek<-h$action[[3]]
        mon_r<-h$action[[2]]
        if ((ctr-bw)<mon_r[[lek]][1]){  #start of bw below range, adjust centre value
          gWidgets2::svalue(h$action[[1]])<-mon_r[[lek]][1]+bw
        }
        if ((ctr+bw)>mon_r[[lek]][2]){   #end of bw over range, adjust centre value
          gWidgets2::svalue(h$action[[1]])<-mon_r[[lek]][2]-bw
        }
      })
    }
    return(gf_byspectra_scaling_data)
  }
  
  
  #***********************************************************************
  build_savgol_widget <-function(sp_list,pf)
    #Build and fills a glayout widget inside the gf_byspectra_scaling_data
    #gframe on the prepro tab for defining wavelength and bandwidth fo spectral
    #data used to compute scaling factor (i.e. mean value over bandwidth).
    #INPUTS:  
    #     sp_list       : the list of selected files (XDatalist)  
    #     pf            : the parent gframe  
    
  {
    #Empty lists for center and range values for the scaling area.
    chk_savgol <- list()
    m <- list()
    w <- list()
    p <- list()
    
    #Array of gspinbutton for mini and maxi
    lyt3<-gWidgets2::glayout(cont=pf)   
    #Cheap trick to control the width : gWidgets2::add extra space
    #Column labels
    lyt3[1,2] <- "     Apply    "
    lyt3[1,3] <- "         Window size      "
    lyt3[1,4] <- "         Deriv. order     "
    lyt3[1,5] <- "         Polyn. order     "
    
    #labels each row. In InSpectoR, get list from XDatalist
    labs<-get_DataType_Names(sp_list)
    
    N<-length(labs)
    
    #Layout gspinbuttons and row labels
    for (k in 1:N){
      lyt3[k+1,1] <- labs[[k]]
      lyt3[k+1,2] <- (chk_savgol[[k]] <- gWidgets2::gcheckbox("",
                                                    handler=function(h,...){ 
                                                      PreProDone <<- FALSE
                                                    }
                      ))
      lyt3[k+1,3] <- (w[[k]] <- gWidgets2::gspinbutton(1,51,2,value=11,cont=lyt3))
      lyt3[k+1,4] <- (m[[k]] <- gWidgets2::gspinbutton(0,2,1,value=0,cont=lyt3))
      lyt3[k+1,5] <- (p[[k]] <- gWidgets2::gspinbutton(1,7,1,value=1,cont=lyt3))
      
    }
    for (k in 1:N){
      addHandlerChanged(w[[k]], action=list(p[[k]]), handler=function(h,..){
        PreProDone<<-FALSE
        p<-gWidgets2::svalue(h$action[[1]])
        w<-gWidgets2::svalue(h$obj)
        if (w<(p+3)){  #force window size to p+3
          gWidgets2::svalue(h$obj)<-p+3
        }
        #Make sur w is odd.
        if (w %% 2 == 0) gWidgets2::svalue(h$obj)<-w+1
      })
      gWidgets2::addHandlerChanged(m[[k]], action=list(p[[k]]), handler=function(h,..){
        PreProDone<<-FALSE
        m<-gWidgets2::svalue(h$obj)
        p<-gWidgets2::svalue(h$action[[1]])
        if (m>=p){  #Force p to m+1
          gWidgets2::svalue(h$action[[1]])<-m+1
        }
      })
      gWidgets2::addHandlerChanged(p[[k]], action=list(w[[k]],m[[k]]), handler=function(h,..){
        PreProDone<<-FALSE
        p<-gWidgets2::svalue(h$obj)
        w<-gWidgets2::svalue(h$action[[1]])
        m<-gWidgets2::svalue(h$action[[2]])
        if (m>=p){  #Force m to (p-1)
          gWidgets2::svalue(h$action[[2]])<-p-1
        }
        if (w<(p+3)){  #Force w to p+3
          gWidgets2::svalue(h$action[[1]])<-p+3
        }
      })
    }
    return(pf)
  }
  
  #***********************************************************************
  Plot_Confusion_Matrix<-function(conf)
  #Creates a plot for a confusion matrix computed with confusionMatrix from caret.
  {
    melted<-reshape2::melt(conf)
    p1<-ggplot2::ggplot(data=melted,ggplot2::aes(x=Prediction,y=Reference,fill=value)) 
    p1<-p1 + ggplot2::geom_tile(colour="grey50")
    p1<-p1 + ggplot2::scale_fill_gradient2(low="blue",high="red",mid="white")
    p1<-p1 + ggplot2::theme(axis.text.x = element_text(angle=60,hjust=1))
    p1<-p1 + ggplot2::theme(text = element_text(size=9))
    p1<-p1 + ggplot2::geom_text(ggplot2::aes(label=value), size=3)
    
    print(p1)
  }
  
  #***********************************************************************
  Apply_PreTreatments <- function(prepro_par=NULL)
  #Apply all pretreatments as defined in the prepro tab.
  #Construct the prepro_param list to store parameters to be saved with models.  
  #   gWidgets2::delete(gf_truncation,gf_truncation$children[[1]])
  # gf_truncation<-build_truncation_widget(XDatalist,gf_truncation,le_r)
  # gWidgets2::delete(gf1,gf1$children[[1]])
  # gf1 <- build_byvalue_scaling_widget(XDatalist,gf1,le_r)
  # gWidgets2::delete(gf_savgol,gf_savgol$children[[1]])
  # gf_savgol <- build_savgol_widget(XDatalist,gf_savgol)
  {
    #Truncation
    if (is.null(prepro_par)){  #retrieve from GUI
      N=length(XData)
      lim_wid <- gf_truncation$children[[1]]
      trunc_limits <- matrix(sapply(lim_wid[-1,-1],gWidgets2::svalue),ncol=2,byrow=FALSE)
      prepro_params$trunc_limits <<- trunc_limits
    }else   #extract from prepro_params
    {
      trunc_limits<-prepro_params$trunc_limits
      #load into gui
      le_r <- lapply(seq_len(nrow(trunc_limits)), function(i) trunc_limits[i,])
      gWidgets2::delete(gf_truncation,gf_truncation$children[[1]])
      gf_truncation<-build_truncation_widget(XDatalist,gf_truncation,le_r)
    }
    ii<-as.list(1:nrow(trunc_limits))
    lapply(ii,function(x){
      wl<-XData[[x]][1,]
      XData_p[[x]] <<- XData[[x]][,((wl>=trunc_limits[x,1]) & (wl<=trunc_limits[x,2]))]
    }) 
    
    #Then per_spectra normalization - value
    if (is.null(prepro_par)){   #retrieve from GUI
      N=length(XData_p)
      if (N>1){
        type <- sapply(gf1$children[[1]][-1,2],gWidgets2::svalue,index=TRUE)
      }else
      {
        type <- gWidgets2::svalue(gf1$children[[1]][-1,2],index=TRUE)
      }
      letest=any(type==2)
      prepro_params$byspectra_scaling_index <<- type
      cntr_n_w_table <- sapply(gf1$children[[1]][-1,-1],gWidgets2::svalue)
      cntr_n_w <- matrix(cntr_n_w_table,ncol=3,byrow=FALSE)
      cntr_n_w <- matrix(as.numeric(cntr_n_w[,-1]),ncol=2,byrow=FALSE)
      prepro_params$cntr_n_w <<- cntr_n_w
    }else   #extract from prepro_params
    {
      type <- prepro_params$byspectra_scaling_index
      letest=any(type==2)
      cntr_n_w <- prepro_par$cntr_n_w
      gWidgets2::delete(gf1,gf1$children[[1]])
      gf1 <- build_byvalue_scaling_widget(XDatalist,gf1,le_r)
      dum <- as.list(seq_along(type))
      sapply(dum, function(ii){
        leGradio <- gf1$children[[1]][-1,2][[ii]]
        gWidgets2::svalue(leGradio,index=TRUE) <- as.numeric(type[ii])
      })
      sapply(dum, function(ii){
        lacase <- gf1$children[[1]][ii+1,3]
        gWidgets2::svalue(lacase) <- cntr_n_w[ii,1]
        lacase <- gf1$children[[1]][ii+1,4]
        gWidgets2::svalue(lacase) <- cntr_n_w[ii,2]
      })
      
    }
    if (letest){
      wl<-lapply(XData_p, function(y) y[1,]) # a list of wl vector, one per spectrum type
      X<-lapply(XData_p, function(y) y[-1,]) #a list of spectral data matrix, no wl.
      ii<-as.list(1:nrow(cntr_n_w))
      XData_p <<- lapply(ii,function(x){
        if (type[x]==2){
          i1<-which(wl[[x]]>=(cntr_n_w[x,1]-cntr_n_w[x,2]))[1]
          i2<-which(wl[[x]]>=(cntr_n_w[x,1]+cntr_n_w[x,2]))[1]
          X[[x]]<-t(apply(X[[x]],1,function(z) z/mean(z[i1:i2])))
        }
        rbind(wl[[x]],X[[x]])  #rebuild matrices with wl
      })
    }
    
    #Then per_spectra normalization - closure
    if (is.null(prepro_par)){   #retrieve from GUI
      N=length(XData)
      if (N>1){
        type <- sapply(gf1$children[[1]][-1,2],gWidgets2::svalue,index=TRUE)
      }else
      {
        type <- gWidgets2::svalue(gf1$children[[1]][-1,2],index=TRUE)
      }
      letest=any(type==3)
    }else  #extract from prepro_params
    {
      type <- prepro_params$byspectra_scaling_index
      letest<-any(type==3)
    }
    if (letest){
      wl<-lapply(XData_p, function(y) y[1,]) # a list of wl vector, one per spectrum type
      X<-lapply(XData_p, function(y) y[-1,]) #a list of spectral data matrix, no wl.
      iis<-as.list(1:length(wl))
      XData_p <<- lapply(iis, function(ii){
        if (type[ii]==3){
          L <- length(wl[[ii]])
          X[[ii]]<-t(apply(X[[ii]],1,function(z) z*L/sum(z))) #normalize matrices
        }
        rbind(wl[[ii]],X[[ii]]) #rebuild matrices with wl
      })
    } 
    
    #Then Savitzky-Golay
    if (is.null(prepro_par)){    #retrieve from GUI
      N=length(XData)
      if (N>1){
        dosavgol <- sapply(gf_savgol$children[[1]][-1,2],gWidgets2::svalue)
        m<-sapply(gf_savgol$children[[1]][-1,4],gWidgets2::svalue)
        p<-sapply(gf_savgol$children[[1]][-1,5],gWidgets2::svalue)
        w<-sapply(gf_savgol$children[[1]][-1,3],gWidgets2::svalue)
      }else
      {
        dosavgol <- gWidgets2::svalue(gf_savgol$children[[1]][-1,2])
        m<-gWidgets2::svalue(gf_savgol$children[[1]][-1,4])
        p<-gWidgets2::svalue(gf_savgol$children[[1]][-1,5])
        w<-gWidgets2::svalue(gf_savgol$children[[1]][-1,3])
      }
      
      letest=any(dosavgol)
     
      prepro_params$do_savgol <<- dosavgol
      prepro_params$m <<- m
      prepro_params$p <<- p
      prepro_params$w <<- w
    }else   #extract from prepro_params
    {
      N=length(XData)
      dum <- as.list(seq_len(N))
      dosavgol<-prepro_params$do_savgol
      letest=any(dosavgol==TRUE)
      m <- prepro_par$m 
      p <- prepro_params$p
      w <- prepro_par$w
      gWidgets2::delete(gf_savgol,gf_savgol$children[[1]])
      gf_savgol <- build_savgol_widget(XDatalist,gf_savgol)
      sapply(dum, function(ii){
        lecheck <- gf_savgol$children[[1]][-1,2][[ii]]
        gWidgets2::svalue(lecheck) <- dosavgol[ii]
        lem <- gf_savgol$children[[1]][-1,4][[ii]]
        gWidgets2::svalue(lem) <- m[ii]
        lep <- gf_savgol$children[[1]][-1,5][[ii]]
        gWidgets2::svalue(lep) <- p[ii]
        lew <- gf_savgol$children[[1]][-1,3][[ii]]
        gWidgets2::svalue(lew) <- w[ii]
      })
    }
    if (letest){
      #retrieve paramaters
      wl<-lapply(XData_p, function(y) y[1,]) # a list of wl vector, one per spectrum type
      #Truncate wl to suit window size (w)
      #Need to adapt length of wl to result of Savitzky-Golay
      iis<-as.list(1:length(wl))
      wl <- lapply(iis,function(k){
          w_2=floor(w[[k]]/2)
          if (dosavgol[k]){
            wl[[k]][(w_2+1):(length(wl[[k]])-w_2)]
          }else
          {
            wl[[k]]
          }
        }
      )
      X<-lapply(XData_p, function(y) y[-1,]) #a list of spectral data matrix, no wl.
      
      
      XData_p <<- lapply(iis, function(k) {
        if (dosavgol[k]){
          X[[k]]<-prospectr::savitzkyGolay(X[[k]],m[k],p[k],w[k])
        }
        rbind(wl[[k]],X[[k]]) #rebuild matrices with wl
      })
    }
    return()
  }
  
  #***********************************************************************
  pr_2_prin <- function(x)
  # Converts objects of class prcomp (Q-mode PCA) to class princomp (R-mode PCA)
  # Bryan Hanson, DePauw Univ, Sept 2009
  # original is modified by adding list elements (these could be removed to save space)
  {
    
    if (!"prcomp" %in% class(x)) stop("The PCA object was not of class prcomp")
    
    # sdev, center and scale for both classes are the same; no change necessary
    
    x$loadings <- x$rotation
    x$scores <- x$x
    x$call <- "No call available (data converted from class prcomp)"
    x$n.obs <- dim(x$x)[1]	
    class(x) <- c("conPCA", class(x))
    x
  }
#Handlers for widgets-------------------------  
  #******************************************************************************
  #Handlers for widgets
  #******************************************************************************
#-------------------------
  #temporary empty widget so it is recognized in handlers
  pick_data_4_PCA <- gWidgets2::gcombobox("")

  
#Widgets on mymain-------------------------  
  #/././././././././././././././././././././././
  #Widgets on mymain
  #/././././././././././././././././././././././
#-------------------------
  #*************************************************************
  nb_handler <- function(h,...)
  # Handler when switching between gnotebook tabs.
  {
    
    #For data_tab initialization
    if (h$page.no==index_data){
      
      if (length(XData)>0){
        #Adjust gspinbuttons for axis limits for plotting raw data
        leslimites <- get_graph_limits(XData_p)
        leslimites$ylimits <-  leslimites$ylimits
        steps <- diff(range(leslimites$ylimits))/20
        #Setup gspinbuttons and reset axis limits to default
        lapply(proplist,gWidgets2::blockHandlers)
        RGtk2::gtkAdjustmentSetUpper(ymaxi$widget$adjustment,leslimites$ylimits[2])
        RGtk2::gtkAdjustmentSetUpper(ymini$widget$adjustment,leslimites$ylimits[2]-steps)
        RGtk2::gtkSpinButtonSetIncrements(ymaxi$widget,steps,5*steps)
        RGtk2::gtkSpinButtonSetIncrements(ymini$widget,steps,5*steps)
        gWidgets2::svalue(proplist[[1]])<-0
        gWidgets2::svalue(proplist[[2]])<-0
        gWidgets2::svalue(proplist[[3]])<-0
        gWidgets2::svalue(proplist[[4]])<-1  #Dummy to make sure plot is updated when setting proplist[[4]] to 0.
        lapply(proplist,gWidgets2::unblockHandlers)
        gWidgets2::svalue(proplist[[4]])<-0
        #update_rawX_plot(Y_selection)
        
      }
    }
    
    #For prepro_tab initialization
    if (h$page.no==index_prepro){
      if (length(XData)>0){
         
      }
    }
    
    #apply selected preprocessing to data in the order
    #on the prepro_tab GUI before proceeding further to PCA, LDA...
    #PreProDone used to control if this is necessary.
    #PreProDone set to FALSE when changing selection of spectrum files.
    if (h$page.no>index_prepro){
      if ((length(XData)>0) & (!PreProDone)){
        Apply_PreTreatments()
        PreProDone<<-TRUE
      }
    }
    
    # For pca_tab initialization:
    #   - Define list of available spectral data files
    #   - Init widgets
    #   - Do intial plot based on defaults.  
    if (h$page.no==index_acp){
      #Compute PCA on all members of XData_p
      N_samples<-nrow(Ys_df_sub)
      #If many samples, display a spinner.
      if (N_samples>500){
        wspin<-MakeSpinner()
      }
      thr=0.001
      lesACPs <<- lapply(XData_p, function(x) prcomp(x[-1,],tol=thr))
      #lesNCPs <<- lapply(lesACPs, function(x) ncol(x$rotation))
      var_exp=lapply(lesACPs,function(x) summary(x)$importance[2,])
      lesNCPs <<- lapply(var_exp, function(x){
        i2<-as.numeric(which(x<=thr)[1])
        if (x[i2]==0) i2=i2-1  #in case where few samples.
        return(2*i2)
      })
      
      if (N_samples>500){
        dispose(wspin)
      }
      
      #Populate pick_data_4_PCA with proper list of spectra types
      #Block handlers
      gWidgets2::blockHandlers(pick_data_4_PCA)
      datatypeIDs<-get_DataType_Names(XDatalist)
      pick_data_4_PCA[] <- datatypeIDs
      gWidgets2::svalue(pick_data_4_PCA, index=TRUE) <- 1
      gWidgets2::unblockHandlers(pick_data_4_PCA)
      
      
      #Initialize to first spectra type: 
      if (length(XData)>0){
        
        gWidgets2::svalue(nCP_label) <- lesNCPs[[1]]
       
        #Populate pickFactor_4_color and pickFactor_4_label
        indi<-unlist(lapply(Ys_df_sub,is.factor))
        pickFactor_4_color[] <- c("None",names(Ys_df_sub)[indi])
        gWidgets2::svalue(pickFactor_4_color,index=TRUE) <- 1
        indi[length(indi)]<-TRUE  #To be able to label with NoSeq (i.e. last column)
        pickFactor_4_label[] <- c("None",names(Ys_df_sub)[indi])
        gWidgets2::svalue(pickFactor_4_label,index=TRUE) <- 1
        #Populate pickPC1 and pickPC2 widgets
        pickPC1[] <- as.character(1:lesNCPs[[1]])
        pickPC2[] <- as.character(1:lesNCPs[[1]])
        gWidgets2::svalue(pickPC1) <- 1
        gWidgets2::svalue(pickPC2) <- 2
        
        #Do plot with defaults
        enabled(updateACPPlotbut)=TRUE
        updateACPPlot(updateACPPlotbut)
      }
    }
    if (h$page.no==index_plsda){
      #Populate Xdat_4_plsda with proper list of spectra types
      #Block handlers
      gWidgets2::blockHandlers(Xdat_4_plsda)
      datatypeIDs<-get_DataType_Names(XDatalist)
      dum<-as.data.frame(datatypeIDs)
      colnames(dum)<-"X Data Files"
      Xdat_4_plsda[] <- dum
      gWidgets2::svalue(Xdat_4_plsda, index=TRUE) <- 1
      gWidgets2::unblockHandlers(Xdat_4_plsda)
      
      #Populate pickFactor_4_plsda 
      indi<-unlist(lapply(Ys_df_sub,is.factor))
      pickFactor_4_plsda[] <- names(Ys_df_sub)[indi]
      gWidgets2::svalue(pickFactor_4_plsda,index=TRUE) <- 1
      
    }
    if (h$page.no==index_pls){
      #Populate Xdat_4_pls with proper list of spectra types
      #Block handlers
      gWidgets2::blockHandlers(Xdat_4_pls)
      datatypeIDs<-get_DataType_Names(XDatalist)
      dum<-as.data.frame(datatypeIDs)
      colnames(dum)<-"X Data Files"
      Xdat_4_pls[] <- dum
      gWidgets2::svalue(Xdat_4_pls, index=TRUE) <- 1
      gWidgets2::unblockHandlers(Xdat_4_pls)
      
      #Populate pickvar_4_pls 
      not_a_factor <- function(x) !(is.factor(x))
      indi<-unlist(lapply(Ys_df_sub,not_a_factor))
      pickvar_4_pls[] <- names(Ys_df_sub)[indi]
      gWidgets2::svalue(pickvar_4_pls,index=TRUE) <- 1
      
      #Populate pickFactor_4_color and pickFactor_4_label
      indi<-unlist(lapply(Ys_df_sub,is.factor))
      pls_pickFactor_4_color[] <- c("None",names(Ys_df_sub)[indi])
      gWidgets2::svalue(pls_pickFactor_4_color,index=TRUE) <- 1
      indi[length(indi)]<-TRUE  #To be able to label with NoSeq (i.e. last column)
      pls_pickFactor_4_label[] <- c("None",names(Ys_df_sub)[indi])
      gWidgets2::svalue(pls_pickFactor_4_label,index=TRUE) <- 1
      
    }
    
  } 
  
  
#Widgets on apply_tab-------------------------  
  #/././././././././././././././././././././././
  #Widgets on apply_tab
  #/././././././././././././././././././././././
#-------------------------
  
  #*************************************************************  
  apply_pca_model<-function(h,...)
  #This:
    # 1. loads pca model (confirm it is a PCA model!)
    # 2. verify data compatibility based on selected data types
    # 3. apply pre-treatments to data
    # 4. apply pca models
    # 5. launch a modal window for viewing results.
  {
    
    MonACP_mod<-function(h,...)
    {
      ncp=lesNCPs[[gWidgets2::svalue(h$obj,index = TRUE)]]
      pickPC1_mod[] <- as.character(1:ncp)
      pickPC2_mod[] <- as.character(1:ncp)
      gWidgets2::svalue(pickPC1_mod) <- 1
      gWidgets2::svalue(pickPC2_mod) <- 2
    }
    
    ff<-""
    dlg = gWidgets2::gbasicdialog(title="Select PCA model", width=500,handler=function(h,...){
      ff<<-lefile
      rm(lefile,envir = .GlobalEnv)
    })
    gWidgets2::size(dlg)<-c(500,500)
    gv=gWidgets2::gvbox(container = dlg,fill=TRUE)
    # btn<-gWidgets2::gbutton("Pick a PCA model file",cont=gv, handler=function(h,...){
    #   dum=getToolkitWidget(dlg)
    #   dum$window$hide()
    #   lefile<<-gfile("Pick a PCA model file",
    #                 filter=list("RData files" = list(patterns=c("*.RData"))))
    #   insert(texte,lefile,do.newline = TRUE)
    #   insert(texte,"\n\n")
    #   load(lefile, envir = calling_enviro)
    #   Encoding(model_descript$description) <- "UTF-8"
    #   insert(texte,model_descript$description)
    #   dum$window$show()
    # })
    
    
    btn<-gWidgets2::gbutton("Pick a PCA model file",cont=gv, handler=function(h,...){
      dum=getToolkitWidget(dlg)
      dum$window$hide()
      lefile<<-gWidgets2::gfile("Pick a PCA model file",
                     filter=list("RData files" = list(patterns=c("*.RData"))))
      gWidgets2::insert(texte,lefile,do.newline = TRUE)
      gWidgets2::insert(texte,"\n")
      load(lefile, envir = .GlobalEnv)
      Encoding(model_descript$description) <- "UTF-8"
      gWidgets2::insert(texte,model_descript$description)
      gWidgets2::insert(texte,"\n")
      gWidgets2::insert(texte,model_descript)
      gWidgets2::insert(texte,"\n")
      gWidgets2::insert(texte,"Truncation limits: ")
      gWidgets2::insert(texte,toString(prepro_params$trunc_limits[,1]))
      gWidgets2::insert(texte,toString(prepro_params$trunc_limits[,2]))
      gWidgets2::insert(texte,"By spectra scaling selection (1-none, ...): ")
      gWidgets2::insert(texte,toString(prepro_params$byspectra_scaling_index))
      gWidgets2::insert(texte,"Center and BW for by value scaling: ")
      gWidgets2::insert(texte,toString(prepro_params$cntr_n_w[,1]))
      gWidgets2::insert(texte,toString(prepro_params$cntr_n_w[,2]))
      gWidgets2::insert(texte,"Flag for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$do_savgol))
      gWidgets2::insert(texte,"Derivative order for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$m))
      gWidgets2::insert(texte,"Polynomial order for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$p))
      gWidgets2::insert(texte,"Window size for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$w))
      dum$window$show()
    })
    
    
    
    texte=gWidgets2::gtext("",container=gv,expand=TRUE)
    dum <- gWidgets2::visible(dlg)
    if (!dum) return()
    
    
    load(ff,envir = .GlobalEnv)
    Encoding(model_descript$description) <- "UTF-8"
    if (!(model_descript$type=="PCA")){
      gWidgets2::gmessage("Not a PCA model file!",title="Warning!",
               icon="warning")
      return()
    }
    
    
    #Check if data type needed by the model are available. If yes, build XData
    # and XData_p as in Make_XDatalist
    
    #Create a list of data types available in the current data set.
    #Then get data types.
    dum=as.list(levels(lesX[,1]))
    xType=get_DataType_Names(dum)  
    #Find indices of data types required by the model
    indix=pmatch(model_descript$datatype,xType)
    if (any(is.na(indix))){  #no match, cannot apply model
      gWidgets2::gmessage(title="ABORTING!",
               msg="Required data type not available in current data set.\nABORTING!")
      return()
    }
    doplot_lesX <<- FALSE
    gWidgets2::svalue(lesX,index=TRUE)<-indix
    doplot_lesX <<- TRUE
    
    
    # if (!all(model_descript$datatype==get_DataType_Names(XDatalist))){
    #   gmessage(msg=paste("Data selection does not match model requirements.\n",
    #                   "Required data types:\n",
    #                   "\t",paste(model_descript$datatype,collapse=", "),sep=""),
    #            title="Warning!",
    #            icon="warning")
    #   return()
    # }
    
 
    #Do pretreatments
    Apply_PreTreatments(prepro_par = prepro_params)
    
    #Apply model
    N_samples<-nrow(Ys_df_sub)
    #If many samples, display a spinner.
    if (N_samples>500){
      wspin<-MakeSpinner()
    }
    iind=as.list(1:length(XData_p))
    acp_pred <- lapply(iind, function(i) predict(lesACPs[[i]],XData_p[[i]][-1,]))
    if (N_samples>500){
      dispose(wspin)
    }
    
    #Show results
    #Build a modal GUI
    dlg <- gWidgets2::gbasicdialog((title="PCA to new data"),
                        parent=mymain,
                        do.buttons = FALSE)
    gv <- gWidgets2::gvbox(container = dlg)
    ggr <- gWidgets2::ggraphics(container=dlg)
    
    
    tmp <- gWidgets2::gframe("Model description",container=gv,expand=FALSE)
    mod_desc_txt <- gWidgets2::gtext(model_descript$description,container = gv)
    
    tmp<-gWidgets2::gframe("Spectral data set for PCA",container=gv,expand=FALSE)
    #Select data
    pick_data_4_PCA_mod <- gWidgets2::gcombobox("")
    datatypeIDs<-get_DataType_Names(XDatalist)
    pick_data_4_PCA_mod[] <- datatypeIDs
    gWidgets2::add(tmp,pick_data_4_PCA_mod,expand=TRUE)
    gWidgets2::svalue(pick_data_4_PCA_mod,index = TRUE)<-1
    
    #To select components
    tmp<-gWidgets2::gframe("Pick PC1 (X-axis for score biplot)",container=gv,expand=FALSE)
    # sélection de PC1 pour graphes d'ACP
    ncp=lesNCPs[[gWidgets2::svalue(pick_data_4_PCA_mod,index = TRUE)]]
    pickPC1_mod<-gWidgets2::gcombobox(as.character(1:ncp))
    gWidgets2::add(tmp,pickPC1_mod)
    tmp <- gWidgets2::gframe("Pick PC2 (Y-axis for score biplot)",container=gv,expand=FALSE)
    # sélection de PC2 pour graphes d'ACP
    pickPC2_mod<-gWidgets2::gcombobox(as.character(1:ncp),selected=2)
    gWidgets2::add(tmp,pickPC2_mod)
    
    
    gWidgets2::addHandlerChanged(pick_data_4_PCA_mod,MonACP_mod)
    
    
    #To select factor for labelling points
    indi<-unlist(lapply(Ys_df_sub,is.factor))
    indi[length(indi)]<-TRUE  #To be able to label with NoSeq (i.e. last column)
    lbl_f_4_label_mod <- gWidgets2::glabel("Labels for score plot")
    pickFactor_4_label_mod<-gWidgets2::gcombobox(c("None",names(Ys_df_sub)[indi]))
    gWidgets2::svalue(pickFactor_4_label_mod,index=TRUE) <- 1
    lbl_c_4_label_mod <- gWidgets2::glabel("Criteria")
    tmp <- gWidgets2::gframe("Point labeling options",container=gv,expand=FALSE,horizontal=FALSE)
    gg <- gWidgets2::ggroup(container=tmp,expand=TRUE)
    gWidgets2::add(gg,lbl_f_4_label_mod)
    gWidgets2::add(gg,pickFactor_4_label_mod,expand=TRUE)
    
    
    tmp <- gWidgets2::gframe("Data ellipse plotting",container=gv,expand=FALSE,horizontal=TRUE)
    lbl_CI_4_ellipse_mod <- gWidgets2::glabel("Pick option",cont=tmp)
    pickCrit_4_ellipse_mod <- gWidgets2::gcombobox(c("None","0.50","0.95"),container=tmp,
                                      expand=TRUE,horizontal=TRUE)
    
    gWidgets2::addSpace(gv,20)
    onebut <- gWidgets2::gbutton("Plot results",container = gv,
                      handler=function(h,...){
                        leX <- gWidgets2::svalue(pick_data_4_PCA_mod,index=TRUE)
                        ncp<-lesNCPs[[leX]]
                        resACP<-lesACPs[[leX]]
                        var_exp=summary(resACP)$importance[2,]
                        #Plot calibration data - NO LABEL ON POINTS!
                        #Retrieve options (factor to apply color or to label)
                        indi_l<-gWidgets2::svalue(pickFactor_4_label_mod,index=TRUE)
                        wh_indi_l <- which(names(Ys_df_sub)==pickFactor_4_label_mod[indi_l])
                        indi_l<-indi_l-1   #first is none
                        if (indi_l>0){
                          pt.label<-Ys_df_sub[,wh_indi_l]
                        }else
                        {
                          pt.label=NULL
                        }
                        #Do plot
                        gWidgets2::visible(ggr)<-TRUE
                        #To define ellipse drawing parameter for BP_Plot_Scores.R
                        ell_option<-gWidgets2::svalue(pickCrit_4_ellipse_mod)
                        #Plot calibration data
                        if (class(colorby)=="data.frame"){  #no color for plotting. Would be a factor for plotting
                          Sc_plot_params <<- BP_Plot_Scores(resACP$x[,gWidgets2::svalue(pickPC1_mod,index=TRUE)],
                                                          resACP$x[,gWidgets2::svalue(pickPC2_mod,index=TRUE)],
                                                          pcs=c(gWidgets2::svalue(pickPC1_mod,index=TRUE),gWidgets2::svalue(pickPC2_mod,index=TRUE)),
                                                          frac=100*c(var_exp[gWidgets2::svalue(pickPC1_mod,index=TRUE)],
                                                                     var_exp[gWidgets2::svalue(pickPC2_mod,index=TRUE)]),
                                                          ell_option=ell_option)
                        }else
                        {
                          Sc_plot_params <<- BP_Plot_Scores(resACP$x[,gWidgets2::svalue(pickPC1_mod,index=TRUE)],
                                                            resACP$x[,gWidgets2::svalue(pickPC2_mod,index=TRUE)],
                                                            pcs=c(gWidgets2::svalue(pickPC1_mod,index=TRUE),gWidgets2::svalue(pickPC2_mod,index=TRUE)),
                                                            frac=100*c(var_exp[gWidgets2::svalue(pickPC1_mod,index=TRUE)],
                                                                       var_exp[gWidgets2::svalue(pickPC2_mod,index=TRUE)]),
                                                            colorby=colorby,
                                                            ell_option=ell_option) 
                        }
                  
                  
                         
                        
                        #Plot new data
                        if (gWidgets2::svalue(pickPC1_mod,index=TRUE)==gWidgets2::svalue(pickPC2_mod,index=TRUE)){
                          #do nothing
                        }else
                        {
                          points(acp_pred[[leX]][,gWidgets2::svalue(pickPC1_mod,index=TRUE)],
                                 acp_pred[[leX]][,gWidgets2::svalue(pickPC2_mod,index=TRUE)],
                                 type="p",
                                 pch=21,col="black",bg="white",cex=2)
                          #add label to points
                          if (indi_l>0) text(acp_pred[[leX]][,gWidgets2::svalue(pickPC1_mod,index=TRUE)],
                               acp_pred[[leX]][,gWidgets2::svalue(pickPC2_mod,index=TRUE)],
                               pt.label,pos=4,cex=0.8)
                        }
                        
                      })
    
    gWidgets2::addSpring(gv)
    
    dum <- gWidgets2::getToolkitWidget(dlg)
    #dum$window$SetDecorations("all")
    dum$window$maximize()
    gWidgets2::visible(dlg)
    
  }
  
  
  #*************************************************************  
  apply_pls_model<-function(h,...)
    #This:
    # 1. loads pls model (confirm it is a pls model!)
    # 2. verify data compatibility based on selected data types
    # 3. apply pre-treatments to data
    # 4. apply pls models
    # 5. launch a modal window for viewing results.
  {
    
    
    ff<-""
    dlg <- gWidgets2::gbasicdialog(title="Select PLS model", width=500,handler=function(h,...){
      ff<<-lefile
      rm(lefile,envir = .GlobalEnv)
    })
    gWidgets2::size(dlg)<-c(500,500)
    gv <- gWidgets2::gvbox(container = dlg,fill=TRUE)
    btn <- gWidgets2::gbutton("Pick a PLS model file",cont=gv, handler=function(h,...){
      dum <- gWidgets2::getToolkitWidget(dlg)
      dum$window$hide()
      lefile <<- gWidgets2::gfile("Pick a PLS model file",
                     filter=list("RData files" = list(patterns=c("*.RData"))))
      gWidgets2::insert(texte,lefile,do.newline = TRUE)
      gWidgets2::insert(texte,"\n")
      load(lefile, envir = .GlobalEnv)
      Encoding(model_descript$description) <- "UTF-8"
      gWidgets2::insert(texte,model_descript$description)
      gWidgets2::insert(texte,"\n")
      gWidgets2::insert(texte,model_descript)
      gWidgets2::insert(texte,"\n")
      gWidgets2::insert(texte,"Truncation limits: ")
      gWidgets2::insert(texte,toString(prepro_params$trunc_limits[,1]))
      gWidgets2::insert(texte,toString(prepro_params$trunc_limits[,2]))
      gWidgets2::insert(texte,"By spectra scaling selection (1-none, ...): ")
      gWidgets2::insert(texte,toString(prepro_params$byspectra_scaling_index))
      gWidgets2::insert(texte,"Center and BW for by value scaling: ")
      gWidgets2::insert(texte,toString(prepro_params$cntr_n_w[,1]))
      gWidgets2::insert(texte,toString(prepro_params$cntr_n_w[,2]))
      gWidgets2::insert(texte,"Flag for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$do_savgol))
      gWidgets2::insert(texte,"Derivative order for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$m))
      gWidgets2::insert(texte,"Polynomial order for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$p))
      gWidgets2::insert(texte,"Window size for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$w))
      dum$window$show()
    })
    
    texte <- gWidgets2::gtext("",container=gv,expand=TRUE)
    dum <- gWidgets2::visible(dlg)
    if (!dum) return()
    
    
    load(ff,envir = .GlobalEnv)
    Encoding(model_descript$description) <- "UTF-8"
    if (!(model_descript$type=="PLS")){
      gWidgets2::gmessage("Not a PLS model file!",title="Warning!",
               icon="warning")
      return()
    }
   
    
    #Check if data type needed by the model are available. If yes, build XData
    # and XData_p as in Make_XDatalist
    
    #Create a list of data types available in the current data set.
    #Then get data types.
    dum=as.list(levels(lesX[,1]))
    xType=get_DataType_Names(dum)  
    #Find indices of data types required by the model
    indix=pmatch(model_descript$datatype,xType)
    if (any(is.na(indix))){  #no match, cannot apply model
      gWidgets2::gmessage(title="ABORTING!",
              msg="Required data type not available in current data set.\nABORTING!")
      return()
    }
    doplot_lesX <<- FALSE
    gWidgets2::svalue(lesX,index=TRUE)<-indix
    doplot_lesX <<- TRUE
    
    
    #Do pretreatments☼
    Apply_PreTreatments(prepro_par = prepro_params)
    
   
     #Apply model
     if (model_descript$aggregation=="concatenate spectra"){
      #ATTN : do not work with XData_p but with the selected items.
      ind_lesX_list<-seq_len(length(model_descript$datatype))
      for (k in ind_lesX_list){
        spdf<-as.data.frame(XData_p[[k]][-1,])
        pre <- model_descript$datatype[k]
        colnames(spdf)<-paste(pre,as.character(XData_p[[k]][1,]),sep="_")
        if (k==1){
          y=spdf
        }else
        {
          y <- cbind(y,spdf)
        }
      }
      pls_set <<- y
    }else
    {
     
    }
    
    
    
    plspreds<-predict(plsFit[[1]],newdata = pls_set,ncomp=pls_ncomp)
    
    
    #Show results
    #Build a modal GUI
    
    dat_tbl=as.data.frame(plspreds[,,1])
    ncol_Ys=ncol(Ys_df)
    dat_tbl=cbind(Ys_df_sub[,-c(ncol_Ys-1,ncol_Ys)],dat_tbl,NoSeq=seq_len(nrow(dat_tbl)))
    colnames(dat_tbl)[ncol_Ys-1]="Predictions"
    
    tmp_w=gWidgets2::gwindow(visible=FALSE,
                  title=paste("PLS on ",
                              paste(model_descript$datatype,collapse = ", "),
                              " with ",pls_ncomp," Lat. var.",sep=""))
    gg=gWidgets2::gvbox(container = tmp_w)
    btn_save_plspred=gWidgets2::gbutton("Save to file",container = gg,
                             handler=function(h,...){
                               outfile=choose.files(default=file.path(basename(gWidgets2::svalue(lefichierY)),"*.txt"),
                                                    caption="Define file name",
                                                    multi=FALSE, 
                                                    filters=Filters[c("txt")])
                               write.table(dat_tbl,
                                           file=outfile,
                                           sep="\t",
                                           dec=".",
                                           row.names = FALSE,
                                           quote=FALSE)
                             })
    gt=gWidgets2::gtable(dat_tbl,container = gg)
    gt$widget$setGridLines('both')
    gWidgets2::visible(tmp_w)<-TRUE
    
  }
  
  #*************************************************************  
  apply_plsda_model<-function(h,...)
    #This:
    # 1. loads plsda model (confirm it is a plsda model!)
    # 2. verify data compatibility based on selected data types
    # 3. apply pre-treatments to data
    # 4. apply plsda model
    # 5. launch a modal window for viewing results.
  {
    
    
    ff<-""
    dlg <- gWidgets2::gbasicdialog(title="Select PLSDA model", width=500,handler=function(h,...){
      ff<<-lefile
      rm(lefile,envir = .GlobalEnv)
    })
    gWidgets2::size(dlg)<-c(500,500)
    gv <- gWidgets2::gvbox(container = dlg,fill=TRUE)
    btn <- gWidgets2::gbutton("Pick a PLSDA model file",cont=gv, handler=function(h,...){
      dum <- gWidgets2::getToolkitWidget(dlg)
      dum$window$hide()
      lefile <<- gWidgets2::gfile("Pick a PLSDA model file",
                     filter=list("RData files" = list(patterns=c("*.RData"))))
      gWidgets2::insert(texte,lefile,do.newline = TRUE)
      gWidgets2::insert(texte,"\n")
      load(lefile, envir = .GlobalEnv)
      Encoding(model_descript$description) <- "UTF-8"
      gWidgets2::insert(texte,model_descript$description)
      gWidgets2::insert(texte,"\n")
      gWidgets2::insert(texte,model_descript)
      gWidgets2::insert(texte,"\n")
      gWidgets2::insert(texte,"Truncation limits: ")
      gWidgets2::insert(texte,toString(prepro_params$trunc_limits[,1]))
      gWidgets2::insert(texte,toString(prepro_params$trunc_limits[,2]))
      gWidgets2::insert(texte,"By spectra scaling selection (1-none, ...): ")
      gWidgets2::insert(texte,toString(prepro_params$byspectra_scaling_index))
      gWidgets2::insert(texte,"Center and BW for by value scaling: ")
      gWidgets2::insert(texte,toString(prepro_params$cntr_n_w[,1]))
      gWidgets2::insert(texte,toString(prepro_params$cntr_n_w[,2]))
      gWidgets2::insert(texte,"Flag for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$do_savgol))
      gWidgets2::insert(texte,"Derivative order for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$m))
      gWidgets2::insert(texte,"Polynomial order for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$p))
      gWidgets2::insert(texte,"Window size for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$w))
      dum$window$show()
    })
    
    texte <- gWidgets2::gtext("",container=gv,expand=TRUE)
    dum <- gWidgets2::visible(dlg)
    if (!dum) return()
    
    
    load(ff,envir = .GlobalEnv)
    Encoding(model_descript$description) <- "UTF-8"
    if (!(model_descript$type=="PLSDA")){
      gWidgets2::gmessage("Not a PLSDA model file!",title="Warning!",
               icon="warning")
      return()
    }
    
    
    #Check if data type needed by the model are available. If yes, build XData
    # and XData_p as in Make_XDatalist
    
    #Create a list of data types available in the current data set.
    #Then get data types.
    dum=as.list(levels(lesX[,1]))
    xType=get_DataType_Names(dum)  
    #Find indices of data types required by the model
    indix=pmatch(model_descript$datatype,xType)
    if (any(is.na(indix))){  #no match, cannot apply model
      gWidgets2::gmessage(title="ABORTING!",
               msg="Required data type not available in current data set.\nABORTING!")
      return()
    }
    doplot_lesX <<- FALSE
    gWidgets2::svalue(lesX,index=TRUE)<-indix
    doplot_lesX <<- TRUE
    
    #Do pretreatments☼
    Apply_PreTreatments(prepro_par = prepro_params)
    
    gWidgets2::blockHandlers(gc_plsda_aggreg)
    gWidgets2::svalue(aggregate_options)$'Aggregation operator: ' <- model_descript$aggregation
    gWidgets2::unblockHandlers(gc_plsda_aggreg)
    
    #Apply model
    if (model_descript$aggregation=="concatenate"){
      #ATTN : do not work with XData_p but with the selected items.
      ind_lesX_list<-seq_len(length(model_descript$datatype))
      for (k in ind_lesX_list){
        spdf<-as.data.frame(XData_p[[k]][-1,])
        pre <- model_descript$datatype[k]
        colnames(spdf)<-paste(pre,as.character(XData_p[[k]][1,]),sep="_")
        if (k==1){
          y=spdf
        }else
        {
          y <- cbind(y,spdf)
        }
      }
      plsda_set <<- list(cbind(Ys_df_sub[,1],y))
    }else
    {
      
      ind_lesX_list<-as.list(seq_len(length(model_descript$datatype)))
      plsda_set <<- lapply(ind_lesX_list, function(ii){
        y<-data.frame(Cl1=Ys_df_sub[,1] ,XData_p[[ii]][-1,])
        colnames(y)[-1]<-as.character(XData_p[[ii]][1,])
        return(y)
      })
    }
      
      plsda_probs <- Predict_plsda(plsdaFit,mydata=plsda_set,probs=TRUE)
      plsda_cl <- Predict_plsda(plsdaFit,mydata=plsda_set,probs=FALSE)
      
  #   #Show results
  #   #Build a modal GUI
  #  
     dat_tbl=as.data.frame(cbind(plsda_cl,plsda_probs))
     ncol_Ys=ncol(Ys_df)
     dat_tbl=cbind(Ys_df_sub[,-c(ncol_Ys-1,ncol_Ys)],dat_tbl,NoSeq=seq_len(nrow(dat_tbl)))
    

    tmp_w <- gWidgets2::gwindow(visible=FALSE,
                  title=paste("PLS on ",
                              paste(model_descript$datatype,collapse = ", "),
                              " with ",pls_ncomp," Lat. var.",sep=""))
    gg <- gWidgets2::gvbox(container = tmp_w)
    btn_save_plspred=gWidgets2::gbutton("Save to file",container = gg,
                             handler=function(h,...){
                               outfile=choose.files(default=file.path(basename(gWidgets2::svalue(lefichierY)),"*.txt"),
                                                    caption="Define file name",
                                                    multi=FALSE,
                                                    filters=Filters[c("txt")])
                               write.table(dat_tbl,
                                           file=outfile,
                                           sep="\t",
                                           dec=".",
                                           row.names = FALSE,
                                           quote=FALSE)
                             })
    gt <- gWidgets2::gtable(dat_tbl,container = gg)
    gt$widget$setGridLines('both')
    gWidgets2::visible(tmp_w)<-TRUE

   }
  
    
#Widgets on data_tab------------------------- 
  #/././././././././././././././././././././././
  #Widgets on data_tab
  #/././././././././././././././././././././././
#-------------------------
  #*************************************************************
  OpenYFile<-function(h,...)
  # Handler for lefichierY widget.
  # Open lefichierY file and load its content into the data_view widget.
  # Identifies associated spectral data file and load their names into the lesX widget.
  {
    
    #Place cursor at the end to make sure the file name is visible to the user
    RGtk2::gtkEditableSetPosition(lefichierY$widget,nchar(gWidgets2::svalue(lefichierY)))
    
    #Load and display Y data file
    if (h$action){ #TRUE when opening the file. FALSE when applying a subsetting
                   #operation or when defining a new factor.
      Ys_df<<-read.table(file=gWidgets2::svalue(lefichierY),header=TRUE,sep="\t",dec=".",
                         na.strings = "")
      #gWidgets2::add a column with sequence number in the data file to allow tracing back when
      #a subset of data is used in the GUI.
      Ys_df<<-cbind(Ys_df,seq_len(nrow(Ys_df)))
      names(Ys_df)[ncol(Ys_df)]<<-"NoSeq_File"
      #gWidgets2::add a column with sequence number in the GUI to facilitate sorting back to initial
      Ys_df<<-cbind(Ys_df,seq_len(nrow(Ys_df)))
      names(Ys_df)[ncol(Ys_df)]<<-"NoSeq"
      #Ys_df_sub is the data frame used in the GUI. It is different from Ys_df when user
      #has selected a subset of data to work on. Ys_df reflects the content of the data
      #files.
      Ys_df_sub <<- Ys_df
      subset_ind <<- rep(TRUE,nrow(Ys_df))
      #Activate buttons that need Ydata
      lapply(btn_needing_Ydata,function(x) enabled(x)<-TRUE)
      gWidgets2::addHandlerChanged(nb,nb_handler)
    }else #Apply subsetting to the data frame
    {
      #Redefine NoSeq
      Ys_df_sub[,ncol(Ys_df)] <<- c(1:nrow(Ys_df_sub))
    }
      
    #Clean up previously loaded tables
    #remove columns from old one
    lescolonnes<-data_view$getColumns()
    sapply(seq_len(length(lescolonnes)), 
          function(i) data_view$removeColumn(lescolonnes[[i]]))
    
    #Start empty
    Y_data <- RGtk2::rGtkDataFrame()
    RGtk2::gSignalConnect(Y_data,"sort-column-changed",Plot_By_Factor)
    data_view$SetModel(Y_data)
    
    #Loading data frame in Y_Data
    Y_data$SetFrame(Ys_df_sub) 
    data_view$SetModel(Y_data)
    
    #Adding columns
    mapply(data_view$insertColumnWithAttributes,
           position=-1,
           title=colnames(Y_data),
           cell=list(RGtk2::gtkCellRendererText( )),
           text=seq_len(ncol(Y_data))-1
    )
    #Make all columns sortable
    sapply(seq_len(ncol(Y_data)),
           function(i) data_view$getColumn(i-1)$setSortColumnId(i-1))
    
    N_samples <- nrow(Ys_df_sub)
    
    if (h$action){  #needed only for a fresh start with a new Y file!
      #Identify associated spectral data files and display names
      lepath<-dirname(gWidgets2::svalue(lefichierY))
      DataDir<<-lepath
      # leY<-basename(gWidgets2::svalue(lefichierY))
      # dum<-dir(lepath,"*_I.txt")
      # ID<-substr(leY, 3, nchar(leY)-4)
      
      ynom <- gWidgets2::svalue(lefichierY)
      dname <- strsplit(tools::file_path_sans_ext(basename(ynom)),"_")[[1]]
      dname=paste(dname[-1],collapse="_")
      #Exclude Y file
      dum<-dir(lepath,paste0("*",dname,"*"))
      dum <- dum[!grepl("Y_",dum)]
      #dum<-dir(lepath,paste0("*",dname,"_I.txt"))
      inds<-grep(dname,dum)
      #For user to select X files to load
      dum = select.list(dum[inds],multiple=T,graphics = T)
      inds <- seq_along(dum)
      #Make sure 1rst column of Ys_df matches 1rst column of XData files.
      #If not, drop the corresponding XData file from the list of lesX and display a warning.
      nn<-length(inds)
      dumw<-gWidgets2::gwindow("Loading spectra files and checking compatibility",
                    width=500,height=30,parent=mymain)
      pg<-gWidgets2::gprogressbar(cont=dumw,value=0)
      
      if (N_samples>10){
        wspin<-MakeSpinner()
      }
      
      All_XData<<-vector("list",nn)
      for (k in 1:nn){
        gWidgets2::svalue(pg)<-floor(100*(k-1)/nn)
        dum1<-read.table(file.path(lepath,dum[inds[k]]),sep="\t", dec=".",header=FALSE)
        test1 <- length(dum1[-1,1])==length(Ys_df[,1])
        test2 <- TRUE
        if (test1) test2 <- any(as.character(dum1[-1,1])!=as.character(Ys_df[,1]))
        #cat("\n",length(dum1[-1,1])," ",length(Ys_df[,1])," ",test1," ",test2)
        if (test2){  #Not matching
          gmessage(paste("File ",dum[inds[k]], " is not valid! Will be ignored.",sep=""),
                   title="WARNING",icon="warning", parent=mymain)
          inds[k]<-0
        }else
          All_XData[[k]]<<-dum1
      }
      dispose(dumw)
      out_ind<-which(inds==0)
      if (length(out_ind)>0) inds<-inds[-out_ind]
      df<-as.data.frame(dum[inds])
      colnames(df)<-"X Data Files"
      assign("df",df,envir = .GlobalEnv)
      #Loading spectra files list in lesX and select first one
      lesX[]<-df
      
      ff <- gWidgets2::svalue(lefichierY)
      yf <- basename(ff)
      pp <- dirname(ff)
      valid_files <<- c(pp,yf,as.character(df[,1]))
      
      indi <- unlist(lapply(All_XData, function(x) !is.null(x)))
      All_XData<<-All_XData[indi]
      
      if (N_samples>10){
        dispose(wspin)
      }
      gWidgets2::svalue(lesX,index=TRUE) <- 1 
     
    }
     
    enabled(file_action$normby) <- TRUE
    enabled(file_action$splitat) <- TRUE
    
    #Clear graphics
    lapply(ggraphs,Clear_graph)
    
    #If there is no "BCK_Files" directory in DataDir,
    #make one and clone DataDir to this new one.
    cwd<-getwd()
    setwd(DataDir)
    dum<-list.dirs()
    if (!(length(grep("BCK_Files",dum))>0)){ #create backups!
      if (gWidgets2::gconfirm("Create backup directory?",icon="question",parent=mymain)){
        toDir<-paste(DataDir,"/BCK_Files_00",sep="")
        dir.create(toDir)
        system(paste('xcopy "', normalizePath(DataDir),'" "', normalizePath(toDir),'"',sep=""))
      }
    }
    setwd(cwd)
  }
  
  
  
  #******************************************************************************
  Find_duplicates <- function(h,...)
    # Handler for btn_find_duplicates
    # Prompt user for keeping 1rst or last items of duplicates
  {
    test=TRUE
    while (test){
      dum=gWidgets2::ginput("Keep first or last in a serie of duplictes? (F or L)", text="F",
                 icon="question",parent=mymain)
      dum=toupper(dum)[1]
      if ((dum=="F") | (dum=="L")) test=FALSE
    }
    fc=data_view$getModel()
    fc=fc[,1]
    if (dum=="F"){
      indd=which(duplicated(fc,fromLast = FALSE))
    }else
    {
      indd=which(duplicated(fc,fromLast = TRUE))
    }
    Y_selection$UnselectAll()
    Nindd=length(indd)
    if (Nindd>1){
      RGtk2::gSignalHandlerBlock(Y_selection,Y_selection_changed_Signal)
      for (k in indd[1:(Nindd-1)]){
        Y_selection$SelectPath(RGtk2::gtkTreePathNewFromIndices(k-1))
      }
      RGtk2::gSignalHandlerUnblock(Y_selection,Y_selection_changed_Signal)
    }
    if (Nindd>0) Y_selection$SelectPath(RGtk2::gtkTreePathNewFromIndices(indd[Nindd]-1)) 
    if (Nindd==0){
      Clear_graph(raw_data_graph)
    }
  }
  
  #******************************************************************************
  Define_Factor <- function(h,...)
    # Handler for btn_newfactor_data
    # Open a window for selecting factors defining interactions.
    # Also, defines the name of the newly created factor.
  {
    
    #Create dialog for user choice of factors to combine
    lesfacs <- names(Filter(is.factor,Ys_df))
    ind_lesfacs <- as.numeric(sapply(lesfacs, function(x) which(names(Ys_df)==x)))
    nfacs <- length(ind_lesfacs)
    dum<-list()
    newf_name<-as.character()
    to_file <- FALSE
    w<-gWidgets2::gbasicdialog(title='Define new factor from interactions', 
                    parent=mymain,
                    handler = function(h,...){
                      to_file <<- gWidgets2::svalue(save_2_file)
                      dum<-gWidgets2::svalue(lesfacteurs,index=TRUE)
                      newf_name<<-gWidgets2::svalue(mynewname)
                      dum<<-list(indices=ind_lesfacs[dum],name=newf_name)
                    },visible=FALSE, horizontal=FALSE)
    
    refer <- gdkWindowGetOrigin(mymain$widget$window)
    RGtk2::gdkWindowMove(w$widget$window,refer$x+50,refer$y+50)
    g <- gWidgets2::ggroup(container = w,horizontal=FALSE,fill=TRUE,expand=TRUE)
    save_2_file <- gWidgets2::gcheckbox("Save new factor to Y file (will erase existing file)",
                                        container = g)
    g1 <- gWidgets2::ggroup(container = g,horizontal=TRUE)
    gWidgets2::glabel('Name of new factor',container = g1)
    mynewname <- gWidgets2::gedit("", container=g1)
    g2 <- gWidgets2::ggroup(container = g, horizontal=FALSE,fill=TRUE,expand=TRUE)
    gWidgets2::glabel('Pick factors', container = g2)
    lesfacteurs<-gWidgets2::gtable(lesfacs,selected=-1,container = g2,fill=TRUE,expand=TRUE,multiple=TRUE)
    gWidgets2::size(lesfacteurs) <- c(200,300)
    selected_button <- gWidgets2::visible(w)
    
    if (selected_button){  #did not cancel, go ahead
      spin<-MakeSpinner()
      
      Ys_df[,dum$name]<<-droplevels(interaction(Ys_df[,dum$indices]))
      cols<-ncol(Ys_df)
      Ys_df <<- Ys_df[c(1:(cols-3),cols,(cols-2),(cols-1))]
      
      #Option to save Y file with new factor
     
      if (to_file) utils::write.table(Ys_df[,1:(cols-2)],file=svalue(lefichierY),sep="\t",row.names=FALSE)
      
      Ys_df_sub<<-Ys_df[subset_ind,]
      #To reload new Ys_df in GUI table
      OpenYFile(list(h=lefichierY,action=FALSE))
      
      # NO NEED TO REDEFINE XData, XDatalist and XData_p as this function changes the Y file only
      # Ys_df_sub<<-Ys_df[]
      # gWidgets2::gmessage("Subset selection has been cancelled!",parent=mymain)
      # subset_ind <<- rep(TRUE,nrow(Ys_df))
      # #To reload new Ys_df in GUI tabl
      # OpenYFile(list(h=lefichierY,action=FALSE))
      # #Apply subsetting to XData and XData_p
      # #Retrieve list of selected X data files.
      # selected<-gWidgets2::svalue(lesX)
      # XDatalist<<-as.list(as.character(selected))
      # inds<-gWidgets2::svalue(lesX,index = TRUE)
      # XData<<-All_XData[inds]
      # #Build XData and XData_p
      # if (length(selected)>0){     #only if a data type is selected
      #   XData <<- lapply(XData,data.matrix)
      #   XData <<- lapply(XData,subset,select=-1)  #Remove sample ID (first column) 
      #   XData <<- lapply(XData,function(x) apply(x,2,as.numeric))  #Make sure data are in numeric format
      #   XData <<- lapply(XData,function(x) x[c(TRUE,subset_ind),]) #TRUE required to preserve 1rst row (wl)
      #   XData_p <<- XData
      #   PreProDone <<- FALSE
      # }
      dispose(spin )
    }
  }
  
  
  #******************************************************************************
  Define_Subset <- function(h,...)
  # Handler for btn_subset_data
  # Open a window with tables for selecting subsets.
  # All factors are presented and selection of levels within all factors is allowed.
  # To be in the subset, a sample must meet selected levels in all factors. 
  # By default, all levels within a factor are selected.
  # If the combination of factors forms an empty subset, all samples are selected and 
  # a warning is issued
  {
    #Use gbasicdialog for a modal window allowing selection of levels in as 
    #many gtable as there are factors in Ys_df. Se Array_Tables_4_Selecting_levels.R
    #and save_pca_scores function in this script for hints.
    #Find factors in Y data: names, column index and total number of factors.
    lesfacs <- names(Filter(is.factor,Ys_df))
    ind_lesfacs <- as.numeric(sapply(lesfacs, function(x) which(names(Ys_df)==x)))
    nfacs <- length(ind_lesfacs)
    
    #Make modal GUI
    lestables <- list()
    dlg<-gWidgets2::gbasicdialog(title="Define subsetting criteria", parent=mymain,
                     handler=function(h,...){
                       sub_index <- rep(TRUE,nrow(Ys_df))
                       all_vals<-lapply(lestables,function (x) as.character(gWidgets2::svalue(x)))
                       for (k in 1:nfacs){
                         lesvals<-all_vals[[k]]
                         indd<-sapply(lesvals, function(x) Ys_df[,ind_lesfacs[k]]==x)
                         indd<-apply(indd,1,sum)
                         sub_index <- sub_index & indd
                       }
                       subset_ind <<- sub_index
                     })
    gWidgets2::size(dlg)=c(800,600)
    refer <- gdkWindowGetOrigin(mymain$widget$window)
    RGtk2::gdkWindowMove(dlg$widget$window,refer$x+50,refer$y+50)
    gg1<-gWidgets2::ggroup(cont=dlg,expand=TRUE,fill=TRUE)
    for (k in 1:nfacs){
      dum<-as.data.frame(levels(Ys_df[,ind_lesfacs[k]]))
      names(dum)=lesfacs[k]
      lestables[[k]] <- gWidgets2::gtable(dum,cont=gg1,fill=TRUE,expand=TRUE,multiple = TRUE)
      gWidgets2::size(lestables[k]) <- c(200,300)
      #Select all by default
      gWidgets2::svalue(lestables[[k]],index=TRUE)<-c(1:nrow(Ys_df))
    }
    selected_button <- gWidgets2::visible(dlg)
    
    if (selected_button){  #did not cancel
      spin<-MakeSpinner()
      
      if (!any(subset_ind)){
        subset_ind <<- rep(TRUE,nrow(Ys_df))
        utils::winDialog("ok","Empty subset resulted from selection. Selecting all samples.")
      }
      
      #Define Ys_df_sub
      Ys_df_sub<<-Ys_df[subset_ind,]
      Ys_df_sub<<-droplevels(Ys_df_sub)
      Ys_df_sub[,ncol(Ys_df_sub)]<<-seq_len(nrow(Ys_df_sub))
      
      #Apply subsetting to XData and XData_p
      #Retrieve list of selected X data files.
      selected<-gWidgets2::svalue(lesX)
      XDatalist<<-as.list(as.character(selected))
      inds<-gWidgets2::svalue(lesX,index = TRUE)
      XData<<-All_XData[inds]
      #Build XData and XData_p
      if (length(selected)>0){     #only if a data type is selected
        XData <<- lapply(XData,data.matrix)
        XData <<- lapply(XData,subset,select=-1)  #Remove sample ID (first column) 
        XData <<- lapply(XData,function(x) apply(x,2,as.numeric))  #Make sure data are in numeric format
        XData <<- lapply(XData,function(x) x[c(TRUE,subset_ind),]) #TRUE required to preserve 1rst row (wl)
        XData_p <<- XData
        Apply_PreTreatments()
        PreProDone <<- TRUE
      }
      
      #Calling OpenYFile with action=FALSE will skip file selection 
      #but will load the appropriate subset of data in the GUI environment.
      OpenYFile(list(h=lefichierY,action=FALSE))
      dispose(spin)  
    }
  }
  
  #******************************************************************************
  Make_XDatalist<-function(selection,...)
  # Handler for the lesX  widget.
  # Create lists for selected spectral data file. 
  # Two lists are created: XData that contains the data and XDatalist that contains  
  # the names of the selected spectral data files.
  # XData_p list is created. Here same as XData. It will be useful to hold preprocessed data.
  # N.B. in XData all data matrices have the same number of columns while in XData_p,
  # the number of columns can be unique to each data matrix.
  {
    #Clear plots of raw data
    Clear_graph(raw_data_graph)
    
    gWidgets2::enabled(selection$obj) <- FALSE  #Make sure handler has proceeded before changing selection in lesX
    
    
    N_samples<-nrow(Ys_df_sub)
    if (N_samples>500){
      wspin<-MakeSpinner()
    }
    #Retrieve list of selected X data files.
    selected<-gWidgets2::svalue(lesX)
    XDatalist<<-as.list(as.character(selected))
    inds<-gWidgets2::svalue(lesX,index = TRUE)
    XData<<-All_XData[inds]
    #Build XData and XData_p
    if (length(selected)>0){     #only if a data type is selected
      XData <<- lapply(XData,data.matrix)
      XData <<- lapply(XData,subset,select=-1)  #Remove sample ID (first column) 
      XData <<- lapply(XData,function(x) apply(x,2,as.numeric))  #Make sure data are in numeric format
      XData <<- lapply(XData,function(x) x[c(TRUE,subset_ind),]) #TRUE required to preserve 1rst row (wl)
      XData_p <<- XData
      
      #As spectral data files selection has changed, must update prepro_tab
      PreProDone <<- FALSE
      le_r<-lapply(XData,function(x) range(x[1,]))
      #Remove previous one
      gWidgets2::delete(gf_truncation,gf_truncation$children[[1]])
      gf_truncation<-build_truncation_widget(XDatalist,gf_truncation,le_r)
      gWidgets2::delete(gf1,gf1$children[[1]])
      gf1 <- build_byvalue_scaling_widget(XDatalist,gf1,le_r)
      gWidgets2::delete(gf_savgol,gf_savgol$children[[1]])
      gf_savgol <- build_savgol_widget(XDatalist,gf_savgol)
      
      #Adjust gspinbuttons for axis limits for plotting raw data
      leslimites <- get_graph_limits(XData_p)
      leslimites$ylimits <-  leslimites$ylimits
      steps <- diff(range(leslimites$ylimits))/20
      #Setup gspinbuttons and reset axis limits to default
      lapply(proplist,gWidgets2::blockHandlers)
      RGtk2::gtkAdjustmentSetUpper(ymaxi$widget$adjustment,leslimites$ylimits[2])
      RGtk2::gtkAdjustmentSetUpper(ymini$widget$adjustment,leslimites$ylimits[2]-steps)
      RGtk2::gtkSpinButtonSetIncrements(ymaxi$widget,steps,5*steps)
      RGtk2::gtkSpinButtonSetIncrements(ymini$widget,steps,5*steps)
      gWidgets2::svalue(proplist[[1]])<-0
      gWidgets2::svalue(proplist[[2]])<-0
      gWidgets2::svalue(proplist[[3]])<-0
      gWidgets2::svalue(proplist[[4]])<-1  #Dummy to make sure plot is updated when setting proplist[[4]] to 0.
      lapply(proplist,gWidgets2::unblockHandlers)
      gWidgets2::svalue(proplist[[4]])<-0
      
      
      #Retrieve row numbers for selected samples in preparation for updating plot on data_tab
      N_ech<<-nrow(XData[[1]])-1L
      selected_rows<-Y_selection$getSelectedRows()
      rows<-list()
      if (length(selected_rows$retval)>0) 
        rows <- sapply ( selected_rows$retval , RGtk2::gtkTreePathGetIndices ) + 1L
      
      #Proceed with plot.
      #N.B. doplot_lesX is a global variable that one can set to block plot updating
      #This can be useful if both new samples and new spectra types are selected at the
      #same time when dragging a rectangle on the plot. In this circumstance, plot don't 
      #need to be updated when data type selection changes.
      if (doplot_lesX){
        if (length(rows)>0){   #rebuild plot if samples are selected.
          subs<-data_view$model[rows,ncol(data_view$GetModel())]
          if (!PreProDone)
            flags <- BP_Spectra_matplot(XData,XDatalist,Ys_df_sub[,1],subs,
                               c(gWidgets2::svalue(xmini),gWidgets2::svalue(xmaxi)),
                               c(gWidgets2::svalue(ymini),gWidgets2::svalue(ymaxi)))
          else
            flags <- BP_Spectra_matplot(XData_p,XDatalist,Ys_df_sub[,1],subs,
                               c(gWidgets2::svalue(xmini),gWidgets2::svalue(xmaxi)),
                               c(gWidgets2::svalue(ymini),gWidgets2::svalue(ymaxi)),TRUE)
          if (flags$xdefault){
            gWidgets2::svalue(xmini)=0
            gWidgets2::svalue(xmaxi)=0
          }
          if (flags$ydefault){
            gWidgets2::svalue(ymini)=0
            gWidgets2::svalue(ymaxi)=0
          }
          PlotIsByFactor <<- FALSE
        }
      }
    }
    
    if (N_samples>500){
      dispose(wspin)
    }
    enabled(selection$obj) <- TRUE  #Allow changing selection in lesX
    if (!have_data[[1]]$widget$visible) lapply(have_data,function(x) x$widget$show())
  }
  
  
  #******************************************************************************
  update_rawX_plot <- function(h,...)
  #Handler for  Y_selection 
  #If used elsewhere for updating plots: in these cases, call is update_rawX_plot(Y_selection)
  {
    gWidgets2::visible(raw_data_graph) <- TRUE
    #Retrieve the row numbers for selected samples.
    selected_rows<-Y_selection$getSelectedRows()
    rows<-list()
    if (length(selected_rows$retval)>0) 
      (rows <- sapply ( selected_rows$retval , RGtk2::gtkTreePathGetIndices ) + 1L)
    if (length(rows)>0){ 
      #retrieve NoSeq (last column) from data_view as indices for parameter subs of BP_Spectra_matplot.R
      subs<-data_view$model[rows,ncol(data_view$GetModel())]
      if (!PreProDone)
        flags <- BP_Spectra_matplot(XData,XDatalist,Ys_df_sub[,1],subs,
                         c(gWidgets2::svalue(xmini),gWidgets2::svalue(xmaxi)),
                         c(gWidgets2::svalue(ymini),gWidgets2::svalue(ymaxi)))
      else
        flags <- BP_Spectra_matplot(XData_p,XDatalist,Ys_df_sub[,1],subs,
                        c(gWidgets2::svalue(xmini),gWidgets2::svalue(xmaxi)),
                        c(gWidgets2::svalue(ymini),gWidgets2::svalue(ymaxi)),TRUE)
      
      #Sets xmini, xmaxi, ymini and ymaxi according to flags (see BP_Spectra_matplot.R)
      if (flags$xdefault){
        gWidgets2::svalue(xmini)=0
        gWidgets2::svalue(xmaxi)=0
      }
      if (flags$ydefault){
        gWidgets2::svalue(ymini)=0
        gWidgets2::svalue(ymaxi)=0
      }
      #Acknowledge this is not a plot by factor graph
      PlotIsByFactor <<- FALSE
      #so user can delete samples
      enabled(btn_remove_samples) <- TRUE
    }else
    {
      #so user cannot delete samples
      enabled(btn_remove_samples) <- FALSE
    }
    return(FALSE)
  }
  
  
  
  #******************************************************************************
  Plot_By_Factor<-function(sorted)
  #Handler for Y_data "sort-column-changed" signal
  #When user selects a column in the Y data file displayed in the
  #Y_selection widget.
  {
    #Find the column
    raw_sort_column<<-sorted
    lacol<-sorted$GetSortColumnId()$sort.column.id+1L
    #If this column is for a factor, proceed
    if (!is.factor(Ys_df_sub[,lacol])){ #option to convert to factor.
      enviro <- environment()
      dlg<-gWidgets2::gbasicdialog(title="Convert to factor", parent=mymain,
                       handler=function(h,...){
                         use_num<-gWidgets2::svalue(gc)
                         nl<-as.numeric(gWidgets2::svalue(ge))
                         assign("use_num",use_num, enviro)
                         assign("nlev",nl, enviro)
                       })
      gv <- gWidgets2::gvbox(cont=dlg)
      gc <- gWidgets2::gcheckbox("Use unique numerical values as factor levels",container = gv)
      gg <- gWidgets2::ggroup(cont=gv)
      gl <- gWidgets2::glabel("  Custom number of levels\n  or 0 for no conversion.",container=gg)
      ge <- gWidgets2::gedit("0",container = gg,width = 2)
      
      out <- gWidgets2::visible(dlg)
      
      if (!out) return()
      
      if (use_num){
        Ys_df_sub[,lacol]<<-factor(Ys_df_sub[,lacol])
        Ys_df[,lacol]<<-factor(Ys_df[,lacol])
      }else
      {  
        if (nlev>0){
          Ys_df_sub[,lacol]<<-cut(Ys_df_sub[,lacol],nlev)
          Ys_df[,lacol]<<-cut(Ys_df[,lacol],nlev)
        }
      }
    }
    
    if (is.factor(Ys_df_sub[,lacol])){
      #Produce a "faked" XData list so can reuse BP_Spectra_matplot
      gWidgets2::visible(raw_data_graph) <- TRUE
      N_sets<-length(XData)
      
      X <- XData
      for (k in 1:N_sets){
        if (PreProDone){
          x<-XData_p[[k]]
          preT=TRUE
        }else
        {
          x<-XData[[k]]
          preT=FALSE
        }
        xp<-as.numeric(x[1,])
        x<-x[-1,]  #remove wavelength/number row
        #Compute mean values by factor level
        moys<-aggregate(x, by=list(Ys_df_sub[,lacol]), FUN=mean)
        #Build the faked XData as X
        X[[k]] <- rbind(xp,as.matrix(moys[,-1]))
        #Get the names of factor levels
        groupID <- as.character(moys[,1])
      }
      #Do plot
      subs<-c(1:length(groupID))
      flags <- BP_Spectra_matplot(X,XDatalist,groupID,subs,
                         c(gWidgets2::svalue(xmini),gWidgets2::svalue(xmaxi)),
                         c(gWidgets2::svalue(ymini),gWidgets2::svalue(ymaxi)),preT)
      #Sets xmini, xmaxi, ymini and ymaxi according to flags (see BP_Spectra_matplot.R)
      if (flags$xdefault){
        gWidgets2::svalue(xmini)=0
        gWidgets2::svalue(xmaxi)=0
      }
      if (flags$ydefault){
        gWidgets2::svalue(ymini)=0
        gWidgets2::svalue(ymaxi)=0
      }
      #Acknowledge that the current plot is by factor.
      PlotIsByFactor <<- TRUE
      isSelectMode <<- FALSE
      isZoomMode <<- TRUE
      #so user cannot delete samples
      enabled(btn_remove_samples) <- FALSE
    }else   #if not a factor, notify the user and do nothing. Column is already sorted
    {
      utils::winDialog("ok","This column is not a factor")
    }
  }
  
  
  #******************************************************************************
  do_scales_raw_plot<-function(h,...)
  # Handler for axis scale properties widgets. Pass plotting update 
  # to update_rawX_plot or Plot_By_Factor according to value of 
  # PlotIsByFactor
  {
    if (PlotIsByFactor){
      Plot_By_Factor(raw_sort_column)  #raw_sort_column is defined if PlotIsByFactor is TRUE
    }else
    {
      update_rawX_plot(Y_selection)
    }
    return(FALSE)
  }
  
  
  #******************************************************************************
  interact_w_spectra_plot <- function(h,...)
  # Manage user interaction with graphics
  # IF isSelectMode is TRUE (isZoomMode is false):  
  #   This finds sample from curve selected by mouse click or samples
  #   isolated with mouse drag over. Then
  #   samples are selected in the table of Ys. Finally, the plot is updated to 
  #   show selected samples for all data types.
  # IF isZoomMode is TRUE (isSelectMode is FALSE):
  #   One click zoom in by a factor of 1.5 around the clicked point.
  #   Dragging a rectangle - zoom to rectangle
  #   Zooming is achieved by interaction with the xlimits, ylimits widgets
  {
    if (isSelectMode){
      #DO NOT FIRE IF THE PLOT IS BY FACTOR LEVELS
      if (!PlotIsByFactor){
        if (identical(h$x[1],h$x[2]) && identical(h$y[1],h$y[2]) ){ #this is a click
          #First find the closest x coordinate.
          idx <- lapply(XData_p, function(x) which.min(abs(x[1,]-h$x[1])))
          #Then scan all XData matrices to find the one where the h$y[1] is the closest
          #Find selected samples
          selected_rows<-Y_selection$getSelectedRows()
          if (length(selected_rows$retval)<2) return()  #only one sample on the plot!
          #Note that in GTK world, indices are zero based!
          rows <- sapply(selected_rows$retval, RGtk2::gtkTreePathGetIndices )+1L
          #ATTN - Have to go back to original row ordering to work on XData
          #Use NoSeq in last column as indices!
          rows_ori <- data_view$model[rows,ncol(data_view$GetModel())]
          #Find which data source in XData_p holds minimum
          indi<-as.list(1:length(XData_p))
          lesminis <- lapply(indi,function(ii) abs(XData_p[[ii]][rows_ori+1,idx[[ii]]]-h$y[1]))
          lesminis <- lapply(lesminis, function(x) min(x))
          ind_lesX <- which.min(unlist(lesminis))
          #Find spectrum in selected member of XData
          lesminis <- pmin(abs(XData_p[[ind_lesX]][rows_ori+1,idx[[ind_lesX]]]-h$y[1]))
          ind_leY <- rows_ori[which.min(lesminis)]
          #Go back to TreeView row numbering
          ind_leY <- rows[rows_ori==ind_leY]
          
          #Change selection in tables.
          #N.B. setting selection in Y_selection triggers the 
          # Y_selection "changed" signal and this updates plot automatically
          #block handler to avoid useless plot update
          RGtk2::gSignalHandlerBlock(Y_selection,Y_selection_changed_Signal)
          doplot_lesX<<-FALSE   #to avoid plotting with changing selection in lesX
         
           #Selecting will update XData, XDatalist, XData_p and reset preprocessing options
          #lesX$set_selected(which(lesX$items[,1]==XDatalist[ind_lesX]) - 1L)
          
          RGtk2::gtkTreeSelectionUnselectAll(Y_selection)
          
          #unblocking handler
          RGtk2::gSignalHandlerUnblock(Y_selection,Y_selection_changed_Signal)
          doplot_lesX<<-TRUE   #to allow plotting with changing selection in lesX
          #Select sample at le_row. This triggers plot update!
          le_row <- RGtk2::gtkTreePathNewFromIndices(ind_leY-1)
          RGtk2::gtkTreeSelectionSelectPath(Y_selection,le_row)
        }else   #this is a rectangle
        {
          #First define a boolean vector to address range within a spectrum
          h$x[1] <- floor(h$x[1])
          h$x[2] <- ceiling(h$x[2])
          if (h$x[1]>=(h$x[2]-2)) h$x[2] <- h$x[1]+2
          #cat("\nx: ",unlist(h$x),"\n")
          #cat("y: ",unlist(h$y),"\n")
          idx<-lapply(XData_p,function(x) findInterval(x[1,],h$x,rightmost.closed = TRUE)==1) 
          #lapply(idx,function(x) cat("idx: ",which(x),"\n"))
          #Find selected samples
          selected_rows<-Y_selection$getSelectedRows()
          if (length(selected_rows$retval)<2) return()  #only one sample on plot!
          
          rows <- sapply(selected_rows$retval, RGtk2::gtkTreePathGetIndices )+1L
          
          #ATTN - Have to go back to original row ordering to work on XData
          rows_ori <- data_view$model[rows,ncol(data_view$GetModel())]
          
          
          indi<-as.list(1:length(XData_p))
          
          #NOT NECESSARY AS CASE OF A SINGLE SAMPLE ON PLOT HANDLED ABOVE!
          # if (length(rows_ori)>1){
          #   in_interval<-lapply(indi, function(ii){
          #     dum=as.matrix(XData_p[[ii]][rows_ori+1,idx[[ii]]])
          #     apply(dum,2,function(y) findInterval(y,h$y))
          #   }) 
          #   ind_lesY <- lapply(in_interval,function(x) apply(x,1,function(y) any(y==1)))
          # }else
          # {
          #   in_interval<-lapply(indi, function(ii) findInterval(XData_p[[ii]][rows_ori+1,idx[[ii]]],h$y))
          #   ind_lesY <- lapply(in_interval,function(x) any(x==1))
          # }
          
          in_interval<-lapply(indi, function(ii){
            dum=as.matrix(XData_p[[ii]][rows_ori+1,idx[[ii]]])
            apply(dum,2,function(y) findInterval(y,h$y))
          })
          ind_lesY <- lapply(in_interval,function(x) apply(x,1,function(y) any(y==1)))
          ind_lesX<-lapply(ind_lesY,"any")
          ind_lesX<-which(as.logical(ind_lesX))
          dum<-colSums(do.call(rbind,ind_lesY)) #do.call to cast list to matrix
          ind_lesY <- which(dum>0)
          if (length(ind_lesY)==0) return()
          
          #Change selection in tables.
          #N.B. setting selection in Y_selection triggers the 
          # Y_selection "changed" signal and this updates plot automatically
          #block handler to avoid useless plot update
          RGtk2::gSignalHandlerBlock(Y_selection,Y_selection_changed_Signal)
          doplot_lesX<<-FALSE   #to avoid plotting with changing selection in lesX
          # dum<-c()
          # for (k in 1:length(ind_lesX)){
          #   dum<-c(dum,which(lesX$items[,1]==XDatalist[ind_lesX[k]]))
          # }
          # #Selecting will update XData, XDatalist, XData_p and reset preprocessing options
          # lesX$set_selected(dum - 1L)
          
          RGtk2::gtkTreeSelectionUnselectAll(Y_selection)
          #Select all but the last ones to keep graph update to last selection
          le_row <<- sapply((rows[ind_lesY]-1),"gtkTreePathNewFromIndices")
          n_rows<-length(le_row)
          dummy<-le_row
          dummy[[n_rows]]<-NULL
          sapply(dummy,function(x) RGtk2::gtkTreeSelectionSelectPath(Y_selection,x))
          ddummy<<-dummy
          #unblocking handler
          RGtk2::gSignalHandlerUnblock(Y_selection,Y_selection_changed_Signal)
          doplot_lesX<<-TRUE   #to allow plotting with changing selection in lesX
          #Selecting last one will launch plot update.
          RGtk2::gtkTreeSelectionSelectPath(Y_selection,le_row[[n_rows]])
        }
      }
    }else   #This is zoom mode
    {
      leslimites<-list(xlimits=NULL,ylimits=NULL)
      #Get what is set in the interface
      leslimites$xlimits<-c(gWidgets2::svalue(xmini),gWidgets2::svalue(xmaxi))
      leslimites$ylimits<-c(gWidgets2::svalue(ymini),gWidgets2::svalue(ymaxi))
      #Get limits from all the data in the selected spectral data files 
      #Find selected samples
      if (!PlotIsByFactor){
        selected_rows<-Y_selection$getSelectedRows()
        rows <- sapply(selected_rows$retval, RGtk2::gtkTreePathGetIndices )+1L
        subs<-data_view$model[rows,ncol(data_view$GetModel())]
      }else
      {
        subs<-data_view$model[nrow(data_view$GetModel()),ncol(data_view$GetModel())]
      }
      XData_tmp<-lapply(XData,function(x) x[c(1,subs+1),])
      dat_limites<-get_graph_limits(XData_tmp)
      zoom_factor<-1.5
      if (identical(h$x[1],h$x[2]) && identical(h$y[1],h$y[2]) ){ #this is a click
        lerange<-diff(leslimites$xlimits)/2
        if (lerange>0){ #xlimits not to default
          leslimites$xlimits<-c(h$x[1]-lerange/zoom_factor, h$x[1]+lerange/zoom_factor)
        }else  #xlimits most likely to default (0,0)
        {
          lerange<-diff(dat_limites$xlimits)/2
          leslimites$xlimits<-c(h$x[1]-lerange/zoom_factor, h$x[1]+lerange/zoom_factor)
        }
        lerange<-diff(leslimites$ylimits)/2
        if (lerange>0){ #ylimits not to default
          leslimites$ylimits<-c(h$y[1]-lerange/zoom_factor, h$y[1]+lerange/zoom_factor)
        }else    #ylimits most likely to default (0,0)
        {
          lerange<-diff(dat_limites$ylimits)/2
          leslimites$ylimits<-c(h$y[1]-lerange/zoom_factor, h$y[1]+lerange/zoom_factor)
        }
        
      }else             #This is a rectangle
      {
        #Get limits from mouse action
        leslimites$xlimits<-c(h$x[1], h$x[2]+(diff(h$x)*0.15))
        leslimites$ylimits<-c(h$y[1], h$y[2])
      }
      lapply(prop_axis_limits_list,gWidgets2::blockHandlers)
      #update widgets giving axis limits
      gWidgets2::svalue(xmini)=leslimites$xlimits[1]
      gWidgets2::svalue(xmaxi)=leslimites$xlimits[2]
      gWidgets2::svalue(ymini)=leslimites$ylimits[1]
      lapply(prop_axis_limits_list,gWidgets2::unblockHandlers)
      gWidgets2::svalue(ymaxi)=leslimites$ylimits[2] #update plot
    }
    return(FALSE)
  }  
  
  
  #******************************************************************************
  remove_raw_samples <- function(h,...)
  #Handler for the btn_remove_samples widget
  {
    if (gWidgets2::gconfirm("Proceed to delete selected samples?",title="Confirm",
                 icon="question",parent=mymain)){
      #First, remove subsetting
      Ys_df_sub<<-Ys_df
      gWidgets2::gmessage("Subset selection has been cancelled!",parent=mymain)
      subset_ind <<- rep(TRUE,nrow(Ys_df))
      #Then ask if an additional backup of current data files is required and create one if desired
      if (gWidgets2::gconfirm("Backup current data files?",title="Confirm",
                   icon="question",parent=mymain)){
        #Retrieve name of last backup (i.e. its number)
        cwd<-getwd()
        setwd(DataDir)
        dum<-list.dirs()
        lesbcks<-grep("BCK_Files",dum)
        nbck<-sapply(lesbcks,function(y){
          l<-nchar(dum[[y]])
          as.numeric(substring(dum[[y]],(l-1)))
        } )
        maxbck<-max(nbck)
        toDir<-paste(DataDir,"/BCK_Files_",sprintf("%02d",maxbck+1),sep="")
        dir.create(toDir)
        system(paste('xcopy "', normalizePath(DataDir),'" "', normalizePath(toDir),'"',sep=""))
        setwd(cwd)
      }
      rows<-list()
      #Create a list of selected samples: selected_rows=h$getSelectedRows()
      selected_rows<-Y_selection$getSelectedRows()
      if (length(selected_rows$retval)>0) 
        (rows <- sapply ( selected_rows$retval , RGtk2::gtkTreePathGetIndices ) + 1L)
      #N.B rows are row numbers in Y_view. To get row number in data file, look at value in previous from 
      #last column where row order in the data files is stored: rows_ori
      rows_ori <- data_view$model[rows,(ncol(data_view$GetModel())-1L)]
      if (length(rows)>0){
        #Create a list of all files to modify
        lepath<-dirname(gWidgets2::svalue(lefichierY))
        ynom <- gWidgets2::svalue(lefichierY)
        dname <- strsplit(tools::file_path_sans_ext(basename(ynom)),"_")[[1]]
        dname=paste(dname[-1],collapse="_")
        #Exclude Y file
        dum<-dir(lepath,paste0("*",dname,"*"))
        dum <- dum[!grepl("Y_",dum)]
        inds<-grep(dname,dum)
        lesSpectraFiles <- dum[inds]
        
        #lesSpectraFiles<-dir(DataDir,"*_I.txt")
        lesFiles <- c(basename(gWidgets2::svalue(lefichierY)),lesSpectraFiles)
        turbifile <- dir(DataDir,"Turbi_*")
        lesFiles <- c(lesFiles,turbifile)
        lesFiles <- as.list(lesFiles)
        #Remove lines in files
        lapply(lesFiles,function(x)
          Remove_line_TextFile(file.path(DataDir,x),rows_ori+1))
       
        #Update All_XData, XData and XData_p
        indices <- rows_ori+1   #First row is the wavelength vector!
        All_XData <<- lapply(All_XData, function(x) x[-indices,])
        XData <<- lapply(XData, function(x) x[-indices,])
        XData_p <<- lapply(XData_p, function(x) x[-indices,])
        
        #Remove rows from Ys_df and Ys_df_sub
        indices <- rows_ori
        Ys_df <<- Ys_df[-indices,]
        #Adjust sequential numbers
        Ys_df[,(ncol(Ys_df)-1)] <<- seq_len(nrow(Ys_df))
        Ys_df[,ncol(Ys_df)] <<- seq_len(nrow(Ys_df))
        Ys_df_sub <<- Ys_df
        #Update subset_ind
        subset_ind <<- rep(TRUE,nrow(Ys_df))
        
        #Update table showing Ys in GUI
       
        Y_data <- RGtk2::rGtkDataFrame(Ys_df_sub)  #create a new rGtkDataFrame with updated data frame
        RGtk2::gSignalConnect(Y_data,"sort-column-changed",Plot_By_Factor)  #Reconnect to handler for sorting
        data_view$setModel(Y_data)  #set the model of the data_view GtkTreeView object.
        
        #Clear plot
        Clear_graph(raw_data_graph)
      }
    }  
    
  }
  
  #******************************************************************************
  #A handler for btn_merge_Ys button
  merge_Ys <- function(h,...)
  #Handler for btn_merge_Ys
  {
    #get current Y data file name
    Current_Yfile=gWidgets2::svalue(lefichierY)
    lepath=dirname(Current_Yfile)
    tomerge_file=utils::choose.files(default="*.txt",caption="Select file to merge",multi=FALSE,filters=Filters[c("txt","All"),])
    merge_df=read.table(file=tomerge_file,header=TRUE,sep="\t",dec=".")
    
    #Find similar column names in both data frame
    is.fact <- sapply(Ys_df, is.factor)
    wh=which(colnames(merge_df) %in% colnames(Ys_df))
    if (length(wh)==0){
      gWidgets2::gmessage("Data frames do not share a common column. ABORTED!",icon="warning",title="CANNOT MERGE")
      return()
    }
    if (length(wh)==1){
      mer_factors=as.list(colnames(merge_df))[[wh]]
    }else
    {
      mer_factors=colnames(merge_df[,wh])
    }
    
    #add dummy to Ys_df to maintain sort order
    
    #Open dialog to select columns used for merging and columns to merge
    ff<-""
    dlg = gWidgets2::gbasicdialog(title="Select data to merge", width=500,handler=function(h,...){
      ff<<-list(key_column=gWidgets2::svalue(gfl)$'Pick key column',tomerge=gWidgets2::svalue(gfl)$'Pick column(s) to add')
      rm(gfl,envir = .GlobalEnv)
    })
    gWidgets2::size(dlg)<-c(500,500)
    gv=gWidgets2::gvbox(container = dlg,fill=TRUE)
    gfl<<-gWidgets2::gformlayout(container = gv, spacing=15)
    gWidgets2::gtable(mer_factors,label="Pick key column",container = gfl,multiple=FALSE)
    gWidgets2::gtable(colnames(merge_df),label="Pick column(s) to add", container = gfl,multiple=TRUE)
    gWidgets2::size(gfl$children$`Pick key column`)=c(20,200)
    gfl$children$`Pick key column`$widget$setGridLines('both')
    gWidgets2::size(gfl$children$`Pick column(s) to gWidgets2::add`)=c(20,200)
    gfl$children$`Pick column(s) to add`$widget$setGridLines('both')
    
    
    
    out=gWidgets2::visible(dlg)
    tomerge_df = merge_df[c(unlist(ff[1]),unlist(ff[2]))]
    merged=merge(Ys_df,tomerge_df,by= unlist(ff[1]))
    
    merged=merged[order(merged$NoSeq_File),]
    drops <- c("NoSeq","NoSeq_File")
    merged=merged[ , !(names(merged) %in% drops)]
    
    #write the new data frame to file
    write.table(merged,file=Current_Yfile,quote=FALSE,sep="\t",dec=".",row.names = FALSE)
    
    #Reopen Y file
    OpenYFile(list(h=Current_Yfile,action=TRUE))
    
    invisible()
    
  }
  
  
  #******************************************************************************
  #A handler for "btn_aggregate_data" button
  # ASK FOR BACKUP
  Aggregate_Spectra <- function(h,...)
    #Handler for btn_aggregate_data
  {
    
    #Define new IDs before aggregating to cover cases of using some samples more than once
    #i.e. several series of repeated measurements on same sample.
    newY <- Ys_df
    
    #Get inputs from user
    lesfacs <- names(Ys_df)
    
    
    gr_var <- 1  #Main key is the sample ID in the first column of Ys_df
    
    
    span_var<-list()
    newf_name<-as.character()
    w<-gWidgets2::gbasicdialog(title='Select variable spanning over sample IDs', #parent=mymain,
                    handler = function(h,...){
                      span_var<<-gWidgets2::svalue(lesfacteurs,index=TRUE)
                    },visible=FALSE, horizontal=FALSE)
    g <- gWidgets2::ggroup(container = w,horizontal=FALSE,fill=TRUE,expand=TRUE)
    g2 <- gWidgets2::ggroup(container = g, horizontal=FALSE,fill=TRUE,expand=TRUE)
    gWidgets2::glabel('Pick a factor', container = g2)
    lesfacteurs<-gWidgets2::gtable(lesfacs,selected=-1,container = g2,fill=TRUE,expand=TRUE,multiple=FALSE)
    gWidgets2::size(lesfacteurs) <- c(200,300)
    gWidgets2::size(w) <- c(500,500)
    gWidgets2::visible(w)
    
    wspin<-MakeSpinner()
    
    if (length(span_var)>0){
      fc = interaction(as.factor(Ys_df[,gr_var]),as.factor(Ys_df[,span_var]),sep = "_")
    }else
    {
      fc = as.factor(Ys_df[,gr_var])
    }
    u_fc = unique(fc)
    dum <- character()
    for (k in (1:length(u_fc))){
      indi=which(fc==u_fc[k])
      for (i in indi)
         dum[i] <- paste0(newY[i,1],paste0("_",which(i==indi)))
    }
    ids <- colnames(newY)[1]
    newY <- newY[,-1]
    newY <- cbind(dum,newY)
    colnames(newY)[1] <- ids
    
    #Aggregate data in XData
    dum=All_XData
    #Update with new IDs in first column so it is consistent with newY
    dum=lapply(dum,function(x){ 
      x[1]=lapply(x[1],as.character)
      x[-1,1] <- as.character(newY[,1])
      return(x)
    })
    
    #Compute means and aggregate
    dum=lapply(dum,function(x) aggregate(x[-1,-1],by=list(x[-1,1]), FUN=mean))
    indi=as.list(1:length(dum))
    
    dum=lapply(indi, function(ii){
      dum[[ii]][1]=lapply(dum[[ii]][1],as.character)
      #Insert column headers
      colnames(dum[[ii]])[1] <- colnames(All_XData[[ii]])[1]
      dum[[ii]] <- rbind(All_XData[[ii]][1,],dum[[ii]])
    })
    
    #Apply proper name to 1rst column
    dum=lapply(dum,function(x){
      x[1]<-lapply(x[1],as.character)
      x[1,1] <- colnames(newY)[1]
      return(x)
    })
    
    #Keep first of a series in the newY file
    newY <- newY[!duplicated(newY[,1]),]
    #aggregate do ordering, so do the same to newY
    newY=newY[order(newY[,1]),]
    #Remove spanning variable
    newY <- newY[-span_var]
    #Remove last 2 columns (sequence numbering)
    nc <- ncol(newY)
    newY=newY[-c((nc-1),nc)]
    
    #save new data files. Keep names but append "_a" to data set name
    #Create subdirectory named "Aggregated"
    if (!dir.exists(file.path(valid_files[1],"Aggregated")))  dir.create(file.path(valid_files[1],"Aggregated"))
    #Yfile first
    Yname <- valid_files[2]
    dname <- strsplit(tools::file_path_sans_ext(Yname),"_")[[1]][2]
    dname <- gsub(dname,paste0(dname,"_a"),dname)
    Yname <- file.path(valid_files[1],"Aggregated",paste0("Y_",dname,".txt"))
    utils::write.table(newY,file=Yname,sep="\t",row.names=FALSE)
    #XDatafiles
    
    dum <- lapply(indi,function(ii){
      Xname <- valid_files[2+ii]
      fparts <- strsplit(tools::file_path_sans_ext(Xname),"_")[[1]]
      fparts[2] <- dname
      fic <- file.path(valid_files[1],"Aggregated",paste0(paste(fparts,collapse="_"),".txt"))
      utils::write.table(dum[[ii]],file=fic,sep="\t",row.names = FALSE, col.names = FALSE, append=FALSE)
    })
    gWidgets2::svalue(lefichierY) <- Yname
    #OpenYFile(list(h=lefichierY,action=TRUE))
    
    dispose(wspin)
  }
  
  #******************************************************************************
  #A handler for "Normalise by factor levels" DataTransform menu item
  Normalize_by <- function(h,...)
    #Handler for Normalize_by
  {
    #Select factor for normalising
    lesfacs <- names(Filter(is.factor,Ys_df))
    le_fac <- select.list(lesfacs,title="Pick a factor",graphics=T)
    norm_fac <- Ys_df[le_fac][,1]
    Nlevels <- nlevels(norm_fac)
    lesnoms <- get_DataType_Names(XDatalist)
    #Normalise
    dum1 <- lapply(seq_along(XData), function(i){
      dum=XData[[i]][-1,]
      
      #Compute mean spectrum by levels of factor
      m_4_cls=by(dum,norm_fac,colMeans)
      
      #Define wavelength range
      i1=ginput(paste0("Lower limit of wavelength range for ", lesnoms[i], " : "),
                title="Defining wavelength range",
                icon="question",
                text=range(XData[[i]][1,])[1])
      i1=as.numeric(i1)
      i1 <- which(XData[[i]][1,]>i1)[1]
      i2=ginput(paste0("Upper limit of wavelength range for", lesnoms[i], " : "),
                title="Defining wavelength range",
                icon="question",
                text=range(XData[[i]][1,])[2])
      i2=as.numeric(i2)
      i2 <- which(XData[[i]][1,]>i2)[1]-1
      
      #Average over wavelength range
      m_4_cls <- lapply(m_4_cls,function(x) mean(x[i1:i2]))
      for (k in 1:Nlevels){
        indi <- norm_fac==levels(norm_fac)[k]
        dum[indi,] <- dum[indi,]/m_4_cls[[k]]  
      }
      XData[[i]][-1,] <<- dum
      return(i)
    })
    XData_p <<- XData
    Apply_PreTreatments()
    #Clear plots of raw data
    Clear_graph(raw_data_graph)
    #Zoom all
    lapply(proplist,"blockHandlers")
    svalue(proplist[[1]])<-0
    svalue(proplist[[2]])<-0
    svalue(proplist[[3]])<-0
    svalue(proplist[[4]])<-1  #Dummy to make sure plot is updated when setting proplist[[4]] to 0.
    lapply(proplist,"unblockHandlers")
    svalue(proplist[[4]])<-0
    
    #make sure a sample is selected for updating plot
    selected_rows<-Y_selection$getSelectedRows()
    if (length(selected_rows$retval)==0){ 
      path = gtkTreePathNewFromIndices(1, -1)  #first row
      gtkTreeSelectionSelectPath(Y_selection, path)
    }
    update_rawX_plot(Y_selection)
  }
  #******************************************************************************
  # A handler for "Split spectra" DataTransform menu item
  # User selects a spectral data type and a wavelength (wavenumber) where to split
  # in two parts.
  # The split is volatile and dissapears when the user select a new Y file or quit
  # the GUI. On leaving the GUI, the splitted spectra are part of the XDatalist, 
  # All_XData, XData and XData_p left in the R global environment.
  Split_at_wv <- function(h,...)
  {
    #Select spectral data type to split
    le_type <- select.list(as.character(lesX[,1]),
                           title="Pick a data type",
                           multiple=F,
                           graphics=T)
    inds <- which(le_type==as.character(lesX[,1]))
    lesdats <<- All_XData[[inds]]
    wls <- lesdats[1,-1]
    wl_range <- range(wls)
    
    splitwl <- numeric()
    dlg <- gWidgets2::gbasicdialog("Select where to cut in 2 pieces", parent=mymain,
                        handler = function(h,...) {splitwl <<- svalue(sl)})
      g <- gWidgets2::ggroup(container = dlg, horizontal = TRUE)
      gWidgets2::glabel("Pick wavelength",cont = g,expang=F)
      sl <- gWidgets2::gslider(from=wl_range[1],
                    to=wl_range[2], by=1, value=mean(wl_range), 
                    cont=g, expand=TRUE, fill= T)
      gWidgets2::size(dlg) <- c(400,80)
      out <- gWidgets2::visible(dlg)
    i2 <- which(wls >= splitwl)[1]
    i3 <- length(wls)+1
    
    #update lesX
    splitchr <- strsplit(as.character(lesX[inds,1]),"_")
    name_lo <- paste0(splitchr[[1]][1],"LO_",paste(splitchr[[1]][-1],collapse = "_"))
    name_hi <- paste0(splitchr[[1]][1],"HI_",paste(splitchr[[1]][-1],collapse = "_"))
    
    lesnoms <- as.character(lesX[,1])
    N_old <- nrow(lesX)
    switch (as.character(N_old),
             "1" = lesnoms <- c(name_lo, name_hi),
             "2" = {  if (inds==1){
                       lesnoms <- c(name_lo, name_hi, lesnoms[2])
                      }else lesnoms <- c(lesnoms[1], name_lo, name_hi)
                  },
             if (inds==1){
               lesnoms <- c(name_lo, name_hi, lesnoms[-1])
             }else if (inds==N_old){
               lesnoms <- c(lesnoms[-N_old], name_lo, name_hi)
             }else
             {
               lesnoms <- c(lesnoms[1:(inds-1)], name_lo, name_hi, lesnoms[(inds+1):length(lesnoms)])
             }
    )
   
    df<-as.data.frame(lesnoms)
    colnames(df)<-"X Data Files"
    assign("df",df,envir = .GlobalEnv)
    #Loading spectra files list in lesX and select first one
    gWidgets2::blockHandlers(lesX)
    lesX[]<-df
    gWidgets2::unblockHandlers(lesX)
    
    #Updata All_XData
    datalo <- lesdats[,(1:i2)]
    datahi <- lesdats[,c(1,(i2+1):i3)]
    data_add <- list(datalo,datahi)
    switch (as.character(N_old),
            "1" = All_XData <<- data_add,
            "2" = { if (inds==1){
                      All_XData <<- c(data_add, All_XData[2])
                    }else All_XData <<- c(All_XData[1], data_add)
            },
            if (inds==1){
              All_XData <<- c(data_add, All_XData[-1])
            }else if (inds==N_old){
              All_XData <<- c(All_XData[-N_old], data_add)
            }else
            {
              All_XData <<- c(All_XData[1:(inds-1)], data_add, All_XData[(inds+1):length(lesnoms)])
            }
    )
    
    #Make selection in lesX to launch Make_XDatalist and this should do it!
    gWidgets2::svalue(lesX, index=TRUE) <- c(inds, inds+1)
    
  }
  
  #******************************************************************************
  addPop_2_raw_ggraphics <- function(g,proplist)
    # Adds popup menu to rigthclick event on ggraphics.
    # To be used in InSpectoR.
    # g is the ggraphics widget and proplist is a list of property widgets for the
    # graphic in g (xmin, ymin, ymin, ymax)
    # Items are:
    #   1. Switch to interactive data select mode.
    #   2. Switch to interactive zoom mode. 
    #   3. Copy to clipboard.
    #   4. Save as pdf. 
    #   5. Save as WMF. 
    #*********************************************************************************************************** 
  # B. Panneton, June 2017
  #*********************************************************************************************************** 
  {
    l <- list() 
    l$selectMode <- gaction("Select sample(s)","Switch to select mode",icon="gtk-find",
                            handler= function(h,...){
                              if (PlotIsByFactor){
                                galert("No select mode when by factor!",
                                       title="ATTENTION!",delay=2,parent=g)
                              }else
                              {
                                isSelectMode<<-TRUE
                                isZoomMode<<-FALSE
                              }
                            })
    l$zoomMode <- gaction("Zoom","Switch to zoom mode",icon="gtk-zoom",
                          handler = function(h,...){
                            isSelectMode <<- FALSE
                            isZoomMode <<- TRUE
                          })
    l$zoomAll <- gaction("Zoom all","Back to full zoom",icon="gtk-zoom-out",
                         handler= function(h,...){
                           #Have to access widgets for setting x and y limits
                           lapply(proplist,"blockHandlers")
                           svalue(proplist[[1]])<-0
                           svalue(proplist[[2]])<-0
                           svalue(proplist[[3]])<-0
                           svalue(proplist[[4]])<-1  #Dummy to make sure plot is updated when setting proplist[[4]] to 0.
                           lapply(proplist,"unblockHandlers")
                           svalue(proplist[[4]])<-0
                         })
    l$copyAction <- gaction("Copy", "Copy current graph to clipboard", icon="copy", 
                            handler=function(h,...) copyToClipboard(g)) 
    l$printActionPDF <- gaction("Save as PDF", "Save current graph", icon="save", 
                                handler=function(h,...) { 
                                  fname <- gfile(gettext("Filename without extension (pdf)"), type="save") 
                                  fname <- paste(fname,".pdf")
                                  if(!file.exists(fname) || gconfirm(gettext("Overwrite file?"))) 
                                    dev.copy2pdf(file=fname)
                                }) 
    l$printActionjpgg <- gaction("Save as WMF", "Save current graph", icon="save", 
                                 handler=function(h,...) { 
                                   fname <- gfile(gettext("Filename without extension (wmf)"), type="save") 
                                   fname <- paste(fname,".wmf")
                                   if(!file.exists(fname) || gconfirm(gettext("Overwrite file?"))) 
                                     dev.copy(win.metafile,file=fname)
                                   dev.off()}) 
    
    addRightclickPopupMenu(g,l)                           
  }
  
#Widgets on prepro_tab-------------------------
  #/././././././././././././././././././././././
  #Widgets on prepro_tab
  #/././././././././././././././././././././././
  
  #******************************************************************************
  save_prepro <- function(h,...)
    #Callback for saving prepro options to a RData file
  {
    if (!PreProDone) {  #needed to define prepro_params.
      Apply_PreTreatments()
      PreProDone <<- TRUE
    }
    ff <- gWidgets2::gfile("Define a file name",type="save")
    exten=tools::file_ext(ff)
    if (!(exten=="RData")) ff=paste0(tools::file_path_sans_ext(ff),".RData")
    #Get a short description from user
    description <- ""
    dlg=gWidgets2::gbasicdialog(title="User input", handler=function(h,...){
      description<<-gWidgets2::svalue(txt)
    })
    gv=gWidgets2::gvbox(container = dlg)
    lab=gWidgets2::glabel("Enter a short description for the model.",cont=gv)
    txt=gWidgets2::gtext(text="",width = 300,container = gv)
    gWidgets2::visible(dlg)
    
    dums <- MakeSpinner()
    
    model_descript=list(type="prepro",
                        description=description,
                        datatype=get_DataType_Names(XDatalist))
    
    save(model_descript, prepro_params, file=ff) 
    
    dispose(dums)
  }
  #******************************************************************************  
  load_prepro<-function(h,...){
    
    ff<-""
    dlg = gWidgets2::gbasicdialog(title="Select PCA model", width=500,handler=function(h,...){
      ff<<-lefile
      rm(lefile,envir = .GlobalEnv)
    })
    gWidgets2::size(dlg)<-c(500,500)
    gv=gWidgets2::gvbox(container = dlg,fill=TRUE)
    
    
    btn<-gWidgets2::gbutton("Pick a PrePro options file",cont=gv, handler=function(h,...){
      dum=getToolkitWidget(dlg)
      dum$window$hide()
      lefile<<-gWidgets2::gfile("Pick a PrePro options file",
                                filter=list("RData files" = list(patterns=c("*.RData"))))
      gWidgets2::insert(texte,lefile,do.newline = TRUE)
      gWidgets2::insert(texte,"\n")
      load(lefile, envir = .GlobalEnv)
      Encoding(model_descript$description) <- "UTF-8"
      gWidgets2::insert(texte,model_descript$description)
      gWidgets2::insert(texte,"\n")
      gWidgets2::insert(texte,model_descript)
      gWidgets2::insert(texte,"\n")
      gWidgets2::insert(texte,"Truncation limits: ")
      gWidgets2::insert(texte,toString(prepro_params$trunc_limits[,1]))
      gWidgets2::insert(texte,toString(prepro_params$trunc_limits[,2]))
      gWidgets2::insert(texte,"By spectra scaling selection (1-none, ...): ")
      gWidgets2::insert(texte,toString(prepro_params$byspectra_scaling_index))
      gWidgets2::insert(texte,"Center and BW for by value scaling: ")
      gWidgets2::insert(texte,toString(prepro_params$cntr_n_w[,1]))
      gWidgets2::insert(texte,toString(prepro_params$cntr_n_w[,2]))
      gWidgets2::insert(texte,"Flag for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$do_savgol))
      gWidgets2::insert(texte,"Derivative order for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$m))
      gWidgets2::insert(texte,"Polynomial order for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$p))
      gWidgets2::insert(texte,"Window size for SavGol: ")
      gWidgets2::insert(texte,toString(prepro_params$w))
      dum$window$show()
    })
    
    
    
    texte=gWidgets2::gtext("",container=gv,expand=TRUE)
    dum <- gWidgets2::visible(dlg)
    if (!dum) return()
    
    
    load(ff,envir = .GlobalEnv)
    Encoding(model_descript$description) <- "UTF-8"
    if (!(model_descript$type=="prepro")){
      gWidgets2::gmessage("Not a prepro options file!",title="Warning!",
                          icon="warning")
      return()
    }
    
    
    #Check if data type needed by the model are available. If yes, build XData
    # and XData_p as in Make_XDatalist
    
    #Create a list of data types available in the current data set.
    #Then get data types.
    dum=as.list(levels(lesX[,1]))
    xType=get_DataType_Names(dum)  
    #Find indices of data types required by the model
    indix=pmatch(model_descript$datatype,xType)
    if (any(is.na(indix))){  #no match, cannot apply model
      gWidgets2::gmessage(title="ABORTING!",
                          msg="Required data type not available in current data set.\nABORTING!")
      return()
    }
    
    #Select required data types
    doplot_lesX <<- FALSE
    gWidgets2::svalue(lesX,index=TRUE)<-indix
    doplot_lesX <<- TRUE
    
    #load Pretreatments options.
    Apply_PreTreatments(prepro_par = prepro_params)
    PreProDone <<- TRUE
    do_scales_raw_plot(xmini)
    
    
  }
  #*************************************************************
    
#-------------------------  
#Widgets on pca_tab-------------------------  
  #/././././././././././././././././././././././
  #Widgets on pca_tab
  #/././././././././././././././././././././././
  
#-------------------------
  
#******************************************************************************
  MonACP <- function(h,...)
  #Called when selecting a spectrum type.
  #Default to first 2 PCs
  {
    #Dummies
    pickPC1[]="1"
    pickPC2[]="1"
    #Maximum number of PCs
    leX <- gWidgets2::svalue(h$obj,index=TRUE)
    ncp<-lesNCPs[[leX]]
    #Populate pickPC1 and pickPC2 widgets and set defaults
    pickPC1[] <- as.character(1:ncp)
    pickPC2[] <- as.character(1:ncp)
    gWidgets2::svalue(pickPC1) <- 1
    gWidgets2::svalue(pickPC2) <- 2
    gWidgets2::svalue(pickACPPlot) <- "Scores"
    #Proceed to plot default score plot
    updateACPPlot(updateACPPlotbut)
    #Enable plot updating
    enabled(updateACPPlotbut)=TRUE
  }
  
  #******************************************************************************
  updateACPPlot <- function(h,...)
  #Callback for updating PC plot
  {
    
    gWidgets2::enabled(subset_from_ACPPlotbut) <- FALSE  #point selection cancelled!
    #ISR_env$SelectedScores <- NULL
    
    #Get spectra type, estimate nb of PCs. Retrieve PCA results from lesNCPs and lesACPS lists.
    leX <- gWidgets2::svalue(pick_data_4_PCA,index=TRUE)
    ncp<-lesNCPs[[leX]]
    gWidgets2::svalue(nCP_label) <- ncp
    resACP<-lesACPs[[leX]]
    var_exp=summary(resACP)$importance[2,]
    #pr_2_prin is for converting prcomp output to princomp output
    pca2<-pr_2_prin(resACP)
    #Compute score and orthogonal distances
    dd<-chemometrics::pcaDiagplot(XData_p[[leX]][-1,],pca2,a=as.numeric(gWidgets2::svalue(nCP_label)),
                                  plot=FALSE,scale=FALSE)
    
    #Retrieve options (factor to apply color or to label)
    indi_c<-gWidgets2::svalue(pickFactor_4_color,index=TRUE)
    wh_indi_c <- which(names(Ys_df_sub)==pickFactor_4_color[indi_c])
    indi_c<-indi_c-1   #first is none!
    indi_l<-gWidgets2::svalue(pickFactor_4_label,index=TRUE)
    wh_indi_l <- which(names(Ys_df_sub)==pickFactor_4_label[indi_l])
    indi_l<-indi_l-1   #first is none!
    
    
    "Plot"
    if (gWidgets2::svalue(pickACPPlot)=="Scores"){                    #Plot scores
     
      if (indi_l>0){                      #label selection
        
        pt.label<-Ys_df_sub[,wh_indi_l]
        labels_select <- gWidgets2::svalue(pickCrit_4_label)
        if (labels_select=="Leverage"){    #Apply label selection rule
           pt.label[dd$SDist<dd$critSD]=NA  #do not show those with small leverage
        }
        if (labels_select=="Orthogonal distance"){
          pt.label[dd$ODist<dd$critOD]=NA #do not show those with small orthogonal distance
        }
        if (labels_select=="Lev. & ortho"){
           pt.label[!((dd$ODist>=dd$critOD) & (dd$SDist>=dd$critSD))]=NA
        }
        
      }
      
      #Do plot
      gWidgets2::blockHandler(nb)
      gWidgets2::svalue(nb)<-index_acp
      gWidgets2::unblockHandler(nb)
      visible(ggacp) <- TRUE
      #To define ellipse drawing parameter for BP_Plot_Scores.R
      ell_option<-gWidgets2::svalue(pickCrit_4_ellipse)
      if (indi_c==0){   #no colors
        if (indi_l==0){  #no labels
          Sc_plot_params <<- BP_Plot_Scores(resACP$x[,gWidgets2::svalue(pickPC1,index=TRUE)],
                                            resACP$x[,gWidgets2::svalue(pickPC2,index=TRUE)],
                                            pcs=c(gWidgets2::svalue(pickPC1,index=TRUE),gWidgets2::svalue(pickPC2,index=TRUE)),
                                            frac=100*c(var_exp[gWidgets2::svalue(pickPC1,index=TRUE)],
                                                   var_exp[gWidgets2::svalue(pickPC2,index=TRUE)]),
                                            ell_option=ell_option)
        }else
        { 
          #no colors with labels
          Sc_plot_params <<- BP_Plot_Scores(resACP$x[,gWidgets2::svalue(pickPC1,index=TRUE)],
                                            resACP$x[,gWidgets2::svalue(pickPC2,index=TRUE)],
                                            pcs=c(gWidgets2::svalue(pickPC1,index=TRUE),gWidgets2::svalue(pickPC2,index=TRUE)),
                                            frac=100*c(var_exp[gWidgets2::svalue(pickPC1,index=TRUE)],
                                                   var_exp[gWidgets2::svalue(pickPC2,index=TRUE)]),
                                            ell_option=ell_option,
                                            pt.label=pt.label)
        }
      }else
      {
        colorby<-Ys_df_sub[,wh_indi_c]
        if (indi_l==0){   #with colors but no labels
          Sc_plot_params <<- BP_Plot_Scores(resACP$x[,gWidgets2::svalue(pickPC1,index=TRUE)],
                                            resACP$x[,gWidgets2::svalue(pickPC2,index=TRUE)],
                                            pcs=c(gWidgets2::svalue(pickPC1,index=TRUE),gWidgets2::svalue(pickPC2,index=TRUE)),
                                            frac=100*c(var_exp[gWidgets2::svalue(pickPC1,index=TRUE)],
                                                   var_exp[gWidgets2::svalue(pickPC2,index=TRUE)]),
                                            ell_option=ell_option,
                                            colorby=colorby)
        }else
        {
          #with colors and labels
          Sc_plot_params <<- BP_Plot_Scores(resACP$x[,gWidgets2::svalue(pickPC1,index=TRUE)],
                                            resACP$x[,gWidgets2::svalue(pickPC2,index=TRUE)],
                                            pcs=c(gWidgets2::svalue(pickPC1,index=TRUE),gWidgets2::svalue(pickPC2,index=TRUE)),
                                            frac=100*c(var_exp[gWidgets2::svalue(pickPC1,index=TRUE)],
                                                   var_exp[gWidgets2::svalue(pickPC2,index=TRUE)]),
                                            ell_option=ell_option,
                                            colorby=colorby, pt.label=pt.label)
        }
      }
    }
    #Loading plot
    if (gWidgets2::svalue(pickACPPlot)=="Loadings"){
      Clear_graph(ggacp)
      #Define colors
      mescols<-c("darkred","blue","green3","salmon","yellow3","black","red3","magenta","gray70","cyan") 
      meslty<-c("solid","88","26","2686","E6","44C4")
      plotbg="gray94"
      legbg="gray99"
      #Get number of first and last loadings
      ind1<-gWidgets2::svalue(pickPC1)
      ind2<-gWidgets2::svalue(pickPC2)
      #Extract required loadings and get the wl
      loads <- resACP$rotation[,ind1:ind2]
      wl <- t(t(XData_p[[leX]][1,]))
      #Compute offset to stack loadings one above the other.
      of <- diff(range(loads))
      of <- ceiling(of*10)/10
      if (ind1!=ind2){
        offs <- matrix((c(1:ncol(loads))-1)*of,ncol=ncol(loads),nrow=nrow(loads),byrow=TRUE)
        #recycling color if necessary
        recycling<-ceiling(ncol(offs)/10)
        mescols<-rep(mescols,recycling)
        nb_curves <- ncol(offs)
      }else
      {
        offs <- matrix(0,ncol=1,nrow=length(loads))
        nb_curves <- 1
      }
      
      #Prepare for plotting
      par(bg=plotbg)
      xrange<-range(XData_p[[leX]][1,])
      xrange[2]=xrange[2]+0.15*diff(xrange) #(room for legend)
      
      #Dummy plot for setting up nicely
      matplot(wl,loads+offs,type="n", bty="n",
                xlab="Wavelength or Wavenumber",ylab="Loadings",
              xlim=xrange)
      #Some decoration
      box(lwd=1,lty=1,col="gray35")
      grid(ny=NA,lwd=3,lty=1,col="white")
      
      
      
      #Plot all curves
      for (k in 1:nb_curves){
        matlines(c(wl[1],xrange[2]),rep(offs[1,k],2), lty=1,
                  col="white",lwd=3)
        matlines(wl,(offs+loads)[,k],type="l",
                   col=mescols[k],lwd=2)
      }
      
      #Plot legend
      letitre<-paste("Offset =",format(of,digits=3,nsmall=2))
      legend("topright", inset=0.01, legend=c(ind1:ind2), 
                bty='o', box.col = legbg, bg=legbg,
                lty=1,lwd=2,seg.len=4, col=mescols,
                cex=0.8,
                title=letitre)
    }
    #plot variance explained
    if (gWidgets2::svalue(pickACPPlot)=="Var. explained")   
    {
      Clear_graph(ggacp)
      i2=which(var_exp<0.005)[1]
      if (ncp>i2) i2=floor(1.2*ncp)
      indices=1:i2
      par(bg="gray90")
      plot(indices,var_exp[indices],type="l",col="Blue",
           lwd=2, panel.first = grid(lty=1 ,col="white"),
           xlab="PC number",
           ylab="Portion of variance explained") 
      abline(v=ncp,col="red")
      legend("topright",legend=c("Var. Expl.","0.001 threshold"),
             lty=c(1,1), col=c("blue","red"),inset=c(0.01,0.01))
    }
    #plot T2 - Q plot
    if (gWidgets2::svalue(pickACPPlot)=="Score and ortho distances")   
    {
      Clear_graph(ggacp)
      if (indi_c==0){
        if (indi_l==0){ #no label
          dum <- BP_Scatter(dd$SDist,dd$ODist,
                            labs=list("Score distance","Orthogonal distance"))
        }else   #labels
        {
          pt.label<-Ys_df_sub[,wh_indi_l]
          dum <- BP_Scatter(dd$SDist,dd$ODist,
                            labs=list("Score distance","Orthogonal distance"),
                            pt.label = pt.label)
        }
        
      }else
      {
        colorby<-Ys_df_sub[,wh_indi_c]
        if (indi_l==0){ #no label
          dum <- BP_Scatter(dd$SDist,dd$ODist,
                            labs=list("Score distance","Orthogonal distance"),
                            colorby=colorby)
        }else  #labels
        {
          pt.label<-Ys_df_sub[,wh_indi_l]
          dum <- BP_Scatter(dd$SDist,dd$ODist,
                            labs=list("Score distance","Orthogonal distance"),
                            colorby=colorby,
                            pt.label=pt.label)
        }
      }
      lines(dum$X_Lim,rep(dd$critOD,2),lty=2,col="red")
      lines(rep(dd$critSD,2),dum$Y_Lim,lty=2,col="red")
    }
  }
  
  #******************************************************************************
  save_pca_scores <- function(h,...)
  #Callback for saving computed PC scores to a text file
  {
    ff <- gWidgets2::gfile("Define a file name - This will be a csv text file",type="save",
                filter=list(".txt"=list(patterns=c("*.txt"))))
    if (nchar(tools::file_ext(ff))==0) ff=paste(ff,".txt",sep="")
    if (length(ff)>0){
      which_ones <- gWidgets2::ginput("Save scores for (A)ll selected spectra or (C)urrent one?",
                         text="A", icon="question")
      if (toupper(which_ones)=="C"){
        leX <- gWidgets2::svalue(pick_data_4_PCA,index=TRUE)
        resACP<-lesACPs[[leX]]
        nbCP <- utils::winDialogString("Select number of PCs",as.character(lesNCPs[[leX]]))
        scores<-resACP$x[,1:as.numeric(nbCP)]
        row.names(scores)=Ys_df_sub[,1]  #Put sample IDs as row names.
        dum<-as.character(XDatalist)
        labs<-substr(dum,start=c(1,1),regexpr("_",dum)-1)
        colnames(scores) <- sapply(colnames(scores),function(y) paste(labs[[1]],"-",y,sep=""))
        utils::write.csv(scores,file=ff)
      }else
      {
        #First inquire about the number of PCs to save for each spectra
        #Use a modal GUI providing estimated N CPs (lesNCPs) as default values.
        #Store the returned values to select the appropriate number of columns
        #in scs below.
        enviro <- environment()
        NCP_toSave<-as.numeric(lesNCPs)
        dlg<-gWidgets2::gbasicdialog(title="Select the number of Principal Components", 
                         handler=function(h,...){
                           valeurs<-numeric()
                           for (k in 1:N){
                             valeurs[k] <- gWidgets2::svalue(val[[k]])
                           }
                           assign("NCP_toSave",valeurs, enviro)
                         })
        dum<-gWidgets2::gframe(cont=dlg,text="SELECT")
        
        lyt<-gWidgets2::glayout(cont=dum)
        lyt[1,1]<-"Data Source"
        lyt[1,2]<-"Nb of PCs"
        
        labs= get_DataType_Names(XDatalist)
        
        N<-length(labs)
        val<-list()
        
      
        #Layout gspinbuttons and row labels
        for (k in 1:N){
          ncp<-min(dim(XData[[k]]))-1  #max number of PCs
          lyt[k+1,1]<-labs[[k]]
          lyt[k+1,2]<-(val[[k]]=gWidgets2::gspinbutton(1,ncp,1,value=NCP_toSave[k],cont=lyt))
        }
        
        gWidgets2::visible(dlg)
        
        
        
        scores<-list()
        for (k in 1:length(XDatalist)){
          resACP<-lesACPs[[k]]
          scs<-resACP$x[,1:NCP_toSave[k]]
          #
          colnames(scs) <- sapply(colnames(scs),function(y) paste(labs[[k]],"-",y,sep=""))
          scores<-cbind(scores,scs)
        }
        row.names(scores)=Ys_df_sub[,1]  #Put sample IDs as row names.
        utils::write.csv(scores,file=ff)
      }
    }
  }
  
  #******************************************************************************
  save_pca_model <- function(h,...)
    #Callback for saving computed PCA model to a RData file
  {
    ff <- gWidgets2::gfile("Define a file name",type="save")
    exten=tools::file_ext(ff)
    if (!(exten=="RData")) ff=paste0(tools::file_path_sans_ext(ff),".RData")
    #Get a short description from user
    description<-""
    dlg=gWidgets2::gbasicdialog(title="User input", handler=function(h,...){
      description<<-gWidgets2::svalue(txt)
    })
    gv=gWidgets2::gvbox(container = dlg)
    lab=gWidgets2::glabel("Enter a short description for the model.",cont=gv)
    txt=gWidgets2::gtext(text="",width = 300,container = gv)
    gWidgets2::visible(dlg)
    
    dums <- MakeSpinner()
    
    model_descript=list(type="PCA",
                        description=description,
                        datatype=get_DataType_Names(XDatalist))
    
    indi_c<-gWidgets2::svalue(pickFactor_4_color,index=TRUE)
    wh_indi_c <- which(names(Ys_df_sub)==pickFactor_4_color[indi_c])
    colorby<-Ys_df_sub[,wh_indi_c]
    
    save(model_descript,prepro_params,lesACPs,lesNCPs,colorby, file=ff) 
    
    dispose(dums)
  }
  
  #******************************************************************************
  interact_w_pca_plot <- function(h,...)
    # Manage user interaction with graphics
    # IF isSelectMode is TRUE (isZoomMode and IsZoomAll are false):  
    #   This finds samples isolated by dragging a rectangle with mouse. 
    #   Just clicking cancel selection.
    # IF isZoomMode is TRUE (isSelectMode and isZoomAll are FALSE):
    #   One click zoom in by a factor of 1.5 around the clicked point.
    #   Dragging a rectangle - zoom to rectangle
    #   Zooming by setting xlim and ylim using par(xaxp=c(min,max,n)) and par(yaxp=c(min,max,n))
    # Is isZoomAll is TRUE (isSelectMode and isZoomMode are false)
    #   Zooms to show all points and swith to isZoomMode=TRUE.
  {
    if (gWidgets2::svalue(pickACPPlot)=="Scores"){
      if (isSelectMode){
        #DO NOT FIRE IF THE PLOT IS BY FACTOR LEVELS
        if (identical(h$x[1],h$x[2]) && identical(h$y[1],h$y[2]) ){ #this is a click
          #cancel selection and redraw
          gWidgets2::enabled(subset_from_ACPPlotbut) <- FALSE
          ISR_env$SelectedScores <- NULL
          with(Sc_plot_params, 
              Sc_plot_params <<- BP_Plot_Scores(Xsc,Ysc,pcs=pcs,frac=frac,ell_option=ell_option,
                                                xlimits=X_Lim, ylimits=Y_Lim,
                                                colorby=colorby, pt.label=pt.label)
          )
        }else   #this is a rectangle
        {
          gWidgets2::enabled(subset_from_ACPPlotbut) <- TRUE
          #Find points in rectangle
          with (Sc_plot_params, {
                  indices <- which(Xsc>h$x[1] & Xsc<h$x[2]
                                   & Ysc>h$y[1] & Ysc<h$y[2])
                  #redraw the points with a small black dot
                  gWidgets2::visible(ggacp) <- TRUE
                  points(Xsc[indices],Ysc[indices],pch=20, cex=0.75, col="white")
                  ISR_env$SelectedScores <- c(indices,ISR_env$SelectedScores)
                  ISR_env$SelectedScores=sort(ISR_env$SelectedScores[!duplicated(ISR_env$SelectedScores)])
                  
                }
          )
          
          
        }
      }else if(isZoomMode)   #This is zoom mode
      {
        zoom_factor<-2
        
        if (identical(h$x[1],h$x[2]) && identical(h$y[1],h$y[2]) ){ #this is a click
          lerange=diff(Sc_plot_params$X_Lim)/2
          xlimits<-c(h$x[1]-lerange/zoom_factor, h$x[1]+lerange/zoom_factor)
          lerange<-diff(Sc_plot_params$Y_Lim)/2
          ylimits<-c(h$y[1]-lerange/zoom_factor, h$y[1]+lerange/zoom_factor)
        }else             #This is a rectangle
        {
          #Get limits from mouse action
          xlimits<-c(h$x[1], h$x[2]+(diff(h$x)*0.15))
          ylimits<-c(h$y[1], h$y[2])
        }
        
        #Set new limits
        gWidgets2::visible(ggacp) <- TRUE
        with(Sc_plot_params, 
             {Sc_plot_params <<- BP_Plot_Scores(Xsc,Ysc,pcs=pcs,frac=frac,ell_option=ell_option,
                            xlimits=xlimits, ylimits=ylimits,
                            colorby=colorby, pt.label=pt.label)
             indices <- ISR_env$SelectedScores
             if (!is.null(indices))
               points(Xsc[indices],Ysc[indices],pch=20, cex=0.75, col="white")
             })
        
      }else  #This is zoom all mode
      {
        visible(ggacp) <- TRUE
        with(Sc_plot_params,
             {Sc_plot_params <<- BP_Plot_Scores(Xsc,Ysc,pcs=pcs,frac=frac,ell_option=ell_option,
                                               colorby=colorby, pt.label=pt.label)
             indices <- ISR_env$SelectedScores
             if (!is.null(indices))
               points(Xsc[indices],Ysc[indices],pch=20, cex=0.75, col="white")
            }
        )
        isZoomAll <<- FALSE
        isZoomMode <<- TRUE
      }
      return(FALSE)
    }
  } 
  
  #******************************************************************************
  addPop_2_pca_ggraphics <- function(g)
    # Adds popup menu to rigthclick event on ggraphics.
    # To be used in InSpectoR - PCA tab.
    # g is the ggraphics widget 
    # Items are: 
    #   1. Copy to clipboard.
    #   2. Save as pdf. 
    #   3. Save as WMF. 
    #*********************************************************************************************************** 
    # B. Panneton, November 2018
    #*********************************************************************************************************** 
  {
    l <- list() 
    
    l$selectMode <- gaction("Select by rectangle (scores only)","Switch to select mode",icon="gtk-find",
                            handler= function(h,...){
                              isSelectMode<<-TRUE
                              isZoomMode<<-FALSE
                              isZoomAll<<-FALSE
                            })
    
    l$zoomMode <- gaction("Zoom (scores only)","Switch to zoom mode",icon="gtk-zoom",
                          handler = function(h,...){
                            isSelectMode <<- FALSE
                            isZoomMode <<- TRUE
                            isZoomAll<<-FALSE
                          })
    
    l$zoomAll <- gaction("Zoom all (scores only)",tooltip = "Right click to complete",icon="gtk-zoom-out",
                         handler= function(h,...){
                           isSelectMode <<- FALSE
                           isZoomMode <<- FALSE
                           isZoomAll<<-TRUE
                           interact_w_pca_plot(ggacp)
                         })
    
    
    l$copyAction <- gaction("Copy", "Copy current graph to clipboard", icon="copy", 
                            handler=function(h,...) copyToClipboard(g)) 
    l$printActionPDF <- gaction("Save as PDF", "Save current graph", icon="save", 
                                handler=function(h,...) { 
                                  fname <- gfile(gettext("Filename without extension (pdf)"), type="save") 
                                  fname <- paste(fname,".pdf")
                                  if(!file.exists(fname) || gconfirm(gettext("Overwrite file?"))) 
                                    dev.copy2pdf(file=fname)
                                }) 
    l$printActionjpgg <- gaction("Save as WMF", "Save current graph", icon="save", 
                                 handler=function(h,...) { 
                                   fname <- gfile(gettext("Filename without extension (wmf)"), type="save") 
                                   fname <- paste(fname,".wmf")
                                   if(!file.exists(fname) || gconfirm(gettext("Overwrite file?"))) 
                                     dev.copy(win.metafile,file=fname)
                                   dev.off()}) 
    
    addRightclickPopupMenu(g,l)                           
  }                
  
#Widgets on plsda_tab-------------------------  
  #/././././././././././././././././././././././
  #Widgets on plsda_tab
  #/././././././././././././././././././././././
  
#-------------------------
  
  #******************************************************************************
  Compute_PLSDA <- function(h,...)
  #Handler for button "Compute model" on PLSDA tab.
  {
    gWidgets2::insert(plsda_outputtext,"\n----------------------------------")
    gWidgets2::insert(plsda_outputtext,"Defining training and test sets.\n")
    
    #Define training and testing sets
    #Retrieve factor to predict
    indi_f<-gWidgets2::svalue(pickFactor_4_plsda,index=TRUE)
    wh_indi_f <- which(names(Ys_df_sub)==pickFactor_4_plsda[indi_f])
    
    #Retrieve proportion of samples for training
    laprop <- as.numeric(gWidgets2::svalue(train_ctrl)$'Prop. of data for training: ')
    #ATTN : do not work with XData_p but with the selected items.
    plsda_inTrain <<- caret::createDataPartition(y=Ys_df_sub[,wh_indi_f],p=laprop, list=FALSE)
    
    if (gWidgets2::svalue(aggregate_options)$'Aggregation operator: '=="concatenate"){
      #ATTN : do not work with XData_p but with the selected items.
      ind_lesX <- gWidgets2::svalue(Xdat_4_plsda,index = TRUE)
      ind_lesX_list<-as.numeric(ind_lesX)
      y <- data.frame(Cl1=Ys_df_sub[,wh_indi_f])
      for (k in ind_lesX_list){
        spdf<-as.data.frame(XData_p[[k]][-1,])
        pre <- get_DataType_Names(XDatalist)[k]
        colnames(spdf)<-paste(pre,as.character(XData_p[[k]][1,]),sep="_")
        y <- cbind(y,spdf)
      }
      plsda_set <<- list(y)
    }else
    {
      ind_lesX <- gWidgets2::svalue(Xdat_4_plsda,index = TRUE)
      ind_lesX_list<-as.list(ind_lesX)
      plsda_set <<- lapply(ind_lesX_list, function(ii){
        y<-data.frame(Cl1=Ys_df_sub[,wh_indi_f] ,XData_p[[ii]][-1,])
        colnames(y)[-1]<-as.character(XData_p[[ii]][1,])
        return(y)
      })
    }
    
    
    
    
    lambdas<-lapply(ind_lesX_list, function(ii){
          return(XData_p[[ii]][1,])
        })
    training <- lapply(plsda_set, function(x) x[plsda_inTrain,])
    testing <- lapply (plsda_set, function(x) x[-plsda_inTrain,])


    
    
    #Training - may take some time
    ctrl<-caret::trainControl(method = gWidgets2::svalue(train_ctrl)$'Resampling method: ',
                       number=as.numeric(gWidgets2::svalue(train_ctrl)$'Number of folds: '),
                       repeats=as.numeric(gWidgets2::svalue(train_ctrl)$'Number of repetitions: '),
                       classProbs = TRUE,
                       returnData = TRUE,
                       allowParallel = TRUE)
    
    
    #Fitting model - This may take some time
    
    gWidgets2::insert(plsda_outputtext,"Training model - takes some time!\n")
    spin<-MakeSpinner()
    if (gWidgets2::svalue(train_options)$'PreProcessing: '=="None"){
      prepro<-NULL
    }else prepro<-gWidgets2::svalue(train_options)$'PreProcessing: ' 

    suppressWarnings(
      plsdaFit <<- lapply(training, function(x){
        caret::train(Cl1~.,
              data=x,
              method="pls",
              tuneLength=as.integer(gWidgets2::svalue(train_options)$'Nb of LVs (max): '),
              trControl=ctrl,
              preProc=prepro,
              probMethod=gWidgets2::svalue(train_options)$'Prediction method: ',
              metric=gWidgets2::svalue(train_options)$'Performance metric: ')
      })
    )
    
    dispose(spin)
    
    #Compute confusion matrices for train and test set for output
    N_modeles<-length(plsdaFit)
    val_pred_cl <-Predict_plsda(plsdaFit,probs=FALSE)
    val_confusionmat<-caret::confusionMatrix(data=val_pred_cl,reference = plsdaFit[[1]]$trainingData[,1])
    
    test_pred_cl <- Predict_plsda(plsdaFit,testing,probs=FALSE)
    test_confusionmat <- caret::confusionMatrix(data=test_pred_cl,reference = testing[[1]][,1])
    
    #Generate output to GUI console.
    lesnoms<-gWidgets2::svalue(Xdat_4_plsda)
    dumind<-as.list(1:length(plsdaFit)) 
    plsda_txt_output<<-lapply(dumind,function(x) utils::capture.output(print(plsdaFit[[x]])))
    if (gWidgets2::svalue(aggregate_options)$'Aggregation operator: '=="concatenate")
      lesnoms=as.list(paste(lesnoms,collapse=" + "))
    for (k in 1:length(plsdaFit)) plsda_txt_output[[k]][1]<<-paste("\n******************\nPLSDA on ",
                                                               lesnoms[[k]], sep = "")
    if (N_modeles>1){
      FUNC<-gWidgets2::svalue(aggregate_options)$'Aggregation operator: ' 
      plsda_txt_output <<- c(plsda_txt_output,
             paste("*******\nConfusion matrix on training data for models aggregated with ",FUNC,
                   ".\n",sep=""),
                           utils::capture.output(print(val_confusionmat)))
      plsda_txt_output <<- c(plsda_txt_output,
                             paste("*******\nConfusion matrix on test data for models aggregated with ",FUNC,
                                   ".\n",sep=""),
                             utils::capture.output(print(test_confusionmat)))
    }else
    {  plsda_txt_output <<- c(plsda_txt_output,
                             paste("*******\nConfusion matrix on training data for model on ",lesnoms[[1]],
                                   ".\n",sep=""),
                             utils::capture.output(print(val_confusionmat)))
        plsda_txt_output <<- c(plsda_txt_output,
                             paste("*******\nConfusion matrix on test data for model on ",lesnoms[[1]],
                                   ".\n",sep=""),
                             utils::capture.output(print(test_confusionmat)))
    }
    
    gWidgets2::insert(plsda_outputtext,"\n----------------------------------")
    lapply(plsda_txt_output,function(x) gWidgets2::insert(plsda_outputtext,x))
   
  }
  
  #******************************************************************************
  Predict_plsda <- function(plsdaFit,mydata=NULL,probs=TRUE)
  #Predicts classes or probabilities on mydata. If mydata is NULL, predicts on training
  #If probs is TRUE, will return a table of probabilities.
  {
    
    N_modeles<-length(plsdaFit)
    if (N_modeles>1){  #Need to aggregate results
      ind_apply<-as.list(1:N_modeles)
      Ps <- lapply(ind_apply, function(ii){
        if (is.null(mydata))
          predict(plsdaFit[[ii]]$finalModel, newdata=plsdaFit[[ii]]$trainingData[,-1], type="prob")
        else
          predict(plsdaFit[[ii]]$finalModel, newdata=mydata[[ii]][,-1], type="prob")
      })
      Ps<-lapply(Ps,drop)  #remove useless third dimension.
      pooled <- Ps[[1]] * NA 
      n <- nrow(pooled) 
      classes <- colnames(pooled) 
      
      FUNC<-gWidgets2::svalue(aggregate_options)$'Aggregation operator: ' 
      for(i in 1:ncol(pooled)) 
      { 
        tmp <- lapply(Ps, function(y, col) y[,col], col = i) 
        tmp <- do.call("rbind", tmp) 
        pooled[,i] <- apply(tmp, 2, function(x) do.call(FUNC,as.list(x))) 
      } 
      pooled <- t(apply(pooled, 1, function(x) x/sum(x))) #the probs combined should be normalized
      classes <- colnames(pooled) 
      val_pred_cl <- pooled
      if (!probs)
        val_pred_cl <- factor(classes[apply(pooled, 1, which.max)], levels = classes) 
    }else  #only one model - no aggregation.
    {
      if (is.null(mydata)){
        val_pred_cl <- predict(plsdaFit[[1]],newdata=plsdaFit[[1]]$trainingData[,-1]
                               , type="prob")
      }else
      {
        val_pred_cl <- predict(plsdaFit[[1]],newdata=mydata[[1]][,-1], type="prob")
      }
      classes<-colnames(val_pred_cl)
      if (!probs)
          val_pred_cl <- factor(classes[apply(val_pred_cl, 1, which.max)], levels = classes) 
    }
    return(val_pred_cl)
  }
  #******************************************************************************
  save_plsda_model <- function(h,...)
    #Callback for saving computed plsda model to a RData file
  {
    ff <- gWidgets2::gfile("Define a file name",type="save")
    exten=tools::file_ext(ff)
    if (!(exten=="RData")) ff=paste0(tools::file_path_sans_ext(ff),".RData")
    #Get a short description from user
    description<-""
    dlg=gbasicdialog(title="User input", handler=function(h,...){
      description<<-gWidgets2::svalue(txt)
    })
    gv=gvbox(container = dlg)
    lab=gWidgets2::glabel("Enter a short description for the model.",cont=gv)
    txt=gtext(text="",width = 300,container = gv)
    visible(dlg)
    
    dums <- MakeSpinner()
    
    model_descript=list(type="PLSDA",
                        description=description,
                        datatype=get_DataType_Names(XDatalist),
                        aggregation=gWidgets2::svalue(aggregate_options)$'Aggregation operator: ')
    pls_ncomp <- lapply(plsdaFit,function(x) x$bestTune$ncomp)
    
    
   
    save(model_descript,prepro_params,plsdaFit,pls_ncomp,file=ff)
    
    dispose(dums)
  }
  
  #******************************************************************************
  clear_plsdafit <- function()
  {
    Clear_graph(plsda_plots)
    plsdaFit <<- NULL
    gWidgets2::insert(plsda_outputtext,"\n***********************************\nMODEL CLEARED!")
    gWidgets2::insert(plsda_outputtext,"\n***********************************\n")
  }
  
  #******************************************************************************
  addPop_2_plsda_ggraphics <- function(g)
    # Adds popup menu to rigthclick event on ggraphics.
    # To be used in InSpectoR - pls tab.
    # g is the ggraphics widget 
    # Items are: 
    #   1. Copy to clipboard.
    #   2. Save as pdf. 
    #   3. Save as WMF. 
    #*********************************************************************************************************** 
    # B. Panneton, November 2018
    #*********************************************************************************************************** 
  {
    l <- list() 
    
    # l$selectMode <- gaction("Select by rectangle (scores only)","Switch to select mode",icon="gtk-find",
    #                         handler= function(h,...){
    #                           isSelectMode<<-TRUE
    #                           isZoomMode<<-FALSE
    #                           isZoomAll<<-FALSE
    #                         })
    
    # l$zoomMode <- gaction("Zoom (scores only)","Switch to zoom mode",icon="gtk-zoom",
    #                       handler = function(h,...){
    #                         isSelectMode <<- FALSE
    #                         isZoomMode <<- TRUE
    #                         isZoomAll<<-FALSE
    #                       })
    # 
    # l$zoomAll <- gaction("Zoom all (scores only)",tooltip = "Right click to complete",icon="gtk-zoom-out",
    #                      handler= function(h,...){
    #                        isSelectMode <<- FALSE
    #                        isZoomMode <<- FALSE
    #                        isZoomAll<<-TRUE
    #                       })
    
    
    l$copyAction <- gaction("Copy", "Copy current graph to clipboard", icon="copy", 
                            handler=function(h,...) copyToClipboard(g)) 
    l$printActionPDF <- gaction("Save as PDF", "Save current graph", icon="save", 
                                handler=function(h,...) { 
                                  fname <- gfile(gettext("Filename without extension (pdf)"), type="save") 
                                  fname <- paste(fname,".pdf")
                                  if(!file.exists(fname) || gconfirm(gettext("Overwrite file?"))) 
                                    dev.copy2pdf(file=fname)
                                }) 
    l$printActionjpgg <- gaction("Save as WMF", "Save current graph", icon="save", 
                                 handler=function(h,...) { 
                                   fname <- gfile(gettext("Filename without extension (wmf)"), type="save") 
                                   fname <- paste(fname,".wmf")
                                   if(!file.exists(fname) || gconfirm(gettext("Overwrite file?"))) 
                                     dev.copy(win.metafile,file=fname)
                                   dev.off()}) 
    
    addRightclickPopupMenu(g,l)                           
  }                
  
#Widgets on pls_tab------------------------- 
  #/././././././././././././././././././././././
  #Widgets on pls_tab
  #/././././././././././././././././././././././
  
#-------------------------
  
  #****************************************************************************
  Compute_PLS <- function(h,...)
    #Handler for button "Compute model" on PLS tab.
  {
    gWidgets2::insert(pls_outputtext,"\n----------------------------------")
    gWidgets2::insert(pls_outputtext,"Defining training and test sets.\n")
    
    #Retrieve variable to predict
    indi_f<-gWidgets2::svalue(pickvar_4_pls,index=TRUE)
    wh_indi_f <- which(names(Ys_df_sub)==pickvar_4_pls[indi_f])
    
    if (gWidgets2::svalue(pls_aggregate_options)$'Aggregation method: '=="concatenate spectra"){
      #ATTN : do not work with XData_p but with the selected items.
      ind_lesX <- gWidgets2::svalue(Xdat_4_pls,index = TRUE)
      ind_lesX_list<-as.numeric(ind_lesX)
      y <- data.frame(V1=Ys_df_sub[,wh_indi_f])
      for (k in ind_lesX_list){
        spdf<-as.data.frame(XData_p[[k]][-1,])
        pre <- get_DataType_Names(XDatalist)[k]
        colnames(spdf)<-paste(pre,as.character(XData_p[[k]][1,]),sep="_")
        y <- cbind(y,spdf)
      }
      pls_set <<- list(y)
    }else
    {
      #ATTN : do not work with XData_p but with the selected items.
      ind_lesX <- gWidgets2::svalue(Xdat_4_pls,index = TRUE)
      ind_lesX_list<-as.list(ind_lesX)
      pls_set <<- lapply(ind_lesX_list, function(ii){
        y<-data.frame(V1=Ys_df_sub[,wh_indi_f] ,XData_p[[ii]][-1,])
        colnames(y)[-1]<-as.character(XData_p[[ii]][1,])
        return(y)
      })
    }
    
    
    lambdas<-lapply(ind_lesX_list, function(ii){
      return(XData_p[[ii]][1,])
    })
    
   
    gWidgets2::insert(pls_outputtext,"Training model - takes some time!\n")
    spin<-MakeSpinner()
    
    #Training - may take some time
    plsFit <<- lapply(pls_set,function(x){
      pls::plsr(formula= V1~., data=x,
           ncomp=as.numeric(gWidgets2::svalue(pls_ctrl)$'Nb of LVs (max): '),
           validation = gWidgets2::svalue(pls_ctrl)$'Resampling method: ',
           scale=gWidgets2::svalue(pls_ctrl)$'Scaling: '=="scale")
    }
    )
    
    if (!(gWidgets2::svalue(pls_ctrl)$'Resampling method: '=="none")){
        plsFit_ncomp <<- lapply(plsFit, function(x) pls::selectNcomp(x,"randomization",plot=FALSE))
    }else
    {
      plsFit_ncomp <<- lapply(plsFit, function(x) NULL)
    }
    dispose(spin)
    N_modeles=length(plsFit)
    
    #Generate output to GUI console.
    lesnoms<-gWidgets2::svalue(Xdat_4_pls)
    dumind<-as.list(1:length(plsFit)) 
    pls_txt_output<<-lapply(dumind,function(x) utils::capture.output(print(plsFit[[x]])))
    
    if (N_modeles>1){
      for (k in 1:length(plsFit)) pls_txt_output[[k]][1]<<-paste("\n******************\nPLS on ",
                                                                 lesnoms[[k]], sep = "")
    }else
    { 
      pls_txt_output[[1]][1]<<-paste("\n******************\nPLS on ",
                                     paste(lesnoms,collapse = ", "), sep = "")
      pls_txt_output <<- c(pls_txt_output,
                             paste("*******\nSummary for pls model on ",paste(lesnoms,collapse = ", "),
                                   ".\n",sep=""),
                             utils::capture.output(summary(plsFit[[1]])),
                             paste("\nSelected number of LVs: ",plsFit_ncomp[[1]],sep=""),
                             "\nRMSEP summary:",
                             utils::capture.output(pls::RMSEP(plsFit[[1]],estimate="all")),
                             "\nR2 summary:",
                             utils::capture.output(pls::R2(plsFit[[1]],estimate="all"))
                            )
    }
    
    gWidgets2::insert(pls_outputtext,"\n----------------------------------")
    lapply(pls_txt_output,function(x) gWidgets2::insert(pls_outputtext,x))
  }
    
    
  #****************************************************************************
  Predict_pls <- function(plsFit,mydata=NULL,probs=TRUE)
    #Predicts classes or probabilities on mydata. If mydata is NULL, predicts on       training dataEllipse()  
    #If probs is TRUE, will return a table of probabilities.
  {
    
    N_modeles<-length(plsFit)
    if (N_modeles>1){  #Need to aggregate results
      
    }else  #only one model - no aggregation.
    {
      
    }
  }
  
  
  #******************************************************************************
  save_pls_model <- function(h,...)
    #Callback for saving computed pls model to a RData file
  {
    ff <- gWidgets2::gfile("Define a file name",type="save")
    exten=tools::file_ext(ff)
    if (!(exten=="RData")) ff=paste0(tools::file_path_sans_ext(ff),".RData")
    #Get a short description from user
    description<-""
    dlg=gWidgets2::gbasicdialog(title="User input", handler=function(h,...){
      description<<-gWidgets2::svalue(txt)
    })
    gv=gWidgets2::gvbox(container = dlg)
    lab=gWidgets2::glabel("Enter a short description for the model.",cont=gv)
    txt=gWidgets2::gtext(text="",width = 300,container = gv)
    gWidgets2::visible(dlg)
    
    dums <- MakeSpinner()
    
    model_descript=list(type="PLS",
                        description=description,
                        datatype=get_DataType_Names(XDatalist),
                        aggregation=gWidgets2::svalue(pls_aggregate_options)$'Aggregation method: ')
    pls_ncomp <- gWidgets2::svalue(pls_user_ncomp)
    
    
    indi_c<-gWidgets2::svalue(pls_pickFactor_4_color,index=TRUE)
    wh_indi_c <- which(names(Ys_df_sub)==pls_pickFactor_4_color[indi_c])
    indi_c<-indi_c-1   #first is none!
    colorby<-Ys_df_sub[,wh_indi_c]
    
    save(model_descript,prepro_params,plsFit,pls_ncomp,colorby,file=ff) 
    
    dispose(dums)
    
  }
  #******************************************************************************
  plot_PLS_preds <- function(h,...)
  {
    #Retrieve options (factor to apply color or to label)
    indi_c<-gWidgets2::svalue(pls_pickFactor_4_color,index=TRUE)
    wh_indi_c <- which(names(Ys_df_sub)==pls_pickFactor_4_color[indi_c])
    indi_c<-indi_c-1   #first is none!
    indi_l<-gWidgets2::svalue(pls_pickFactor_4_label,index=TRUE)
    wh_indi_l <- which(names(Ys_df_sub)==pls_pickFactor_4_label[indi_l])
    indi_l<-indi_l-1   #first is none!
    
    lewhich=gWidgets2::svalue(plspred_which)
    if ((gWidgets2::svalue(pls_ctrl)$'Resampling method: '=="none") & (lewhich=="validation"))
    {
      gmessage("Invalid option when no validation!\nWill plot for training.")
      lewhich="train"
    }
    pls_plots$widget$window$hide()
    pl<-plot(plsFit[[1]],plottype = "prediction",
             ncomp=gWidgets2::svalue(pls_user_ncomp),
             which=lewhich, type="n")
    pls_plots$widget$window$show()
    gWidgets2::visible(pls_plots) <- TRUE
    pl=as.data.frame(pl)
    
    letitre=paste(lewhich," - ",gWidgets2::svalue(pls_user_ncomp)," L.V.",sep="")
    
    if ((length(wh_indi_c)==0) & (length(wh_indi_l)==0)){
      PLS_scatter_params <<- BP_Scatter(pl$measured,pl$predicted, 
                                        titre=letitre, labs=list("Measured","Predicted"))
    }
    if (!(length(wh_indi_c)==0) & (length(wh_indi_l)==0)){
        colorby<-Ys_df_sub[,wh_indi_c]
        PLS_scatter_params <<- BP_Scatter(pl$measured,pl$predicted,
                                        colorby=colorby, titre=letitre,
                                        labs=list("Measured","Predicted"))
    }
    if (!(length(wh_indi_c)==0) & (!length(wh_indi_l)==0)){  
        colorby<-Ys_df_sub[,wh_indi_c]
        PLS_scatter_params <<- BP_Scatter(pl$measured,pl$predicted,
                                          colorby=colorby,
                                          pt.label=Ys_df_sub[,wh_indi_l],
                                          titre=letitre,
                                          labs=list("Measured","Predicted"))
    }
    if ((length(wh_indi_c)==0) & (!length(wh_indi_l)==0)){  
      PLS_scatter_params <<- BP_Scatter(pl$measured,pl$predicted,
                                        pt.label=Ys_df_sub[,wh_indi_l],
                                        titre=letitre,
                                        labs=list("Measured","Predicted"))
    
    } 
    
    gWidgets2::svalue(pls_plot_type) <- "Prediction"
  }
  
  #******************************************************************************
  interact_w_pls_plot <- function(h,...)
    # Manage user interaction with graphics
    # IF isZoomMode is TRUE (isSelectMode and isZoomAll are FALSE):
    #   One click zoom in by a factor of 1.5 around the clicked point.
    #   Dragging a rectangle - zoom to rectangle
    #   Zooming by setting xlim and ylim using par(xaxp=c(min,max,n)) and par(yaxp=c(min,max,n))
    # Is isZoomAll is TRUE (isSelectMode and isZoomMode are false)
    #   Zooms to show all points and swith to isZoomMode=TRUE.
  {
    
    if (gWidgets2::svalue(pls_plot_type) == "Prediction"){
      if(isZoomMode)   #This is zoom mode
      {
        zoom_factor<-2
        
        if (identical(h$x[1],h$x[2]) && identical(h$y[1],h$y[2]) ){ #this is a click
          lerange=diff(PLS_scatter_params$X_Lim)/2
          xlimits<-c(h$x[1]-lerange/zoom_factor, h$x[1]+lerange/zoom_factor)
          lerange<-diff(PLS_scatter_params$Y_Lim)/2
          ylimits<-c(h$y[1]-lerange/zoom_factor, h$y[1]+lerange/zoom_factor)
        }else             #This is a rectangle
        {
          #Get limits from mouse action
          xlimits<-c(h$x[1], h$x[2]+(diff(h$x)*0.15))
          ylimits<-c(h$y[1], h$y[2])
        }
        
        #Set new limits
        gWidgets2::visible(pls_plots) <- TRUE
        with(PLS_scatter_params, 
             PLS_scatter_params <<- BP_Scatter(Xsc,Ysc,labs,titre,
                                               xlimits=xlimits, ylimits=ylimits,
                                               colorby=colorby, pt.label=pt.label)
        )
        
      }else  #This is zoom all mode
      {
        gWidgets2::visible(pls_plots) <- TRUE
        with(PLS_scatter_params,
             PLS_scatter_params <<- BP_Scatter(Xsc,Ysc,labs,titre,
                                               colorby=colorby, pt.label=pt.label)
        )
        isZoomAll <<- FALSE
        isZoomMode <<- TRUE
      }
      return(FALSE)
    }
  }  
  
  #******************************************************************************
  clear_plsfit <- function()
  {
    Clear_graph(pls_plots)
    plsFit <<- NULL
    gWidgets2::svalue(pls_user_ncomp) <- 1
    gWidgets2::insert(pls_outputtext,"\n***********************************\nMODEL CLEARED!")
    gWidgets2::insert(pls_outputtext,"\n***********************************\n")
  }
  
  #******************************************************************************
  addPop_2_pls_ggraphics <- function(g)
    # Adds popup menu to rigthclick event on ggraphics.
    # To be used in InSpectoR - pls tab.
    # g is the ggraphics widget 
    # Items are: 
    #   1. Copy to clipboard.
    #   2. Save as pdf. 
    #   3. Save as WMF. 
    #*********************************************************************************************************** 
    # B. Panneton, November 2018
    #*********************************************************************************************************** 
  {
    l <- list() 
    
    # l$selectMode <- gaction("Select by rectangle (scores only)","Switch to select mode",icon="gtk-find",
    #                         handler= function(h,...){
    #                           isSelectMode<<-TRUE
    #                           isZoomMode<<-FALSE
    #                           isZoomAll<<-FALSE
    #                         })
    
    # l$zoomMode <- gaction("Zoom (scores only)","Switch to zoom mode",icon="gtk-zoom",
    #                       handler = function(h,...){
    #                         isSelectMode <<- FALSE
    #                         isZoomMode <<- TRUE
    #                         isZoomAll<<-FALSE
    #                       })
    # 
    # l$zoomAll <- gaction("Zoom all (scores only)",tooltip = "Right click to complete",icon="gtk-zoom-out",
    #                      handler= function(h,...){
    #                        isSelectMode <<- FALSE
    #                        isZoomMode <<- FALSE
    #                        isZoomAll<<-TRUE
    #                      })
    
    
    l$copyAction <- gaction("Copy", "Copy current graph to clipboard", icon="copy", 
                            handler=function(h,...) copyToClipboard(g)) 
    l$printActionPDF <- gaction("Save as PDF", "Save current graph", icon="save", 
                                handler=function(h,...) { 
                                  fname <- gfile(gettext("Filename without extension (pdf)"), type="save") 
                                  fname <- paste(fname,".pdf")
                                  if(!file.exists(fname) || gconfirm(gettext("Overwrite file?"))) 
                                    dev.copy2pdf(file=fname)
                                }) 
    l$printActionjpgg <- gaction("Save as WMF", "Save current graph", icon="save", 
                                 handler=function(h,...) { 
                                   fname <- gfile(gettext("Filename without extension (wmf)"), type="save") 
                                   fname <- paste(fname,".wmf")
                                   if(!file.exists(fname) || gconfirm(gettext("Overwrite file?"))) 
                                     dev.copy(win.metafile,file=fname)
                                   dev.off()}) 
    
    addRightclickPopupMenu(g,l)                           
  }                
  
  
#GUI building-------------------------
  #******************************************************************************
  #GUI building
  #******************************************************************************
  
  #Base window tab-------------------------
  #Base window with logo and an empty notebook
  #-------------------------
  main_title<-paste("InSpectoR",ISR_env$laversion," - ",ISR_env$ladate)
  mymain <- gWidgets2::gwindow(main_title, visible=FALSE,
                  width=MainWidth, height=MainHeight,
                 handler=function(h,...){
                   options(warn=-1)
                   if (parcomp) parallel::stopCluster(mycluster)
                   # remove variables from .GlobalEnv
                   var_2_keep <- c("All_XData",
                                   "XDatalist",
                                   "XData",
                                   "XData_p",
                                   "lesACPs",
                                   "lesNCPs",
                                   "prepro_params",
                                   "Ys_df",
                                   "Ys_df_sub",
                                   "subset_ind",
                                   "DataDir",
                                   "plsdaFit",
                                   "plsFit")
                   to_remove_indi <- !(ls(.GlobalEnv) %in% var_2_keep)
                   rm(list=ls(.GlobalEnv)[to_remove_indi],envir = .GlobalEnv)
                   options(warn=0)})
  main_group <- gWidgets2::ggroup(cont=mymain)
  logo <- gWidgets2::gimage(container=main_group)
  nb <- gWidgets2::gnotebook(container=main_group,expand=TRUE,index=TRUE)
  
  #Data selection tab-------------------------
  #Tab for data selection and viewing
  #-------------------------
  data_tab <- gWidgets2::ggroup(container=nb,label="Data selection",horizontal=FALSE)
  nb$add_tab_tooltip(index_data, "Select data for subsequent processing and modeling. View raw data")
  
  file_action = list(merge=gaction("Merge Data sets", icon="convert", handler=Match_Dataset_Multiple),
                     normby=gaction("Normalise by factor levels", icon="spike", handler=Normalize_by),
                     splitat=gaction("Split spectra", icon="split", handler=Split_at_wv))

  menubarlist <- list(DataTransform=file_action) 
  f_menu <- gmenu(menubarlist, cont = mymain) 
  
  
  all_data_tab <- gWidgets2::gvbox(container=data_tab,expand=TRUE)
  
  lefichierY <- gWidgets2::gfilebrowse("Pick a Y file:",container=all_data_tab,expand=FALSE,
                        quote=FALSE,
                        filter=list("Y files" = list(patterns=c("Y*.*"))),
                        handler=OpenYFile, action=TRUE)    #action=TRUE when opening a new file
                                                           #use action=FALSE to apply custom subset to data    
  RGtk2::gtkButtonSetLabel(lefichierY$button,"OPEN a Y file")
  RGtk2::gtkWidgetSetTooltipText(lefichierY$button,
                    "Opening a new Y file will reset the entire application. It is a fresh start!")
  
  xdata_select_group<-gWidgets2::ggroup(container=all_data_tab,expand=FALSE)
  tmp <- gWidgets2::gframe("Spectral data files - Select",container=xdata_select_group,expand=FALSE)
  lesX <- gWidgets2::gtable(c(""),container=tmp,multiple=TRUE, 
                 expand=FALSE)
  gWidgets2::addHandlerSelectionChanged(lesX,Make_XDatalist) #,doplot=TRUE)
  gWidgets2::size(lesX) <- c(280,100)
  gWidgets2::tooltip(lesX)<-paste("A list of available spectrum files. ",
                       "Pick from the list the ones you want to use in viewing and subsequent processing/modeling.",
                       sep="")
  
  btn_xdata_select_group <- gWidgets2::ggroup(container=xdata_select_group,horizontal=FALSE,expand=FALSE)
  gWidgets2::addSpring(btn_xdata_select_group)  #to stack buttons at the bottom of the group
  
  btn_newfactor_data <- gWidgets2::gbutton("Define new factor",container=btn_xdata_select_group,
                             handler=Define_Factor)
  gWidgets2::tooltip(btn_newfactor_data) <- paste("Open a modal dialog to define a new factor",
                                       "based on interactions between existing factors.",
                                       sep=" ")
  
  btn_aggregate_data <- gWidgets2::gbutton("Aggregate spectra",container=btn_xdata_select_group,
                                handler=Aggregate_Spectra)
  gWidgets2::tooltip(btn_aggregate_data) <- paste("Groups measurements by sample ID (first column of Ys_df)",
                                       "Users can select a variable to span over the grouping variable.",
                                       "Aggregation is performed by the interaction ",
                                       "between the grouping and the spanning variables.",
                                       "When no spanning variable, measurements with the same sample ID",
                                       "are aggregated. Otherwise, aggregation is over the range of the",
                                       "spanning variable.",
                                       "The aggregating function is mean.",
                                       sep=" ")
  
  
  btn_subset_data <- gWidgets2::gbutton("Define subset",container=btn_xdata_select_group,
                             handler=Define_Subset)
  gWidgets2::tooltip(btn_subset_data) <- paste("Open a modal dialog to define a subset of data. ",
                                    "Subsets are defined from factors identified in the Y file. ",
                                    "For each factor, one can choose levels to construct subsets. ",
                                    "Selection from all factors are ANDed to construct a subset. ",
                                    "All samples are selected by default or when the defined subset is empty.",
                                    sep="")
  
  btn_find_duplicates <- gWidgets2::gbutton("Find duplicate samples from Col. 1",container=btn_xdata_select_group,
                                 handler=Find_duplicates)
  gWidgets2:: tooltip(btn_find_duplicates) <- paste("Will find duplicates in YData table (bottom left). ",
                                       "Duplicates are selected. ",
                                       "User will be prompted to keep first or last duplicate",
                                        sep="")
  gWidgets2::enabled(btn_find_duplicates) <- FALSE
                             
  btn_remove_samples <- gWidgets2::gbutton("Remove selected samples",container=btn_xdata_select_group,
                                handler=remove_raw_samples)
  gWidgets2::tooltip(btn_remove_samples) <- paste("Will remove the samples selected in the YData table (bottom left). ",
                                       "Samples are removed from the data files. ",
                                       "User will be prompted to accept/reject the creation on an incremental backup.",
                                       sep="")
  gWidgets2::enabled(btn_remove_samples) <- FALSE
  
  btn_merge_Ys <- gWidgets2::gbutton("Merge Y data files",container=btn_xdata_select_group,
                                handler=merge_Ys)
  gWidgets2::tooltip(btn_merge_Ys) <- paste("Will merge data from 2 Y files according to user choices.",sep="")
  gWidgets2::enabled(btn_merge_Ys) <- FALSE
  
  prop_xdata_select_group <- gWidgets2::ggroup(container=xdata_select_group,horizontal=FALSE)
  gWidgets2::addSpring(prop_xdata_select_group)
  gf2 <- gWidgets2::gframe("Properties",container=prop_xdata_select_group,horizontal=FALSE)
  gWidgets2::tooltip(gf2)<-"Manually define limits of x and y-axis. To reset to default, put all these to 0."
  gWidgets2::glabel("X-axis min",container=gf2,anchor=c(-1,0))
  xmini <- gWidgets2::gspinbutton(0,2000,10,value=0, container=gf2)
  gWidgets2::addHandlerChanged(xmini,do_scales_raw_plot)
  gWidgets2::glabel("X-axis max",container=gf2,anchor=c(-1,0))
  xmaxi <- gWidgets2::gspinbutton(0,2000,10,value=0, container=gf2)
  gWidgets2::addHandlerChanged(xmaxi,do_scales_raw_plot)
  gWidgets2::glabel("Y-axis min",container=gf2,anchor=c(-1,0))
  ymini <- gWidgets2::gspinbutton(-10000,100000,1000,value=0, digits=1, container=gf2)
  gWidgets2::addHandlerChanged(ymini,do_scales_raw_plot)
  gWidgets2::glabel("Y-axis max",container=gf2,anchor=c(-1,0))
  ymaxi <- gWidgets2::gspinbutton(0.0,100000.0,1000.0,value=0, digits=1, container=gf2)
  gWidgets2::addHandlerChanged(ymaxi,do_scales_raw_plot)
  prop_axis_limits_list=c(xmini,xmaxi,ymini,ymaxi)
  
  hints_xdata_select_group <-gWidgets2::ggroup(container=xdata_select_group,horizontal=FALSE,expand=TRUE)
  gf3 <- gWidgets2::gframe("Hints",container=hints_xdata_select_group,horizontal=FALSE,expand=TRUE)
  hints_gtext <- gWidgets2::gtext("----- PLOT AREA -----
Right click to select 'Select', 'Zoom' modes and 'Copy' or 'Save'
Select mode: click for nearest or drag rectangle.
Zoom mode: click to zoom around or drag rectangle.
Right click and select Zoom all.
Right click to Copy/Save

----- TABLE AREA -----
Click/CtrlClick/ShiftClick in table to select samples
Click on column for sorting and plot by factor. If column 
    does not contain a factor, user is prompted to convert 
    to a factor. Cancel when prompted to avoid conversion.
                       
----- PROPERTIES -----
Enter X/Y limits for plotting.
If min=0 and max=0 -> reset to full scale.",
      font.attr=list(size="small", foreground="blue", style="oblique"),
                            container=gf3,expand=TRUE)
  
  #Bottom with data table to the left and plots to the right
  bottom_data_tab <- gWidgets2::ggroup(container=all_data_tab,
                            horizontal=TRUE, expand=TRUE)
  gf<-gWidgets2::gframe("Y data",container=bottom_data_tab,expand=TRUE,fill=TRUE)
  gWidgets2::tooltip(gf)<-"See Hints section to the TOP RIGHT for possible interactions."
  #The Y data table. Use gtkTreeView instead of gtable for more convenient column
  #sorting behaviour.  For now, just an invisible empty table. 
  Y_data<-RGtk2::rGtkDataFrame()
  data_view<-RGtk2::gtkTreeView(Y_data)
  Y_selection<-data_view$getSelection()
  Y_selection$setMode("multiple")
  sw<-RGtk2::gtkScrolledWindow()
  gWidgets2::add(gf,sw,expand=TRUE,fill=TRUE)
  sw$add(data_view)
  #add two handlers for plotting spectra and mean spectra by factor level.
  RGtk2::gSignalConnect(Y_data,"sort-column-changed",Plot_By_Factor)
  Y_selection_changed_Signal <- RGtk2::gSignalConnect(Y_selection,"changed",update_rawX_plot)
  
  #The plotting area
  raw_data_graph <- gWidgets2::ggraphics(container=bottom_data_tab,fill=TRUE,
                              width=100, height=100,
                              expand=TRUE,no_popup=TRUE)
  gWidgets2::tooltip(raw_data_graph)<-"See Hints section to the TOP RIGHT for possible interactions."
  #This adds options associated with a right click on the figure (Copy, Save as PDF...)
  addPop_2_raw_ggraphics(raw_data_graph,prop_axis_limits_list)
  #Handle to find sample associated with mouse click on graph
  gWidgets2::addHandlerChanged(raw_data_graph,interact_w_spectra_plot)
  
  
  #Apply tab-------------------------
  #Tab for application of stored models
  #-------------------------
  apply_tab <- gWidgets2::ggroup(container=nb,label="Apply models",horizontal=FALSE)
  nb$add_tab_tooltip(index_applymods, "Apply stored models to current data.")
  
  apply_models_btn_ggroup <- gWidgets2::gvbox(container=apply_tab,expand=TRUE)
  
  btn_apply_PCA <- gWidgets2::gbutton("Apply PCA model",
                           container = apply_models_btn_ggroup,
                           handler = apply_pca_model)
  btn_apply_PLS <- gWidgets2::gbutton("Apply PLS model",
                           container = apply_models_btn_ggroup,
                           handler = apply_pls_model)
  btn_apply_PLSDA <- gWidgets2::gbutton("Apply PLSDA model",
                           container = apply_models_btn_ggroup,
                           handler = apply_plsda_model)
  
  
  
  #Prepro tab -------------------------
  #Tab for preprocessing
  #-------------------------
  # NOTE : PreProDone must be set to FALSE when any option is changed.
  prepro_tab <- gWidgets2::ggroup(cont=nb,label="PrePro", use.scrollwindow = T)
  nb$add_tab_tooltip(index_prepro, "Define data preprocessing steps to selected spectrum types. They will be applied from top to bottom.")
  
  pp_xlim_select_group=gWidgets2::ggroup(cont=prepro_tab,horizontal=FALSE,expand=TRUE)
  gf_truncation<-gWidgets2::gframe("Set wavelenght/wavenumber limits", 
                        container=pp_xlim_select_group)
  #pp_truncate_limits <- list()
  #lyt is just a dummy. A useful one is created
  #on the fly as the selection of spectral data files changes (lesX)
  lyt <- gWidgets2::glayout(container=gf_truncation)
  
  gf1 <- gWidgets2::gframe(cont=pp_xlim_select_group,"Scaling on a per spectrum basis",horizontal = FALSE)
  
 
  
  
  lyt2<-gWidgets2::glayout(container=gf1)
  
  gf_savgol <- gWidgets2::gframe("Parameters for Savitzky-Golay filtering",container=pp_xlim_select_group,horizontal=FALSE)
  
  lyt3<-gWidgets2::glayout(container=gf_savgol)
  
 
  gWidgets2::addSpace(pp_xlim_select_group,10)
  
  btn_apply_pretreatment <- gWidgets2::gbutton("Apply",container = pp_xlim_select_group,
                                    handler = function(h,...){
                                      Apply_PreTreatments()
                                      PreProDone<<-TRUE
                                    })
  
  btn_save_pretreatment <- gWidgets2::gbutton("Save preprocessing options",
                                              container = pp_xlim_select_group,
                                              handler = save_prepro)
  
  btn_load_pretreatment <- gWidgets2::gbutton("Load preprocessing options",
                                              container = pp_xlim_select_group,
                                              handler = load_prepro)
  
  gWidgets2::addSpring(pp_xlim_select_group)
  
  
  
  #PCA tab-------------------------
  #Tab for PCA
  #-------------------------
  pca_tab <- gWidgets2::ggroup(cont=nb,label="PCA")
  nb$add_tab_tooltip(index_acp, "Perform PCA on preprocessed data and graphically display results.")
  
  group_acp<-gWidgets2::ggroup(horizontal=FALSE, container=pca_tab)
  
  #Frame for PCA computation setup
  tmp<-gWidgets2::gframe("Spectral data set for PCA",container=group_acp,expand=FALSE)
  #Select data
  pick_data_4_PCA <- gWidgets2::gcombobox("",handler=MonACP)
  gWidgets2::add(tmp,pick_data_4_PCA,expand=TRUE)
  
  
  #Show number of PCs for 99,5% of variance explained.
  tmp<-gWidgets2::gframe("Estimated number of significant PCs",container=group_acp,horizontal=FALSE)
  nCP_label <- gWidgets2::glabel("Wait for data",expand=TRUE)
  gWidgets2::font(nCP_label) <- list(weight="bold")
  gWidgets2::add(tmp,nCP_label)
  
  gWidgets2::addSpace(group_acp,10)
  
  # To select type of plot
  pickACPPlot=gWidgets2::gcombobox(c("Scores","Var. explained","Score and ortho distances","Loadings"))
  tmp<-gWidgets2::gframe("Choose plot type",container=group_acp,expand=FALSE)
  gWidgets2::add(tmp,pickACPPlot,expand=TRUE)
  
  #To select components
  tmp<-gWidgets2::gframe("Pick PC1 (X-axis for score biplot or first loading)",container=group_acp,expand=FALSE)
  # sélection de PC1 pour graphes d'ACP
  pickPC1<-gWidgets2::gcombobox(c("1","2"))
  gWidgets2::add(tmp,pickPC1)
  tmp<-gWidgets2::gframe("Pick PC2 (Y-axis for score biplot or last loading)",container=group_acp,expand=FALSE)
  # sélection de PC2 pour graphes d'ACP
  pickPC2<-gWidgets2::gcombobox(c("1","2"),selected=2)
  gWidgets2::add(tmp,pickPC2)
  
  #To select factor to use to define point colors
  lbl_f_4_col <- gWidgets2::glabel("Grouping for score plot")
  pickFactor_4_color<-gWidgets2::gcombobox(c(""))
  tmp<-gWidgets2::gframe("Point coloring options",container=group_acp,expand=FALSE)
  gg<-gWidgets2::ggroup(container=tmp,expand=TRUE)
  gWidgets2::add(gg,lbl_f_4_col)
  gWidgets2::add(gg,pickFactor_4_color,expand=TRUE)
  
  #To select factor for labelling points
  lbl_f_4_label <- gWidgets2::glabel("Labels for score plot")
  pickFactor_4_label<-gWidgets2::gcombobox(c(""))
  lbl_c_4_label <- gWidgets2::glabel("Criteria")
  pickCrit_4_label<-gWidgets2::gcombobox(c("All","Leverage","Orthogonal distance","Lev. & ortho"))
  tmp<-gWidgets2::gframe("Point labeling options",container=group_acp,expand=FALSE,horizontal=FALSE)
  gg<-gWidgets2::ggroup(container=tmp,expand=TRUE)
  gWidgets2::add(gg,lbl_f_4_label)
  gWidgets2::add(gg,pickFactor_4_label,expand=TRUE)
  gg2<-gWidgets2::ggroup(container=tmp,expand=TRUE)
  gWidgets2::add(gg2,lbl_c_4_label)
  gWidgets2::add(gg2,pickCrit_4_label,expand=TRUE)
  gWidgets2::svalue(pickCrit_4_label,index=TRUE) <- 1
  
  tmp<-gWidgets2::gframe("Data ellipse plotting",container=group_acp,expand=FALSE,horizontal=TRUE)
  lbl_CI_4_ellipse<-gWidgets2::glabel("Pick option",cont=tmp)
  pickCrit_4_ellipse<-gWidgets2::gcombobox(c("None","0.50","0.95"),container=tmp,expand=TRUE,horizontal=TRUE)
  
  #Button to update plot
  updateACPPlotbut<-gWidgets2::gbutton("Update plot", handler=updateACPPlot)
  gWidgets2::add(group_acp,updateACPPlotbut)
  gWidgets2::enabled(updateACPPlotbut)=TRUE
  
  #Button to remove select points on score plots from subset
  #Should only be active when points are selected
  subset_from_ACPPlotbut<-gWidgets2::gbutton("Remove selected from subset",
                                  handler=function(h,...){
                                    #Subsets ********************
                                    #List of selected is wrt NoSeq but subset_ind is wrt to order in XData
                                    #Update subset_ind with selected scores
                                    ind_subset <- Ys_df_sub[ISR_env$SelectedScores,ncol(Ys_df_sub)-1]
                                    ledum <<- rep(TRUE,length(subset_ind))
                                    ledum[ind_subset] <<- FALSE
                                    subset_ind <<- subset_ind & ledum
                                    
                                    spin<-MakeSpinner()
                                    
                                    #Define Ys_df_sub
                                    Ys_df_sub<<-Ys_df[subset_ind,]
                                    Ys_df_sub<<-droplevels(Ys_df_sub)
                                    Ys_df_sub[,ncol(Ys_df_sub)]<<-seq_len(nrow(Ys_df_sub))
                                    
                                    #Apply subsetting to XData and XData_p
                                    #Retrieve list of selected X data files.
                                    selected<-gWidgets2::svalue(lesX)
                                    XDatalist<<-as.list(as.character(selected))
                                    inds<-gWidgets2::svalue(lesX,index = TRUE)
                                    XData<<-All_XData[inds]
                                    #Build XData and XData_p
                                    if (length(selected)>0){     #only if a data type is selected
                                      XData <<- lapply(XData,data.matrix)
                                      XData <<- lapply(XData,subset,select=-1)  #Remove sample ID (first column) 
                                      XData <<- lapply(XData,function(x) apply(x,2,as.numeric))  #Make sure data are in numeric format
                                      XData <<- lapply(XData,function(x) x[c(TRUE,subset_ind),]) #TRUE required to preserve 1rst row (wl)
                                      XData_p <<- XData
                                      PreProDone <<- FALSE
                                    }
                                    
                                    #Calling OpenYFile with action=FALSE will skip file selection 
                                    #but will load the appropriate subset of data in the GUI environment.
                                    OpenYFile(list(h=lefichierY,action=FALSE))
                                    #***************************
                                    
                                    #Preprocess*****************
                                    Apply_PreTreatments()
                                    PreProDone <<- TRUE
                                    #***************************
                                    
                                    #Recompute PCA**************
                                    N_samples<-nrow(Ys_df_sub)
                                    thr=0.001
                                    lesACPs <<- lapply(XData_p, function(x) prcomp(x[-1,],tol=thr))
                                    #lesNCPs <<- lapply(lesACPs, function(x) ncol(x$rotation))
                                    var_exp=lapply(lesACPs,function(x) summary(x)$importance[2,])
                                    lesNCPs <<- lapply(var_exp, function(x){
                                      i2<-as.numeric(which(x<=thr)[1])
                                      if (x[i2]==0) i2=i2-1  #in case where few samples.
                                      return(i2)
                                    })
                                    leX <- gWidgets2::svalue(pick_data_4_PCA,index=TRUE)
                                    gWidgets2::svalue(nCP_label) <- lesNCPs[[leX]]
                                    #***************************
                                    
                                    #update current plot on PCA tab****
                                    updateACPPlot(updateACPPlotbut)
                                    #***************************
                                    
                                    ISR_env$SelectedScores <- NULL
                                    
                                    dispose(spin)
                                  }
  )
  gWidgets2::add(group_acp,subset_from_ACPPlotbut)
  gWidgets2::enabled(updateACPPlotbut)=FALSE
  
  
  save_pca_scores_btn <- gWidgets2::gbutton("Save scores",cont=group_acp,handler=save_pca_scores)
  
  save_pca_model_btn <- gWidgets2::gbutton("Save PCA model",cont=group_acp,handler=save_pca_model)
  
  gWidgets2::addSpring(group_acp)
  
  ggacp <- gWidgets2::ggraphics(cont=pca_tab,
                     width=100, height=100,
                     expand=TRUE,no_popup=TRUE)
  addPop_2_pca_ggraphics(ggacp)
  #Handle to find sample associated with mouse click on graph
  gWidgets2::addHandlerChanged(ggacp,interact_w_pca_plot)
  
  #PLSDA tab-------------------------
  #Tab for PLSDA
  #-------------------------
  plsda_tab <- gWidgets2::ggroup(cont=nb,label="PLSDA")
  nb$add_tab_tooltip(index_plsda, "Perform PLSDA on preprocessed data and graphically display results. Text outputs are displayed in the Console frame.")

  group_plsda<-gWidgets2::ggroup(horizontal=FALSE, container=plsda_tab, expand=TRUE )
  #Top container for model defs and graphics
  topgroup<-gWidgets2::ggroup(cont=group_plsda,expand=TRUE)

  #Modeling definitions
  choice_frm <- gWidgets2::gframe("Modeling definitions", cont=topgroup)
  choicegroup<-gWidgets2::gvbox(cont=choice_frm,expand=FALSE,spacing=10)
  tmp<-gWidgets2::gframe("Spectral data set(s) for PLSDA",container=choicegroup,expand=FALSE)
  gWidgets2::size(tmp) <- c(200,74)
  #Select data
  Xdat_4_plsda <- gWidgets2::gtable(c(""),container=tmp,multiple=TRUE, 
                 expand=TRUE)
  
  
  gWidgets2::tooltip(Xdat_4_plsda)<-paste("A list of available spectrum files. ",
                       "Pick one from the list to compute plsda on this data only.",
                       "Pick more than one to compute an ensemble of models ",
                       "that will be aggregated.",
                       sep="")
  
  #To select factor to predict
  lbl_f_4_plsda <- gWidgets2::glabel("Pick a factor: ")
  pickFactor_4_plsda<-gWidgets2::gcombobox(c(""))
  tmp<-gWidgets2::gframe("Factor to predict",container=choicegroup,expand=FALSE)
  gWidgets2::add(tmp,lbl_f_4_plsda)
  gWidgets2::add(tmp,pickFactor_4_plsda,expand=TRUE)
  
  
  #Widgets for selecting train control options
  tmp<-gWidgets2::gframe("Train control options",container = choicegroup,expand=FALSE)
  train_ctrl <- gWidgets2::gformlayout(container = tmp)
  #Below the ellipsize argument to make sure gWidgets2::gbutton sizes to show full content.
  gc_plsda_prop <- gWidgets2::gcombobox(seq(0.05,0.95,0.05),label="Prop. of data for training: ",
            container = train_ctrl,ellipsize="none")
  gWidgets2::svalue(train_ctrl)$'Prop. of data for training: ' <- 0.6
  gc_plsda_resamp <- gWidgets2::gcombobox(c("cv","repeatedcv","LOOCV"),label="Resampling method: ",
            container = train_ctrl,ellipsize="none")
  gc_plsda_nfolds <- gWidgets2::gcombobox(seq(3,20,1),label="Number of folds: ",
            container = train_ctrl,ellipsize="none")
  gc_plsda_nreps <- gWidgets2::gcombobox(seq(1,10,1),label="Number of repetitions: ",
            container = train_ctrl,ellipsize="none")
  
  #Widgets for training parameters
  tmp<-gWidgets2::gframe("Training options",container = choicegroup,expand=FALSE)
  train_options<-gWidgets2::gformlayout(container=tmp)
  gc_plsda_nlvs <- gWidgets2::gcombobox(seq(1,50,1),label="Nb of LVs (max): ",
            container = train_options,ellipsize="none")
  gc_plsda_prepro <- gWidgets2::gcombobox(c("None","center"),label="PreProcessing: ",
            container = train_options,ellipsize="none")
  gc_plsda_predmethod <- gWidgets2::gcombobox(c("softmax"),label="Prediction method: ",
            container = train_options,ellipsize="none")
  gc_plsda_perfmetric <- gWidgets2::gcombobox(c("Kappa","Accuracy"),label="Performance metric: ",
            container = train_options,ellipsize="none")
  
  #How to combine models
  tmp<-gWidgets2::gframe("Aggregation options",cont=choicegroup,expand=FALSE)
  aggregate_options<-gWidgets2::gformlayout(container=tmp,expand=TRUE)
  gc_plsda_aggreg <- gWidgets2::gcombobox(c("concatenate","median","max","prod","mean"),label="Aggregation operator: ",
            container = aggregate_options,ellipsize="none")
  
  dum1<-gWidgets2::gbutton("Compute model",container=choicegroup,
               handler=function(h,...){
                 visible(plsda_plots)<-TRUE  #clear plot area
                plot(0,type="n",axes=FALSE,xlab="",ylab="")
                 Compute_PLSDA(h,...)
                 })
  
  btn_save_plsda <-gWidgets2::gbutton("Save model",container=choicegroup,handler=save_plsda_model)
  
  #Plotting area
  output_frm<-gWidgets2::gframe("Output area",cont=topgroup,horizontal = FALSE,expand=TRUE)
  top_output_group<-gWidgets2::ggroup(container = output_frm,horizontal = TRUE,expand=TRUE)
  plot_param_frm<-gWidgets2::gframe("Output controls",cont=top_output_group,expand=FALSE,horizontal=FALSE)
  
  gWidgets2::addSpace(top_output_group,10)
  #plotgroup<-ggroup(cont=plot_frm,expand=TRUE,horizontal=FALSE)
  plsda_plots<-gWidgets2::ggraphics(cont=top_output_group,
                         width=100, height=450,
                         expand=TRUE,no_popup=TRUE)
  
  addPop_2_plsda_ggraphics(plsda_plots)
  
  plot_Acc_vs_VLs <- gWidgets2::gbutton("Accuracy - VL plot", container = plot_param_frm,
                             handler = function(h,...){
                               nom_lesX <- gWidgets2::svalue(Xdat_4_plsda,index=FALSE)
                               accs <- lapply(plsdaFit,function(p) p$results$Accuracy)
                               accsSD <- lapply(plsdaFit,function(p) p$results$AccuracySD)
                               N <- length(accs[[1]])
                               grp <- rep(unlist(nom_lesX),each=N)
                               if (gWidgets2::svalue(aggregate_options)$'Aggregation operator: '=="concatenate")
                                 grp <- rep(paste(nom_lesX,collapse=" + "),N*length(accs))
                               df=data.frame(GRP=grp,
                                             VLs=rep(1:N,times=length(accs)),
                                             accs=unlist(accs),
                                             accsSD=unlist(accsSD))
                               p <- ggplot2::ggplot(df,ggplot2::aes(x=VLs,y=accs)) + 
                                 ggplot2::geom_errorbar(ggplot2::aes(ymin=accs-accsSD,ymax=accs+accsSD),width=0.2,color="mediumblue") + 
                                 ggplot2::geom_line(color="blue",lwd=1)
                               if (length(levels(df$GRP)) > 1){
                                 p <- p +   ggplot2::facet_grid(df$GRP~.)
                               }else
                               {
                                 p <- p + ggplot2::ggtitle(levels(df$GRP))
                               }
                               gWidgets2::visible(plsda_plots)<-TRUE
                               print(p)
                             })
  
  confmat_frm <- gWidgets2::gframe("Confusion matrix", container=plot_param_frm)
  conf_radio <- gWidgets2::gradio(items=c("Validation","Test"), selected=1, container = confmat_frm)
  
  
  plot_confmat <- gWidgets2::gbutton("Plot", container=confmat_frm,
                          expand=FALSE,
                          handler=function(h,...)
                            {
                              if (gWidgets2::svalue(conf_radio)=="Validation"){
                                pred_cl <- Predict_plsda(plsdaFit,probs=FALSE)
                                confusionmat<-caret::confusionMatrix(data=pred_cl,reference = plsdaFit[[1]]$trainingData[,1])
                              }else
                              {
                                testing <- lapply (plsda_set, function(x) x[-plsda_inTrain,])
                                pred_cl <- Predict_plsda(plsdaFit,testing,probs=FALSE)
                                confusionmat<-caret::confusionMatrix(data=pred_cl,reference = testing[[1]][,1])
                              }
                          
                              gWidgets2::visible(plsda_plots)<-TRUE
                              Plot_Confusion_Matrix(confusionmat$table)
                          })
  
  
  
  
  prob_boxp_frm <- gWidgets2::gframe("Probs boxplot", container=plot_param_frm)
  prob_boxp_radio <- gWidgets2::gradio(items=c("Validation","Test"), selected=1, container = prob_boxp_frm)
  plot_prob_boxp <- gWidgets2::gbutton("Plot", container=prob_boxp_frm,
                          expand=FALSE,
                          handler=function(h,...)
                          {
                            if (gWidgets2::svalue(prob_boxp_radio)=="Validation"){
                              pred_prob <- Predict_plsda(plsdaFit,probs=TRUE)
                              dum1<-data.frame(cl=plsdaFit[[1]]$trainingData[,1],pred_prob)
                              
                            }else
                            {
                              testing <- lapply (plsda_set, function(x) x[-plsda_inTrain,])
                              pred_prob <- Predict_plsda(plsdaFit,testing,probs=TRUE)
                              dum1<-data.frame(cl=testing[[1]][,1],pred_prob)
                              
                            }
                            dum2<-tidyr::gather(dum1,Pred,Prob,-cl,factor_key = TRUE)
                            levels(dum2$cl)=paste("True: ",levels(dum2$cl),sep="")
                            gWidgets2::visible(plsda_plots)<-TRUE
                            p<-ggplot2::ggplot(dum2,ggplot2::aes(Pred,Prob))
                            p<-p+ggplot2::geom_boxplot()+ggplot2::facet_wrap(~cl)
                            p<-p + ggplot2::theme(text = element_text(size=9))
                            p<-p + ggplot2::theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1))
                            gWidgets2::visible(plsda_plots)<-TRUE
                            print(p)
                          })
  
  
  
  
  
  
  show_probs_but <- gWidgets2::gbutton("Show probs",container=plot_param_frm,
                            handler=function(h,...){
                              pr<-Predict_plsda(plsdaFit,plsda_set,probs=TRUE)
                              pr<-round(pr,2)
                              lesdiffs<-t(apply(pr,1,sort,decreasing=TRUE))
                              lesdiffs<-lesdiffs[,1]-lesdiffs[,2]
                              lescl<-Predict_plsda(plsdaFit,plsda_set,probs=FALSE)
                              letest <- lescl==plsda_set[[1]][,1]
                              pr<-data.frame(NoSeq=Ys_df_sub$NoSeq,
                                            True_Cl=plsda_set[[1]][,1],
                                            Pred_Cl=lescl,
                                            Test=letest,
                                            minDiff=lesdiffs,pr)
                              dumw<-gWidgets2::gwindow("Probabilities", visible=FALSE)
                              dumg <- gWidgets2::ggroup(cont=dumw)
                              dumtbl <- gWidgets2::gtable(pr, cont=dumg,expand=TRUE, fill=TRUE)
                              dumf<-gWidgets2::gfilter(dumtbl,
                                                       initial.vars = data.frame(colnames(pr)[c(2,3)],
                                                                                 colnames(pr)[c(2,3)],
                                                                                 c("single","multiple"),
                                                                                 stringsAsFactors=FALSE),
                                            allow.edit = TRUE,
                                            container = dumg,
                                            handler=function(h,...) {
                                              gWidgets2::visible(dumtbl)<-h$obj$get_value()
                                            }
                              )
                              gWidgets2::visible(dumw) <- TRUE
                            })
 
  
  
  plot_b_coeff_but <-gWidgets2::gbutton("Plot B-coeff", container=plot_param_frm,
                             handler=function(h,...){
                               if (gWidgets2::svalue(aggregate_options)$'Aggregation operator: '=="concatenate"){
                                 gWidgets2::visible(plsda_plots) <- TRUE 
                                 x_id=plsdaFit[[1]]$coefnames
                                 x_id=strsplit(x_id,"_")
                                 N_xlabels <- length(x_id)
                                 x_id=unlist(x_id)
                                 x_cls=x_id[seq(1,2*N_xlabels,2)]
                                 xlabels=as.numeric(x_id[seq(2,2*N_xlabels,2)])
                                 nc=plsdaFit[[1]]$bestTune$ncomp
                                 # df=data.frame(y=plsdaFit[[1]]$finalModel$coefficients[,1,nc],
                                 #               Cl=x_cls, 
                                 #               xlabs=as.numeric(xlabels))
                                 
                                 
                                 df=data.frame(Cl=x_cls, 
                                               xlabs=as.numeric(xlabels))
                                 for (k in 1:length(plsdaFit[[1]]$levels)){
                                   eval(parse(text=paste0("tt=data.frame(",plsdaFit[[1]]$levels[k],"=plsdaFit[[1]]$finalModel$coefficients[,k,nc])")))
                                   df = cbind(df,tt)
                                 }
                                 
                                 df1 <- tidyr::gather(data=df,key="Class",value,-Cl,-xlabs)
                                 
                                 p <- ggplot2::ggplot(df1, ggplot2::aes(x=xlabs,y=value,colour=Class))+
                                   #geom_line(size=1.2) + 
                                   ggplot2::geom_smooth() + ggplot2::facet_wrap(~Cl, ncol=1,scales="free_y")
                                   #facet_grid(Cl ~ .)
                                 gWidgets2::visible(plsda_plots) <- TRUE
                                 print(p)
                               }else
                               {
                                 coeffs <- lapply(plsdaFit, function(x) coef(x$finalModel))
                                 ind_lesX <- gWidgets2::svalue(Xdat_4_plsda,index = TRUE)
                                 nom_lesX <- gWidgets2::svalue(Xdat_4_plsda,index=FALSE)
                                 ind_lesX_list<-as.list(ind_lesX)
                                 wl <- lapply(ind_lesX_list, function(ii) XData_p[[ii]][1,])
                                 plotframe<-cbind(coeffs[[1]],
                                                 data.frame(wl=wl[[1]],
                                                            source=rep(nom_lesX[[1]],length(wl[[1]]))))
                                 N_Classes<-ncol(coeffs[[1]])
                                 colnames(plotframe)[1:N_Classes]<-colnames(coeffs[[1]])
                                 if (length(wl)>1){
                                   for (k in (2:length(wl))){
                                      dum<-cbind(coeffs[[k]],
                                                data.frame(wl=wl[[k]],
                                                           source=rep(nom_lesX[[k]],length(wl[[k]]))))
                                      colnames(dum)[1:N_Classes]<-colnames(plotframe)[1:N_Classes]
                                      plotframe<-rbind(plotframe,dum)
                                   }
                                 }
                                 dum<-reshape2::melt(plotframe,id.vars=(N_Classes+c(1,2)),
                                          variable.name="Class",
                                          value.name="B_coeffs")
                                 p <- ggplot2::ggplot(dum,ggplot2::aes(x=wl,y=B_coeffs,colour=Class))
                                 p <- p + ggplot2::geom_smooth() + ggplot2::facet_wrap(~source, ncol=1,scales="free_y")
                                 gWidgets2::visible(plsda_plots)<-TRUE
                                 print(p)
                               }
                             })
  
  
  #Bottom container for console output
  bottomgroup<-gWidgets2::gvbox(cont=group_plsda,expand=TRUE)
  console_frm <- gWidgets2::gframe("Console output",cont=output_frm,expand=TRUE)
  
  plsda_outputtext<-gWidgets2::gtext("",cont=console_frm,
                    font.attr=list(family="monospace",foreground="blue",size="x-small")
                    ,expand=TRUE,fill=TRUE)
  #Controls for saving console output to a text file
  outputgroup<-gWidgets2::ggroup(cont=bottomgroup)
  plsda_2_txt_button<-gWidgets2::gbutton(text="Console to txt",cont=outputgroup,expand=FALSE,fill=FALSE,
                     handler=function(h,...){
                       dum<-unlist(plsda_txt_output)
                       lelabel<-gWidgets2::ginput("Enter an output ID:","",icon="question",
                                      parent=mymain)
                       write(c("----------------------",
                               lelabel,"\n",dum,"\n"),
                             file=gWidgets2::svalue(plsda_txt_filename),append=TRUE)
                     })
  plsda_txt_filename_definebutton<-
    gWidgets2::gbutton(text="Pick an output file name",cont=outputgroup,
       expand=FALSE,fill=FALSE, handler=function(h,...){
         initfname<-gWidgets2::svalue(plsda_txt_filename)
         gWidgets2::svalue(plsda_txt_filename)<-gWidgets2::gfile("Pick a file name for output",
                                    type="open",intial.filename=basename(initfname),
                                    initial.dir=dirname(initfname))
       })
  filelabel<-gWidgets2::glabel("Output file name: ",cont=outputgroup)
  plsda_txt_filename<-gWidgets2::gedit(file.path(getwd(),"Default_PLSDA_Ouptut.txt"),
                    cont=outputgroup,expand=TRUE,handler=function(h,...){
                      #Place cursor at the end to make sure the file name is visible to the user
                      RGtk2::gtkEditableSetPosition(h$obj$widget,
                                             nchar(gWidgets2::svalue(h$obj))) 
                    })
  
  #Place cursor at the end to make sure the file name is visible to the user
  RGtk2::gtkEditableSetPosition(plsda_txt_filename$widget, nchar(gWidgets2::svalue(plsda_txt_filename))) 
  
  
  w_4_plsdareset = list(pickFactor_4_plsda,Xdat_4_plsda,
                        gc_plsda_prop,
                        gc_plsda_resamp,
                        gc_plsda_nfolds,
                        gc_plsda_nreps,
                        gc_plsda_nlvs,
                        gc_plsda_prepro,
                        gc_plsda_predmethod,
                        gc_plsda_perfmetric,
                        gc_plsda_aggreg)
  lapply(w_4_plsdareset, function(x){
    gWidgets2::addHandlerClicked(x,handler = function(h,...){
      if (!is.null(plsdaFit)) clear_plsdafit()
    })
  })
  
  
  #PLS tab-------------------------
  #Tab for PLS
  #-------------------------
  pls_tab <- gWidgets2::ggroup(cont=nb,label="PLS")
  nb$add_tab_tooltip(index_pls, "Perform PLS on preprocessed data and graphically display results. Text outputs are displayed in the Console frame.")
  
  group_pls<-gWidgets2::ggroup(horizontal=FALSE, container=pls_tab, expand=TRUE )
  #Top container for model defs and graphics
  topgroup<-gWidgets2::ggroup(cont=group_pls,expand=TRUE)
  
  #Modeling definitions
  choice_frm <- gWidgets2::gframe("Modeling definitions", cont=topgroup)
  choicegroup<-gWidgets2::gvbox(cont=choice_frm,expand=FALSE,spacing=10)
  tmp<-gWidgets2::gframe("Spectral data set(s) for PLS",container=choicegroup,expand=FALSE)
  gWidgets2::size(tmp) <- c(200,74)
  #Select data
  Xdat_4_pls <- gWidgets2::gtable(c(""),container=tmp,multiple=TRUE, 
                         expand=TRUE)
  
  gWidgets2::tooltip(Xdat_4_pls)<-paste("A list of available spectrum files. ",
                               "Pick one from the list to compute pls on this data only.",
                               "Pick more than one to compute a model ",
                               "on concatenated spectra.",
                               sep="")
  
  #To select variable to predict
  var_f_4_pls <- gWidgets2::glabel("Pick a variable: ")
  pickvar_4_pls<-gWidgets2::gcombobox(c(""))
  tmp<-gWidgets2::gframe("Variable to predict",container=choicegroup,expand=FALSE)
  gWidgets2::add(tmp,var_f_4_pls)
  gWidgets2::add(tmp,pickvar_4_pls,expand=TRUE)
  
  
  #Widgets for control options
  tmp<-gWidgets2::gframe("Train control options",container = choicegroup,expand=FALSE)
  pls_ctrl <- gWidgets2::gformlayout(container = tmp)
  #Below the ellipsize argument to make sure gWidgets2::gbutton sizes to show full content.
  gc_pls_resamp <- gWidgets2::gcombobox(c("none","CV","LOO"),label="Resampling method: ",
            container = pls_ctrl,ellipsize="none")
  gc_pls_nlvs <- gWidgets2::gcombobox(seq(1,50,1),label="Nb of LVs (max): ",
            container = pls_ctrl,ellipsize="none",
            handler = function(h,...) {
              pls_user_ncomp$widget$SetRange(1,gWidgets2::svalue(h$obj))
              plsscore_1$widget$SetRange(1,gWidgets2::svalue(h$obj))
              plsscore_2$widget$SetRange(1,gWidgets2::svalue(h$obj))
              clear_plsfit()
            })
  gc_pls_scaling <- gWidgets2::gcombobox(c("none","scale"),label="Scaling: ",
            container = pls_ctrl,ellipsize="none")
  
  #How to combine models
  tmp<-gWidgets2::gframe("Aggregation options",cont=choicegroup,expand=FALSE)
  pls_aggregate_options<-gWidgets2::gformlayout(container=tmp,expand=TRUE)
  gc_pls_aggreg <- gWidgets2::gcombobox(c("concatenate spectra"),label="Aggregation method: ",
            container = pls_aggregate_options,ellipsize="none")
  
  dum1<-gWidgets2::gbutton("Compute model",container=choicegroup,
                handler=function(h,...){
                  gWidgets2::visible(pls_plots)<-TRUE  #clear plot area
                  plot(0,type="n",axes=FALSE,xlab="",ylab="")
                  Compute_PLS(h,...)
                })
  
  dum2 <- gWidgets2::gbutton("Save model",container=choicegroup,
                  handler = save_pls_model)
  
  
  #Plotting area
  output_frm<-gWidgets2::gframe("Output area",cont=topgroup,horizontal = FALSE,expand=TRUE)
  top_output_group<-gWidgets2::ggroup(container = output_frm,horizontal = TRUE,expand=TRUE)
  plot_param_frm<-gWidgets2::gframe("Output controls",cont=top_output_group,expand=FALSE,horizontal=FALSE)
  
  gWidgets2::addSpace(top_output_group,10)
  
  pls_plots<-gWidgets2::ggraphics(cont=top_output_group,
                         width=100, height=450,
                         expand=TRUE,no_popup=TRUE)
  addPop_2_pls_ggraphics(pls_plots)
  gWidgets2::addHandlerChanged(pls_plots,interact_w_pls_plot)
  
  pls_plot_type <- gWidgets2::glabel("None",container = plot_param_frm)
  
  plsvalidation_frm <- gWidgets2::gframe("Validation plot", container=plot_param_frm)
  
  plot_plsvalidation <- gWidgets2::gbutton("Plot", container=plsvalidation_frm,
                                expand=FALSE,
                                handler=function(h,...)
                                {
                                  # validationplot(plsFit[[1]],estimate="all",
                                  #                legendpos="bottomleft",
                                  #                lwd=2)
                                  dats = pls::RMSEP(plsFit[[1]], estimate="all")
                                  df=as.data.frame(t(dats$val[,1,]))
                                  df$x=seq_len(dim(df)[1])-1
                                  plot_data <- reshape2::melt(df,id.var="x")
                                  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x=x,y=value,group=variable,colour=variable))+
                                    ggplot2::geom_point(size=4, shape=21, stroke=1.2)+
                                    ggplot2::geom_line(ggplot2::aes(lty=variable),lwd=1.2)+
                                    ggplot2::xlab("Nb of latent variables") + ggplot2::ylab("RMSE")
                                  gWidgets2::visible(pls_plots) <- TRUE
                                  print(p)
                                  gWidgets2::svalue(pls_plot_type) <- "Validation"
                                })
  
  dum <- gWidgets2::ggroup(container = plsvalidation_frm,horizontal=TRUE)
  gWidgets2::addSpace(dum,10)
  dum2 <- gWidgets2::glabel("Nb of latent var.",container = dum)
  pls_user_ncomp <- gWidgets2::gspinbutton(from = 1, to = 1, by = 1, container = dum)
  
  
  plspred_frm <- gWidgets2::gframe("Prediction plot", container=plot_param_frm)
  dum0 <- gWidgets2::gvbox(container = plspred_frm)
  # dum <- ggroup(container = dum0,horizontal=TRUE)
  # dum2 <- gWidgets2::glabel("Nb of latent var.",container = dum)
  # plspred_ncomp <- gspinbutton(from = 1, to = 1, by = 1, container = dum)
  plspred_which <- gWidgets2::gradio(items=c("train","validation"), selected=2, container = dum0)
  dumf <- gframe("Color by", container = dum0)
  pls_pickFactor_4_color<-gWidgets2::gcombobox(c(""),container = dumf)
  dumf2 <- gframe("Label with", container = dum0)
  pls_pickFactor_4_label <- gWidgets2::gcombobox(c(""),container = dumf2)
  
  
  
  plot_plspred <- gWidgets2::gbutton("Plot", container=plspred_frm,
                          expand=FALSE,
                          handler= plot_PLS_preds)
  
  plscoeff_frm <- gWidgets2::gframe("B-coeff plot", container=plot_param_frm)
  # dum0 <- gvbox(container = plscoeff_frm)
  # dum <- ggroup(container = dum0,horizontal=TRUE)
  # dum2 <- gWidgets2::glabel("Nb of latent var.",container = dum)
  # plscoeff_ncomp <- gspinbutton(from = 1, to = 1, by = 1, container = dum)
  plot_pls_b_coeff_but <-gWidgets2::gbutton("Plot", container=plscoeff_frm,
                             handler=function(h,...){
                                x_id=rownames(plsFit[[1]]$coefficients)
                                x_id=strsplit(x_id,"_")
                                N_xlabels <- length(x_id)
                                x_id=unlist(x_id)
                                x_cls=x_id[seq(1,2*N_xlabels,2)]
                                xlabels=as.numeric(x_id[seq(2,2*N_xlabels,2)])
                                nc=gWidgets2::svalue(pls_user_ncomp)
                                df=data.frame(y=plsFit[[1]]$coefficients[,1,nc],
                                                 Cl=x_cls, 
                                                 xlabs=as.numeric(xlabels))
                                
                                p <<- ggplot2::ggplot(df, ggplot2::aes(x=xlabs,y=y,colour=Cl))+
                                  ggplot2::geom_line() + ggplot2::xlab("Wavelength or Wavenumber") + ggplot2::ylab("B-coefficients") +
                                  ggplot2::facet_wrap( ~ Cl, ncol=1, scales="free_y") + ggplot2::theme(legend.position="none") +
                                  ggplot2::theme(strip.background = ggplot2::element_rect(fill="grey80"))
                                  gWidgets2::visible(pls_plots) <- TRUE
                                  
                                  
                                visible(pls_plots) <- TRUE
                                print(p)
                                
                                gWidgets2::svalue(pls_plot_type) <- "B-coeff"
                                
                             })

  plsscores_frm <- gWidgets2::gframe("Score plot", container=plot_param_frm)
  dum0 <- gWidgets2::gvbox(container = plsscores_frm)
  dum <- gWidgets2::ggroup(container = dum0,horizontal=TRUE)
  dum2 <- gWidgets2::glabel("First latent var.",container = dum)
  plsscore_1 <- gWidgets2::gspinbutton(from = 1, to = 1, by = 1, container = dum)
  dum <- gWidgets2::ggroup(container = dum0,horizontal=TRUE)
  dum2 <- gWidgets2::glabel("Second latent var.",container = dum)
  plsscore_2 <- gWidgets2::gspinbutton(from = 1, to = 1, by = 1, container = dum)
  plot_pls_scores_but <-gWidgets2::gbutton("Plot", container=plsscores_frm,
                                 handler=function(h,...){
                                   visible(pls_plots) <- TRUE
                                   par(bg="gray90")
                                   pls::scoreplot(plsFit[[1]],
                                             comps=c(gWidgets2::svalue(plsscore_1),gWidgets2::svalue(plsscore_2)),
                                             col="blue", pch=20, cex=1.2,
                                             panel.first = grid(col="white",lty=1))
                                   
                                   gWidgets2::svalue(pls_plot_type) <- "Score"
                                  })
  
  plstable_frm <- gWidgets2::gframe("Pred. table", container=plot_param_frm)
  dum0 <- gWidgets2::gvbox(container = plstable_frm)
  dum <- gWidgets2::ggroup(container = dum0,horizontal=TRUE)
  # dum2 <- gWidgets2::glabel("Nb of latent var.",container = dum)
  # plstable_ncomp <- gspinbutton(from = 1, to = 1, by = 1, container = dum)
  nc<-gWidgets2::svalue(pls_user_ncomp)
  btn_show_predict <- gWidgets2::gbutton("Show pred. table",container = dum,
                              handler=function(h,...){
                                pl<-plot(plsFit[[1]],plottype = "prediction",
                                         ncomp=gWidgets2::svalue(pls_user_ncomp),
                                         which="train")
                                dat_tbl=as.data.frame(pl)
                                pl<-plot(plsFit[[1]],plottype = "prediction",
                                         ncomp=gWidgets2::svalue(pls_user_ncomp),
                                         which="validation")
                                dat_tbl=cbind(dat_tbl,as.data.frame(pl)[,2])
                                dat_tbl=cbind(Ys_df_sub[,1],dat_tbl,NoSeq=seq_len(nrow(dat_tbl)))
                                colnames(dat_tbl)[c(1,3,4)]=c(colnames(Ys_df_sub)[1],
                                                              "Training","Validation")
                                
                                plot(dat_tbl[,c(-1,-5)])
                                
                                tmp_w=gWidgets2::gwindow()
                                gg=gWidgets2::gvbox(container = tmp_w)
                                btn_save_plspred=gWidgets2::gbutton("Save to file",container = gg,
                                                         handler=function(h,...){
                                                           outfile=choose.files(default=file.path(basename(gWidgets2::svalue(lefichierY)),"*.txt"),
                                                                                caption="Define file name",
                                                                                multi=FALSE, 
                                                                                filters=Filters[c("txt")])
                                                           utils::write.table(dat_tbl,
                                                                       file=outfile,
                                                                       sep="\t",
                                                                       dec=".",
                                                                       row.names = FALSE,
                                                                       quote=FALSE)
                                                         })
                                gt=gWidgets2::gtable(dat_tbl,container = gg)
                                gt$widget$setGridLines('both')
                              })
  
  #Bottom container for console output
  bottomgroup<-gWidgets2::gvbox(cont=group_pls,expand=TRUE)
  console_frm <- gWidgets2::gframe("Console output",cont=output_frm,expand=TRUE)
  
  pls_outputtext<-gWidgets2::gtext("",cont=console_frm,
                          font.attr=list(family="monospace",foreground="blue",size="x-small")
                          ,expand=TRUE,fill=TRUE)
  #Controls for saving console output to a text file
  outputgroup<-gWidgets2::ggroup(cont=bottomgroup)
  pls_2_txt_button<-gWidgets2::gbutton(text="Console to txt",cont=outputgroup,expand=FALSE,fill=FALSE,
                              handler=function(h,...){
                                dum<-unlist(pls_txt_output)
                                lelabel<-gWidgets2::ginput("Enter an output ID:","",icon="question",
                                                parent=mymain)
                                write(c("----------------------",
                                        lelabel,"\n",dum,"\n"),
                                      file=gWidgets2::svalue(pls_txt_filename),append=TRUE)
                              })
  pls_txt_filename_definebutton<-
    gWidgets2::gbutton(text="Pick an output file name",cont=outputgroup,
            expand=FALSE,fill=FALSE, handler=function(h,...){
              initfname<-gWidgets2::svalue(pls_txt_filename)
              gWidgets2::svalue(pls_txt_filename)<-gWidgets2::gfile("Pick a file name for output",
                                                type="open",intial.filename=basename(initfname),
                                                initial.dir=dirname(initfname))
            })
  filelabel<-gWidgets2::glabel("Output file name: ",cont=outputgroup)
  pls_txt_filename<-gWidgets2::gedit(file.path(getwd(),"Default_PLS_Ouptut.txt"),
                            cont=outputgroup,expand=TRUE,handler=function(h,...){
                              #Place cursor at the end to make sure the file name is visible to the user
                              RGtk2::gtkEditableSetPosition(h$obj$widget,
                                                     nchar(gWidgets2::svalue(h$obj))) 
                            })
  
  #Place cursor at the end to make sure the file name is visible to the user
  RGtk2::gtkEditableSetPosition(pls_txt_filename$widget, nchar(gWidgets2::svalue(pls_txt_filename))) 
  
  w_4_plsreset = list(Xdat_4_pls,
                      pickvar_4_pls,
                      gc_pls_resamp,
                      gc_pls_scaling,
                      gc_pls_aggreg)
  lapply(w_4_plsreset, function(x){
    gWidgets2::addHandlerClicked(x,handler = function(h,...){
      if (!is.null(plsFit)) clear_plsfit()
    })
  })
  

  #******************************************************************************
  #Initialise-------------------------
  #Initialise
  #-------------------------
 
  # #Initialise l'interface
  gWidgets2::svalue(logo) <- logoimg 
  gWidgets2::visible(mymain)<-TRUE
  btn_needing_Ydata<-list(btn_subset_data,btn_newfactor_data, btn_find_duplicates,
                          btn_merge_Ys,
                          btn_aggregate_data)
  lapply(btn_needing_Ydata,function(x) enabled(x)<-FALSE)
  have_data<-list(nb[index_applymods],nb[index_prepro],nb[index_acp],nb[index_plsda],nb[index_pls])
  lapply(have_data, function(x) x$widget$hide())
  #List all ggraphics objects
  ggraphs <- list(raw_data_graph,ggacp,plsda_plots)
  #List of axis limits properties on data selection tab.
  proplist <- list(xmini,xmaxi,ymini,ymaxi)
  #This will set the delay for showing tooltip to 1000 msec.
  dum<-RGtk2::gtkWidgetGetSettings(mymain$widget)
  RGtk2::gtkSettingsSetLongProperty(dum,'gtk-tooltip-timeout',1000,0)
  RGtk2::gtkSettingsSetLongProperty(dum,'gtk-tooltip-browse-timeout',60,0)
  
  if (MainHeight<701){
    RGtk2::gtkWindowResize(mymain$widget,MainWidth,700)
  }else
  {
    RGtk2::gtkWindowResize(mymain$widget,MainWidth,MainHeight)
  }
  gWidgets2::svalue(nb) <- index_data
  
  if (!is.null(yfile)){
    gWidgets2::svalue(lefichierY) <- yfile
    #OpenYFile(list(h=lefichierY,action=TRUE))
  }
  
  ## disable normalis
  enabled(file_action$normby) <- FALSE
  enabled(file_action$splitat) <- FALSE
  
  # dum='A'
  # while (toupper(dum) != 'Q') dum=readline("\nq to quit: ")
  
  #list(win=mymain)  #returned so can be used to keep GUI active.
                    #See runInSpectoR.m
}                                    #end Main
