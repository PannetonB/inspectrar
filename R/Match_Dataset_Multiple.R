#' Match_Dataset_Multiple
#' 
#' Undocumented function to merge data in the special case where several 
#' instruments are used to acquire data but the data types are not the
#' same for each of the instrument.
#'
#' @return
#'  Nothing. Creates a new set of data files (Y and X files) of the
#'  merged data
#'  
#' @export
#'
#' @examples
#' Match_Dataset_Multiple()
#'
Match_Dataset_Multiple <- function(...){
  #Choisir un fichier Y
  #Sommaire par instrument indiquant pr?sence/absence de fluo ou reflect
  #Pour un instrument o? on a les deux mesures, InSpector permet de choisir
  #un sous-ensemble pour cet instrument pour constituer un jeu de donn?es conforme
  #et sans NA.
  #Si on a un ou des instruments avec de la fluo seulement (Fs) et un ou des instruments
  #avec de la r?flectance seulement (Rs) faire choisir un instrument parmi les Fs et
  #un instrument parmi les Rs. ?videmment si il n'y a qu'un Fs ou qu'un Rs, l'usager
  #n'a pas ? choisir.
  
  #Ouvre le fichier Y
  library(tcltk)
  lefichierY <- choose.files(caption="Choisir un fichier Y",
                             default=file.path(getwd(),"Y*.txt"),
                             filters = Filters[c("txt"),])
  Ys_df <- read.table(lefichierY,header=TRUE,sep="\t",dec=".")
  
  
  
  #Sommaire par instrument indiquant pr?sence/absence de fluo ou reflect
  lesinstrus <- as.character(unique(Ys_df$Instrument))
  N_instrus <- length(lesinstrus)
  if (N_instrus==1){
    cat("Seulement 1 instrument!")
    return("Quitte")
  }
  lepath<-dirname(lefichierY)
  DataDir <-lepath
  ynom <- lefichierY
  dname <- strsplit(tools::file_path_sans_ext(basename(ynom)),"_")[[1]]
  dname=paste(dname[-1],collapse="_")
  #Exclude Y file
  dum<-dir(lepath,paste0("*",dname,"*"))
  dum <- dum[!grepl("Y_",dum)]
  #dum<-dir(lepath,paste0("*",dname,"_I.txt"))
  inds<-grep(dname,dum)
  #Make sure 1rst column of Ys_df matches 1rst column of XData files.
  #If not, drop the corresponding XData file from the list of lesX and display a warning.
  nn<-length(inds)
  All_XData <- vector("list",nn)
  for (k in 1:nn){
    dum1<-read.table(file.path(lepath,dum[inds[k]]),sep="\t", dec=".",header=FALSE)
    test1 <- length(dum1[-1,1])==length(Ys_df[,1])
    test2 <- TRUE
    if (test1) test2 <- any(as.character(dum1[-1,1])!=as.character(Ys_df[,1]))
    #cat("\n",length(dum1[-1,1])," ",length(Ys_df[,1])," ",test1," ",test2)
    if (test2){  #Not matching
      gWidgets2::gmessage(paste("File ",dum[inds[k]], " is not valid! Will be ignored.",sep=""),
                          title="WARNING",icon="warning")
      inds[k]<-0
    }else
      All_XData[[k]]<-dum1
  }
  hasFluo=rep(TRUE,N_instrus)
  hasReflect=rep(TRUE,N_instrus)
  for (k in 1:N_instrus){
    filtre <- Ys_df$Instrument==lesinstrus[k]
    unEx <- which(startsWith(dum,"EX"))[1]
    if (all(is.na(All_XData[[unEx]][c(FALSE,filtre),-1]))) hasFluo[k] <- FALSE
    unTr <- which(startsWith(dum,"Tr"))[1]
    if (all(is.na(All_XData[[unTr]][c(FALSE,filtre),-1]))) hasReflect[k] <- FALSE
  }
  cat("\nDonn?es de fluorescence pour: ", lesinstrus[hasFluo],"\n")
  cat("\nDonn?es de transmittance/r?flectance pour: ", lesinstrus[hasReflect])
  
  
  instr_4_fluo=tk_select.list(as.character(lesinstrus[hasFluo]),multiple=T,
                              title="Choisir les instruments pour des donn?es de fluorescence")
  instr_4_reflect=tk_select.list(as.character(lesinstrus[hasReflect]),multiple=F,
                                 title="Choisir un instrument pour des donn?es de r?flectance")
  
 
   Ys_df_new <- data.frame()
   N_instr_4_fluo <- length(instr_4_fluo)
   All_XData_new <- All_XData
   for (kk in seq_len(N_instr_4_fluo)){
    #Cr?e un nouveau Ys_df en enlevant toutes les lignes qui ne correspondent pas ? l'instrument pour la fluo
    indi <- Ys_df$Instrument==instr_4_fluo[kk]
    Ys_df_new <- rbind(Ys_df_new,Ys_df[indi,])
    
    #Cr?e les fichier de fluo correspondant pour les donn?es
    indj <- which(startsWith(dum,"EX"))
    if (kk==1){
      for (jj in indj) All_XData_new[[jj]] <- All_XData[[jj]][c(1,(which(indi)+1)),]
    }else
    {
      for (jj in indj) All_XData_new[[jj]] <- rbind(All_XData_new[[jj]],All_XData[[jj]][which(indi)+1,])
    }
  }
  
  #Pour les donn?es de r?flectance, enlever toutes les lignes avec les NAs et rbind pour chaque inst de fluo
  ind_tr <- which(startsWith(dum,"Tr"))
  df=All_XData_new[[ind_tr]]
  ind_na <- rowSums(is.na(df)) != ncol(df[,-1])
  All_XData_new[[ind_tr]] <- All_XData_new[[ind_tr]][ind_na,]
  df_tr <- All_XData_new[[ind_tr]][-1,]
  for (kk in 2:N_instr_4_fluo) All_XData_new[[ind_tr]] <- rbind(All_XData_new[[ind_tr]], df_tr)
  
  #Enregistre les donn?es dans un nouveau jeu de donn?es
  
  fname <- winDialogString('Nom du jeu de donn?es', 'Test')
  newpath <- tk_choose.dir(default = lepath, caption = "Choisir un r?pertoire pour stocker les donn?es")
  
  Yname <- paste0("Y_",fname,".txt")
  Xnames <- gsub(dname,fname,dum)
  
  write.table(Ys_df_new, file=file.path(newpath,Yname), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  indi <- as.list(seq(1,length(All_XData_new)))
  lapply(indi,function(ii) 
    write.table(All_XData_new[[ii]], file=file.path(newpath,Xnames[[ii]]), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE))
}
