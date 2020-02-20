showStockIcons <- function()
{  
  les_icons <- getStockIcons()
  N <- length(les_icons)
  nrangee <- 25
  ncolonne <- N %/% nrangee
  reste <- N %% nrangee
  w <- gwindow()
  gtop <- ggroup(cont=w,horizontal = T)
  for (i in 1:ncolonne){
    g <- ggroup(cont=gtop,horizontal = F)
    for (j in 1:nrangee){
      unicon <- les_icons[[(i-1)*nrangee+j]]
      b <- gbutton(unicon,cont=g)
      
    }
  }
  g <- ggroup(cont=gtop,horizontal = F)
  for (j in 1:reste){
    unicon <- les_icons[[25*ncolonne+j]]
    b <- gbutton(unicon,cont=g)
  }
}  