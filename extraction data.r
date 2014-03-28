###### Ce fichier contient les fonctions permettant d'extraire nos données contenues dans des dossiers adjacents
# Il permet aussi de calculer les rendements, outils de travail au même titre que les cours des actions eux mêmes (d'où leur présence dans ce fichier)

# la gestion des dates présentes dans les tableaux des cours nécessite le package suivant
library("date")

periode <- function(data){
  return(rbind(date.mmddyy(min(data[,"Date"])),date.mmddyy(max(data[,"Date"]))))
}


# 1) Renvoie un tableau dont les colonnes sont les valeurs de clôture de différents titres

# titres = notre portefeuille
# date1 et date2 = les intervalles choisis, par défaut toute la période commune aux titres
# type = J, H ou M (journalier, hebdomadaire ou mensuel)

load_data <- function(deb, fin, titres = NULL, type = "J"){
  # le path est à changer selon chaque utilisateur
  path_data="D:/cours/2A/Statap/Data/Cours/"
  path_data=paste(path_data,type,"/",sep="")
  
  # on sélectionne les fichiers contenus dans le dossier indiqué
  files <- list.files(path=path_data)
  date1 = as.numeric(as.date(deb,order = "dmy"))
  date2 = as.numeric(as.date(fin,order = "dmy"))
  
  # Fonction pour sélectionner les cours entre les dates que l'on veut
  selection <- function(data){
    if (!is.null(date1) & !is.null(date2)){
      subset(data, subset=((Date >= date1) & (Date <= date2))) -> data
    }
    else {if (!is.null(date1)) {
      subset(data, subset=(Date >= date1)) -> data
    }
    else {if (!is.null(date2)){
      subset(data, subset=(Date <= date2)) -> data
    }}}
    return(data)
  }
  
  file=paste(path_data,files[1],sep="")
  
  # on lit les tableaux via read
  data <- read.csv(file,sep=";",header=TRUE)
  data = data[,c("Date","Close")]
  data[,"Date"] = sapply(data[,"Date"], function(x){as.date(as.character(x),order = "dmy")})
  data = selection(data)
  colnames(data)[2] <- substring(files[1],1,nchar(files[1])-6)
  
  periode = array(dim = c(3,46))
  periode[1,1] = colnames(data)[2]
  periode[2:3,1] = periode(data)
  
  for (i in 2:length(files)){
    file=paste(path_data,files[i],sep="")
    temp <-read.csv(file,sep=";",header=TRUE)
    temp = temp[,c("Date","Close")]
    temp[,"Date"] = sapply(temp[,"Date"], function(x){as.date(as.character(x),order = "dmy")})
    colnames(temp)[2] <- substring(files[i],1,nchar(files[i])-6)
    
    periode[1,i] = colnames(temp)[2]
    periode[2:3,i] = periode(temp)
    data = merge(data, temp, by = "Date", all = TRUE, sort = FALSE)
  }
  
  data = data[,2:dim(data)[2]]
  data = data[,colSums(is.na(data))<nrow(data)] # On supprime les colonnes entièrement vides
  
  #Je donne des infos sur la taille finale de la df:
  print(paste("La taille de la table finale est : ", as.character(dim(data))))
  print(na.omit(periode))
  if (!is.null(titres)) {
    data = data[,titres]
  }
  return(na.omit(data))
}

# 2) Calcul des log-rendements (ou des rendements standards)
rendements <- function(data) {
  d1 = dim(data)[1]
  d2 = dim(data)[2]
  Y=numeric(d1)
  for (i in 1:d2){
    
    # Rendements standards
    #Y0=(data[1:(d1-1),i]-data[2:d1,i])/data[2:d1,i]
    
    # Log-rendements
    Y0=log(data[1:(d1-1),i]/data[2:d1,i])
    Y=cbind(Y,Y0)
  }
  Y=Y[,2:ncol(Y)]
  colnames(Y)=colnames(data)
  return(Y)
}