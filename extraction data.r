###### Ce fichier contient les fonctions permettant d'extraire nos donn?es contenues dans des dossiers adjacents
# Il permet aussi de calculer les rendements, outils de travail au m?me titre que les cours des actions eux m?mes (d'o? leur pr?sence dans ce fichier)

# la gestion des dates pr?sentes dans les tableaux des cours necessite le package suivant
library("date")

periode <- function(data){
  return(rbind(date.mmddyy(min(data[,"Date"])),date.mmddyy(max(data[,"Date"]))))
}


# 1) Renvoie un tableau dont les colonnes sont les valeurs de cloture de differents titres

# titres = notre portefeuille
# date1 et date2 = les intervalles choisis, par d?faut toute la p?riode commune aux titres
# type = J, H ou M (journalier, hebdomadaire ou mensuel)

load_data <- function(deb, fin, titres = NULL, type = "J"){
  # le path est ? changer selon chaque utilisateur
  path_data="/Users/Bouxtehouve/Documents/ENSAE/2A/Projet Statapp/Cours/"
  path_data=paste(path_data,type,"/",sep="")
  
  # on s?lectionne les fichiers contenus dans le dossier indiqu?
  titres <- sort(titres)
  files <- list.files(path=path_data)
  files <- sort(files)
  files2 <- c()
  for (i in 1:length(files)){
    files2[i] <- substring(files[i],1,nchar(files[i])-6)
  }
  date1 = as.numeric(as.date(deb,order = "dmy"))
  date2 = as.numeric(as.date(fin,order = "dmy"))
  
  # Fonction pour s?lectionner les cours entre les dates que l'on veut
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
    if (files2[i] %in% titres){
      file=paste(path_data,files[i],sep="")
      temp <-read.csv(file,sep=";",header=TRUE)
      temp = temp[,c("Date","Close")]
      temp[,"Date"] = sapply(temp[,"Date"], function(x){as.date(as.character(x),order = "dmy")})
      temp <- selection(temp)
      colnames(temp)[2] <- substring(files[i],1,nchar(files[i])-6)
      
      periode[1,i] = colnames(temp)[2]
      periode[2:3,i] = periode(temp)
      data = merge(data, temp, by = "Date", all = TRUE, sort = FALSE)
    }
  }
  
  data <- data[order(data$Date),]
  rownames(data) <- data$Date
  #data = data[,2:dim(data)[2]]
  data = data[,colSums(is.na(data))<nrow(data)] # On supprime les colonnes enti?rement vides
  
  #Je donne des infos sur la taille finale de la df:
  print("La taille de la table finale (avant suppression des lignes manquantes) est : ")
  print(dim(data))
  print(periode[,colSums(is.na(periode))<nrow(periode)])
  data$Date <- NULL
  #if (!is.null(titres)) {
  #  data = data[,titres]
  #}
  return(na.omit(data))
}

# 2) Calcul des log-rendements (ou des rendements standards)
rendements <- function(data) {
  d1 = dim(data)[1]
  d2 = dim(data)[2]
  Y=numeric(d1)
  for (i in 1:d2){
    
    # Rendements standards
    #Y0=(data[2:d1,i]-data[1:(d1-1),i])/data[1:(d1-1),i]
    
    # Log-rendements
    Y0=log(data[2:d1,i]/data[1:(d1-1),i])
    Y=cbind(Y,Y0)
  }
  Y=Y[,2:ncol(Y)]
  colnames(Y)=colnames(data)
  return(Y)
}

# 3) Calcul des rendements journaliers actualis?s des OAT et agr?gat avec les rendements classiques de ci-dessus
global_return<-function(data,type="J"){
  actu=365
  if (type=="M"){actu=12}
  d=(1+data[,"oat"]/100)^(1/(10*actu))-1
  Y <- rendements(data)
  Y[,'oat'] <- d
  return(Y)
}