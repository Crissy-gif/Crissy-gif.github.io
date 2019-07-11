rm(list = ls(all=TRUE))
setwd("~/Documents/Pelago_Project/Multiple Datasets/UCP_data")
source("~/Documents/Functions:Templates/Functions.R")

library(tidyr)
library(dplyr)
library(stringr)
library(tidyverse)

#------------------load data--------------------------------

#file <- read.csv("2DDR_AZ_1808_UCP_proteins.csv", sep = ";")
file <- read.csv("2DDR_AZ_1808_UCP_proteins.csv", sep = ";", colClasses = c(T1='character', T2='character', T3='character', T4='character', T5='character', 
                                                                         T6='character', T7='character', T8='character', T9='character', T10='character'))

#-------------------variable selection----------------------
# unite() paste together multple columns into one. T1:T10 are individual columns and we wnat them in the one column.

vari_sel <- file

vari_sel$T_combo <- vari_sel %>%
  unite("T_combo", T1:T10, sep=",") %>%
  .$T_combo %>%
  str_split(",") %>%
  lapply(FUN = function(val){
      val %>%
        str_trim %>%
        as.numeric
    })

vari_sel2 <- vari_sel %>%
  select(Acc.Number, T_combo)
  
####################################################################################################
#R modified

protein_dist_3d <- function(protein1, protein2){
  amp_length <- length(protein1$T_combo)
  
  amp1 <- protein1$T_combo[[1]]
  amp2 <- protein2$T_combo[[1]]
  amp1[which(is.nan(amp1))] <- 0
  amp2[which(is.nan(amp2))] <- 0
  rmse <- sqrt(mean((amp1 - amp2)^2))
  
  #rmse_ctl_diff <- sqrt(mean(((protein1$amp - protein1$amp_ctl) - (protein2$amp - protein2$amp_ctl))^2))
  #total <- rmse + rmse_ctl_diff
  
  total <- rmse
  total
}

protein_dist_3d_na2 <- function(protein1, protein2){
  amp_length <- length(protein1$T_combo)
  amp1 <- protein1$T_combo[[1]]
  amp2 <- protein2$T_combo[[1]]
  nans <- is.nan(amp1) | is.nan(amp2) # All indices with nans in either protein1 or protein 2
  rmse <- sqrt(mean((amp1[!nans] - amp2[!nans])^2))
  rmse
}

#######################################################################################################
library(doParallel)
library(foreach)
n_threads <- 10
cl <- makeCluster(n_threads) # 10 threads
doParallel::registerDoParallel(cl)
foreach::getDoParRegistered()
source("~/Documents/Functions:Templates/Functions.R")
#a_cool_dist_matey <- dist_matrix(proteins = head(vari_sel2, 50), protein_dist_function = protein_dist_3d, n_batches = 10)
a_cool_dist_mateyer <- dist_matrix(proteins = vari_sel2, protein_dist_function = protein_dist_3d_na2, n_batches = 50)
#a_cool_dist_mateyer2 <- dist_matrix(proteins = vari_sel2, protein_dist_function = protein_dist_3d_na2, n_batches = 50)

xy <- my_mds(a_cool_dist_mateyer)
which(is.na(a_cool_dist_mateyer))
sum(is.na(a_cool_dist_mateyer))
hupp <- apply(a_cool_dist_mateyer, MARGIN=1, FUN=function(x) sum(is.na(x)))
length(hupp[hupp == 0])
length(hupp)
a_not_cool_dist_mateyer <- a_cool_dist_mateyer[hupp == 0, hupp == 0]
xy <- my_mds(a_not_cool_dist_mateyer)
vari_sel2[hupp==0, 'mds_x'] <- xy$x
vari_sel2[hupp==0, 'mds_y'] <- xy$y
plot(xy$x,xy$y)

write.csv(vari_sel2[,c('Acc.Number', 'mds_x', 'mds_y')], "a_non_cool_thing_stuff1.csv")

#######################################################################################################
#######################################################################################################
#######################################################################################################



#Do string split on each temp
#Modify dist function
 #Use the mean in the each temp
#convert them in an array

# dist_T

protein_dist <- function(protein1, protein2){
  amp_length <- length(protein1$Amplitude)
  
  amp1 <- protein1$Amplitude2[[1]]
  amp2 <- protein2$Amplitude2[[1]]
  amp1[which(is.nan(amp1))] <- 0
  amp2[which(is.nan(amp2))] <- 0
  rmse <- sqrt(mean((amp1 - amp2)^2))
  
  #rmse_ctl_diff <- sqrt(mean(((protein1$amp - protein1$amp_ctl) - (protein2$amp - protein2$amp_ctl))^2))
  #total <- rmse + rmse_ctl_diff
  
  total <- rmse
  total
}

#############################testing#################################

x1 <- vari_sel_clean[1, 1:11]
x2 <- vari_sel_clean[2, 1:11]


##################################################################

m <- x1$T1_1[[1]]
u <- x2$T1_1[[1]]
rmse <- sqrt(mean((m - u)^2))
##################################################################
c <- function(y, y2){
 
  for (i in y){
    for (j in y2){

      k1 <- y$i[[1]]
      k2 <- y2$j[[1]]
      rmse <- sqrt(mean((k1 - k2)^2))
      total <- rmse
     
    }
total
  }
  
}

#library(magicfor)
#magic_for()


for(i in x1[2:8]){
  for(j in x2[2:8]){
    k1 <- i[[1]]
    k2 <- j[[1]] 
    rmse <- sqrt(mean((k1 - k2)^2))
   results <- rmse
    results
  }
}


######################################################

#Try out tSNE  on the matrix
#Find out how to get rid of the zeros for the cmds scale.

########################################################################################################
######################################################################################################
install.packages("tsne")
library(tsne)

tsne_out <- tsne(a_not_cool_dist_mateyer, perplexity = 1000, max_iter = 5000)
 write.csv (tsne_out, "tsne_out_1000.csv")
######################################################################################################
z <- 0#I dont know what the heck is this!
epc <- function(z) {
  z <<- z + 1
  filename <- paste("d:\\plot", x, "jpg", sep=".")
  cat("> Plotting TSNE to ", filename, " ")
  
  # plot to d:\\plot.x.jpg file of 2400x1800 dimension
  jpeg(filename, width=2400, height=1800)
  
  plot(z, t='n', main="T-SNE")
  text(z, labels=rownames(mydata))
  dev.off()
}

tsne_data <- tsne(tsne_out, k=5, epoch_callback=epc, max_iter=1000, epoch=100)



########################################################################################################

