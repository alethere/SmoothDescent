## ----setup, include = FALSE---------------------------------------------------
# knitr::opts_chunk$set(
#   collapse = TRUE,
#   comment = "#>"
# )
library(SmoothDescent)
library(knitr)

## ----eval = F-----------------------------------------------------------------
#  devtools::install_github("https://github.com/Alethere/SmoothDescent")
#  #to install the vignette as well use (takes some time)
#  devtools::install_github("https://github.com/Alethere/SmoothDescent", build_vignettes = T)

## -----------------------------------------------------------------------------
sd1 <- smooth_descent(geno = geno, homologue = hom, map = map,
                      ploidy = 2, p1name = "P1", p2name = "P2")


#The observed IBD matrix
graphical_genotype(sd1$obsIBD)
#The predicted IBD matrix
graphical_genotype(sd1$predIBD)
#The error matrix
graphical_genotype(sd1$error)
#The new IBD matrix
graphical_genotype(sd1$newIBD)

knitr::kable(sd1$newgeno[1:10,1:10])

## ---- eval = F----------------------------------------------------------------
#  #There's no preliminary map made in this example, so it's calculated at the beggining
#  #We do not evaluate this chunk to build the vignette faster
#  #Results equivalent to sd_iter$iter1
#  sm1 <- smooth_map(geno = geno, homologue = hom ,
#                    ploidy = 2, p1name = "P1", p2name = "P2",
#                    mapping_ndim = 3)

## -----------------------------------------------------------------------------
#Using the stored object to build the vignette faster
sm1 <- sd_iter$iter1
#Same outputs as before, plus a newmap, recdist and r2 parameters
graphical_genotype(sm1$predIBD)

## -----------------------------------------------------------------------------
#This indicates the relationship between recombination frequency and distance of the new map
recdist_plot(sm1$recdist)

#This compares the old and new maps obtained
#tau indicates order correlation between maps
iterplot(list(premap = sm1$oldmap,iter1 = sm1$newmap),
         main = "Two map comparison")

## ----eval = F-----------------------------------------------------------------
#  #We do not run this code to build the vignette faster
#  #results are equivalent to sd_iter object
#  sm5 <- smooth_map_iter(geno = geno, homologue = hom, iters = 5,
#                      ploidy = 2, p1name = "P1", p2name = "P2",
#                      mapping_ndim = 3, verbose = F)
#  

## -----------------------------------------------------------------------------
#Use stored object for a faster vignette building
sm5 <- sd_iter

#The results of each iteration are equivalent to the output of smooth_map()
graphical_genotype(sm5$iter1$predIBD)

#Use extract() to access results from all iterations
#To access the r2 of each iteration
r2 <- extract(sm5,"r2")
r2

#Or the newmaps of each iteration
newmaps <- extract(sm5,"newmap", simplify = F)
iterplot(newmaps)

## -----------------------------------------------------------------------------
#A genotype matrix with column names and row names
knitr::kable(geno[1:15,1:15])

#A homologue matrix with at least row names
knitr::kable(hom[1:15,])

#A map data.frame with at least two columns: marker and position
knitr::kable(map[1:15,])

## -----------------------------------------------------------------------------
sd1 <- smooth_descent(geno,hom,map, 
                      ploidy = 2, p1name = "P1", p2name = "P2")

## ---- eval = F----------------------------------------------------------------
#  #Result equivalent to sd_iter$iter1
#  sm1 <- smooth_map(geno,hom,map,
#                    ploidy = 2, p1name = "P1", p2name = "P2")

## ---- eval = F----------------------------------------------------------------
#  #Result equivalent to sd_iter
#  sm5 <- smooth_map_iter(geno,hom,map, iters = 5,
#                         ploidy = 2, p1name = "P1", p2name = "P2",
#                         verbose = F)

## -----------------------------------------------------------------------------
#Outputs of smooth descent
names(sd1)

## -----------------------------------------------------------------------------
kable(sd1$obsIBD[[1]][1:15,1:10])
kable(sd1$predIBD[[1]][1:15,1:10],digits = 3)

## ---- fig.width= 3, fig.height=5----------------------------------------------
predicted06 <- extract(sd1$predIBD,"Ind06")
kable(predicted06[1:15,])
#This individual inherited homologues 2 and 4

## -----------------------------------------------------------------------------
all(sd1$error[[1]] == sd1$error[[2]])

## -----------------------------------------------------------------------------
sum(sd1$oldgeno != sd1$newgeno)
mean(sd1$oldgeno != sd1$newgeno)

## -----------------------------------------------------------------------------
#Predicted errors
sum(sd1$error[[1]])
#Corrected errors
sum(sd1$oldgeno != sd1$newgeno)

## -----------------------------------------------------------------------------
names(sm1)
names(sm5$iter1)

## -----------------------------------------------------------------------------
#Shows difference between the names of the two elements
setdiff(names(sm1),names(sd1))

knitr::kable(sm1$newmap[1:15,])
knitr::kable(sm1$recdist[1:15,])
sm1$r2
sm1$tau
sm1$eliminated #In this case no markers have been eliminated

## -----------------------------------------------------------------------------
r2 <- extract(sm5,"r2")

#I usually round the r2, since many iterations might have very similar r2. 
#We normally want a relatively high r2 but a low iteration number.
which.max(round(r2,2))

## -----------------------------------------------------------------------------
#results from  smooth_map()
names(sm1$rec) #One element per type of IBD matrix
names(sm1$rec$pred) #One element for each parent
names(sm1$rec$pred$p1) #Two outputs

#Let's look at the predicted recombinations
p1rec <- sm1$rec$pred$p1
p1rec$individual
kable(p1rec$on_map)

#Let's compare with observed recombinations
p1rec_obs <- sm1$rec$obs$p1
#Clearly, some individuals have a high number of genotyping errors
#look at Ind19 for example. 
kable(data.frame(predicted = p1rec$individual[1:20],
          observed = p1rec_obs$individual[1:20]))


## -----------------------------------------------------------------------------
#Here we see the IBD matrices for all homologues and all individuals
graphical_genotype(sd1$predIBD)
#Or for all individuals but a single homologue
graphical_genotype(sd1$predIBD$H1_P1, main = "Only for homologue 1")

## ---- fig.width= 3------------------------------------------------------------
ibd <- extract(sd1$pred,c("Ind02"))
#Title and colours can also be easily changed.
graphical_genotype(ibd, main = "Individual 02", 
                   col = hcl.colors(300,"Dark Mint"))

## -----------------------------------------------------------------------------
ibd <- extract(sd1$pred,c("Ind02","Ind14","Ind07"))
graphical_genotype(ibd, col = hcl.colors(300,"Sunset"))

## -----------------------------------------------------------------------------
#We can create a list of maps from the output of smooth_map_iter()
maplist <- extract(sm5,"newmap",simplify = F)
#With this we don't see the first map (iteration 0)
iterplot(maplist)

maplist <- c(sm5$iter1["oldmap"],maplist)
iterplot(maplist,main = "Comparison of all maps")

## -----------------------------------------------------------------------------
par(mfrow = c(1,2))
recdist_plot(sm5$iter1$recdist)
recdist_plot(sm5$iter5$recdist)

## -----------------------------------------------------------------------------
#Here we will test the effect of a different prediction intervals
sd_default <- smooth_descent(geno = geno, homologue = hom, map = map,
                             ploidy = 2, p1name = "P1", p2name = "P2",
                             verbose = F)
sd_int_large <- smooth_descent(geno = geno, homologue = hom, map = map,
                             ploidy = 2, p1name = "P1", p2name = "P2",
                             prediction_interval = 30,verbose = F)
sd_int_small <- smooth_descent(geno = geno, homologue = hom, map = map,
                             ploidy = 2, p1name = "P1", p2name = "P2",
                             prediction_interval = 2,verbose = F)
H1 <- list(default = sd_default$predIBD$H1_P1[,1:5],
           larger_interval = sd_int_large$predIBD$H1_P1[,1:5],
           smaller_interval = sd_int_small$predIBD$H1_P1[,1:5])

graphical_genotype(H1, main = "IBD of homologue 1 for 5 individuals")

## -----------------------------------------------------------------------------
#Here we will test the effect of different prediction thresholds.
sd_default <- smooth_descent(geno = geno, homologue = hom, map = map,
                             ploidy = 2, p1name = "P1", p2name = "P2",
                             verbose = F)
#We will only accept imputations if IBD probabilities are 1
sd_th_large <- smooth_descent(geno = geno, homologue = hom, map = map,
                             ploidy = 2, p1name = "P1", p2name = "P2",
                             prediction_threshold = 1,verbose = F)
#We accept any imputation
sd_th_small<- smooth_descent(geno = geno, homologue = hom, map = map,
                             ploidy = 2, p1name = "P1", p2name = "P2",
                             prediction_threshold = 0,verbose = F)

#The proportion of changes
def_change <- mean(sd_default$newgeno != sd_default$oldgeno)
large_change <- mean(sd_th_large$newgeno != sd_default$oldgeno)
small_change <- mean(sd_th_small$newgeno != sd_default$oldgeno)

print(paste0("By default we impute ",round(def_change*100,2),"% of genotypes"))
print(paste0("With increased stringency we impute ",
             round(large_change*100,2),"% of genotypes"))
print(paste0("With decreased stringency we impute ",
             round(small_change*100,2),"% of genotypes"))

## -----------------------------------------------------------------------------
#Here we will test the effect of error thresholds
sd_default <- smooth_descent(geno = geno, homologue = hom, map = map,
                             ploidy = 2, p1name = "P1", p2name = "P2",
                             verbose = F)
#Only if the error probability is above 0.95 we consider it an error (very lax)
sd_th_large <- smooth_descent(geno = geno, homologue = hom, map = map,
                             ploidy = 2, p1name = "P1", p2name = "P2",
                             error_threshold = 0.95,verbose = F)
#Any marker with an error probability above 0.5 is considered an error (very stringent)
sd_th_small<- smooth_descent(geno = geno, homologue = hom, map = map,
                             ploidy = 2, p1name = "P1", p2name = "P2",
                             error_threshold = 0.5,verbose = F)

#The proportion of errors
#All error matrices have the same markers identified as errors
def_change <- mean(sd_default$error$H1_P1)
large_change <- mean(sd_th_large$error$H1_P1)
small_change <- mean(sd_th_small$error$H1_P1)

print(paste0("By default we consider ",round(def_change*100,2),"% of genotypes as errors"))
print(paste0("With decreased stringency we cosider ",
             round(large_change*100,2),"% of genotypes as errors"))
print(paste0("With increased stringency we consider ",
             round(small_change*100,2),"% of genotypes as errors"))

## -----------------------------------------------------------------------------
#For a single marker
calc_IBD(geno = 1,p1hom = c(h1 = 0, h2 = 1), p2hom = c(h3 = 0,h4 = 0), ploidy = 2)
#Tetraploid probabilities
calc_IBD(geno =2,p1hom = c(0,0,1,0),p2hom = c(0,1,0,1), ploidy = 4)
#Impossible genotype
calc_IBD(geno = 2,p1hom = c(0,1),p2hom = c(0,0), ploidy = 2)

#For a single individual the genotype input must be a one column matrix
calc_IBD(geno[1:10,3,drop = F],
         p1hom = hom[1:10,1:2],
         p2hom = hom[1:10,3:4],ploidy = 2)

## -----------------------------------------------------------------------------
#We calculate IBDs removing the genotypes of the parents
IBD <- calc_IBD(geno[,-1:-2],p1hom = hom[,1:2],p2hom = hom[,3:4])

#One parental homologue for one individual: a named vector
ibd <- predict_IBD(IBD$H1_P1[,3],map = map)
graphical_genotype(matrix(ibd,ncol = 1))

#Many individuals, a single homologue: a matrix
ibd <- predict_IBD(IBD$H1_P1,map = map)
graphical_genotype(ibd, main = "A single homologue")

#All homologues of one individual
ibd <- predict_IBD(extract(IBD,"Ind02"),map = map)
graphical_genotype(ibd, main = "A single individual")

#All homologues of all individuals: a list of IBD matrices
ibd <- predict_IBD(IBD,map = map)
graphical_genotype(ibd)

#Multiple individuals
ibd <- predict_IBD(extract(IBD,c("Ind02","Ind14","Ind07")),map = map)
graphical_genotype(ibd)

## -----------------------------------------------------------------------------
#We calculate IBDs removing the genotypes of the parents
IBD <- calc_IBD(geno[,-1:-2],p1hom = hom[,1:2],p2hom = hom[,3:4])
#In this case there's no change to the IBD matrix so there's no change in genotypes
newgeno <- genotype(IBD,homologue = hom)
all(newgeno == geno[,-1:-2])

#For one individual
ibd <- extract(IBD,"Ind02")
gen <- genotype(ibd,hom)
#Again all genotypes are the same as the originals
all(gen == geno[,"Ind02"])

