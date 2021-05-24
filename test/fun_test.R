##### Test file #####
#This file uses all functions and checks the output in search for errors.

for(file in list.files("R/",full.names = T)){ source(file)}

#Calc IBD tests -------
source("R/Utils.R")
dos <- as.matrix(read.table("test/test_geno.txt",header = T))
dos <- dos[,-1:-2] #we take out parental genotypes
hom <- as.matrix(read.table("test/test_hom.txt",header = T))
#Numeric
p2_test <- calc_IBD(dos[1,1],hom[1,1:2],hom[1,3:4],ploidy = 2)

if(!is.matrix(p2_test)){
  stop("Output of calc_IBD.numeric is not matrix")
}
if(!ncol(p2_test) == 2*2){
  stop("Length of calc_IBD.numeric output is not ploidy*2")
}
if(!all(p2_test == c(0.5,0.5,1,0))){
  stop("Output of calc_IBD.numeric does not coincide with saved output")
}

#Matrix
ploidy = 2
matrix_test <- calc_IBD(dos,p1hom = hom[,1:2],p2hom = hom[,3:4],ploidy = 2)
#saveRDS(matrix_test,"test/calc_IBD.matrix.RDS")


if(!is.list(matrix_test)){
  stop("Output of calc_IBD.matrix is not list")
}
if(!length(matrix_test) == ploidy*2){
  stop("Length of calc_IBD.matrix output is not ploidy*2")
}
if(!all(sapply(matrix_test,ncol) == ncol(dos))){
  stop("Not all elements from output of calc_IBD.matrix have ncol == ncol(dos)")
}
saved_mtest <- readRDS("test/calc_IBD.matrix.RDS")
if(!identical(matrix_test,saved_mtest)){
  stop("Output of calc_IBD.matrix does not coincide with saved output")
}


#Diagnostics ----------
source("R/Utils.R")
maplist <- list(pre = read.table("test/test_map.txt"),
                flat = readRDS("test/test_flat_map.RDS")$locimap,
                sphere = readRDS("test/test_sphere_map.RDS")$locimap)
maplist$flat$marker <- maplist$flat$locus
maplist$sphere$marker <- maplist$sphere$locus

#Rec count
ibdlist <- readRDS("test/calc_IBD.matrix.RDS")

#Matrix method
rc_mat <- rec_count(ibd = ibdlist[[1]],map = maplist$pre)
#saveRDS(rc_mat,"test/rec_count.matrix.RDS")
rc_test <- readRDS("test/rec_count.matrix.RDS")
if(!identical(rc_mat,rc_test)){
  stop("Rec_count.matrix output does not coincide with stored output")
}
if(sum(rc_mat$individual) != sum(rc_mat$on_map$count)){
  stop("Recombination counts per individual and per window do not coincide")
}

#List method
rc_list <- rec_count(ibd = ibdlist,map = maplist$pre)
#saveRDS(rc_list,"test/rec_count.list.RDS")
rc_test <- readRDS("test/rec_count.list.RDS")
if(!identical(rc_list,rc_test)){
  stop("Rec_count.mlist output does not coincide with stored output")
}

#Genotype prediction -------------
source("R/Utils.R")
ibd <- readRDS("test/pred_IBD.list.RDS")
hom <- read.table("test/test_hom.txt",header = T)

#Test of list method
geno <- genotype(ibd,hom,ploidy = 2)
#saveRDS(geno,"test/test_geno.list.RDS")
test_geno <- readRDS("test/test_geno.list.RDS")
if(!identical(geno,test_geno)) stop("Output from genotype.list does not coincide with stored result")
if(!all(ncol(geno) == sapply(ibd,ncol))) stop("Output from genotype.list does not contain as many columns as IBD input")

#Test NA production
ib <- lapply(ibd,function(i) i > 0.8)
na_test <- all(is.na(geno) == (Reduce('+',ib) < 2))
if(!na_test) stop("Output from genotype.list does not contain NAs in the correct places")

#Matrix method
ind_ibd <- sapply(ibd,'[',,1)
geno <- genotype(ind_ibd,homologue = hom)
#saveRDS(geno,"test/test_geno.matrix.RDS")
test_geno <- readRDS("test/test_geno.matrix.RDS")
if(!identical(geno,test_geno)) stop("Output from genotype.matrix does not coincide with stored result")
if(!length(geno) == nrow(ind_ibd)) stop("Output from genotype.matrix does not contain as many genotypes as IBD input")

#Test NA production
should_be_na <- !rowSums(ind_ibd > 0.8) == 2
is_na <- is.na(geno)
if(!all(should_be_na == is_na)){
  stop("Output from matrix.list does not contain NAs in the correct palces")
}

#Mapping ---------------
geno <- read.table("test/test_geno.txt",header = T)

#This tests both linkpar and linkdf_shortcut function
linkdf <- linkdf_shortcut(geno,p1name = "P1",p2name = "P2",ploidy = 2, ncores = 1)
#saveRDS(linkdf,"test/test_linkdf.RDS")
linkdf_test <- readRDS("test/test_linkdf.RDS")
if(!identical(linkdf,linkdf_test)){
  stop("Output of linkdf_shortcut does not match the stored output")
}

#Calculating number of informative pair-wise markers
type <- sapply(1:nrow(geno),function(i) paste0(geno$P1[i],geno$P2[i]))
type_count <- table(type)
self_links <- type_count * (type_count -1 )/2
#This is only valid for diploids
other_link <- (type_count[3]*type_count[1:2]  )
if(sum(self_links,other_link) != nrow(linkdf)){
  stop("Not all informative pair-wise combinations have been calculated in linkdf_shortcut function.")
}

#Mapping test
flat_map <- mdsmap(linkdf,ndim = 2)
#saveRDS(flat_map,"test/test_flat_map.RDS")
sphere_map <- mdsmap(linkdf,ndim = 3)
#saveRDS(sphere_map,"test/test_sphere_map.RDS")

test_map <- readRDS("test/test_flat_map.RDS")
if(!identical(flat_map,test_map)){
  warning("2D map output does not coincide with stored map.\n This might be due to random variation in the EM algorithm. Check map output.")
}

test_map <- readRDS("test/test_sphere_map.RDS")
if(!identical(sphere_map,test_map)){
  warning("3D map output does not coincide with stored map.\n This might be due to random variation in the EM algorithm. Check map output.")
}

#Prediction of IBD -------------
source("R/Utils.R")
map <- read.table("test/test_map.txt",header = T)
IBD <- readRDS("test/calc_IBD.matrix.RDS")

#List
pred <- predict_IBD(IBD,map,interval = 10)
#saveRDS(pred,"test/pred_IBD.list.RDS")
if(!is.list(pred)) stop("Output of predict_IBD.list is not list")
if(any(sapply(IBD,dim) != sapply(pred,dim))) stop("Output of predict_IBD.list does not match input")
pred_test <- readRDS("test/pred_IBD.list.RDS")
if(!identical(pred,pred_test)) stop("Output of predict_IBD.list does not match stored output")

#Matrix
ib <- IBD$H1_P1[,-1:-2]
pred <- predict_IBD(ib,map,interval = 10)
#saveRDS(pred,"test/pred_IBD.matrix.RDS")
if(!is.matrix(pred)) stop("Output of predict_IBD.matrix is not matrix")
if(any(dim(ib) != dim(pred))) stop("Output of predict_IBD.matrix does not match input")
pred_test <- readRDS("test/pred_IBD.matrix.RDS")
if(!identical(pred,pred_test)) stop("Output of predict_IBD.matrix does not match stored output")

#Numeric
ib <- IBD$H1_P1[map$marker,3]
pred <- predict_IBD(ib,map,interval = 10)
#saveRDS(pred,"test/pred_IBD.numeric.RDS")
if(!is.numeric(pred)) stop("Output of predict_IBD.numeric is not numeric")
if(length(pred) != nrow(map)) stop("Output of predict_IBD.numeric has wrong length")
pred_test <- readRDS("test/pred_IBD.numeric.RDS")
if(!identical(pred,pred_test)) stop("Output of predict_IBD.numeric does not match stored output")

#Smooth wrappers ----------------
geno <- read.table("test/test_geno.txt",header = T)
map <- read.table("test/test_map.txt",header = T)
hom <- read.table("test/test_hom.txt",header = T)

smooth_test <- smooth_descent(geno,map,homologue = hom,ploidy = 2,
                              p1name = "P1",p2name = "P2",verbose = T)

smoothmap_test <- smooth_map(geno,map,homologue = hom,ploidy = 2,
                             p1name = "P1",p2name = "P2",verbose = T,
                             mapping_ndim = 3)


iterplot(smoothmap_test[c("oldmap","newmap")])



