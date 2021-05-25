#Data making
geno <- readRDS("../Strawberry Genotyping/Data/genotyped/processed/final_binlist_1.RDS")

physical <- sapply(geno$name,function(b){
  pos <- as.numeric(extract(strsplit(b,":"),2))
  mean(pos,trim = 0.1)
})
physical <- physical[geno$seq == "1A"]
geno <- subset(geno$geno,geno$seq == "1A")
hom <- readRDS("../Strawberry Genotyping/Data/map/all_chrom/homologue_matrix_chrom1.RDS")
hom <- hom$`1A`
map <- readRDS("../Strawberry Genotyping/Data/map/all_chrom/preliminary_maps_chrom1.RDS")
map <- map$`1A`

map <- map[order(map$position),]
map <- map[,c("marker","position")]
map$physical <- physical[map$marker]
geno <- geno[map$marker,]
hom <- hom[map$marker,]

inds <- extract(strsplit(colnames(geno),"_"),1)[-1:-2]
numcode <- as.numeric(factor(inds))
inds <- sprintf("Ind%02d",numcode)
inds[duplicated(inds)] <- paste0(inds[duplicated(inds)],"_2")
colnames(geno) <- c(c("P1","P2"),inds)
geno <- geno[, c("P1","P2",sort(inds))]


map$marker <- paste0("m",1:nrow(map))
rownames(geno) <- map$marker
rownames(hom) <- map$marker

colnames(hom) <- c("H1_P1","H2_P1","H3_P2","H4_P2")

write.table(map,file = "test/test_map.txt",quote = F)
write.table(geno,file = "test/test_geno.txt",quote = F)
write.table(hom,file = "test/test_hom.txt",quote = F)
