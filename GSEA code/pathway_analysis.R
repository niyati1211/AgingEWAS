# loading packages
library(missMethyl)
library(data.table)
library(dplyr)

# identifying file names for the age-DNAm association results to be analyzed 
file_list <- c("age_EWAS_results_breast_manifest_sv5updated.txt",
               "age_EWAS_results_colon_manifest_sv10updated.txt",
               "age_EWAS_results_lung_manifest_sv10updated.txt",
               "age_EWAS_results_ovary_manifest_sv5updated.txt",
               "age_EWAS_results_prostate_manifest_sv10updated.txt",
               "age_EWAS_results_testis_manifest_sv5updated.txt",
               "age_EWAS_results_wb_manifest_sv5updated.txt")

# importing hallmark geneset database
hallmark <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds"))

for (i in 1:length(file_list)) { 
  print("################")
  print(i)
  print(file_list[i])
  results_final<-read.table(file_list[i],header=T,sep=" ",row.names=1)
  sig.cpg <- (results_final %>% filter(adj.P.Val < 0.05))$Name
  all.cpg <- results_final$Name
  pathway_GO_FDR <- data.frame(gometh(sig.cpg, all.cpg, collection = ("GO"), array.type = c("EPIC")))
  pathway_KEGG_FDR <- data.frame(gometh(sig.cpg, all.cpg, collection = ("KEGG"), array.type = c("EPIC")))
  hallmark_gsa <- gsameth(sig.cpg, all.cpg, collection=hallmark, array.type = c("EPIC"))
  hallmark_gsa <- data.frame(hallmark_gsa)
  hallmark_gsa$Description <- rownames(hallmark_gsa)
  
  write.table(pathway_GO_FDR,file=paste0("/gpfs/data/phs/groups/Projects/GTEx/Age_TL_EWAS/Lin_analysis/EWAS_results/SV_updated_results/age/missMethyl_output/randi_pathway_GO_FDR_",file_list[i]),row.names=F,sep="\t")
  write.table(pathway_KEGG_FDR,file=paste0("/gpfs/data/phs/groups/Projects/GTEx/Age_TL_EWAS/Lin_analysis/EWAS_results/SV_updated_results/age/missMethyl_output/randi_pathway_KEGG_FDR_",file_list[i]),row.names=F,sep="\t")
  write.table(hallmark_gsa,file=paste0("/gpfs/data/phs/groups/Projects/GTEx/Age_TL_EWAS/Lin_analysis/EWAS_results/SV_updated_results/age/missMethyl_output/randi_hallmark_gsa_FDR_",file_list[i]),row.names=F,sep="\t")
  
  print(paste("Number of Significant GO terms (FDR):",nrow(pathway_GO_FDR%>%filter(FDR < 0.05))))
  print(paste("Number of Significant KEGG terms (FDR):",nrow(pathway_KEGG_FDR%>%filter(FDR < 0.05))))
  print(paste("Number of Significant Hallmark GS (FDR):",nrow(hallmark_gsa%>%filter(FDR < 0.05))))
}

# plotting venn diagram of KEGG results
cpg_GS_list <- vector(length=5,mode="list")
for (i in 1:5) {
  tmp_gs_df <- fread(paste0("/Volumes/phs/groups/Projects/GTEx/Age_TL_EWAS/Lin_analysis/EWAS_results/SV_updated_results/age/missMethyl_output/randi_pathway_KEGG_FDR_",file_list[i]),header=T)
  cpg_GS_list[[i]] <- (tmp_gs_df %>% filter(FDR<0.05))$Description
}
png("/Volumes/phs/groups/Projects/GTEx/Age_TL_EWAS/Lin_analysis/EWAS_results/SV_updated_results/age/missMethyl_output/venn_KEGG.png",width = 1000, height = 1000, units = "px", pointsize = 25)
print(venn(cpg_GS_list,ilab=TRUE, zcolor = "style",snames=c("breast","colon","lung","ovary","prostate")))
dev.off()
all_gs<- c(as.vector(cpg_GS_list[1]),
           as.vector(cpg_GS_list[2]),
           as.vector(cpg_GS_list[3]),
           as.vector(cpg_GS_list[4]),
           as.vector(cpg_GS_list[5])
)
all_gs <- unlist(cpg_GS_list)
data.table(table(all_gs) %>% sort(decreasing = T)) %>% filter(N>=2)

# plotting venn diagram of Hallmark gene set results 
cpg_GS_list <- vector(length=5,mode="list")
for (i in 1:5) {
  tmp_gs_df <- fread(paste0("/Volumes/phs/groups/Projects/GTEx/Age_TL_EWAS/Lin_analysis/EWAS_results/SV_updated_results/age/missMethyl_output/randi_hallmark_gsa_FDR_",file_list[i]),header=T)
  cpg_GS_list[[i]] <- (tmp_gs_df %>% filter(FDR<0.05))$Description
}
png("/Volumes/phs/groups/Projects/GTEx/Age_TL_EWAS/Lin_analysis/EWAS_results/SV_updated_results/age/missMethyl_output/venn_HALLMARK.png",width = 1000, height = 1000, units = "px", pointsize = 25)
print(venn(cpg_GS_list,ilab=TRUE, zcolor = "style",snames=c("breast","colon","lung","ovary","prostate")))
dev.off()

all_gs<- c(as.vector(cpg_GS_list[1]),
           as.vector(cpg_GS_list[2]),
           as.vector(cpg_GS_list[3]),
           as.vector(cpg_GS_list[4]),
           as.vector(cpg_GS_list[5])
)
all_gs <- unlist(cpg_GS_list)
data.table(table(all_gs) %>% sort(decreasing = T)) %>% filter(N>=2)

######################################################################
# plotting venn diagram of GO gene set results 
cpg_GS_list <- vector(length=5,mode="list")
for (i in 1:5) {
  tmp_gs_df <- fread(paste0("/Volumes/groups/Projects/GTEx/Age_TL_EWAS/Lin_analysis/EWAS_results/SV_updated_results/age/missMethyl_output/randi_GO_FDR_",file_list[i]),header=T)
  cpg_GS_list[[i]] <- (tmp_gs_df %>% filter(FDR<0.05))$Description
}
png("/Volumes/phs/groups/Projects/GTEx/Age_TL_EWAS/Lin_analysis/EWAS_results/SV_updated_results/age/missMethyl_output/venn_HALLMARK.png",width = 1000, height = 1000, units = "px", pointsize = 25)
print(venn(cpg_GS_list,ilab=TRUE, zcolor = "style",snames=c("breast","colon","lung","ovary","prostate")))
dev.off()

breast <- fread("/Volumes/groups/Projects/GTEx/Age_TL_EWAS/Lin_analysis/EWAS_results/SV_updated_results/age/missMethyl_output/randi_pathway_GO_FDR_age_EWAS_results_breast_manifest_sv5updated.txt")
colon <- 
all_gs<- c(as.vector(cpg_GS_list[1]),
           as.vector(cpg_GS_list[2]),
           as.vector(cpg_GS_list[3]),
           as.vector(cpg_GS_list[4]),
           as.vector(cpg_GS_list[5])
)
all_gs <- unlist(cpg_GS_list)
data.table(table(all_gs) %>% sort(decreasing = T)) %>% filter(N>=2)

