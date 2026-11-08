#!/usr/bin/env Rscript
# Recompute the frozen26patient ECM covariate from the original1717gene panel.
# Run after merging Code_resources_update into the analysis repository.
args <- commandArgs(trailingOnly=TRUE)
repo <- if(length(args)>=1L)normalizePath(args[[1]],winslash="/",mustWork=TRUE)else getwd()
out <- if(length(args)>=2L)args[[2]]else file.path(repo,"validation_patient_ecm")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(edgeR)})
root <- file.path(repo,"resources","ecm_direction")
d <- fread(file.path(root,"GSE183904_primary_target_panel_counts.tsv.gz"))
m <- fread(file.path(root,"GSE183904_primary_sample_patient_map.tsv"))
frozen <- fread(file.path(root,"GSE183904_primary_tumor_frozen_ECM_scores.tsv"))
sets <- fread(file.path(root,"frozen_ECM_highlow50.tsv"))
stopifnot(nrow(m)==26L,uniqueN(d$gene_symbol)==1717L,all(d$sample_id %in% m$sample_id))
counts <- d[,.(count=sum(count)),by=.(sample_id,gene_symbol)]
wide <- dcast(counts,gene_symbol~sample_id,value.var="count",fill=0)
mat <- as.matrix(wide[,-"gene_symbol"]);rownames(mat)<-wide$gene_symbol;storage.mode(mat)<-"integer"
y <- calcNormFactors(DGEList(mat),method="TMM")
lcpm <- cpm(y,log=TRUE,prior.count=1)
# Same sample-wise arrangement and gene-wise sample SD as the original script16.
x <- t(lcpm)
z <- matrix(NA_real_,nrow(x),ncol(x),dimnames=dimnames(x))
for(j in seq_len(ncol(x))){v<-x[,j];mu<-mean(v,na.rm=TRUE);s<-sd(v,na.rm=TRUE)
 if(is.finite(s)&&s>0)z[,j]<-(v-mu)/s else z[is.finite(v),j]<-0}
groups <- split(sets$gene_symbol,sets$ecm_group)
scores <- lapply(groups,function(gs){g<-intersect(unique(gs),colnames(z));stopifnot(length(g)>=2L);rowMeans(z[,g,drop=FALSE])})
result <- data.table(sample_id=rownames(z),ECM_C1_HIGH50=scores$C1_HIGH50,ECM_C2_HIGH50=scores$C2_HIGH50)
result[,ECM_direction:=ECM_C2_HIGH50-ECM_C1_HIGH50]
check <- merge(result,frozen[,.(sample_id,ECM_direction_frozen=ECM_direction)],by="sample_id")
stopifnot(nrow(check)==26L)
err <- max(abs(check$ECM_direction-check$ECM_direction_frozen))
stopifnot(err<1e-10)
fwrite(check,file.path(out,"patient_ECM_recomputed_comparison.tsv"),sep="\t")
fwrite(data.table(sample_id=colnames(y$counts),target_panel_library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors,effective_library_size=y$samples$lib.size*y$samples$norm.factors),file.path(out,"patient_ECM_TMM_factors.tsv"),sep="\t")
fwrite(as.data.table(t(lcpm[unique(sets$gene_symbol),,drop=FALSE]),keep.rownames="sample_id"),file.path(out,"patient_ECM_highlow100_logCPM.tsv"),sep="\t")
writeLines(c(paste0("patients\t",nrow(check)),paste0("target_panel_genes\t",nrow(mat)),paste0("max_ECM_absolute_difference\t",format(err,digits=17)),paste0("R\t",R.version.string),paste0("edgeR\t",as.character(packageVersion("edgeR")))),file.path(out,"verification.tsv"))
print(data.table(patients=nrow(check),target_panel_genes=nrow(mat),max_ECM_absolute_difference=err))
