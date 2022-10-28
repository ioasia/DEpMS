### Benchmark (adapted from https://github.com/ibludau/ProteoformAnanlysis/blob/main/PerformanceEvaluation/InterlabBenchmark_final_paper.R)

#' ## Load CCprofiler and setup the work environment
# library(devtools)
# options(warn=-1)

# install_github("CCprofiler/CCprofiler", ref =  "proteoformLocationMapping")
#devtools::load_all("~/Desktop/projects/ProteoformProject/CCprofiler/")
library('CCprofiler')
library('data.table')
library('ggplot2')
library('UpSetR')

setwd("~/Downloads/DEpMS_project/depms_benchmark/")

#####################################
### Create dataset ##################
#####################################

dt <- fread("site02_global_q_0.01_applied_to_local_global.txt") # PXD004886

# rename protein
dt[,ProteinName:=gsub("1/","",ProteinName)]
dt <- subset(dt, ProteinName != "iRT_protein") # remove iRT peptides
dt <- subset(dt, ProteinName != "sp|AQUA30|AQUA30") # remove AQUA peptides

# rename runs
dt[,run:=gsub("Site2_AQUA_HEK_","",run_id)]
dt[,run:=gsub("_180714.mzXML.gz","",run)]

# determine day
dt[,day:=gsub("S.*_SW_","",run)]
dt[,day:=gsub("_.*","",day)]

# subset to important columns
dt_sub <- subset(dt, select=c("run","day","ProteinName","FullPeptideName","peptide_group_label", "Intensity"))

# aggregate charge states
dt_sub[,pep_int:=sum(Intensity), by=c("ProteinName","FullPeptideName","run")]
dt_sub <- unique(subset(dt_sub, select=c("run","day","ProteinName","FullPeptideName","pep_int")))

# subset to sufficient data per condition
dt_sub[,n_day1 := sum(day=="day1"), by=c("FullPeptideName")]
dt_sub[,n_day3 := sum(day=="day3"), by=c("FullPeptideName")]
dt_sub[,n_day5 := sum(day=="day5"), by=c("FullPeptideName")]
dt_sub[,min_n_per_day := min(n_day1,n_day3,n_day5),by=c("FullPeptideName")]
dt_sub <- subset(dt_sub, min_n_per_day==7)


# subset to proteins with > 4 peptides:
dt_sub[,n_pep:=length(unique(FullPeptideName)),by=c("ProteinName")]
dt_sub <- dt_sub[n_pep>=4]

# median normalization
dt_sub[,log2_int:=log2(pep_int)]
dt_sub[,median_perRun:=median(log2_int), by="run"]
dt_sub[,median_median:=mean(median_perRun)]
dt_sub[,diff_median:=median_median-median_perRun]
dt_sub[,norm_log2_int := log2_int+diff_median]
dt_sub[,norm_int := 2^norm_log2_int]

# introduce variation
set.seed(1)
dt_sub[,diff_fac_3 := runif(1, min = 1, max = 6), by="ProteinName"]
set.seed(2)
dt_sub[,diff_fac_5 := runif(1, min = 1, max = 6), by="ProteinName"]
#dt_sub[,diff_norm_int:=ifelse(day %in% c("day3","day5"), diff_fac*norm_int, norm_int)]
dt_sub[,diff_norm_int:=ifelse(day == "day3", diff_fac_3*norm_int, norm_int)]
dt_sub[,diff_norm_int:=ifelse(day == "day5", diff_fac_5*norm_int, diff_norm_int)]

# randomly select 1000 proteins to perturb
set.seed(22)
proteins_to_perturb = sample(unique(dt_sub$ProteinName),1000)
dt_sub[,perturbed_protein := ifelse(ProteinName %in% proteins_to_perturb, TRUE, FALSE), by="ProteinName"]

# determine reduction factor for each protein
# sample from uniform distribution
set.seed(44)
dt_sub[,red_fac:=ifelse(perturbed_protein, runif(1, min = 0.01, max = 0.90), 1), by="ProteinName"]

interlab_benchmark_data <- copy(dt_sub)

saveRDS(interlab_benchmark_data, "interlab_benchmark_data.rds")


##################################
### Introduce pep. perturbation ##
##################################

rm(list=ls())
gc()
interlab_benchmark_data <- readRDS("interlab_benchmark_data.rds")

generatePerturbedProfiles <- function(input_data, nf_peptides_to_perturb="random"){
  dt_input <- copy(input_data)
  
  set.seed(66)
  
  if (nf_peptides_to_perturb == "random") {
    dt_input[,max_perturbed_peptides := ceiling(n_pep*0.5), by="ProteinName"]
    dt_input[,n_perturbed_peptides := max(sample(seq(2,max(2,max_perturbed_peptides),1), 1),2), by="ProteinName"]
  } else if (nf_peptides_to_perturb >= 1) {
    dt_input[,n_perturbed_peptides := nf_peptides_to_perturb, by="ProteinName"]
  } else if (nf_peptides_to_perturb < 1) {
    dt_input[,n_perturbed_peptides := ceiling(n_pep*nf_peptides_to_perturb), by="ProteinName"]
    dt_input[,n_perturbed_peptides := ifelse(n_perturbed_peptides<2,2,n_perturbed_peptides), by="ProteinName"]
  }
  
  dt_input[,perturbed_peptides := paste(unique(FullPeptideName)[c(1:n_perturbed_peptides)],collapse=";"), by = c("ProteinName")]
  
  # reduce day 5 by reduction factor for peptides in perturbed peptides
  dt_input[,perturbed_peptide:=((FullPeptideName %in% unlist(strsplit(perturbed_peptides,";"))) & (perturbed_protein==TRUE)), by=c("FullPeptideName","ProteinName")]
  
  dt_input[,mod_pep_int:=ifelse(((day=="day5") & (perturbed_peptide) & (perturbed_protein)), diff_norm_int*red_fac, diff_norm_int)]
  
  setnames(dt_input, c("FullPeptideName","ProteinName","run","mod_pep_int"), c("peptide_id","protein_id","filename","intensity"))
  
  # format for CCprofiler
  input_data <- subset(dt_input, select=c("peptide_id","protein_id","filename","intensity"))
  
  input_data_annotation <- unique(subset(dt_input,select=c("filename","day")))
  setorderv(input_data_annotation, c("day","filename"))
  input_data_annotation[,fraction_number:=.I]
  input_data_annotation <- subset(input_data_annotation,select=c("filename","fraction_number"))
  
  # Plot
  pdf(paste0("plot_example_perturbed_proteins_",nf_peptides_to_perturb,".pdf"), width = 5, height = 3)
  col <- c("TRUE"="magenta4", "FALSE"="#6699cc")
  perturbed_proteins <- unique(dt_input[perturbed_protein==TRUE]$protein_id)
  for (i in c(1:100)){
    x <- dt_input[protein_id==perturbed_proteins[i]]
    
    x$filename <- factor(x$filename, levels=input_data_annotation$filename)
    
    filename_label <- c(rep("day1", 7),rep("day3", 7),rep("day5", 7))
    rep_label <- rep(seq(1,7,1),3)
    filename_label <- paste(filename_label, rep_label, sep="_")
    
    p <- ggplot(x, aes(x=filename, y=intensity, color=perturbed_peptide, group=peptide_id)) + 
      geom_point(alpha=0.9) +
      geom_line(alpha=0.7) +
      theme_classic() +
      theme(legend.position = "bottom") +
      scale_color_manual(values = col, name="perturbed peptide") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_x_discrete(labels=filename_label) +
      xlab("sample") +
      ggtitle(paste0(perturbed_proteins[i]," \nperturbation factor = ",round(unique(x$red_fac),digits = 3)))
    print(p)
  }
  dev.off()
  
  pdf(paste0("plot_example_non_perturbed_proteins_",nf_peptides_to_perturb,".pdf"), width = 5, height = 3)
  proteins_not_perturbed = unique(dt_input[perturbed_protein == FALSE]$protein_id)
  col <- c("TRUE"="magenta4", "FALSE"="#6699cc")
  for (i in c(1:50)){
    x <- dt_input[protein_id==proteins_not_perturbed[i]]
    
    x$filename <- factor(x$filename, levels=input_data_annotation$filename)
    
    filename_label <- c(rep("day1", 7),rep("day3", 7),rep("day5", 7))
    rep_label <- rep(seq(1,7,1),3)
    filename_label <- paste(filename_label, rep_label, sep="_")
    
    p <- ggplot(x, aes(x=filename, y=intensity, color=perturbed_peptide, group=peptide_id)) + 
      geom_point(alpha=0.9) +
      geom_line(alpha=0.7) +
      theme_classic() +
      theme(legend.position = "bottom") +
      scale_color_manual(values = col, name="perturbed peptide") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_x_discrete(labels=filename_label) +
      xlab("sample") +
      ggtitle(paste0(perturbed_proteins[i]," \nperturbation factor = ",round(unique(x$red_fac),digits = 3)))
    print(p)
  }
  dev.off()
  
  
  return(dt_input)
}

random_perm_data <- generatePerturbedProfiles(input_data=interlab_benchmark_data, nf_peptides_to_perturb="random")
perm_1pep_data <- generatePerturbedProfiles(input_data=interlab_benchmark_data, nf_peptides_to_perturb=1)
perm_2pep_data <- generatePerturbedProfiles(input_data=interlab_benchmark_data, nf_peptides_to_perturb=2)
perm_025pep_data <- generatePerturbedProfiles(input_data=interlab_benchmark_data, nf_peptides_to_perturb=0.25)
perm_050pep_data <- generatePerturbedProfiles(input_data=interlab_benchmark_data, nf_peptides_to_perturb=0.5)


# Combine all data
all_data <- list('random' = random_perm_data, 
                 '1pep' = perm_1pep_data, 
                 '2pep' = perm_2pep_data, 
                 '025pep' = perm_025pep_data, 
                 '050pep' = perm_050pep_data)

#####################################
### CCprofiler analysis #############
#####################################

library(ggpubr)
library(gridExtra)
library(grid)

performCCprofilerAnalysis <- function(data, name, score_cutoff = 0.1, adj_pval_cutoff =  0.1){
  
  dt_input <- copy(data)
  
  input_data_annotation <- unique(subset(dt_input,select=c("filename","day")))
  setorderv(input_data_annotation, c("day","filename"))
  input_data_annotation[,fraction_number:=.I]
  input_data_annotation <- subset(input_data_annotation,select=c("filename","fraction_number"))
  
  traces <- importPCPdata(input_data = dt_input,
                          fraction_annotation = input_data_annotation)
  
  trace_annotation <- unique(subset(dt_input, select=c("peptide_id",
                                                       "n_pep","n_perturbed_peptides",
                                                       "perturbed_protein","perturbed_peptide","red_fac")))
  
  traces <- annotateTraces(traces,
                           trace_annotation,
                           traces_id_column = "id",
                           trace_annotation_id_column = "peptide_id")
  
  ## Remove traces with 0 standard deviation (can't be clustered)
  zerovar <- apply(getIntensityMatrix(traces), 1, var) > 0
  traces_zerovar <- subset(traces,
                           trace_subset_ids = names(zerovar[zerovar]))
  
  #' ## Remove single peptide genes
  traces_multiPep <- filterSinglePeptideHits(traces_zerovar)
  
  #' ## Estimate proteoform scores by correlation clustering
  traces_corr <- calculateGeneCorrMatrices(traces_multiPep)
  
  traces_clustered <- clusterPeptides(
    traces_corr,
    method = "average", plot = F, PDF=F,
    name=paste0("ProteoformClusters_interlab_",name))
  
  traces_clusteredInN <- cutClustersInNreal(traces_clustered, clusterN = 2,
                                            min_peptides_per_cluster = 2)
  
  traces_scored <- calculateProteoformScore(traces_clusteredInN)
  
  scores <- unique(traces_scored$trace_annotation[,.(proteoform_score,proteoform_score_z,proteoform_score_pval_adj,perturbed_protein), by=.(protein_id)])
  scores <- scores[!is.na(scores$proteoform_score)]
  
  
  
  library("ggExtra")
  scatter <- ggplot(scores, 
                    aes(x=proteoform_score,y=-log10(proteoform_score_pval_adj), 
                        color=perturbed_protein)) + 
    scale_color_manual(values=c("#ff7700", "green4"), name= "protein with proteoform") +
    geom_point(alpha=0.5) +
    geom_vline(xintercept = score_cutoff, colour='darkgrey', lty='dashed') +
    geom_hline(yintercept = -log10(adj_pval_cutoff), colour='darkgrey', lty='dashed') +
    theme_classic() +
    ylab("-log10 ( adj. pvalue )") +
    xlab("proteoform score") +
    theme(legend.position="bottom") 
  
  scatter_m <- ggMarginal(scatter, type = "histogram", groupFill = TRUE, groupColour = TRUE, 
                          bins=75, position = "identity", size = 4)
  
  pdf(paste0("pseudo_volcano_hist_",name,".pdf"), height=5, width=5)
  print(scatter_m)
  dev.off()
  
  traces_proteoforms <- annotateTracesWithProteoforms(
    traces_scored, score_cutoff = score_cutoff, adj_pval_cutoff =  adj_pval_cutoff)
  
  return(traces_proteoforms)
}

performFDRbenchmark <- function(traces, name="traces", score_cutoff = 0.1, adj_pval_cutoff =  0.1){
  score_thresholds = seq(0,0.35,0.05)
  pval_thresholds = sort(unique(c(c(1 %o% 10^(-12:1)),seq(0,1,0.05))))
  
  n_row = nrow(expand.grid(score_thresholds,pval_thresholds))
  
  res <- data.table(score_thresholds=rep(0,n_row),
                    pval_thresholds=rep(0,n_row),
                    P = length(unique(traces$trace_annotation[perturbed_protein==TRUE]$protein_id)),
                    N = length(unique(traces$trace_annotation[perturbed_protein==FALSE]$protein_id)),
                    n_prot_proteoforms = 0,
                    TP = 0,
                    FP = 0,
                    N_prot_noMistakes = 0,
                    N_prot_withMistakes = 0,
                    TP_noMistakes = 0,
                    TP_withMistakes = 0
  )
  
  idx = 0
  for (i in seq(1,length(score_thresholds),1)) {
    for (j in seq(1,length(pval_thresholds),1)) {
      idx = idx+1
      res$score_thresholds[idx] <- score_thresholds[i]
      res$pval_thresholds[idx] <- pval_thresholds[j]
      
      traces_proteoforms <- annotateTracesWithProteoforms(
        traces, score_cutoff = score_thresholds[i],  adj_pval_cutoff = pval_thresholds[j])
      
      res$n_prot_proteoforms[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1)]$protein_id))
      res$TP[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1) & (perturbed_protein)]$protein_id))
      res$FP[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1) & (perturbed_protein==FALSE)]$protein_id))
      
      traces_proteoforms$trace_annotation[, n_proteoforms_per_perturbed_group := length(unique(proteoform_id)), by=c("protein_id","perturbed_peptide")]
      traces_proteoforms$trace_annotation[, n_perturbed_groups_per_proteoform := length(unique(perturbed_peptide)), by=c("protein_id","proteoform_id")]
      
      traces_proteoforms$trace_annotation[, correct_peptide := ((n_proteoforms_per_perturbed_group==1) & (n_perturbed_groups_per_proteoform==1))]
      
      traces_proteoforms$trace_annotation[, protein_without_mistakes := all(correct_peptide), by="protein_id"]
      
      res$N_prot_noMistakes[idx] <- length(unique(traces_proteoforms$trace_annotation[protein_without_mistakes==TRUE]$protein_id))
      res$N_prot_withMistakes[idx] <- length(unique(traces_proteoforms$trace_annotation[protein_without_mistakes==FALSE]$protein_id))
      
      res$TP_noMistakes[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1) & (perturbed_protein) & (protein_without_mistakes==TRUE)]$protein_id))
      res$TP_withMistakes[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1) & (perturbed_protein) & (protein_without_mistakes==FALSE)]$protein_id))
    }
  }
  
  res[,FDR := FP/(FP+TP)]
  res[,FN := P-TP]
  res[,F1 := TP/(TP+(0.5*(FP+FN)))]
  
  res[,TPR := TP/P]
  res[,FPR := FP/N]
  
  res[,percent_TP_perfect := TP_noMistakes/TP]
  res[,percent_perfect := N_prot_noMistakes/(P+N)]
  
  res$score_thresholds <- round(res$score_thresholds, digits = 2)
  
  pdf(paste0("FDR_benchmark_",name,".pdf"), width=6.5, height=5)
  okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "darkgrey")
  res$score_thresholds_fac = as.factor(res$score_thresholds)
  p <- ggplot(res, aes(x=pval_thresholds,y=FDR, group=score_thresholds, colour=score_thresholds_fac)) + 
    geom_point() + geom_line() + 
    xlim(0,0.5) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = okabe) +
    xlab("adj. p-value threshold") + 
    labs(colour = "proteoform score \n threshold") +
    theme_classic()
  print(p)
  p <- ggplot(res, aes(x=pval_thresholds,y=FDR, group=score_thresholds, colour=score_thresholds_fac)) + 
    geom_point() + geom_line() + 
    xlim(0,0.17) +
    ylim(0,0.17) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = okabe) +
    xlab("adj. p-value threshold") + 
    labs(colour = "proteoform score \n threshold") +
    theme_classic()
  print(p)
  dev.off()
  
  res[,selected_pval := ifelse(pval_thresholds==adj_pval_cutoff,"TRUE","FALSE")]
  pval_col <- c("TRUE"="black","FALSE"="#00ff0000")
  
  pdf(paste0("ROC_percentPerfect_",name,"_pval_",adj_pval_cutoff,".pdf"), width=6.5, height=5)
  p <- ggplot(res[score_thresholds==score_cutoff], aes(x=FPR,y=TPR)) + 
    geom_line(colour="grey",alpha=0.5) +
    geom_abline(intercept = 0, slope = 1, color='grey') +
    scale_fill_gradientn(limits = c(0,1),
                         colours=c("navyblue", "darkmagenta", "darkorange1"),
                         breaks=seq(0,1,0.25)) +
    geom_point(aes(fill=percent_perfect, color=selected_pval,stroke=selected_pval), shape = 21, size = 3, stroke=1) +
    scale_color_manual(values = pval_col) +
    geom_abline(intercept = 0, slope = 1) +
    xlim(0,1) +
    ylim(0,1) +
    theme_classic()
  print(p)
  dev.off()
  
  return(res)
}

random_perm_scored <- performCCprofilerAnalysis(data=random_perm_data, name="random", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
# plotProteoformPvalHist(random_perm_scored, name="random_perm_scored_pval_hist", PDF=T)
# plotProteoformVolcano(random_perm_scored, name="random_perm_scored_volcano", PDF=T)

perm_1pep_scored <- performCCprofilerAnalysis(data=perm_1pep_data, name="perm_1pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
perm_2pep_scored <- performCCprofilerAnalysis(data=perm_2pep_data, name="perm_2pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
perm_025pep_scored <- performCCprofilerAnalysis(data=perm_025pep_data, name="perm_025pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
perm_050pep_scored <- performCCprofilerAnalysis(data=perm_050pep_data, name="perm_050pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)

random_perm_scored_annot <- random_perm_scored$trace_annotation
random_perm_scored_annot$annot <- 'random'
perm_1pep_scored_annot <- perm_1pep_scored$trace_annotation
perm_1pep_scored_annot$annot <- '1pep'
perm_2pep_scored_annot <- perm_2pep_scored$trace_annotation
perm_2pep_scored_annot$annot <- '2pep'
perm_025pep_scored_annot <- perm_025pep_scored$trace_annotation
perm_025pep_scored_annot$annot <- '025pep'
perm_050pep_scored_annot <- perm_050pep_scored$trace_annotation
perm_050pep_scored_annot$annot <- '050pep'


scores_copf <- rbind(random_perm_scored_annot, 
                     perm_1pep_scored_annot,
                     perm_2pep_scored_annot,
                     perm_025pep_scored_annot,
                     perm_050pep_scored_annot)

saveRDS(scores_copf, "scores_copf.rds")

fdr_bench_random_perm <- performFDRbenchmark(traces = random_perm_scored, name="random_perm", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
fdr_bench_perm_1pep <- performFDRbenchmark(traces = perm_1pep_scored, name="perm_1pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
fdr_bench_perm_2pep <- performFDRbenchmark(traces = perm_2pep_scored, name="perm_2pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
fdr_bench_perm_025pep <- performFDRbenchmark(traces = perm_025pep_scored, name="perm_025pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
fdr_bench_perm_050pep <- performFDRbenchmark(traces = perm_050pep_scored, name="perm_050pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)


res_copf <- rbind(fdr_bench_random_perm[,nPerm := "random"],
                  fdr_bench_perm_1pep[,nPerm := "1"],
                  fdr_bench_perm_2pep[,nPerm := "2"],
                  fdr_bench_perm_025pep[,nPerm := "25%"],
                  fdr_bench_perm_050pep[,nPerm := "50%"])

saveRDS(res_copf, "res_copf.rds")

#####################################
### PeCorA analysis #################
#####################################
library(devtools)
# install_github("jessegmeyerlab/PeCorA")
library(PeCorA)

performPecoraAnalysis <- function(data){
  input_data_pecora <- copy(data)
  
  setnames(input_data_pecora,c("peptide_id","protein_id","intensity","day"), c("Peptide","Protein","Normalized.Area","Condition"))
  input_data_pecora[,"Peptide.Modified.Sequence":=Peptide]
  
  ann <- unique(subset(input_data_pecora, select=c("filename","Condition")))
  setorderv(ann,c("Condition","filename"))
  ann[,BioReplicate:=c(1:7), by="Condition"]
  
  input_data_pecora <- merge(input_data_pecora,ann,by=c("filename","Condition"))
  
  input_data_pecora <- subset(input_data_pecora, select=c("Protein","Peptide","Peptide.Modified.Sequence","Condition","BioReplicate","Normalized.Area"))
  
  pecora_df <- setDF(input_data_pecora)
  
  scaled_peptides <- PeCorA_preprocessing(pecora_df,
                                          area_column_name=6,
                                          threshold_to_filter=min(pecora_df$Normalized.Area),
                                          control_name="day1")
  
  disagree_peptides <- PeCorA(scaled_peptides)
  
  pecora_res <- data.table(disagree_peptides)
  pecora_res[,"peptide_id":=gsub("_all","",peptide)]
  
  trace_annotation <- unique(subset(data, select=c("peptide_id","n_pep","n_perturbed_peptides",
                                                   "perturbed_protein","perturbed_peptide","red_fac")))
  
  pecora_res <- merge(pecora_res,trace_annotation,by="peptide_id")
  
  pecora_res[, adj_adj_pval := p.adjust(adj_pval, method = "BH")]
  
  return(pecora_res)
}


pecoraFDRbenchmark <- function(data, name="pecora", pval_col = "adj_pval", pval_cutoff =  0.05){
  
  pecora_res <- copy(data)
  
  all_pos_prot <- length(unique(pecora_res[(perturbed_protein==TRUE)]$protein))
  all_neg_prot <- length(unique(pecora_res[(perturbed_protein==FALSE)]$protein))
  
  all_pos_pep <- length(unique(pecora_res[perturbed_peptide==TRUE]$peptide))
  all_neg_pep <- length(unique(pecora_res[perturbed_peptide==FALSE]$peptide))
  
  score_thresholds = sort(unique(c(c(1 %o% 10^(-12:1)),seq(0,1,0.05))))
  
  pecora_res_df <- data.table(score_thresholds=score_thresholds,
                              P = all_pos_prot,
                              N = all_neg_prot,
                              n_prot_proteoforms = 0,
                              TP = 0,
                              FP = 0,
                              N_prot_noMistakes = 0,
                              N_prot_withMistakes = 0,
                              TP_noMistakes = 0,
                              TP_withMistakes = 0
  )
  
  for (i in seq(1,length(score_thresholds),1)) {
    pecora_res_sig <- subset(pecora_res, get(pval_col)<=score_thresholds[i])
    
    pecora_res_df$n_prot_proteoforms[i] <- length(unique(pecora_res_sig$peptide))
    
    TP_pep = length(unique(pecora_res_sig[perturbed_peptide==TRUE]$peptide))
    FP_pep = length(unique(pecora_res_sig[perturbed_peptide==FALSE]$peptide))
    TP_prot = length(unique(pecora_res_sig[perturbed_protein==TRUE]$protein))
    FP_prot = length(unique(pecora_res_sig[perturbed_protein==FALSE]$protein))
    
    pecora_res_df$TP[i] <- TP_prot
    pecora_res_df$FP[i] <- FP_prot
    
    pecora_res[, TP_perturbed := (perturbed_peptide==TRUE) & (get(pval_col)<=score_thresholds[i]), by= "peptide"]
    pecora_res[, FP_perturbed := (perturbed_peptide==FALSE) & (get(pval_col)<=score_thresholds[i]), by= "peptide"]
    
    pecora_res[, correct_assignment := (((perturbed_peptide==TRUE) & (get(pval_col)<=score_thresholds[i])) | ((perturbed_peptide==FALSE) & (get(pval_col)>score_thresholds[i]))) , by= "peptide"]
    
    pecora_res[, perfect_prot := all(correct_assignment), by="protein"] 
    pecora_res[, min_adj_pval := min(get(pval_col)), by="protein"]
    
    pecors_score <- unique(subset(pecora_res, select=c("protein","perturbed_protein","perfect_prot","min_adj_pval")))
    pecora_res_df$N_prot_noMistakes[i] <- sum(pecors_score$perfect_prot)
    pecora_res_df$N_prot_withMistakes[i] <- nrow(pecors_score[perfect_prot==FALSE])
    
    pecora_res_df$TP_noMistakes[i] <- sum(pecors_score[(perturbed_protein==TRUE) & (min_adj_pval<=score_thresholds[i])]$perfect_prot)
    pecora_res_df$TP_withMistakes[i] <- nrow(pecors_score[(perturbed_protein==TRUE) & (min_adj_pval<=score_thresholds[i]) & (perfect_prot==FALSE)])
    
  }
  
  pecora_res_df[,TPR:=TP/P]
  pecora_res_df[,FPR:=FP/N]
  pecora_res_df[,percent_TP_perfect := TP_noMistakes/TP]
  pecora_res_df[,percent_perfect := N_prot_noMistakes/(P+N)]
  
  pecora_res_df[,FDR := FP/(FP+TP)]
  pecora_res_df[,FN := P-TP]
  pecora_res_df[,F1 := TP/(TP+(0.5*(FP+FN)))]
  
  pecora_res_df[,pval_thresholds:=score_thresholds]
  pecora_res_df[,selected_pval := ifelse(pval_thresholds==pval_cutoff,TRUE,FALSE)]
  pval_col <- c("TRUE"="black","FALSE"="#00ff0000")
  
  pdf(paste0("pecora_ROC_percentPerfect_",name,"_pval_",pval_cutoff,".pdf"), width=6.5, height=5)
  p <- ggplot(pecora_res_df, aes(x=FPR,y=TPR)) + 
    geom_line(colour="grey",alpha=0.5) +
    geom_abline(intercept = 0, slope = 1, color='grey') +
    scale_fill_gradientn(limits = c(0,1),
                         colours=c("navyblue", "darkmagenta", "darkorange1"),
                         breaks=seq(0,1,0.25)) +
    geom_point(aes(fill=percent_perfect, color=selected_pval,stroke=selected_pval), shape = 21, size = 3, stroke=1) +
    scale_color_manual(values = pval_col) +
    geom_abline(intercept = 0, slope = 1) +
    xlim(0,1) +
    ylim(0,1) +
    theme_classic()
  print(p)
  dev.off()
  
  return(pecora_res_df)
}


random_perm_pecora <- performPecoraAnalysis(data=random_perm_data)
fdr_bench_random_perm_pecora <- pecoraFDRbenchmark(data=random_perm_pecora, name="random_perm", pval_col = "adj_pval", pval_cutoff =  0.1)
fdr_bench_random_perm_pecora_adj <- pecoraFDRbenchmark(data=random_perm_pecora, name="random_perm_adj", pval_col = "adj_adj_pval", pval_cutoff =  0.1)

perm_1pep_pecora <- performPecoraAnalysis(data=perm_1pep_data)
fdr_bench_perm_1pep_pecora <- pecoraFDRbenchmark(data=perm_1pep_pecora, name="perm_1pep", pval_col = "adj_pval", pval_cutoff =  0.1)
fdr_bench_perm_1pep_pecora_adj <- pecoraFDRbenchmark(data=perm_1pep_pecora, name="perm_1pep_adj", pval_col = "adj_adj_pval", pval_cutoff =  0.1)

perm_2pep_pecora <- performPecoraAnalysis(data=perm_2pep_data)
fdr_bench_perm_2pep_pecora <- pecoraFDRbenchmark(data=perm_2pep_pecora, name="perm_2pep", pval_col = "adj_pval", pval_cutoff =  0.1)
fdr_bench_perm_2pep_pecora_adj <- pecoraFDRbenchmark(data=perm_2pep_pecora, name="perm_2pep_adj", pval_col = "adj_adj_pval", pval_cutoff =  0.1)

perm_025pep_pecora <- performPecoraAnalysis(data=perm_025pep_data)
fdr_bench_perm_025pep_pecora <- pecoraFDRbenchmark(data=perm_025pep_pecora, name="perm_025pep", pval_col = "adj_pval", pval_cutoff =  0.1)
fdr_bench_perm_025pep_pecora_adj <- pecoraFDRbenchmark(data=perm_025pep_pecora, name="perm_025pep_adj", pval_col = "adj_adj_pval", pval_cutoff =  0.1)

perm_050pep_pecora <- performPecoraAnalysis(data=perm_050pep_data)
fdr_bench_perm_050pep_pecora <- pecoraFDRbenchmark(data=perm_050pep_pecora, name="perm_050pep", pval_col = "adj_pval", pval_cutoff =  0.1)
fdr_bench_perm_050pep_pecora_adj <- pecoraFDRbenchmark(data=perm_050pep_pecora, name="perm_050pep_adj", pval_col = "adj_adj_pval", pval_cutoff =  0.1)


random_perm_pecora$annot <- 'random'
perm_1pep_pecora$annot <- '1pep'
perm_2pep_pecora$annot <- '2pep'
perm_025pep_pecora$annot <- '025pep'
perm_050pep_pecora$annot <- '050pep'


scores_pecora <- rbind(random_perm_pecora, 
                       perm_1pep_pecora,
                       perm_2pep_pecora,
                       perm_025pep_pecora,
                       perm_050pep_pecora)

saveRDS(scores_pecora, "scores_pecora.rds")


res_pecora <- rbind(fdr_bench_random_perm_pecora[,c("nPerm","pval_col") := list("random","adj_pval")],
                    fdr_bench_perm_1pep_pecora[,c("nPerm","pval_col") := list("1","adj_pval")],
                    fdr_bench_perm_2pep_pecora[,c("nPerm","pval_col") := list("2","adj_pval")],
                    fdr_bench_perm_025pep_pecora[,c("nPerm","pval_col") := list("25%","adj_pval")],
                    fdr_bench_perm_050pep_pecora[,c("nPerm","pval_col") := list("50%","adj_pval")],
                    fdr_bench_random_perm_pecora_adj[,c("nPerm","pval_col") := list("random","adj_adj_pval")],
                    fdr_bench_perm_1pep_pecora_adj[,c("nPerm","pval_col") := list("1","adj_adj_pval")],
                    fdr_bench_perm_2pep_pecora_adj[,c("nPerm","pval_col") := list("2","adj_adj_pval")],
                    fdr_bench_perm_025pep_pecora_adj[,c("nPerm","pval_col") := list("25%","adj_adj_pval")],
                    fdr_bench_perm_050pep_pecora_adj[,c("nPerm","pval_col") := list("50%","adj_adj_pval")])

saveRDS(res_pecora, "res_pecora.rds")

#####################################
### DEpMS analysis  #####
#####################################
# Functions
source('~/Downloads/DEpMS_project/depms_functions.R')

# Directory
path_code <- '~/Downloads/DEpMS_project/depms_benchmark/'
path_data <- '~/Downloads/DEpMS_project/depms_benchmark/'
path_figures <- '~/Downloads/DEpMS_project/depms_benchmark/'


# Parameters
ioi_id <- 'protein_id'
cor_method <-  'pearson'


performDepmsAnalysis <- function(dt_input,  fitting) {
  
  input_data_annotation <- unique(subset(dt_input,select=c("filename","day")))
  setorderv(input_data_annotation, c("day","filename"))
  input_data_annotation[,fraction_number:=.I]
  input_data_annotation <- subset(input_data_annotation,select=c("filename","fraction_number"))
  
  
  traces <- importPCPdata(input_data = dt_input,
                          fraction_annotation = input_data_annotation)
  
  trace_annotation <- unique(subset(dt_input, select=c("peptide_id",
                                                       "n_pep","n_perturbed_peptides",
                                                       "perturbed_protein","perturbed_peptide","red_fac")))
  
  traces <- annotateTraces(traces,
                           trace_annotation,
                           traces_id_column = "id",
                           trace_annotation_id_column = "peptide_id")
  
  ## Remove traces with 0 standard deviation (can't be clustered)
  zerovar <- apply(getIntensityMatrix(traces), 1, var) > 0
  traces_zerovar <- subset(traces,
                           trace_subset_ids = names(zerovar[zerovar]))
  
  #' ## Remove single peptide genes
  traces_multiPep <- filterSinglePeptideHits(traces_zerovar)
  
  
  # Peptide qunatification
  peptide_quant <- as.data.frame(getIntensityMatrix(traces_multiPep))
  
  # log2 transformation
  peptide_quant <- log2(peptide_quant)
  
  # Relative quantification
  peptide_quant <- sweep(peptide_quant, 1, apply(peptide_quant, 1, median), '-')
  
  
  # Protein - peptide annotation
  protein_peptide_annotation <- as.data.frame(traces_multiPep$trace_annotation)
  colnames(protein_peptide_annotation)[1] <- 'Peptide_sequence'
  
  
  # Derive protein quantification from peptide table
  protein_quant <- peptide_quant
  protein_quant$ioi <- protein_peptide_annotation[[ioi_id]]
  protein_quant <- protein_quant[!grepl(';', protein_quant$ioi), ]
  protein_quant <- protein_quant[!is.na(protein_quant$ioi), ]
  protein_quant <- protein_quant[protein_quant$ioi != '', ]
  protein_quant <- as.data.frame(protein_quant %>% group_by(ioi) %>%
                                   summarise_all(median, na.rm = TRUE))
  rownames(protein_quant) <- protein_quant$ioi
  protein_quant$ioi <- NULL
  

  # # check proteins with more than 1 peptide
  freq_ioi <- table(protein_peptide_annotation$protein_id)
  tocheck_ioi <- names(freq_ioi)[freq_ioi > 1]
  ####################################################################################################################################
  # Median, IQR modeling 
  global_sim <- peptide.protein.similarity(peptide_quant = peptide_quant, 
                                           protein_quant = protein_quant, 
                                           protein_peptide_annotation = protein_peptide_annotation, 
                                           tocheck_ioi = tocheck_ioi)
  
  ####################################################################################################################################
  # Profile construction
  collapsed_profiles <- profile.generation(peptide_quant = peptide_quant,
                                           protein_quant = protein_quant, 
                                           protein_peptide_annotation = protein_peptide_annotation, 
                                           global_sim = global_sim,
                                           tocheck_ioi = tocheck_ioi)
  

  collapsed_quant_case <- collapsed_profiles$collapsed_quant_case
  support_case <- collapsed_profiles$support_case
  collapsed_quant_background <- collapsed_profiles$collapsed_quant_background
  support_background <- collapsed_profiles$support_background
  ####################################################################################################################################
  ## PC analysis
  total_n_pcs <- 2
  dat <- list('Protein_level' = protein_quant,
              'Collapsed_profile_level' = collapsed_quant_case)

  day_factor <- factor(rep(paste0('day', c(1,3,5)), each = 7),
         levels = unique(rep(paste0('day', c(1,3,5)), each = 7)))
  days_cols <- c('red', 'blue', 'green')
  names(days_cols) <- unique(day_factor)
  
  lapply(names(dat), function(i) {

    input_dat <- dat[[i]]

    if(i == 'Collapsed_profile_level') {
      tokeepforpca <- unique(support_case$id[support_case$weight == 1])
      input_dat <- input_dat[tokeepforpca, ]
    }
    pca_dt <- prcomp(t(input_dat[complete.cases(input_dat), ]))

    p_var <- plot.expl.variance(pca_dt, i)

    pdf(paste0(path_figures, 'PCA_variance_',i,'.pdf'), width = 6 ,height = 5)
    plot(p_var)
    dev.off()

    lapply(c('Days'), function(j) {
      
      cols_var <- days_cols
      
      p_proj <- plot.pcs(pca_res = pca_dt,
                         n_pcs =  total_n_pcs,
                         var_to_plot = day_factor,
                         cols_var = cols_var,
                         plot_title = paste0(i, ' - ', j))

      pdf(paste0(path_figures, 'PCA_', i, '_', j, '.pdf'), width = 6 ,height = 5)
      plot(p_proj)
      dev.off()
    })

    # p_assoc <- plot.pc.associations(pca_res = pca_dt,
    #                                 tocheckmetadata = as.data.frame(day_factor),
    #                                 n_pcs = total_n_pcs,
    #                                 plot_title = i)
    # 
    # pdf(paste0(path_figures, 'PCA_', i, '_associations.pdf'), width = 6 ,height = 8)
    # plot(p_assoc)
    # dev.off()

  })
  
  ####################################################################################################################################
  # Fitting
  if(fitting == 'group') {
    
    clusters <- factor(rep(paste0('day', c(1,3,5)), each = 7),
                     levels = unique(rep(paste0('day', c(1,3,5)), each = 7)))
    design_mat <- model.matrix(~0 + clusters)
    colnames(design_mat) <- levels(clusters)
    
    # All comparisons
    tocompare <- combn(colnames(design_mat), 2)
    tocompare <- apply(tocompare, 2, paste, collapse = '-')
    
    depms.res.all <-  do.call(rbind, lapply(colnames(design_mat), function(design_var) {
      
      print(design_var)
      
      design_idx <- grep(design_var, tocompare)
      
      depms.res <- do.call(rbind, lapply(design_idx, function(n_comp) {
        
        print(tocompare[n_comp])
        design_vars <- unlist(strsplit(tocompare[n_comp], '-'))
        samples_id <- which(clusters %in% design_vars)
        
        # Background
        tofit_null <- collapsed_quant_background
        
        tokeep_var_null <- do.call(cbind, lapply(design_vars, function(k) {
          dat_sub <- tofit_null[, which(design_mat[, k] == 1)]
          idx <- apply(dat_sub, 1, function(i) var(i, na.rm = TRUE) !=0)
          }))
        tokeep_var_null <- apply(tokeep_var_null, 1, all)
        tofit_null <- tofit_null[tokeep_var_null, ]
    

        
        print(paste('background dataset:', nrow(tofit_null), ncol(tofit_null)))
        # boxplot(tofit_null)
        
        
        fit_null <- lmFit(tofit_null, design_mat)
        contrasts <- makeContrasts(contrasts = tocompare[n_comp], levels =  colnames(design_mat))
        fit2_null <- contrasts.fit(fit_null, contrasts)
        fit3_null <- eBayes(fit2_null)
        
        fit3_null$count <- quant.for.variance(fit = fit3_null, 
                                              support_dt = support_background, 
                                              protein_peptide_annotation = protein_peptide_annotation,
                                              psm_vec = NULL)
        fit4_null <- spectraCounteBayesGAM(fit = fit3_null, coef_col = 1)
        # VarianceBoxplotGAM(fit4_null,n=30,main="",xlab="Count")
        # VarianceScatterplotGAM(fit4_null,main="")
        # ResidualplotGAM(fit4_null)
        
        
        DEpMS.results_null <- outputResultGAM(fit = fit4_null, coef_col = 1)
      
        ## Remove NAs
        DEpMS.results_null <- DEpMS.results_null[complete.cases(DEpMS.results_null), ]
        
        
        
        # Case 
        tofit <- collapsed_quant_case
        
        print(paste('case dataset:', nrow(tofit), ncol(tofit)))
        
        fit <- lmFit(tofit, design_mat)
        fit2 <- contrasts.fit(fit, contrasts)
        fit3 <- eBayes(fit2)
        
        
        ## Quantification for variance model
        fit3$count <- quant.for.variance(fit = fit3, 
                                         support_dt = support_case, 
                                         protein_peptide_annotation = protein_peptide_annotation,
                                         psm_vec = NULL)
        fit4 <- spectraCounteBayesGAM(fit3, coef_col = 1)
        
        # VarianceBoxplotGAM(fit4,n=30,main="",xlab="Count")
        # VarianceScatterplotGAM(fit4,main="")
        # ResidualplotGAM(fit4)
        
        DEpMS.results <- outputResultGAM(fit = fit4, coef_col = 1)
        DEpMS.results <- DEpMS.results[complete.cases(DEpMS.results), ]
        
        # Q-value calculation
        res <- q.value.calculation(support_case = support_case,
                                   # support_background = support_background,
                                   DEpMS.results_null = DEpMS.results_null, 
                                   DEpMS.results = DEpMS.results, 
                                   protein_peptide_annotation = protein_peptide_annotation, 
                                   design_var = tocompare[n_comp],
                                   plot_density = TRUE)
        
        return(res)
      }))
      
      depms.res$day <- design_var
      return(depms.res)
    }))
    
    write.table(depms.res.all, paste0(path_data, 'depms_results_groups.txt'), sep = '\t')
    # depms.res <- read.delim(paste0(path_data, 'depms_results_groups.txt'), sep = '\t')
  
    logFC <- depms.res.all$logFC
    
  
    ## Step-wise summarization
    # Summarize within groups
    top_n <- 2
    
    depms.res.all_within  <- do.call(rbind, lapply(colnames(design_mat), function(i) {
      print(i)
      depms.res.all.sub <- depms.res.all[depms.res.all$day == i, ]
      
      res <- depms.res.all.sub %>% group_by(proteoform, day) %>%
        arrange(q_value, desc(abs(logFC))) %>%
        dplyr::slice(1:top_n) %>%
        group_by(proteoform, day) %>%
        arrange(desc(q_value), desc(abs(logFC))) %>%
        dplyr::slice(1)
      
      # Change sign 
      direction_change_idx <- sapply(strsplit(res$design_var, '-'), function(j) grep(i,j) == 2)
      res[direction_change_idx, c("logFC", "t",  "sca.t", "sca.t_weighted")] <- -res[direction_change_idx, c("logFC", "t",  "sca.t", "sca.t_weighted")]
      
      return(res)
      
    }))
    
    # Summarize across groups by min
    depms.res.all_across <- depms.res.all_within %>% group_by(proteoform) %>%
      arrange(desc(-q_value), desc(abs(sca.t_weighted)), .by_group = TRUE) %>%
      dplyr::slice(1)
    
    
    depms.res.all_across$design_var <-  depms.res.all_across$day
    depms.res.all_across$day <- NULL
    depms.res.all_across$ioi <- depms.res.all$ioi[match(depms.res.all_across$proteoform, depms.res.all$proteoform)]
    

    
  } else if (fitting == 'pca') {
    
    n_pcs <- 1
    # tokeepforpca <- unique(support_case$id)
    input_dat <- collapsed_quant_case
    tokeepforpca <- support_case$id[support_case$weight == 1]
    input_dat <- input_dat[tokeepforpca, ]
    
    pca_dt <- prcomp(t(input_dat[complete.cases(input_dat), ]))
    pc_vals <- pca_dt$x[,n_pcs, drop = FALSE]
    design_mat <- model.matrix(~ pc_vals)
    colnames(design_mat)[2] <- paste0('PC', n_pcs)
    
    depms.res.all_across <-  do.call(rbind, lapply(colnames(design_mat)[2], function(design_var) {
      
      # Background
      tofit_null <- collapsed_quant_background
      
      tokeep_var_null <- apply(tofit_null, 1, function(i) var(i, na.rm = TRUE) !=0)
      tofit_null <- tofit_null[tokeep_var_null, ]
      
      print(paste('background dataset:', nrow(tofit_null), ncol(tofit_null)))
      
      
      fit_null <- lmFit(tofit_null, design_mat)
      fit2_null <- fit_null
      fit3_null <- eBayes(fit2_null)
      fit3_null$count <- quant.for.variance(fit = fit3_null, 
                                            support_dt = support_background, 
                                            protein_peptide_annotation = protein_peptide_annotation,
                                            psm_vec = NULL)
      
      fit4_null <- spectraCounteBayesGAM(fit = fit3_null, coef_col = design_var)
      # VarianceBoxplotGAM(fit4_null,n=30,main="",xlab= "Count")
      # VarianceScatterplotGAM(fit4_null,main="")
      # ResidualplotGAM(fit4_null)
      
      DEpMS.results_null = outputResultGAM(fit = fit4_null, coef_col = design_var)
      
      ## Remove NAs
      DEpMS.results_null <- DEpMS.results_null[complete.cases(DEpMS.results_null), ]
      
      
      # Case 
      tofit <- collapsed_quant_case
      
      tokeep_var <- apply(tofit, 1, function(i) var(i, na.rm = TRUE) !=0)
      tofit <- tofit[tokeep_var, ]
      
      print(paste('case dataset:', nrow(tofit), ncol(tofit)))
      fit <- lmFit(tofit, design_mat)
      fit2 <- fit
      fit3 <- eBayes(fit2)
      
      ## Quantification for variance model
      fit3$count <- quant.for.variance(fit = fit3, 
                                       support_dt = support_case, 
                                       protein_peptide_annotation = protein_peptide_annotation,
                                       psm_vec = NULL)
      fit4 <- spectraCounteBayesGAM(fit3, coef_col = design_var)
      # VarianceBoxplotGAM(fit4,n=30,main="",xlab="Count")
      # VarianceScatterplotGAM(fit4,main="")
      # ResidualplotGAM(fit4)
      
      
      DEpMS.results = outputResultGAM(fit = fit4, coef_col = design_var)
      DEpMS.results <- DEpMS.results[complete.cases(DEpMS.results), ]
      
      # Q-value calculation
      res <- q.value.calculation(support_case = support_case,
                                 # support_background c= support_background,
                                 DEpMS.results_null = DEpMS.results_null, 
                                 DEpMS.results = DEpMS.results, 
                                 protein_peptide_annotation = protein_peptide_annotation, 
                                 design_var = design_var,
                                 plot_density = TRUE)
      
      return(res)
    }))
    
    logFC <- depms.res.all_across$logFC
    
        }
  ####################################################################################################################################
  
  
  # Summarize proteoforms
  summarized_proteoforms <- summarize.proteoforms(depms.res = depms.res.all_across,
                                                  collapsed_profiles = collapsed_profiles, 
                                                  protein_peptide_annotation = protein_peptide_annotation)
  
  # Add annotation
  annot_cols <-  c("n_pep", "n_perturbed_peptides", "perturbed_protein", "perturbed_peptide", "red_fac", "n_peptides")
  summarized_proteoforms[, annot_cols] <- protein_peptide_annotation[match(summarized_proteoforms$peptide, protein_peptide_annotation$Peptide_sequence), annot_cols]
  
  # return(list('summarized_proteoforms' = summarized_proteoforms,
  #             'logFC' = logFC))
  
  return(summarized_proteoforms)
  
}


depmsFDRbenchmark <- function(data, name="DEpMS", pval_col = "q_value_protein", pval_cutoff, fc_thres){
  
  
  depms_res <- copy(data)
  
  dt_input <- all_data[[unique(depms_res$annot)]]
  
  input_data_annotation <- unique(subset(dt_input,select=c("filename","day")))
  setorderv(input_data_annotation, c("day","filename"))
  input_data_annotation[,fraction_number:=.I]
  input_data_annotation <- subset(input_data_annotation,select=c("filename","fraction_number"))
  
  
  traces <- importPCPdata(input_data = dt_input,
                          fraction_annotation = input_data_annotation)
  
  trace_annotation <- unique(subset(dt_input, select=c("peptide_id",
                                                       "n_pep","n_perturbed_peptides",
                                                       "perturbed_protein","perturbed_peptide","red_fac")))
  
  traces <- annotateTraces(traces,
                           trace_annotation,
                           traces_id_column = "id",
                           trace_annotation_id_column = "peptide_id")
  
  ## Remove traces with 0 standard deviation (can't be clustered)
  zerovar <- apply(getIntensityMatrix(traces), 1, var) > 0
  traces_zerovar <- subset(traces,
                           trace_subset_ids = names(zerovar[zerovar]))
  
  #' ## Remove single peptide genes
  traces_multiPep <- filterSinglePeptideHits(traces_zerovar)
  
  score_thresholds = sort(unique(c(c(1 %o% 10^(-12:1)),seq(0,1,0.05))))
  
  depms_res_df <- data.table(score_thresholds=score_thresholds,
                             P = length(unique(traces_multiPep$trace_annotation[perturbed_protein==TRUE]$protein_id)),
                             N = length(unique(traces_multiPep$trace_annotation[perturbed_protein==FALSE]$protein_id)),
                             n_prot_proteoforms = 0,
                             TP = 0,
                             FP = 0,
                             N_prot_noMistakes = 0,
                             N_prot_withMistakes = 0,
                             TP_noMistakes = 0,
                             TP_withMistakes = 0
  )
  
  
  for (i in seq(1,length(score_thresholds),1)) {
    
    depms_res_sig <- subset(depms_res, get(pval_col)<=score_thresholds[i] & abs(depms_res$fc_protein) >= fc_thres) 
    
    depms_res_df$n_prot_proteoforms[i] <- length(unique(depms_res_sig$ioi))
    # res$n_prot_proteoforms[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1)]$protein_id))
    
    depms_res_df$TP[i] <- length(unique(depms_res_sig$ioi[depms_res_sig$perturbed_protein==TRUE]))
    depms_res_df$FP[i] <-  length(unique(depms_res_sig$ioi[depms_res_sig$perturbed_protein==FALSE]))
    
    depms_res_sig <- depms_res_sig %>% group_by(ioi, peptide) %>%
      mutate(n_proteoforms_per_perturbed_group = length(unique(id)))
    
    depms_res_sig <- depms_res_sig %>% group_by(ioi, id) %>%
      mutate(n_perturbed_groups_per_proteoform = length(unique(perturbed_peptide)))
    
    
    depms_res_sig$correct_peptide <- depms_res_sig$n_proteoforms_per_perturbed_group == 1 & depms_res_sig$n_perturbed_groups_per_proteoform == 1
    
    depms_res_sig <-  depms_res_sig %>% group_by(ioi) %>%
      mutate(protein_without_mistakes = all(correct_peptide))
    
    depms_res_df$N_prot_noMistakes[i] <- length(unique(depms_res_sig$ioi[depms_res_sig$protein_without_mistakes==TRUE]))
    depms_res_df$N_prot_withMistakes[i] <- length(unique(depms_res_sig$ioi[depms_res_sig$protein_without_mistakes==FALSE]))
    
    depms_res_df$TP_noMistakes[i] <- length(unique(depms_res_sig$ioi[depms_res_sig$perturbed_protein & depms_res_sig$protein_without_mistakes]))
    depms_res_df$TP_withMistakes[i] <- length(unique(depms_res_sig$ioi[depms_res_sig$perturbed_protein & !depms_res_sig$protein_without_mistakes]))
    
  }
  
  depms_res_df$FDR <- depms_res_df$FP/(depms_res_df$FP + depms_res_df$TP)
  depms_res_df$FN <- depms_res_df$P - depms_res_df$TP
  depms_res_df$F1 <- depms_res_df$TP/(depms_res_df$TP + 0.5*(depms_res_df$FP + depms_res_df$FN))
  depms_res_df$TPR <- depms_res_df$TP/depms_res_df$P
  depms_res_df$FPR <- depms_res_df$FP/depms_res_df$N
  depms_res_df$percent_TP_perfect <- depms_res_df$TP_noMistakes/depms_res_df$TP
  depms_res_df$percent_perfect <- depms_res_df$N_prot_noMistakes/(depms_res_df$P + depms_res_df$N)
  
  depms_res_df$pval_thresholds <- depms_res_df$score_thresholds
  
  
  
  
  
  depms_res_df$selected_pval <-ifelse(depms_res_df$pval_thresholds==pval_cutoff,TRUE,FALSE)
  pval_col <- c("TRUE"="red","FALSE"="#00ff0000")
  
  # pdf(paste0("ROC_percentPerfect_",name,"_pval_", pval_cutoff, ".pdf"), width=6.5, height=5)
  p <- ggplot(depms_res_df, aes(x=FPR,y=TPR)) + 
    geom_line(colour="grey",alpha=0.5) +
    geom_abline(intercept = 0, slope = 1, color='grey') +
    scale_fill_gradientn(limits = c(0,1),
                         colours=c("navyblue", "darkmagenta", "darkorange1"),
                         breaks=seq(0,1,0.25)) +
    geom_point(aes(fill=percent_TP_perfect, color=selected_pval,stroke=selected_pval), shape = 21, size = 3, stroke=2) +
    scale_color_manual(values = pval_col) +
    geom_abline(intercept = 0, slope = 1) +
    xlim(0,1) +
    ylim(0,1) +
    theme_classic()
  print(p)
  # dev.off()
  
  return(depms_res_df)
}


# Group-wise
random_perm_depms_group <- performDepmsAnalysis(dt_input = random_perm_data, fitting = 'group')
perm_1pep_depms_group <- performDepmsAnalysis(dt_input = perm_1pep_data, fitting = 'group')
perm_2pep_depms_group <- performDepmsAnalysis(dt_input = perm_2pep_data, fitting = 'group')
perm_025pep_depms_group <- performDepmsAnalysis(dt_input = perm_025pep_data, fitting = 'group')
perm_050pep_depms_group <- performDepmsAnalysis(dt_input = perm_050pep_data, fitting = 'group')

random_perm_depms_group$annot <- 'random'
perm_1pep_depms_group$annot <- '1pep'
perm_2pep_depms_group$annot <- '2pep'
perm_025pep_depms_group$annot <- '025pep'
perm_050pep_depms_group$annot <- '050pep'


scores_depms_group <- rbind(random_perm_depms_group, 
                     perm_1pep_depms_group,
                     perm_2pep_depms_group,
                     perm_025pep_depms_group,
                     perm_050pep_depms_group)

saveRDS(scores_depms_group, "scores_depms_group.rds")


fdr_bench_random_perm_depms_group <- depmsFDRbenchmark(data = random_perm_depms_group, name="random_perm", pval_col = "q_value_protein", pval_cutoff = 0.1, fc_thres = 0)
fdr_bench_perm_1pep_depms_group <- depmsFDRbenchmark(data = perm_1pep_depms_group, name="perm_1pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0)
fdr_bench_perm_2pep_depms_group <- depmsFDRbenchmark(data = perm_2pep_depms_group, name="perm_2pep",pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0)
fdr_bench_perm_025pep_depms_group <- depmsFDRbenchmark(data = perm_025pep_depms_group, name="perm_025pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0)
fdr_bench_perm_050pep_depms_group <- depmsFDRbenchmark(data = perm_050pep_depms_group, name="perm_050pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0)


res_depms_group <- rbind(fdr_bench_random_perm_depms_group[,nPerm := "random"],
                   fdr_bench_perm_1pep_depms_group[,nPerm := "1"],
                   fdr_bench_perm_2pep_depms_group[,nPerm := "2"],
                   fdr_bench_perm_025pep_depms_group[,nPerm := "25%"],
                   fdr_bench_perm_050pep_depms_group[,nPerm := "50%"])

saveRDS(res_depms_group, paste0(path_data, "res_depms_group.rds"))


# PCA
random_perm_depms_pca <- performDepmsAnalysis(dt_input = random_perm_data, fitting = 'pca')
perm_1pep_depms_pca <- performDepmsAnalysis(dt_input = perm_1pep_data, fitting = 'pca')
perm_2pep_depms_pca <- performDepmsAnalysis(dt_input = perm_2pep_data, fitting = 'pca')
perm_025pep_depms_pca <- performDepmsAnalysis(dt_input = perm_025pep_data, fitting = 'pca')
perm_050pep_depms_pca <- performDepmsAnalysis(dt_input = perm_050pep_data, fitting = 'pca')

random_perm_depms_pca$annot <- 'random'
perm_1pep_depms_pca$annot <- '1pep'
perm_2pep_depms_pca$annot <- '2pep'
perm_025pep_depms_pca$annot <- '025pep'
perm_050pep_depms_pca$annot <- '050pep'


scores_depms_pca <- rbind(random_perm_depms_pca, 
                         perm_1pep_depms_pca,
                         perm_2pep_depms_pca,
                         perm_025pep_depms_pca,
                         perm_050pep_depms_pca)

saveRDS(scores_depms_pca, "scores_depms_pca.rds")

fdr_bench_random_perm_depms_pca <- depmsFDRbenchmark(data = random_perm_depms_pca, name="random_perm", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0)
fdr_bench_perm_1pep_depms_pca <- depmsFDRbenchmark(data = perm_1pep_depms_pca, name="perm_1pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0)
fdr_bench_perm_2pep_depms_pca <- depmsFDRbenchmark(data = perm_2pep_depms_pca, name="perm_2pep",pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0 )
fdr_bench_perm_025pep_depms_pca <- depmsFDRbenchmark(data = perm_025pep_depms_pca, name="perm_025pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0)
fdr_bench_perm_050pep_depms_pca <- depmsFDRbenchmark(data = perm_050pep_depms_pca, name="perm_050pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0)



res_depms_pca <- rbind(fdr_bench_random_perm_depms_pca[,nPerm := "random"],
                       fdr_bench_perm_1pep_depms_pca[,nPerm := "1"],
                       fdr_bench_perm_2pep_depms_pca[,nPerm := "2"],
                       fdr_bench_perm_025pep_depms_pca[,nPerm := "25%"],
                       fdr_bench_perm_050pep_depms_pca[,nPerm := "50%"])

saveRDS(res_depms_pca, paste0(path_data, "res_depms_pca.rds"))


# Threshold group-wise 1
fdr_bench_random_perm_depms_group_thres <- depmsFDRbenchmark(data = random_perm_depms_group, name="random_perm", pval_col = "q_value_protein", pval_cutoff = 0.1, fc_thres = 1)
fdr_bench_perm_1pep_depms_group_thres <- depmsFDRbenchmark(data = perm_1pep_depms_group, name="perm_1pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 1)
fdr_bench_perm_2pep_depms_group_thres <- depmsFDRbenchmark(data = perm_2pep_depms_group, name="perm_2pep",pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 1)
fdr_bench_perm_025pep_depms_group_thres <- depmsFDRbenchmark(data = perm_025pep_depms_group, name="perm_025pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 1)
fdr_bench_perm_050pep_depms_group_thres <- depmsFDRbenchmark(data = perm_050pep_depms_group, name="perm_050pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 1)


res_depms_group_thres <- rbind(fdr_bench_random_perm_depms_group_thres[,nPerm := "random"],
                         fdr_bench_perm_1pep_depms_group_thres[,nPerm := "1"],
                         fdr_bench_perm_2pep_depms_group_thres[,nPerm := "2"],
                         fdr_bench_perm_025pep_depms_group_thres[,nPerm := "25%"],
                         fdr_bench_perm_050pep_depms_group_thres[,nPerm := "50%"])

saveRDS(res_depms_group_thres, paste0(path_data, "res_depms_group_thres.rds"))


# Threshold PCA 0.01
fdr_bench_random_perm_depms_pca_thres <- depmsFDRbenchmark(data = random_perm_depms_pca, name="random_perm", pval_col = "q_value_protein", pval_cutoff = 0.1, fc_thres = 0.01)
fdr_bench_perm_1pep_depms_pca_thres <- depmsFDRbenchmark(data = perm_1pep_depms_pca, name="perm_1pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0.01)
fdr_bench_perm_2pep_depms_pca_thres <- depmsFDRbenchmark(data = perm_2pep_depms_pca, name="perm_2pep",pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0.01)
fdr_bench_perm_025pep_depms_pca_thres <- depmsFDRbenchmark(data = perm_025pep_depms_pca, name="perm_025pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0.01)
fdr_bench_perm_050pep_depms_pca_thres <- depmsFDRbenchmark(data = perm_050pep_depms_pca, name="perm_050pep", pval_col = "q_value_protein", pval_cutoff =  0.1, fc_thres = 0.01)


res_depms_pca_thres <- rbind(fdr_bench_random_perm_depms_pca_thres[,nPerm := "random"],
                             fdr_bench_perm_1pep_depms_pca_thres[,nPerm := "1"],
                             fdr_bench_perm_2pep_depms_pca_thres[,nPerm := "2"],
                             fdr_bench_perm_025pep_depms_pca_thres[,nPerm := "25%"],
                             fdr_bench_perm_050pep_depms_pca_thres[,nPerm := "50%"])

saveRDS(res_depms_pca_thres, paste0(path_data, "res_depms_pca_thres.rds"))
#####################################
### Compare DEpMS & PeCorA & CCprofiler #####
#####################################
performBenchmarkComparison3 <- function(data1, data2, data3, method1, method2, method3, name="depms_copf_pecora_benchmark", ccprofiler_score_cutoff = 0.1){
  
  res1 <- copy(data1)
  res2 <- copy(data2)
  res3 <- copy(data3)

  
  res1[, method:=method1]
  res2[, method:=method2]
  res2[,score_thresholds_fac:=NULL]
  res3[, method:=method3]
  
  
  res2 <- setcolorder(res2, names(res1))
  res3 <- setcolorder(res3, names(res1))
  
  res3[,pval_col:=NULL]
 
  res_combi <- rbind(res1, res2, res3)
  
  
  
  cols <- c("PeCorA" = "blue4", 
            "COPF" = "green4", 
            "PeCorA_adj" = "blue4",
            "COPF_0.1" = "green4",
            'DEpMS-Group' = 'orange4', 
            "DEpMS-Group-thres" = 'orange4',
            'DEpMS-PCA' = 'orange1',
            'DEpMS-PCA-thres' = 'orange1'
            )
  
  pshapes <- c("PeCorA" = 15, 
               "PeCorA_adj" = 15, 
               "COPF" = 19,  
               "COPF_0.1" = 19,
               'DEpMS-Group' = 17,
               'DEpMS-Group-thres' = 17,
               'DEpMS-PCA' = 4,
               'DEpMS-PCA-thres' = 4
               )
  
  res_combi$method <- factor(res_combi$method, levels = unique(res_combi$method))
  
  pdf(paste0("ROC_benchmark_",name, ".pdf"), height=4, width=9)
  p <- ggplot(res_combi, aes(x=FPR, y= TPR, shape=method, color=percent_TP_perfect)) + 
    geom_line(alpha=0.5) +
    geom_point() +
    scale_color_gradientn(limits = c(0,1),
                          colours=c("navyblue", "darkmagenta", "darkorange1"),
                          breaks=seq(0,1,0.25)) +
    geom_point(data=res_combi[(selected_pval==TRUE)], color="red", shape=1, size=3) +
    scale_shape_manual(values = pshapes,  limits = force) +
    xlim(0,1) +
    ylim(0,1) +
    facet_wrap(. ~ nPerm) +
    theme_classic() +
    theme(legend.position="bottom") 
  print(p)
  dev.off()
  
  pdf(paste0("FDR_benchmark_",name,".pdf"), height=4, width=9)
  p <- ggplot(res_combi, aes(x=pval_thresholds, y= FDR, colour=method, shape=method, group=method)) + 
    geom_abline(intercept = 0, slope = 1, color='grey') +
    xlab("adj. p-value") +
    ylab("empirical FDR") +
    geom_line()+
    geom_point()+
    scale_color_manual(values = cols, limits = force) +
    scale_shape_manual(values = pshapes,  limits = force) +
    geom_point(data=res_combi[(selected_pval==TRUE)], color="red", shape=1, size=3) +
    xlim(0,1)+
    ylim(0,1)+
    facet_wrap(. ~ nPerm) +
    theme_classic() +
    theme(legend.position="bottom") +
    guides(shape = 'none')
  print(p)
  dev.off()
  
  return(res_combi)
}

res_combi_group <- performBenchmarkComparison3(
                                         data1 = res_depms_group,  
                                         data2 = res_copf[score_thresholds == 0], 
                                         data3 = res_pecora[pval_col=="adj_pval"], 
                                         method1 = 'DEpMS-Group',
                                         method2 = 'COPF',
                                         method3  = 'PeCorA',
                                         name="depms_copf0_pecora_benchmark")

saveRDS(res_combi_group, "res_combi_group.rds")

res_combi_pca <- performBenchmarkComparison3(data1 = res_depms_pca,  
                                         data2= res_copf[score_thresholds == 0], 
                                         data3 = res_pecora[pval_col=="adj_pval"], 
                                         method1 = 'DEpMS-PCA',
                                         method2 = 'COPF',
                                         method3  = 'PeCorA',
                                         name="DEpMS_PCA_copf0_pecora_benchmark")

saveRDS(res_combi_pca, "res_combi_pca.rds")



res_combi_group_thres <- performBenchmarkComparison3(data1 = res_depms_group_thres,  
                                                   data2= res_copf[score_thresholds == 0.1], 
                                                   data3 = res_pecora[pval_col=="adj_adj_pval"], 
                                                   method1 = 'DEpMS-Group-thres',
                                                   method2 = 'COPF_0.1',
                                                   method3  = 'PeCorA_adj',
                                                   name="DEpMS_Group_thres_copf01_pecoraAdj_benchmark")
saveRDS(res_combi_group_thres, "res_combi_group_thres.rds")



res_combi_pca_thres <- performBenchmarkComparison3(data1 = res_depms_pca_thres,  
                                             data2= res_copf[score_thresholds == 0.1], 
                                             data3 = res_pecora[pval_col=="adj_adj_pval"], 
                                             method1 = 'DEpMS-PCA-thres',
                                             method2 = 'COPF_0.1',
                                             method3  = 'PeCorA_adj',
                                             name="DEpMS_PCA_thres_copf01_pecoraAdj_benchmark")
saveRDS(res_combi_pca_thres, "res_combi_pca_thres.rds")





paper_sel <- c("1", "2", "50%")
performBenchmarkComparison3(data1 = res_depms_group[nPerm %in% paper_sel],  
                            data2 = res_copf[score_thresholds == 0][nPerm %in% paper_sel], 
                            data3 = res_pecora[nPerm %in% paper_sel][pval_col=="adj_pval"], 
                            method1 = 'DEpMS-Group',
                            method2 = 'COPF',
                            method3  = 'PeCorA',
                            name="depms_Group_copf_pecora_benchmark_PAPER")


performBenchmarkComparison3(data1 = res_depms_pca[nPerm %in% paper_sel],  
                            data2 = res_copf[score_thresholds == 0][nPerm %in% paper_sel], 
                            data3 = res_pecora[nPerm %in% paper_sel][pval_col=="adj_pval"], 
                            method1 = 'DEpMS-PCA',
                            method2 = 'COPF',
                            method3  = 'PeCorA',
                            name="depms_PCA_copf_pecora_benchmark_SUPPLEMENT")


performBenchmarkComparison3(data1 = res_depms_group_thres[nPerm %in% paper_sel],  
                            data2 = res_copf[score_thresholds == 0.1][nPerm %in% paper_sel], 
                            data3 = res_pecora[nPerm %in% paper_sel][pval_col=="adj_adj_pval"], 
                            method1 = 'DEpMS-Group-thres',
                            method2 = 'COPF_0.1',
                            method3  = 'PeCorA_adj',
                            name="depms_Group-thres_copf01_pecoraAdj_benchmark__PAPER1")


performBenchmarkComparison3(data1 = res_depms_pca_thres[nPerm %in% paper_sel],  
                            data2 = res_copf[score_thresholds == 0.1][nPerm %in% paper_sel], 
                            data3 = res_pecora[nPerm %in% paper_sel][pval_col=="adj_adj_pval"], 
                            method1 = 'DEpMS-PCA-thres',
                            method2 = 'COPF_0.1',
                            method3  = 'PeCorA_adj',
                            name="depms_PCA-thres_copf01_pecoraAdj_benchmark_SUPPLEMENT1")



####################################################################################################################################
## Overlap
thres_qvalue <- 0.1

# scores_pecora <- readRDS('scores_pecora.rds')
# scores_copf <- readRDS('scores_copf.rds')
# scores_depms <- readRDS('scores_depms_group.rds')
# scores_depms_pca <- readRDS('scores_depms_pca.rds')


lapply(unique(scores_depms_group$annot), function(j) {
  
  copf_ids <- unique(scores_copf$protein_id[Reduce(intersect, list(which(scores_copf$annot == j),
                                                          which(scores_copf$proteoform_score_pval_adj <= thres_qvalue),
                                                          which(scores_copf$perturbed_protein == TRUE)))])

  pecora_ids <- unique(scores_pecora$protein[Reduce(intersect, list(which(scores_pecora$annot == j),
                                                                   which(scores_pecora$adj_pval <= thres_qvalue),
                                                                   which(scores_pecora$perturbed_protein == TRUE)))])
  
  depms_group_ids <-  unique(scores_depms_group$ioi[Reduce(intersect, list(which(scores_depms_group$annot == j), 
                                                                      which(scores_depms_group$q_value_proteoform <= thres_qvalue), 
                                                                      which(scores_depms_group$perturbed_protein == TRUE)))])

  depms_pca_ids <- unique(scores_depms_pca$ioi[Reduce(intersect, list(which(scores_depms_pca$annot == j), 
                                                               which(scores_depms_pca$q_value_proteoform <= thres_qvalue), 
                                                               which(scores_depms_pca$perturbed_protein == TRUE)))])
  
 
  depms_group_thres_ids <- unique(scores_depms_group$ioi[Reduce(intersect, list(which(scores_depms_group$annot == j),
                                                                        which(scores_depms_group$q_value_proteoform <= thres_qvalue),
                                                                        which(abs(scores_depms_group$fc_protein) >= 1),
                                                                        which(scores_depms_group$perturbed_protein == TRUE)))])

  
  depms_pca_thres_ids <- unique(scores_depms_pca$ioi[Reduce(intersect, list(which(scores_depms_pca$annot == j),
                                                                                which(scores_depms_pca$q_value_proteoform < thres_qvalue),
                                                                                which(abs(scores_depms_pca$fc_protein) >= 0.01),
                                                                                which(scores_depms_pca$perturbed_protein == TRUE)))])
  # correct_ids <- union(scores_depms_pca$ioi[intersect(which(scores_depms_pca$perturbed_peptide == TRUE),
  #                                                     which(scores_depms_pca$annot == j))],
  #                      scores_depms$ioi[intersect(which(scores_depms$peptide == TRUE),
  #                                                 which(scores_depms_pca$annot == j))])
  
  depms_scores_list = list(
    copf_ids,
    pecora_ids,
    depms_group_ids,
    depms_pca_ids, 
    depms_group_thres_ids,
    depms_pca_thres_ids )
  names(depms_scores_list) <- c(
    'COPF',
    'PeCorA',
    'DEpMS-Group',
    'DEpMS-PCA', 
    'DEpMS-Group-thres',
    'DEpMS-PCA-thres')
  
  
  p <-  upset(data = fromList(depms_scores_list),  mb.ratio = c(0.7, 0.3),
              order.by = "freq", nsets = length(depms_scores_list), text.scale = 1.3)
  pdf(file = paste0(path_figures, 'Upset_plot_', j, '.pdf'), width = 6, height = 5.5, onefile = FALSE)
 
  print(p)
  dev.off()
  
  })

####################################################################################################################################
# pdf("ROC_benchmark_depms_copf_pecora_nPermComp.pdf", width=6, height=4)
# ggplot(res_combi_pca[nPerm != 'random'], aes(x=FPR, y= TPR, color=nPerm, group=nPerm)) + 
#   geom_abline(intercept = 0, slope = 1, color='grey') +
#   geom_line(alpha=0.5) +
#   geom_point() +
#   scale_color_manual(values =c("#A9C9A4","#699864","#308014","#2F4F2F")) +
#   #geom_point(data=res_combi[(selected_pval==TRUE)], color="red", shape=1, size=3) +
#   xlim(0,1) +
#   ylim(0,1) +
#   facet_wrap(. ~ method) +
#   theme_classic()
# dev.off()

###