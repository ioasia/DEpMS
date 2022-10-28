## Functions 

# Libraries
library(data.table)
library(dplyr)
library(Biostrings)
library(GenomicRanges)
library(igraph)
library(Hmisc)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(mgcv)
library(reshape2)
library(DEqMS)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(rtracklayer)
library(Gviz)
library(biomaRt)
library(httr)


## Peptide filtering
filter.peptides <- function(peptide_quant,
                            protein_quant, 
                            protein_peptide_annotation, 
                            ioi_id,
                            ioi_to_db_match,
                            db_match_fasta,
                            PSM_vec,
                            MS1_vec) {
  
  
  
  # Coordinates
  unique_iois <- unique(protein_peptide_annotation[[ioi_id]])
  
  pb  <- txtProgressBar(1, length(tocheck_ioi), style=3)
  
  misscleaved_modified <- do.call(rbind, lapply(1:length(unique_iois), function(i) {
    
    setTxtProgressBar(pb, i)
    ioi = unique_iois[i]
    
    # print(ioi)
    peptide_dt <- protein_peptide_annotation[which(protein_peptide_annotation[[ioi_id]] %in% ioi), ,drop = FALSE]
    peptide_dt$aa_seq <-  gsub('[^[A-Za-z]', '', peptide_dt$Peptide_sequence)
    unique_aa_seq <- unique(peptide_dt$aa_seq)
    
    
    db_match <- ioi_to_db_match[[ioi]]
    protein_seq <- db_match_fasta[intersect(db_match, names(db_match_fasta))]
    # protein_seq <- as.character(fasta_file[grep(k, names(fasta_file))])
    
    
    coords_total <- do.call(rbind, lapply(1:length(unique_aa_seq), function(m) {
      
      coords <- as.data.frame(unlist(vmatchPattern(pattern = unique_aa_seq[m], subject = protein_seq)))
      coords$aa_seq <-  unique_aa_seq[m]
      coords$protein <- sapply(strsplit(coords$names, ' '), function(o) o[1])
      coords
      
    }))
    
    coords_total <- do.call(rbind, rep(split(coords_total, coords_total$aa_seq),  table(peptide_dt$aa_seq)))
    
    peptide_coords <- GRanges(
      seqnames = coords_total$protein,
      ranges = IRanges(start = coords_total$start, end = coords_total$end),)
    
    
    overlap_pairs <- as.data.frame(GenomicRanges::findOverlaps(peptide_coords, peptide_coords, type = 'within'))
    overlap_pairs$queryHits <- coords_total$aa_seq[overlap_pairs$queryHits]
    overlap_pairs$subjectHits <- coords_total$aa_seq[overlap_pairs$subjectHits]
    
    
    # Remove single hits
    freq_pairs <- table(overlap_pairs$queryHits)
    overlap_pairs <- overlap_pairs[overlap_pairs$queryHits %in% names(freq_pairs)[freq_pairs > 1], ]
    
    if(nrow(overlap_pairs) == 0) {
      return(NULL)
    }
    
    graph_pep <- graph_from_data_frame(overlap_pairs, directed = FALSE)
    
    # pdf(paste0(path_figures,'peptide_graphs/', k, '_graph.pdf'), width = 10, height = 10)
    # plot(simplify(graph_pep))
    # dev.off()
    
    
    dg <- decompose.graph(simplify(graph_pep))
    
    graph_pep_overlap <- lapply(dg, function(g) {
      names(V(g))
    })
    names(graph_pep_overlap) <-  paste0(ioi, '_', 1:length(graph_pep_overlap))
    overlap_pairs <- reshape2::melt(graph_pep_overlap)
    
    peptide_dt$group <- overlap_pairs$L1[match(peptide_dt$aa_seq,overlap_pairs$value)]
    peptide_dt <- peptide_dt[!is.na(peptide_dt$group), ]
    
    peptide_dt$sim <-  as.numeric(cor(t(peptide_quant[peptide_dt$Peptide_sequence, ]),  
                                      t(protein_quant[ioi, ]), method = cor_method, use = 'pairwise.complete.obs'))
    
    
    return(peptide_dt)
    
  }
  ))
  
  
  ## add MS1
  if(!is.null(MS1_vec)) {
    
    # Boxplot
    boxplot_dt <-  misscleaved_modified[, c("Peptide_sequence","sim"), with = FALSE]
    boxplot_dt$MS1 <- MS1_vec[misscleaved_modified$Peptide_sequence]
    
    boxplot_dt$group_MS1 <- cut2(boxplot_dt$MS1, m = nrow(boxplot_dt)/10)
    p <- ggplot(boxplot_dt, aes(x = group_MS1 , y = sim)) +
      geom_violin() +
      geom_boxplot(width = 0.1) +
      labs(x = 'Average MS1 intensity' , y = 'Peptide-protein similarity') +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(hjust = 1, angle = 45))
    
    pdf(paste0(path_figures, 'protein_peptide_sim_MS1.pdf'), width = 8, height = 5)
    plot(p)
    dev.off()
  }
  
  if(!is.null(PSM_vec)) {
    boxplot_dt <-  misscleaved_modified[,  c("Peptide_sequence","sim"),  with = FALSE]
    boxplot_dt$min_PSM <- min_PSM_vec[misscleaved_modified$Peptide_sequence]
    
    boxplot_dt$group_min_PSM <- cut2(boxplot_dt$min_PSM, m = nrow(boxplot_dt)/10)
    p <- ggplot(boxplot_dt, aes(x = group_min_PSM , y = sim)) +
      geom_violin() +
      geom_boxplot(width = 0.1) +
      labs(x = 'min(PSMs)' , y = 'Peptide-protein similarity') +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    
    pdf(paste0(path_data, 'protein_peptide_sim_min_PSM.pdf'), width = 8, height = 5)
    plot(p)
    dev.off()
    
  }
  
  
  # Summarize modifications
  sorted_misscleaved_modified <- misscleaved_modified %>% 
    group_by(group) %>%
    arrange(desc(sim), .by_group = TRUE) %>%
    dplyr::slice(1)
  
  protein_peptide_annotation_misscleaved_modified <- protein_peptide_annotation[protein_peptide_annotation$Peptide_sequence %in% sorted_misscleaved_modified$Peptide_sequence, ]
  protein_peptide_annotation_rest <- protein_peptide_annotation[!protein_peptide_annotation$Peptide_sequence %in% misscleaved_modified$Peptide_sequence, ]
  
  final_tab <- rbind.data.frame(protein_peptide_annotation_misscleaved_modified, protein_peptide_annotation_rest)
  
  
  return(final_tab)
  
}
###################################################################################################################################################
## peptide - protein similarity 
peptide.protein.similarity <- function(peptide_quant, protein_quant, protein_peptide_annotation, tocheck_ioi) {
  
  pb  <- txtProgressBar(1, length(tocheck_ioi), style=3)
  
  protein_peptide_sim <- do.call(rbind, lapply(1:length(tocheck_ioi), function(i) {
    
    setTxtProgressBar(pb, i)
    
    ioi = tocheck_ioi[i]
    protein_quant_ioi <- protein_quant[ioi, ]
    id_peps <- protein_peptide_annotation$Peptide_sequence[protein_peptide_annotation[[ioi_id]] %in% ioi]
    peptide_quant_ioi <- peptide_quant[id_peps, ]
    
    peptide_sim_ioi <- cor(t(peptide_quant_ioi),t(protein_quant_ioi), method = cor_method, use = 'pairwise.complete.obs')
    colnames(peptide_sim_ioi) <- NULL
    
    res <- data.frame(ioi = ioi,sim = peptide_sim_ioi,n_pep = length(id_peps))
  }))
  
  # Median similarity
  protein_peptide_sim <- protein_peptide_sim[protein_peptide_annotation$Peptide_sequence, ]
  protein_peptide_sim_median <- as.data.frame(protein_peptide_sim %>% group_by(ioi) %>%
                                                summarise(m = median(sim, na.rm = TRUE) ))
  rownames(protein_peptide_sim_median) <- protein_peptide_sim_median$ioi
  protein_peptide_sim_median$ioi <- NULL
  protein_peptide_sim_median$n_pep <- protein_peptide_sim$n_pep[match(rownames(protein_peptide_sim_median), protein_peptide_sim$ioi)]
  
  
  n_ioi <- nrow(protein_peptide_sim_median)
  p_sim_median <- ggplot(protein_peptide_sim_median, mapping = aes(x = log2(n_pep), y = m)) +
    geom_point() +
    geom_smooth(method = 'gam') + 
    annotate(geom = 'text', x = Inf, y = Inf, label = bquote(n[ioi] == .(n_ioi)), hjust = 1.2, vjust = 1.2) +
    labs(x = 'Number of peptides per protein - log2', y = 'Median peptide-protein similarity') +
    theme_classic() 
  
  pdf(paste0(path_figures, 'protein_peptide_median_sim.pdf'), width = 6, height= 4.5)
  plot(p_sim_median)
  dev.off()
  
  
  # IQR
  protein_peptide_sim_iqr <- as.data.frame(protein_peptide_sim %>% group_by(ioi) %>%
                                             summarise(m = IQR(sim, na.rm = TRUE) ))
  rownames(protein_peptide_sim_iqr) <- protein_peptide_sim_iqr$ioi
  protein_peptide_sim_iqr$ioi <- NULL
  protein_peptide_sim_iqr$n_pep <- protein_peptide_sim$n_pep[match(rownames(protein_peptide_sim_iqr), protein_peptide_sim$ioi)]
  
  
  p_sim_iqr <- ggplot(protein_peptide_sim_iqr, mapping = aes(x = log2(n_pep), y = m)) +
    geom_point() +
    geom_smooth(method = 'gam') + 
    annotate(geom = 'text', x = Inf, y = Inf, label = bquote(n[ioi] == .(n_ioi)), hjust = 1.2, vjust = 1.2) +
    labs(x = 'Number of peptides per protein - log2', y = 'Similarity interquartile range') +
    theme_classic() 
  
  pdf(paste0(path_figures, 'protein_peptide_median_iqr.pdf'), width = 6, height= 4.5)
  plot(p_sim_iqr)
  dev.off()
  
  # Scaling factor
  min_peptide_support <- min(protein_peptide_sim_median$n_pep)
  scaling_factor_dt <- cbind(protein_peptide_sim_median,  'iqr' = protein_peptide_sim_iqr$m)
  scaling_factor_dt$group <- ifelse(scaling_factor_dt$n_pep > min_peptide_support, paste0('> ',min_peptide_support, ' peptides'),  paste0('= ',min_peptide_support, ' peptides'))
  
  
  scaling_factor <- (median(scaling_factor_dt$m[scaling_factor_dt$n_pep == min_peptide_support]) -  median(scaling_factor_dt$m[scaling_factor_dt$n_pep > min_peptide_support]))/ 
    (sum(scaling_factor_dt$n_pep == min_peptide_support)/nrow(scaling_factor_dt) * median(scaling_factor_dt$iqr[scaling_factor_dt$n_pep == min_peptide_support]) + sum(scaling_factor_dt$n_pep > min_peptide_support)/nrow(scaling_factor_dt) * median(scaling_factor_dt$iqr[scaling_factor_dt$n_pep > min_peptide_support]))
  
  scaling_factor <- max(0, scaling_factor)
  
  
  p_scaling_factor <- ggplot(scaling_factor_dt,aes(x = group, y = m)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.color = NA) +
    geom_hline(yintercept = c(median(scaling_factor_dt$m[scaling_factor_dt$n_pep == min_peptide_support]),
                              median(scaling_factor_dt$m[scaling_factor_dt$n_pep > min_peptide_support])),col ='red', lty = 2) +
    labs(x = '', y = 'Median peptide-protein similarity') +
    annotate(geom = 'text', x = Inf, y = Inf, label = paste0('scaling factor = ', round(scaling_factor, 2)),
             hjust = 1.2, vjust = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 1, angle = 45))
  
  
  pdf(paste0(path_figures, 'protein_peptide_iqr_scaling_factor.pdf'), width = 5, height= 5)
  plot(p_scaling_factor)
  dev.off()
  
  
  
  
  list('median_sim' = protein_peptide_sim_median,
       'iqr_sim' = protein_peptide_sim_iqr,
       'scaling_factor' = scaling_factor)
  
}
###################################################################################################################################################
## Putative proteoform (profile) generation
profile.generation <- function(peptide_quant, protein_quant, protein_peptide_annotation, global_sim, tocheck_ioi) {
  
  scaling_factor <- global_sim$scaling_factor
  
  ## Prediction models
  global_sim_models <- list('median_sim_mod' = gam(m ~ s(log2(n_pep),  bs = "cs"),
                                                   data = global_sim$median_sim, method="REML"),
                            'iqr_sim_mod' = gam(m ~ s(log2(n_pep),  bs = "cs"),
                                                data = global_sim$iqr_sim, method="REML") )
  
  pb  <- txtProgressBar(1, length(tocheck_ioi), style=3)
  
  collapsed_quant <- lapply(1:length(tocheck_ioi), function(i) {

    # print(i)
    setTxtProgressBar(pb, i)
    ioi = tocheck_ioi[i]
    
    # Protein profile
    protein_quant_ioi <- protein_quant[ioi, ]
    
    # Peptide profile
    id_peps <- protein_peptide_annotation$Peptide_sequence[grep(paste0('^', ioi, '$'), protein_peptide_annotation[[ioi_id]])]
    peptide_quant_ioi <- peptide_quant[id_peps,  ]
    
    
    ##########################################################################################
    # Check network modularity for protein with n_pep >= 4. If modularity >= 0.3, cluster according to community detection 
    if(length(id_peps) >=4) {
      peptide_quant_sim <-  reshape2::melt(cor(t(peptide_quant_ioi[id_peps,]), method = cor_method, use = 'pairwise.complete.obs'))
      
      peptide_quant_sim <- peptide_quant_sim[as.character(peptide_quant_sim$Var1) < as.character(peptide_quant_sim$Var2), ]
      peptide_quant_sim$value[peptide_quant_sim$value < 0] <- 0
      g <- graph_from_data_frame(peptide_quant_sim, directed = FALSE)
      cluster_res <- cluster_louvain(g, weights = E(g)$value)
      
      modularity_val <- cluster_res$modularity[1]
    } else {
      modularity_val = 0
    }
    
    # Case of good clusters
    if(modularity_val >= 0.3) {
      # print(ioi)
      membership_g <- sort(membership(cluster_res))
      collapsed_profiles <- do.call(rbind, lapply(unique(membership_g), function(j) {
        dt <- peptide_quant_ioi[names(membership_g)[membership_g == j], , drop = FALSE]
        apply(dt, 2, median, na.rm = TRUE)
      }))
      
      rownames(collapsed_profiles) <- paste0('P', unique(membership_g), '_mod_', ioi)
      collapsed_profiles_dtp <- t(t(collapsed_profiles) - as.numeric(protein_quant_ioi))
      
      weight <- apply(collapsed_profiles_dtp, 1, function(i) min(1,sd(as.numeric(i), na.rm = TRUE)/
                                                                   sd(as.numeric(protein_quant_ioi), na.rm = TRUE)))
      
      membership_dtp <- data.frame('peptide' = names(membership_g),
                                   'membership' = membership_g,
                                   'id' =  paste0('P', membership_g, '_mod_', ioi),
                                   'ioi' = ioi)
      membership_dtp$weight <- weight[membership_dtp$id]
      
      return(list('profile' = collapsed_profiles_dtp,
                  'membership' = membership_dtp)
      )
    }
    
    #########################################################################################
    # Peptide - protein similarity
    res_sim <- as.numeric(cor(t(peptide_quant_ioi),
                              t(protein_quant_ioi), method = cor_method, use = 'pairwise.complete.obs'))
    
    # Similar-to-protein (stp) peptides
    thres_sim <- as.numeric(predict(global_sim_models$median_sim_mod, newdata = global_sim$median_sim[ioi, ])) 
    peptide_quant_stp <- peptide_quant_ioi[res_sim >= thres_sim & !is.na(res_sim), , drop = FALSE]
    stp_id_peps <- rownames(peptide_quant_stp)
    
    memberiship_stp <- data.frame('peptide' = stp_id_peps,
                                  'membership'= rep('sim',length(stp_id_peps)),
                                  'id' = rep(paste0('P_sim_', ioi),length(stp_id_peps)),
                                  'ioi' = rep(ioi, length(stp_id_peps)),
                                  'weight' = rep(1, length(stp_id_peps)))
    
    
    # Dissimilar-to-protein (dtp) peptides
    thres_dtp <- thres_sim - as.numeric(predict(global_sim_models$iqr_sim_mod, newdata = global_sim$iqr_sim[ioi, ])) * scaling_factor
    peptide_quant_dtp <- peptide_quant_ioi[res_sim < thres_dtp & !is.na(res_sim), , drop = FALSE]
    dtp_id_peps <- rownames(peptide_quant_dtp)
    
    
    # Case of zero (dis)similar-to-protein peptides 
    if(length(stp_id_peps) + length(dtp_id_peps) == 0) {
      return(NULL)
    }
    
    
    # Case of no similar-to-protein peptides
    if(length(stp_id_peps) == 0) {
      
      collapsed_profile_all <- apply(peptide_quant_ioi, 2, median, na.rm = TRUE)
      
    }  else {
      
      collapsed_profile_all <- apply(peptide_quant_stp, 2, median, na.rm = TRUE)
      
    }
    
    #  Case of no dissimilar-to-protein peptides
    if(length(dtp_id_peps) == 0) {
      
      collapsed_profile_main <- t(collapsed_profile_all) - as.numeric(protein_quant_ioi)
      rownames(collapsed_profile_main) <- paste0('main_', ioi)
      
      cor_p <- cor.test(as.numeric(collapsed_profile_main),  as.numeric(protein_quant_ioi), 
                        use = 'pairwise.complete.obs', method = cor_method)
      
      if(is.na(cor_p$p.value)) {
        return(NULL)
      } else if(cor_p$p.value > 0.05) {
        return(list('profile' = collapsed_profile_main,
                    'membership' = data.frame('peptide' = stp_id_peps, 
                                              'membership' = 'main', 
                                              'id' =  rownames(collapsed_profile_main),
                                              'ioi' = ioi,
                                              'weight' = 1)
        )
        )
      } else {
        return(NULL)
      }
    } 
    
    
    # Consider <4 peptides cases
    # Case of single dissimilar-to-protein peptide
    if(length(dtp_id_peps) == 1) {
      
      
      
      # Difference from main profile (median)
      collapsed_profiles_dtp <- t(t(peptide_quant_dtp) - as.numeric(collapsed_profile_all))
      
      
      rownames(collapsed_profiles_dtp) <- paste0('P1_', ioi)
      weight <- min(1,sd(as.numeric(peptide_quant_dtp), na.rm = TRUE)/
                      sd(as.numeric(collapsed_profile_all), na.rm = TRUE))
      memberiship_dtp <- data.frame('peptide' = dtp_id_peps,
                                    'membership' = '1',
                                    'id' =  rownames(collapsed_profiles_dtp),
                                    'ioi' = ioi,
                                    'weight' = weight)
      
      if(length(stp_id_peps) == 0) {
        collapsed_profile_main <- NULL
        membership_main <- NULL
      } else {
        collapsed_profile_main <- t(collapsed_profile_all) - as.numeric(protein_quant_ioi)
        rownames(collapsed_profile_main) <- paste0('main_', ioi)
        
        cor_p <- cor.test(as.numeric(collapsed_profile_main),  as.numeric(protein_quant_ioi), 
                          use = 'pairwise.complete.obs', method = cor_method)
        
        if(is.na(cor_p$p.value)) {
          collapsed_profile_main <- NULL
          membership_main <- NULL
        } else if(cor_p$p.value < 0.05) {
          collapsed_profile_main <- NULL
          membership_main <- NULL
        } else {
          membership_main <- data.frame('peptide' = stp_id_peps, 
                                        'membership' = 'main', 
                                        'id' =  rownames(collapsed_profile_main),
                                        'ioi' = ioi,
                                        'weight' = 1)
        }
      }
      
      
      return(list('profile' = rbind(collapsed_profile_main, collapsed_profiles_dtp),
                  'membership' = rbind(membership_main, memberiship_dtp, memberiship_stp))
      )
    }
    
    ## Graph construction
    peptide_quant_sim <- reshape2::melt(cor(t(peptide_quant_ioi[dtp_id_peps, ]), method = cor_method, use = 'pairwise.complete.obs'))
    peptide_quant_sim <- peptide_quant_sim[as.character(peptide_quant_sim$Var1) > as.character(peptide_quant_sim$Var2), ]
    
    peptide_quant_sim$value[(peptide_quant_sim$value) < thres_sim] <- NA
    peptide_quant_sim <- peptide_quant_sim[!is.na(peptide_quant_sim$value), ]
    
    
    # graph components
    sim_graph <-  graph_from_data_frame(peptide_quant_sim, directed = FALSE)
    toadd <-  setdiff(dtp_id_peps, names(V(sim_graph)))
    
    # Add singletons
    if(length(toadd) > 0) {
      sim_graph <- add_vertices(graph = sim_graph, nv = length(toadd))
      V(sim_graph)$name[rev(1:length(V(sim_graph)))[1:length(toadd)]] <- toadd
    }
    
    # plot(sim_graph,vertex.label=NA)
    # components(sim_graph)
    
    # Check graph coreness for graphs > 2 nodes
    dg <- decompose.graph(sim_graph)
    dg_size <- unlist(lapply(dg, function(i) length(V(i))))
    graphstocheck <- which(dg_size > 2)
    dg <- decompose.graph(sim_graph)
    dg_size <- unlist(lapply(dg, function(i) length(V(i))))
    graphstocheck <- which(dg_size > 3)
    if(length(graphstocheck) != 0) {
      
      sim_graph_sub <- lapply(graphstocheck, function(i) {
        g_sub <- dg[[i]]
        # plot(g_sub, vertex.label=NA)
        # tokeepsimg <- coreness(g_sub) >= round(log2(length(V(g_sub))))
        tokeepsimg <- igraph::degree(g_sub) > 1
        g_sub_con <- induced_subgraph(g_sub, names(tokeepsimg)[tokeepsimg])
        # plot(g_sub_con, vertex.label=NA)
        if(length(V(g_sub_con)) > 1) {
          # Create data.frame
          g_sub_con_peps <- names(V(g_sub_con))
          graph_dt <- reshape2::melt(cor(t(peptide_quant_ioi[c(g_sub_con_peps,stp_id_peps), ]), method = cor_method, use = 'pairwise.complete.obs'))
          
          graph_dt <- graph_dt[as.character(graph_dt$Var1) > as.character(graph_dt$Var2), ]
          graph_dt <- graph_dt[graph_dt$value > 0 & !is.na(graph_dt$value), ]
          g_sub_con <- graph_from_data_frame(graph_dt, directed = FALSE)
          # plot(g_sub_con, vertex.label=NA)
          
          g_sub_con_clus <- cluster_louvain(g_sub_con, weights = E(g_sub_con)$value)
          g_sub_con_clus_member <- membership(g_sub_con_clus)[g_sub_con_peps] 
          n_clust <- unique(g_sub_con_clus_member)
          g_clust <- lapply(n_clust, function(j) {
            g_clust_subg <- induced_subgraph(g_sub_con, membership(g_sub_con_clus) == j)
            g_clust_subg <- induced_subgraph(g_clust_subg, names(V(g_clust_subg)) %in% g_sub_con_peps)
            
            # plot(g_clust_subg)
          })
        } else {
          NULL
        }
      })
      
      sim_graph_sub <- unlist(sim_graph_sub, recursive = FALSE)
      
      # Remove 0 edge graphs
      if(is.null(sim_graph_sub)) {
        sim_graph <- Reduce(igraph::union, dg[-graphstocheck])
      } else {
        sim_graph_subs <- do.call(igraph::union, sim_graph_sub)
        # plot(sim_graph_subs, vertex.label=NA)
        if(length(V(sim_graph_subs)) > 0) {
          sim_graph <- igraph::union(sim_graph_subs, dg[-graphstocheck])
          # plot(sim_graph)
        } else {
          sim_graph <- Reduce(igraph::union, dg[-graphstocheck])
        }
      }
    }
    # plot(sim_graph, vertex.label=NA)
    
    if(is.null(sim_graph)) {
      return(NULL)
    }
    
    # Collapse connected components
    sim_graph_members <- components(sim_graph)$membership
    sim_graph_members <- data.frame('peptide' = names(sim_graph_members),
                                    'membership' = sim_graph_members,
                                    'id' =  paste0('P',sim_graph_members, '_', ioi),
                                    'ioi' = ioi)
    
    collapsed_profiles <- do.call(rbind, lapply(unique(sim_graph_members$membership), function(j) {
      dt <- peptide_quant_ioi[sim_graph_members$peptide[sim_graph_members$membership == j], , drop = FALSE]
      apply(dt, 2, median, na.rm = TRUE)
    }))
    
    collapsed_profiles <- rbind(collapsed_profile_all,
                                collapsed_profiles)
    
    collapsed_profiles_names <- c(paste0('main', '_', ioi),
                                  paste0('P', unique(sim_graph_members$membership), '_', ioi))
    
    rownames(collapsed_profiles) <- collapsed_profiles_names
    
    # Peptide support
    pep_support <- sort(table(sim_graph_members$membership), decreasing = TRUE)
    names(pep_support) <- collapsed_profiles_names[as.numeric(names(pep_support)) + 1]
    
    
    # Re-factor according to peptide support
    collapsed_profiles_names <- c(collapsed_profiles_names[1],names(pep_support))
    
    # Filter collapsed patterns that correlate
    sim_collapsed_profiles <- reshape2::melt(cor(t(collapsed_profiles[collapsed_profiles_names, ]),method = cor_method, use = 'pairwise.complete.obs'))
    sim_collapsed_profiles <- sim_collapsed_profiles[as.character(sim_collapsed_profiles$Var1) < as.character(sim_collapsed_profiles$Var2), ]
    sim_collapsed_profiles <- sim_collapsed_profiles[sim_collapsed_profiles$value >= thres_dtp & !is.na(sim_collapsed_profiles$value), ]
    
    if(nrow(sim_collapsed_profiles) > 0) {
      sim_collapsed_profiles[, 1:2] <- do.call(rbind, lapply(1:nrow(sim_collapsed_profiles), function(i) {
        res <- as.character(sort(factor(unlist(sim_collapsed_profiles[i,c(1:2)]), levels = collapsed_profiles_names)))
        
      }))}
    
    tokeep_final <- setdiff(rownames(collapsed_profiles),sim_collapsed_profiles$Var2)
    
    
    if(length(tokeep_final) == 0) {
      return(NULL)
    }
    
    collapsed_profiles <- collapsed_profiles[c(paste0('main', '_', ioi), 
                                               tokeep_final), , drop = FALSE]
    
    # Weight by standard deviation ratio
    weight <- apply(collapsed_profiles[-1, , drop = FALSE], 1, function(i) min(1,sd(as.numeric(i), na.rm = TRUE)/
                                                                                 sd(as.numeric(collapsed_profile_all), na.rm = TRUE)))
    sim_graph_members$weight <- weight[sim_graph_members$id]
    
    collapsed_profiles <- t(t(collapsed_profiles[-1, , drop = FALSE]) - as.numeric(collapsed_profiles[1, ]))
    memberiship_dtp <- sim_graph_members[sim_graph_members$id %in% tokeep_final,]
    
    
    if(length(stp_id_peps) > 0) {
      
      
      collapsed_profile_main <- t(collapsed_profile_all) - as.numeric(protein_quant_ioi)
      rownames(collapsed_profile_main) <- paste0('main_', ioi)
      cor_p <- cor.test(as.numeric(collapsed_profile_main),  as.numeric(protein_quant_ioi), 
                        use = 'pairwise.complete.obs', method = cor_method)
      if(is.na(cor_p$p.value)) {
        collapsed_profile_main <- NULL
        membership_main <- NULL
      } else if(cor_p$p.value < 0.05) {
        collapsed_profile_main <- NULL
        membership_main <- NULL
      } else {
        membership_main = data.frame('peptide' = stp_id_peps, 
                                     'membership' = 'main', 
                                     'id' =  paste0('main_', ioi),
                                     'ioi' = ioi,
                                     'weight' = 1)
      }
    } else {
      collapsed_profile_main <- NULL
      membership_main <- NULL
    }
    
    return(list(
      # 'profile' = collapsed_profiles,
      'profile' = rbind(collapsed_profile_main, collapsed_profiles),
      'membership' = rbind(membership_main, memberiship_dtp, memberiship_stp))
    )
    # }
  }
  )
  
  # Collapsed data sets
  collapsed_quant_profiles <- do.call(rbind, lapply(collapsed_quant, function(i) {
    i[[1]]
  }))
  
  main_profiles_idx <- grep('main', rownames(collapsed_quant_profiles))
  collapsed_quant_background <- collapsed_quant_profiles[main_profiles_idx, ]
  collapsed_quant_case <- collapsed_quant_profiles[-main_profiles_idx, ]
  
  
  # Peptide support
  collapsed_quant_support <- do.call(rbind, lapply(collapsed_quant, function(i) {
    i[[2]]
  }))
  main_support_idx <- grep('main',collapsed_quant_support$id)
  support_background <- collapsed_quant_support[main_support_idx, ]
  
  support_stp_idx <- grep('sim',collapsed_quant_support$id)
  support_stp <- collapsed_quant_support[support_stp_idx, ]
  
  
  support_case <- collapsed_quant_support[-c(main_support_idx,support_stp_idx), ]
  
  
  res <- list('collapsed_quant_case' = collapsed_quant_case, 
              'support_case' = support_case,
              'collapsed_quant_background' = collapsed_quant_background,
              'support_background' = support_background,
              'support_stp' = support_stp)
  
  return(res)
  
}
###################################################################################################################################################
## Choose number of PCs
plot.expl.variance <- function(pca_res, plot_title) {
  
  # variance explained
  pca_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
  
  cumsum_var <- cumsum(pca_explained)
  
  ## Plot explained variance significance
  expl_var <- data.frame(
    'pcs' = 1:length(pca_explained),
    'Variance'= 100 * pca_explained,
    'Cumulative_variance' = 100 * cumsum_var)
  expl_var <- reshape2::melt(expl_var, id.vars = c('pcs'))
  
 
  
  p_var_expl <- ggplot(expl_var, aes(x = pcs, y = value, group = variable)) +
    geom_point() +
    geom_line() + 
    facet_wrap(.~variable, ncol = 1, scales = 'free') +
    labs(x = 'PCs', y = 'Variance explained (%)', title =  plot_title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}


plot.pc.associations <- function(pca_res, tocheckmetadata, n_pcs, plot_title) {
  
  ## Correlation of PCs with other variables
  pcs_vec <-  pca_res$x[, c(1:n_pcs), drop = FALSE]
  
  pc_variable_pvals <- do.call(cbind, lapply(1:ncol(pcs_vec), function(i) {
    
    pc_values <- pcs_vec[,i]
    res <- do.call(rbind, lapply(1:ncol(tocheckmetadata), function(j) {
      
      var_vec <- tocheckmetadata[,j]
      
      if(class(var_vec) %in% c('character', 'factor')) {
        
        # Kruskal-Wallis test
        kw_test <- kruskal.test(pc_values, var_vec)
        
        pval = kw_test$p.value
        
      } else if(class(var_vec)  == 'numeric') {
        
        spearman_test <- cor.test(pc_values, var_vec, method = 'spearman')
        pval <- spearman_test$p.value
        
      } else {return(NULL)}
      
    }))
  }))
  colnames(pc_variable_pvals) <- colnames(pcs_vec)
  rownames(pc_variable_pvals) <- colnames(tocheckmetadata)
  
  
  # Plot associations
  pc_variable_pvals_dt <- reshape2::melt(pc_variable_pvals)
  
  
  pc_variable_pvals_plot <- ggplot(pc_variable_pvals_dt, aes(x = Var1, y = -log10(value))) +
    geom_point(size = 3) +
    geom_point(data = pc_variable_pvals_dt[pc_variable_pvals_dt$value < 0.05, ], col = 'red') +
    labs(x = '', y = 'nominal p-value (-log10)', title = plot_title, subtitle =  'PCs-variables corrrelation') + 
    facet_grid(Var2~.) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
}  


plot.pcs <- function(pca_res, n_pcs, var_to_plot, cols_var, plot_title)  {
  
  # variance explained
  pr_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
  
  ## Plot explained variance significance
  expl_var <- data.frame(
    'pcs' = 1:length(pr_explained),
    'var'= pr_explained)
  
  ## Plots
  pca_plot_dt <- data.frame('pca' = pca_res$x[, c(1:n_pcs), drop = FALSE],
                            # 'variable' = as.character(tocheck$IDH1.y)
                            'variable' = var_to_plot)
  
  
  
  combination_pcs <- combn(1:n_pcs, 2)
  pca_plots <- apply(combination_pcs, 2, function(i) {
    
    pca_protein_plot <- ggplot(data = pca_plot_dt, aes(x = pca_plot_dt[, i[1]], y = pca_plot_dt[, i[2]], 
                                                       fill = variable)) + 
      geom_point(pch = 21, size = 4, col = 'black') + 
      labs(x = paste0('PCA-', i[1], ' (',paste0(round(pr_explained[i[1]], 2) * 100, '% of variance explained)')), 
           y = paste0('PCA-', i[2], ' (',paste0(round(pr_explained[i[2]],2) * 100, '% of variance explained)')),
           title = '') + 
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    if(class(pca_plot_dt$variable) == 'numeric') {
      pca_protein_plot <- pca_protein_plot + scale_fill_gradient(low = 'white',high = 'red')
    } 
    if(!is.null(cols_var)) {
      pca_protein_plot <- pca_protein_plot + scale_fill_manual(name = '', values = cols_var) 
    }
    
  })
  
  p <- wrap_plots(pca_plots) +  plot_annotation(title = plot_title) + plot_layout(guides = 'collect')  & theme(plot.title = element_text(hjust = 0.5))
  
}

###################################################################################################################################################

## add quantification information for variance modeling
quant.for.variance <- function(fit, support_dt, protein_peptide_annotation, psm_vec) {
  
  # PSM level
  if(!is.null(psm_vec)) {
    # PSM level
    support_dt$min_PSM <- psm_vec[support_dt$peptide]
    psm_count <- support_dt %>% group_by(id) %>%
      summarise(s = sum(min_PSM))
    
    # Add 1 as pseudo-count
    res <- psm_count$s[match(rownames(fit$coefficients), psm_count$id)] + 1
    return(res)
  } else {
    # Peptide level
    pep_per_ioi <- table(support_dt$id)
    pep_per_ioi <- pep_per_ioi[rownames(fit$coefficients)]
    res <- as.numeric(pep_per_ioi)
    return(res)
  } 
}

###################################################################################################################################################
spectraCounteBayesGAM <-function(fit,coef_col) {
  
  ################################################
  #  The function was adapted from:
  #  Function for IBMT (Intensity-based Moderated 
  #  T-statistic) Written by Maureen Sartor
  #  University of Cincinnati, 2006
  ################################################
  logVAR<-log(fit$sigma^2)
  df<-fit$df.residual
  numgenes<-length(logVAR[df>0])
  df[df==0]<-NA
  eg<-logVAR-digamma(df/2)+log(df/2)
  names(fit$count) <- rownames(fit$coefficients)
  output<-fit
  output$fit.method <- 'gam'
  
  x<-log2(fit$count)
  gam.model <- gam(logVAR ~ s(x,  bs = "cs"), method="REML")
  y.pred <- fitted(gam.model)
  output$model <- gam.model
  
  egpred<-y.pred-digamma(df/2)+log(df/2)
  
  myfct<- (eg-egpred)^2 - trigamma(df/2)
  
  mean.myfct<-mean(myfct,na.rm=TRUE)
  
  priordf<-vector()
  testd0<-vector()
  
  for (i in seq(1,numgenes*10)) {
    testd0[i]<-i/10
    priordf[i]= abs(mean.myfct-trigamma(testd0[i]/2))
    if (i>2) {
      if (priordf[i-2]<priordf[i-1]) { break }
    }
  }
  
  d0<-testd0[match(min(priordf),priordf)] # prior degree found
  
  s02<-exp(egpred + digamma(d0/2) - log(d0/2)) # calculate prior variance
  
  post.var <- (d0*s02 + df*fit$sigma^2)/(d0+df)
  post.df <- d0+df
  
  # fit
  # # sca.t and scc.p stands for spectra count adjusted t and p values.  
  sca.t<-as.matrix(fit$coefficients[,coef_col]/(fit$stdev.unscaled[,coef_col]
                                                *sqrt(post.var)))
  sca.p<-as.matrix(2*pt(abs(sca.t),post.df,lower.tail = FALSE))
  
  output$sca.t<-sca.t
  output$sca.p<-sca.p
  output$sca.postvar<-post.var
  output$sca.priorvar<-s02
  output$sca.dfprior<-d0
  
  return (output)
}

outputResultGAM <-function(fit,coef_col){
  results.table = limma::topTable(fit,coef = coef_col,n= Inf)
  
  results.table$gene = rownames(results.table)
  results.table$count = fit$count[results.table$gene]
  
  results.table$sca.t = fit$sca.t[results.table$gene, ]
  results.table$sca.P.Value = fit$sca.p[results.table$gene, ]
  results.table$sca.adj.pval = p.adjust(results.table$sca.P.Value,
                                        method = "BH")
  results.table = results.table[order(results.table$sca.P.Value), ]
}


VarianceBoxplotGAM <- function (fit, n=20, xlab="count",
                                ylab = "log(Variance)", main=""){
  x <- fit$count
  y <- fit$sigma^2
  
  df.temp <- data.frame(pep_count =x, variance = y )
  df.temp.filter <- df.temp[df.temp$pep_count<=n,]
  y.pred <- predict(fit$model,data.frame(x=log2(seq(1,n))))
  
  boxplot(log(variance)~pep_count,df.temp.filter, xlab=xlab,
          ylab = ylab, main=main)
  lines(seq(1,n),y.pred,col='red',lwd=3)
  
}

VarianceScatterplotGAM <- function(fit, xlab="log2(count)",  ylab = "log(Variance)", main=""){
  x <- fit$count
  y <- fit$sigma^2
  
  y.pred <- fitted(fit$model)
  
  plot(log2(x),log(y),ylab=ylab,xlab=xlab,main= main)
  k = order(x)
  lines(log2(x[k]),y.pred[k],col='red',lwd=3)
}

ResidualplotGAM <- function (fit, xlab="log2(count)",
                             ylab="Variance(fitted - observed)", main=""){
  x <- fit$count
  y <- residuals(fit$model)
  
  plot(log2(x),y, pch=20, cex=0.5,ylab=ylab,xlab=xlab,main= main)
}

###################################################################################################################################################
## Proteoform ordering based on peptide support

# Ordering of proteoforms
proteoform.order <- function(support_case) {
  
  res <- support_case %>% group_by(ioi, id) %>%
    mutate(count_peptides = n()) %>%
    group_by(ioi) %>%
    arrange(desc(count_peptides), .by_group = TRUE) %>%
    distinct(id, .keep_all= TRUE) %>%
    mutate(o = dense_rank(-count_peptides)) %>%
    dplyr::select(c(ioi, id, o))
  
  return(res)
}

###################################################################################################################################################
## q-value calculation
q.value.calculation <- function(support_case, DEpMS.results_null, DEpMS.results, protein_peptide_annotation,design_var, plot_density) {
  
  
  # Add standard deviation weight to the case distribution
  DEpMS.results$weight <- support_case$weight[match(DEpMS.results$gene, support_case$id)]
  DEpMS.results$sca.t_var <- DEpMS.results$sca.t *  DEpMS.results$weight
  DEpMS.results$ioi <- sapply(strsplit( DEpMS.results$gene, '_'), function(i) i[[length(i)]])
  
  
  ## Penalize multiple ids
  # Peptide level
  pep_per_ioi <- table(support_case$id)
  DEpMS.results$n_peptides <- as.numeric(pep_per_ioi[DEpMS.results$gene])
  
  DEpMS.results <- DEpMS.results %>% 
    group_by(ioi) %>%
    arrange(desc(n_peptides), desc(abs(sca.t_var)), .by_group = TRUE) %>%
    mutate(o = order(n_peptides, abs(sca.t_var), decreasing = TRUE))
  
  # Proteoform order
  proteoform_ord <- proteoform.order(support_case)
  DEpMS.results$pep_ord <- proteoform_ord$o[match(DEpMS.results$gene, proteoform_ord$id)]
  starting_pep_ord <- DEpMS.results$pep_ord[DEpMS.results$o == 1]
  names(starting_pep_ord) <- DEpMS.results$ioi[DEpMS.results$o == 1]
  
  DEpMS.results$pep_ord <- starting_pep_ord[DEpMS.results$ioi]
  DEpMS.results$o_cor <- DEpMS.results$o +  DEpMS.results$pep_ord -1
  
  
  DEpMS.results$sca.t_weighted <- sapply(1:nrow(DEpMS.results), function(i) DEpMS.results$sca.t_var[i] * 
                                           min(1, DEpMS.results$n_peptides[i]/DEpMS.results$o_cor[i]) )
  DEpMS.results <- DEpMS.results[order(abs(DEpMS.results$sca.t_weighted), decreasing = TRUE), ]
  
  
  # Calculate q-values
  t_null <- abs(DEpMS.results_null$sca.t)
  t_values <- abs(DEpMS.results$sca.t_weighted)
  
  
  t_steps <- sort(abs(t_values), decreasing = TRUE)
  e_fdr <- sapply(t_steps, function(i) {
    n_null <- (sum(t_null >= i)/length(t_null))/ ((sum(t_values >= i)/length(t_values)) )
  })
  
  q_values <- c()
  q_values[1] <-  min(e_fdr)
  q_values[2:nrow(DEpMS.results)] <- sapply(2:nrow(DEpMS.results), function(i) {
    min(e_fdr[-c(1:(i-1))])
  })
  
  DEpMS.results$q_value <- q_values
  
  DEpMS.results$design_var <- design_var
  
  
  if(plot_density) {
    
    density_plot <- t.stat.density.plot(t_null, t_values)
    
    pdf(paste0(path_figures, 'density_plot_', design_var, '.pdf'), width = 5, height = 4.5)
    plot(density_plot)
    dev.off()
  }
  
  
  # Re-order columns
  colstokeep <- c('gene', "ioi", 'logFC', 'AveExpr', 'B', 't', 'sca.t', 'sca.t_weighted', "sca.P.Value" , "sca.adj.pval", 
                  'q_value', 'count', 'n_peptides', 'design_var')
  res <- DEpMS.results[, colstokeep]
  colnames(res)[1] <- 'proteoform'
  return(res)
  
  
}

###################################################################################################################################################
## T-statistic density plot
t.stat.density.plot <- function(t_null, t_values) {
  
  distr_plot <- data.frame('t_stat' = c(t_null, t_values),
                           'hypothesis'= c(rep('background', length(t_null)),
                                           rep('case', length(t_values))))

  
  label_background <- paste0('n[background] == ', sum(distr_plot$hypothesis == 'background'))
  label_case <- paste0('n[case] == ', sum(distr_plot$hypothesis == 'case'))
  
  density_plot <- ggplot(distr_plot, aes(x = t_stat, fill = hypothesis)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(name = 'Hypothesis', values = c('case' = 'indianred3', 'background' = 'steelblue3'),
                      breaks = c('background', 'case')) +
    scale_x_continuous(expand = c(0,0,0.1,0.1)) +
    scale_y_continuous(expand = c(0,0,0.1,0.1)) +
    labs(x = 'weighted t-stat', y = 'Density') +
    annotate(geom = 'text', x = Inf, y = Inf, label =  bquote(.(label_background)), hjust = 1, vjust = 1.2, col= 'red', parse = TRUE) +
    annotate(geom = 'text', x = Inf, y = Inf, label =   bquote(.(label_case)),  hjust = 1, vjust = 3.2, col= 'red', parse = TRUE) +
    theme_classic()
  
  
}
###################################################################################################################################################
# DEpMS Volcano plot
 depms.volcano.plot <- function(depms.res, design_var, thres_qvalue, thres_fc) {
  
   depms.res[['design_var']] <-  design_var
   depms.res$thres_q_logical <-  depms.res$q_value <= thres_qvalue
   depms.res$thres_fc_logical <- abs( depms.res$logFC) >= thres_fc
  
   depms.res$sig <- depms.res$thres_q_logical &  depms.res$thres_fc_logical
  
  plot_text <- depms.res %>% group_by(design_var) %>%
    summarise(s = sum(sig))
  
  
  
  res <- ggplot(depms.res, aes(x = logFC, y = -log10(sca.P.Value), size = log2(n_peptides))) +
    geom_point(col = 'grey', alpha = 0.5) +
    geom_point(pch = 21, data = depms.res[depms.res$sig, ],fill ='red') +
    scale_y_continuous(expand = c(0, 0, 0.1, 0.1)) + 
    labs(x = 'Fold change (log2)', y = 'Nominal p-value (-log10)', parse = TRUE) +
    geom_vline(xintercept = 0, col = 'black') +
    geom_text(inherit.aes = FALSE, data = plot_text, mapping = aes(x = Inf, y = Inf, label = paste0('n = ', s)), col = 'red', vjust = 1.2, hjust = 1.2) +
    facet_grid(.~design_var) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(size = guide_legend(title = '# Peptides - log2'))
  
}

###################################################################################################################################################
## Summarize proteoforms
## Summarize proteoforms
summarize.proteoforms <-  function(depms.res, collapsed_profiles, protein_peptide_annotation) {
  
  iois <- unique(depms.res$ioi)
  
  pb  <- txtProgressBar(1, length(iois), style=3)
  
  proteoform_annotation <- do.call(rbind, lapply(1:length(iois), function(i) {
    
    setTxtProgressBar(pb, i)
    
    ioi <- iois[i]
    depms.res.ioi <- depms.res[depms.res$ioi %in% ioi, ]
    sig_prots <- gsub('P', '', sapply(strsplit(depms.res.ioi$proteoform, '_'), function(i) i[[1]]))
    
    dt_dis <- collapsed_profiles$support_case
    dt_dis <- dt_dis[dt_dis$membership %in% sig_prots, ] 
    
    
    if(all(grepl('_mod_', depms.res.ioi$proteoform))) {
      
      dt_dis <- dt_dis[grepl(paste0('^', ioi, '$'), dt_dis$ioi), ]
      dt_dis <- dt_dis[, c('peptide', 'ioi','membership', 'id')]
      
      # Use minimun q-value across proteoforms and comparisons
      q_value_fc_proteoform <- depms.res.ioi %>% group_by(proteoform, design_var) %>% arrange(q_value, desc(abs(logFC)))
      dt_dis[, c('q_value_proteoform', 'fc_proteoform')] <-  q_value_fc_proteoform[match(dt_dis$id, q_value_fc_proteoform$proteoform),
                                                                                   c('q_value', 'logFC')]
      
      # Minimum proteoform q-value
      dt_dis$q_value_protein <- min(dt_dis$q_value_proteoform, na.rm = TRUE)
      
      # Fold-change corresponding to the proteoform with the minimum q-value
      dt_dis$fc_protein <- dt_dis$fc_proteoform[which.min(dt_dis$q_value_proteoform)]
      
      rownames(dt_dis) <- NULL
      
      return(dt_dis)
      
    } else {
      
      dt_sim <- collapsed_profiles$support_stp
      
      dt_sim <- dt_sim[dt_sim$ioi == ioi, ]
      dt_sim  <- dt_sim[,  c('peptide', 'ioi','membership', 'id')]
      
      dt_dis <- dt_dis[dt_dis$ioi == ioi, ]
      dt_dis  <- dt_dis[,  c('peptide', 'ioi','membership', 'id')]
      
      
      dt_res <- rbind(dt_dis, dt_sim)
      
      q_value_fc_proteoform <- depms.res.ioi %>% group_by(proteoform, design_var) %>% arrange(q_value,  desc(abs(logFC)))
      
      # Minimum q-value across variables
      dt_res[, c('q_value_proteoform', 'fc_proteoform')] <-  q_value_fc_proteoform[match(dt_res$id, q_value_fc_proteoform$proteoform),
                                                                                   c('q_value', 'logFC')]
      
      # Minimum proteoform q-value
      dt_res$q_value_protein <- min(dt_res$q_value_proteoform, na.rm = TRUE)
      
      # Fold-change corresponding to the proteoform with the minimum q-value
      min_value <- min(dt_res$q_value_proteoform, na.rm = TRUE)
      max_fc <- dt_res$fc_proteoform[which(dt_res$q_value_proteoform == min_value)]
      dt_res$fc_protein <- max_fc[which.max(abs(max_fc))]
      
      rownames(dt_res) <- NULL
      
      return(dt_res)
    }
  }))
  
  return(proteoform_annotation)
}
###################################################################################################################################################
# Peptide support
peptide.support.plot <- function(summarized_proteoforms, thres_qvalue, thres_fc) {
  
  
  peptides_per_proteoform <- summarized_proteoforms %>% group_by(ioi, id) %>%
    mutate(n = n())
  
  
  
  peptides_per_proteoform$n_group <- cut2(peptides_per_proteoform$n, m = nrow(peptides_per_proteoform)/10)
  peptides_per_proteoform$sig <- ifelse(peptides_per_proteoform$q_value_proteoform <= thres_qvalue & abs(peptides_per_proteoform$fc_proteoform) >= thres_fc, 'dtp_profiles_significant',
                                        'dtp_profiles_non_significant')
  peptides_per_proteoform$sig[is.na(peptides_per_proteoform$sig)] <- 'stp_profiles'
  peptides_per_proteoform$sig <- factor(peptides_per_proteoform$sig, levels = c('stp_profiles',
                                                                                'dtp_profiles_non_significant',
                                                                                'dtp_profiles_significant'))
  
  peptides_per_proteoform <- peptides_per_proteoform[!duplicated(peptides_per_proteoform$id),]
  plot_text <- as.data.frame(table(peptides_per_proteoform$sig))
  
  plot_text$label <-  apply(plot_text, 1, function(i)  {
    toconvert <- paste0('n[',i[1],']', ' == ', i[2])
  })
  
  
  # Calculate percentages
  summarized_proteoforms_perc <- peptides_per_proteoform %>% group_by(sig, n_group) %>% 
    summarise(n = n()) %>%
    group_by(sig) %>%
    mutate(perc = 100 *n/sum(n))
  
  
  # Remove NAs
  res <- ggplot(summarized_proteoforms_perc, aes(x = n_group, y = perc, fill = sig)) +
    geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single"), col = 'black') +
    scale_fill_manual(values = c('stp_profiles' = 'grey', 'dtp_profiles_non_significant' = 'deepskyblue', 
                                 'dtp_profiles_significant' = 'red')) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0.1)) + 
    geom_text(inherit.aes = FALSE, data = plot_text, mapping = aes(x = Inf, y = Inf, label = bquote(.(label))), col = 'black', vjust = c(1.2, 3.2, 5.2), hjust = 1, parse = TRUE) +
    labs(x = '# peptides', y = 'Percentage (%)') +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) + 
    guides(fill = guide_legend(title = ''))
  
  
}

###################################################################################################################################################
## Quality plot
depms.feature.plot <- function(summarized_proteoforms, thres_qvalue, thres_fc, quality_vec, var_id) {
  
  quality_vec <- quality_vec[summarized_proteoforms$peptide]
  summarized_proteoforms$quality_vec <- cut2(quality_vec, m = length(quality_vec)/10)
  summarized_proteoforms$sig <- ifelse(summarized_proteoforms$q_value_proteoform < thres_qvalue & abs(summarized_proteoforms$fc_proteoform) >= thres_fc, 'dtp_significant',
                                       'dtp_non_significant')
  summarized_proteoforms$sig[is.na( summarized_proteoforms$sig)] <- 'stp'
  summarized_proteoforms$sig <- factor(summarized_proteoforms$sig, levels = c('stp',
                                                                              'dtp_non_significant',
                                                                              'dtp_significant'))
  
  plot_text <- as.data.frame(table(summarized_proteoforms$sig))
  
  plot_text$label <-  apply(plot_text, 1, function(i)  {
    toconvert <- paste0('n[',i[1],']', ' == ', i[2])
  })
  
  
  # Calculate percentages
  summarized_proteoforms_perc <- summarized_proteoforms %>% group_by(sig, quality_vec) %>% 
    summarise(n = n()) %>%
    group_by(sig) %>%
    mutate(perc = 100 *n/sum(n))
  
  
  # Remove NAs
  res <- ggplot(summarized_proteoforms_perc, aes(x = quality_vec, y = perc, fill = sig)) +
    geom_bar(stat ="identity", position = position_dodge2(width = 0.9, preserve = "single"), col = 'black') +
    scale_fill_manual(values = c('stp' = 'grey', 'dtp_non_significant' = 'deepskyblue', 
                                 'dtp_significant' = 'red')) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0.1)) + 
    geom_text(inherit.aes = FALSE, data = plot_text, mapping = aes(x = Inf, y = Inf, label = bquote(.(label))), col = 'black', vjust = c(1.2, 3.2, 5.2), hjust = 1, parse = TRUE) +
    labs(x = var_id, y = 'Percentage (%)') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(title = ''))
  
}

###################################################################################################################################################
peptide.feature.plot <- function(ioi, summarized_proteoforms, protein_peptide_annotation, thres_qvalue, thres_fc, quality_vec, var_id, proteoform_col) {
  
  
  summarized_proteoforms_ioi <- summarized_proteoforms[summarized_proteoforms$ioi %in% ioi, ]
  summarized_proteoforms_ioi_sig <- summarized_proteoforms_ioi[which(summarized_proteoforms_ioi$q_value_proteoform <= thres_qvalue &
                                                                       abs(summarized_proteoforms_ioi$fc_proteoform) >= thres_fc), ]
  
  # Similar-to-protein peptides
  stp_id_peps <- summarized_proteoforms_ioi$peptide[grep('sim', summarized_proteoforms_ioi$id)]
  
  # if no significant putative proteoforms, plot all peptides 
  if(nrow(summarized_proteoforms_ioi_sig) == 0) {
    
    
    summarized_proteoforms_ioi <- summarized_proteoforms_ioi[!summarized_proteoforms_ioi$peptide %in% stp_id_peps, ]
    
    peptide_seqs <- c(summarized_proteoforms_ioi$peptide, stp_id_peps)
    peptide_quant_ioi <- peptide_quant[peptide_seqs, ]
    
    row_split <- c(paste0('P', summarized_proteoforms_ioi$membership), rep('similar-to-protein', length(stp_id_peps)))
    row_split_cols <- rep('grey', length(unique(row_split)))
    names(row_split_cols) <- unique(row_split)
    row_split_cols[paste0('P', unique(summarized_proteoforms_ioi$membership))]  <- proteoform_col[order(unique(summarized_proteoforms_ioi$membership))]
    
    
    plot_dt <- data.frame('peptide' = peptide_seqs, 
                          'membership' = row_split, 
                          'variable' = quality_vec[peptide_seqs])
    
  } else {
    
    sig_id_peps <- summarized_proteoforms_ioi[which(summarized_proteoforms_ioi$q_value_proteoform <= thres_qvalue & 
                                                      abs(summarized_proteoforms_ioi$fc_proteoform) >= thres_fc),]
    
    row_split <- c(paste0('P', sig_id_peps$membership), rep('similar-to-protein', length(stp_id_peps)))
    
    row_split_cols <- rep('grey', length(unique(row_split)))
    names(row_split_cols) <- unique(row_split)
    row_split_cols[paste0('P', unique(sig_id_peps$membership))]  <- proteoform_col[order(unique(sig_id_peps$membership))]
    
    
    row_split <- factor(row_split, levels = unique(row_split))
    
    
    peptide_seqs <- c(sig_id_peps$peptide, stp_id_peps)
    peptide_quant_ioi <- peptide_quant[peptide_seqs, ]
    # peptide_var <- apply(peptide_quant_ioi, 1 , var, na.rm = TRUE)
    # peptide_na <- apply(peptide_quant_ioi, 1 ,function(i) sum(is.na(i)))
    
    plot_dt <- data.frame('peptide' = peptide_seqs, 
                          'membership' = row_split, 
                          'variable' = quality_vec[peptide_seqs])
  }
  
  plot_dt$peptide <- factor(plot_dt$peptide, levels = rev(plot_dt$peptide))
  p <- ggplot(plot_dt, aes(x = peptide, y = variable)) +
    geom_point(aes(col = membership), size = 2) + 
    scale_color_manual(values = row_split_cols) +
    facet_grid(membership ~ ., space = 'free', scales = 'free') + 
    labs(x = '', y = var_id, title = ioi) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.spacing = unit(0, "lines")) +
    coord_flip() +
    guides(color = 'none')
  
  return(p) 
}
###################################################################################################################################################
## Proteoform composition
proteoform.composition.plot <- function(depms.res.all, summarized_proteoforms, thres_qvalue, thres_fc) {
  
  proteins_sig_dt <- summarized_proteoforms[which(summarized_proteoforms$q_value_protein <= thres_qvalue &
                                                    abs(summarized_proteoforms$fc_protein) >= thres_fc), ]
  idx_proteoforms_sig <- which(proteins_sig_dt$q_value_proteoform <= thres_qvalue &  abs(proteins_sig_dt$fc_proteoform) >= thres_fc)
  
  proteoforms_sig <- unique(proteins_sig_dt$id[idx_proteoforms_sig])
  proteins_sig <- unique(proteins_sig_dt$ioi)
  
  n_proteoforms <- length(proteoforms_sig)
  n_proteins <- length(proteins_sig)
  
  text_labels <- c(paste0('n[proteoforms] == ', n_proteoforms),
                   paste0('n[proteins] == ', n_proteins))
  
  depms.res.all_sig <- depms.res.all[depms.res.all$proteoform %in% proteoforms_sig & depms.res.all$q_value <= thres_qvalue & abs(depms.res.all$logFC) >= thres_fc, ]
  freq_proteoform <- as.data.frame(table(depms.res.all_sig$design_var))
  
  p_freq_proteoform <- ggplot(freq_proteoform, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = 'identity', col = 'black') +
    geom_text(aes(label = Freq), vjust = -1) +
    annotate(geom = 'text', x = Inf, y = Inf, label = bquote(.(text_labels)), 
             vjust = c(1,3), hjust = c(1,1),  parse = TRUE) + 
    scale_y_continuous(expand = c(0,0, 1, 1)) +
    labs(x = '', y = paste0('# significant proteoforms (q-value <= ', thres_qvalue,', |log2FC| >= ', round(thres_fc, 3), ')')) +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    guides(fill = 'none')
  
  
  # Keep outlying, similar, and modular proteoforms
  proteins_proteoforms_sim_sig_dt <- proteins_sig_dt[c(which(proteins_sig_dt$q_value_proteoform <= thres_qvalue &
                                                               abs(proteins_sig_dt$fc_proteoform) >= thres_fc),
                                                       grep(paste(c('sim', 'mod'), collapse = '|'), proteins_sig_dt$id)), ]
  
  proteins_proteoforms_sim_sig_dt <- proteins_proteoforms_sim_sig_dt[!duplicated(proteins_proteoforms_sim_sig_dt$id), ]
  proteoforms_per_protein <- table(proteins_proteoforms_sim_sig_dt$ioi)
  
  # Exclude single-case proteoforms
  proteoforms_per_protein <- proteoforms_per_protein[proteoforms_per_protein != 1]
  proteoforms_per_protein <- as.data.frame(table(proteoforms_per_protein))
  proteoforms_per_protein$perc <- proteoforms_per_protein$Freq/sum(proteoforms_per_protein$Freq)
  
  proteoforms_per_protein <- proteoforms_per_protein %>% 
    mutate(perc = Freq/sum(Freq),
           csum = rev(cumsum(rev(perc))), 
           pos = perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), perc/2, pos))
  
  p_proteoforms_per_protein <- ggplot(proteoforms_per_protein, aes(x = "" , y = perc, fill = proteoforms_per_protein)) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Greens") +
    geom_label_repel(data = proteoforms_per_protein,
                     aes(y = pos, label = paste0(Freq, ' (', 100*round(perc, 3), "%", ')')),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    theme_void() +
    guides(fill = guide_legend(title = '# proteoforms per protein'))
  
  
  
  # Proteoform uniqueness per comparison
  freq_unique <- table(table(depms.res.all_sig$proteoform))
  freq_unique <- as.data.frame(freq_unique)
  
  freq_unique <- freq_unique %>% 
    mutate(perc = Freq/sum(Freq),
           csum = rev(cumsum(rev(perc))), 
           pos = perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), perc/2, pos))
  
  
  p_freq_unique <- ggplot(freq_unique, aes(x = "" , y = perc, fill = Var1)) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Greens") +
    geom_label_repel(data = freq_unique,
                     aes(y = pos, label = paste0(Freq, ' (', 100*round(perc, 3), "%", ')')),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    theme_void() +
    guides(fill = guide_legend(title = '# Comparisons'))
  
  p_proteoform_description <- p_freq_proteoform/p_proteoforms_per_protein/p_freq_unique
  
}
###################################################################################################################################################
## Heatmap
heatmap.plot <- function(ioi, protein_quant, peptide_quant, summarized_proteoforms, thres_qvalue, thres_fc, proteoform_col, sort_samples_idx, heatmap_title, top_annotation = ha) {
  
  
  # Protein profile
  protein_quant_ioi <- protein_quant[ioi, ]
  
  # Significant putative proteoforms
  summarized_proteoforms_ioi <- summarized_proteoforms[summarized_proteoforms$ioi %in% ioi, ]
  summarized_proteoforms_ioi_sig <- summarized_proteoforms_ioi[which(summarized_proteoforms_ioi$q_value_proteoform <= thres_qvalue 
                                                                     & abs(summarized_proteoforms_ioi$fc_proteoform) >= thres_fc), ]
  
  # Similar-to-protein peptides
  stp_id_peps <- summarized_proteoforms_ioi$peptide[grep('sim', summarized_proteoforms_ioi$id)]
  
  # if modular proteoform, plot all peptides
  if(all(grepl('_mod_', summarized_proteoforms_ioi$id))) {
    
    summarized_proteoforms_ioi <- summarized_proteoforms_ioi[!summarized_proteoforms_ioi$peptide %in% summarized_proteoforms_ioi_sig$peptide, ]
    
    peptide_quant_ioi <- peptide_quant[c(summarized_proteoforms_ioi_sig$peptide, summarized_proteoforms_ioi$peptide), ]
    
    if(nrow(summarized_proteoforms_ioi) == 0) {
      row_split <- c(paste0('P', summarized_proteoforms_ioi_sig$membership), 'protein')
    } else {
      row_split <- c(paste0('P', summarized_proteoforms_ioi_sig$membership), paste0('P', summarized_proteoforms_ioi$membership), 'protein')
      
    }
    
    row_split_cols <- rep('grey', length(unique(row_split)))
    names(row_split_cols) <- unique(row_split)
    
    row_split_cols[paste0('P', unique(summarized_proteoforms_ioi_sig$membership))]  <- proteoform_col[order(unique(summarized_proteoforms_ioi_sig$membership))]
    row_split_cols['protein']  <- 'black'
    
    if(nrow(summarized_proteoforms_ioi_sig) > 0 ) {
      row_split_cols[paste0('P', unique(summarized_proteoforms_ioi_sig$membership))]  <- proteoform_col[order(unique(summarized_proteoforms_ioi_sig$membership))]
    }
    row_split_cols['protein']  <- 'black'
    
    
    ra <- rowAnnotation(df = data.frame('Profile' = row_split),
                        col = list('Profile' = row_split_cols),
                        gp = gpar(col = "black"),
                        show_annotation_name = FALSE)
    
    
  } else {
    
    # if no significant putative proteoforms, plot all peptides 
    if(nrow(summarized_proteoforms_ioi_sig) == 0) {
      
      summarized_proteoforms_ioi <- summarized_proteoforms_ioi[!summarized_proteoforms_ioi$peptide %in% stp_id_peps, ]
      
      peptide_seqs <- c(summarized_proteoforms_ioi$peptide, stp_id_peps)
      peptide_quant_ioi <- peptide_quant[peptide_seqs, ]
      
      row_split <- c(paste0('P', summarized_proteoforms_ioi$membership), rep('similar-to-protein', length(stp_id_peps)), 'protein')
      row_split_cols <- rep('grey', length(unique(row_split)))
      names(row_split_cols) <- unique(row_split)
      row_split_cols[paste0('P', unique(summarized_proteoforms_ioi$membership))]  <- proteoform_col[order(unique(summarized_proteoforms_ioi$membership))]
      row_split_cols['similar-to-protein'] = 'grey'
      row_split_cols['protein']  <- 'black'
      
      
      ra <- rowAnnotation(df = data.frame('Profile' = row_split),
                          col = list('Profile' = row_split_cols),
                          gp = gpar(col = "black"),
                          show_annotation_name = FALSE)
    } else {
      
      # if significant putative proteoforms exist, plot proteoform-specific peptides
      peptide_quant_ioi <- peptide_quant[c(summarized_proteoforms_ioi_sig$peptide, stp_id_peps), ]
      
      row_split <- c(paste0('P', summarized_proteoforms_ioi_sig$membership), rep('similar-to-protein', length(stp_id_peps)), 'protein')
      
      row_split_cols <- rep('grey', length(unique(row_split)))
      names(row_split_cols) <- unique(row_split)
      row_split_cols[paste0('P', unique(summarized_proteoforms_ioi_sig$membership))]  <- proteoform_col[order(unique(summarized_proteoforms_ioi_sig$membership))]
      row_split_cols['protein']  <- 'black'
      
      
      
      
      ra <- rowAnnotation(df = data.frame('Profile' = row_split),
                          col = list('Profile' = row_split_cols),
                          gp = gpar(col = "black"),
                          show_annotation_name = FALSE)
    }
  }
  
  row_split <- factor(row_split, levels = unique(row_split))
  
  ht <- Heatmap(
    as.matrix(rbind(peptide_quant_ioi[, sort_samples_idx],
                    protein_quant_ioi[, sort_samples_idx])),
    row_title_rot = 0,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(title = 'Log2 ratio'),
    row_split = row_split,
    row_title = NULL,
    cluster_column_slices = FALSE,
    column_title = heatmap_title,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = 'left',
    right_annotation = ra,  
    top_annotation = top_annotation,
    use_raster = FALSE,
    row_gap = unit(0, "mm"), 
    column_gap = unit(0, "mm"),
    border = TRUE)
  
  
  return(ht)
}
####################################################################################################################################
## Usage barplot
usage.barplot <- function(ioi, collapsed_quant_case, summarized_proteoforms, thres_qvalue, thres_fc, annotation_vec, annotation_cols, sort_samples_idx,  plot_title) {
  
  proteoform_sig_ioi <- summarized_proteoforms[summarized_proteoforms$ioi %in% ioi, ]
    
  proteoform_sig_ioi <-  proteoform_sig_ioi[which(proteoform_sig_ioi$q_value_proteoform <= thres_qvalue & 
                                                    abs(proteoform_sig_ioi$fc_proteoform) >= thres_fc), ]
  proteoform_sig_ioi <-  unique(proteoform_sig_ioi$id)
  
  
  dt <- reshape2::melt(collapsed_quant_case[proteoform_sig_ioi, , drop = FALSE])
  
  # Re-order samples
  dt$Var2 <-factor(dt$Var2, levels = levels(dt$Var2)[sort_samples_idx])
  
  if(!is.null(annotation_vec)) {
    
    dt$group <- annotation_vec[as.character(dt$Var2)]
    
    p <- ggplot(dt, aes(x = Var2, y = value, fill = group)) +
      geom_bar(stat = 'identity')  +
      geom_hline(yintercept = 0) +
      scale_fill_manual(name = '', values = annotation_cols) +
      facet_wrap(.~Var1, ncol = 1) +
      labs(x =  'samples', y = 'proteoform usage', title = plot_title) +
      theme_classic() +
      theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            plot.title = element_text(hjust = 0.5)) 
  } else {
    
    p <- ggplot(dt, aes(x = Var2, y = value)) +
      geom_hline(yintercept = 0) +
      geom_bar(stat = 'identity')  +
      facet_wrap(.~Var1, ncol = 1) +
      labs(x =  'samples', y = 'proteoform usage', title = plot_title) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            axis.ticks.x = element_blank(),
            plot.title = element_text(hjust = 0.5)) 
  }
  
  return(p)
  
}
####################################################################################################################################
## Peptide gene coordinates
peptide.coordinates <- function(peptide, protein_peptide_annotation, ioi_id, db_match_id, ioi_to_db_match, db_match_fasta, gtf_file, gtf_ioi_id, gtf_db_match_id, positional) {
  
  # print(peptide)
  idx <- match(peptide, protein_peptide_annotation$Peptide_sequence)
  ioi <- protein_peptide_annotation[[ioi_id]][idx]
  db_match_ioi <- protein_peptide_annotation[[db_match_id]][idx]
  
  db_match_ioi <- unlist(strsplit(db_match_ioi,';'))
  
  exons <- gtf_file[values(gtf_file)[, gtf_ioi_id] %in% ioi & gtf_file$type == 'exon', ]
  cds <- gtf_file[values(gtf_file)[, gtf_db_match_id] %in% db_match_ioi & gtf_file$type == 'CDS', ]
  
  if(length(cds) == 0) {
    paste0('no cds found for ', db_match_ioi, '-', ioi)
    return(NULL)
  }
  
  cds <- cds[values(cds)[, gtf_ioi_id] == ioi, ]
  db_match_ioi <- unique(values(cds)[, gtf_db_match_id])[1]
  cds <- gtf_file[values(gtf_file)[, gtf_db_match_id] %in% db_match_ioi & gtf_file$type == 'CDS', ]
  
  cds_length <- width(cds)
  sum_cds_length <-  c(0, cumsum(cds_length)[-length(cds)])
  start_cds <- sum_cds_length + 1
  end_cds <- start_cds + cds_length - 1
  
  new_cds <- GRanges(
    seqnames = seqnames(cds),
    ranges = IRanges(start = start_cds, end = end_cds),
    strand = strand(cds))
  
  peptide_seq <- gsub('[^A-Z]', '', peptide)
  protein_seq <- db_match_fasta[db_match_ioi]
  
  coords <- as.data.frame(unlist(vmatchPattern(peptide_seq, protein_seq)))
  
  
  if(all(is.na(coords))) {
    print(paste0('peptide ', peptide, ' not mapping to the protein sequence'))
    return(NULL)
  }
  
  if(!is.null(positional)) {
    
    coords[1,'start'] <- coords[1,'start'] + positional - 1
    coords[1,'end'] <- coords[1,'start'] 
    
  }
  
  
  trans_coords_start <- (coords[1,'start'] - 1) * 3  + 1
  trans_coords_end <- coords[1,'end'] * 3
  
  trans_coords <-  GRanges(
    seqnames = runValue(seqnames(cds)),
    ranges = IRanges(start = trans_coords_start, end = trans_coords_end),
    strand = runValue(strand(cds))
  )
  
  overlapping_regions <- findOverlaps(trans_coords, new_cds)
  cds_idx <- to(overlapping_regions)
  new_cds <- new_cds[cds_idx]
  old_cds <- cds[cds_idx] 
  
  if(length(cds_idx) == 1) {
    if(runValue(strand(trans_coords)) == '-') {
      peptide_start <- end(old_cds) - (end(trans_coords) - start(new_cds))
      peptide_end <- end(old_cds)[1] - (start(trans_coords) - start(new_cds)[1])
      chr_coord <- paste(peptide_start, peptide_end, sep = '-')
    } else {
      peptide_start <- start(old_cds) + (start(trans_coords) - start(new_cds))
      peptide_end <- start(old_cds) + (end(trans_coords) - start(new_cds))
      chr_coord <- paste(peptide_start, peptide_end, sep = '-')
    }
  } else {
    if(runValue(strand(trans_coords)) == '-') {
      peptide_start1 <- start(old_cds)[1]
      peptide_start2 <- start(old_cds)[1] + (end(new_cds)[1] - start(trans_coords)) 
      chr_coord1 <- paste(peptide_start1, peptide_start2, sep= '-')
      
      peptide_end1 <-  end(old_cds)[length(old_cds)] -  (end(trans_coords) - start(new_cds)[length(old_cds)]) 
      peptide_end2 <- end(old_cds)[length(old_cds)] 
      chr_coord2 <- paste(peptide_end1, peptide_end2, sep= '-')
      
      if(length(cds_idx) > 2) {
        inter_chr_coord <- old_cds[-c(1, length(old_cds))]
        inter_chr_coord <- paste(start(inter_chr_coord), end(inter_chr_coord), sep = '-')
        chr_coord <- paste(chr_coord1, inter_chr_coord, chr_coord2, sep = ':')
      } else {
        chr_coord <- paste(chr_coord1, chr_coord2, sep = ':')
      }
      
    } else {
      peptide_start1 <- start(old_cds)[1] + (start(trans_coords) - start(new_cds)[1])
      peptide_start2 <- end(old_cds)[1]
      chr_coord1 <- paste(peptide_start1, peptide_start2, sep= '-')
      
      peptide_end1 <- start(old_cds)[length(old_cds)]
      peptide_end2 <- start(old_cds)[length(old_cds)] + (end(trans_coords) - start(new_cds)[length(old_cds)])
      chr_coord2 <- paste(peptide_end1, peptide_end2, sep= '-')
      
      if(length(cds_idx) > 2) {
        inter_chr_coord <- old_cds[-c(1, length(old_cds))]
        inter_chr_coord <- paste(start(inter_chr_coord), end(inter_chr_coord), sep = '-')
        chr_coord <- paste(chr_coord1, inter_chr_coord, chr_coord2, sep = ':')
      } else {
        chr_coord <- paste(chr_coord1, chr_coord2, sep = ':')
      }
      
    }
  }
  
  exon_number <- unique(cds$exon_number[cds_idx])
  res_map <- data.frame(
    'peptide' = peptide,
    'ioi' = ioi,
    'chr' = unique(seqnames(old_cds)),
    'strand' = unique(strand(old_cds)),
    'chr_coords' = chr_coord
    # 'exon' = paste(exon_number,collapse = ':'),
    # 'exon_id' = paste(unique(exons$exon_id[match(exon_number,exons$exon_number)]),collapse = ':')
  )
  
  all_db_match_ids <- ioi_to_db_match[[ioi]]
  protein_seqs <- db_match_fasta[all_db_match_ids]
  all_db_match_ids <- names(unlist(vmatchPattern(peptide_seq, protein_seqs)))
  
  res_map$db_match_ids <- paste(all_db_match_ids, collapse = ',')

  return(res_map)
}
###################################################################################################################################################
## Track plot
track.plot <- function(ioi, 
                       summarized_proteoforms,
                       protein_peptide_annotation, 
                       ioi_id,
                       db_match_id,
                       ioi_to_db_match, 
                       db_match_fasta,
                       gtf_file, 
                       gtf_ioi_id, 
                       gtf_db_match_id,
                       thres_qvalue,
                       thres_fc,
                       genome_version, 
                       biomart_obj,
                       biomart_filters,
                       proteoform_col,
                       positional) {
  
  
  # Putative proteoforms
  summarized_proteoforms_ioi <- summarized_proteoforms[summarized_proteoforms$ioi %in% ioi, ]
  
  # Significant proteoforms
  summarized_proteoforms_ioi_sig <- summarized_proteoforms_ioi[which(summarized_proteoforms_ioi$q_value_proteoform <= thres_qvalue &
                                                                       abs(summarized_proteoforms_ioi$fc_proteoform) >= thres_fc), ]
  
  peptidestoplot <- summarized_proteoforms_ioi$peptide
  
  peptide_coordinates_dt <- do.call(rbind, lapply(peptidestoplot, function(i) {
    res <-  peptide.coordinates(peptide = i, 
                                protein_peptide_annotation = protein_peptide_annotation,
                                ioi_id = ioi_id,
                                gtf_db_match_id = gtf_db_match_id,
                                db_match_id <- db_match_id,
                                ioi_to_db_match = ioi_to_db_match,
                                db_match_fasta = db_match_fasta,
                                gtf_file = gtf_file,
                                gtf_ioi_id = gtf_ioi_id,
                                positional = positional)
  }))
  
  
  # Case of no peptide coordinates
  if(is.null(peptide_coordinates_dt)) {
    return(NULL)
  }
  
  chr <- unique(peptide_coordinates_dt$chr)
  # Capture patches
  if(!chr %in%  c(1:22, 'X', 'Y')) {
    print(paste0(ioi, ' mapping to patched region ', chr))
    return(NULL)
  }
  
  # Add membership, significance
  peptide_coordinates_dt$proteof_membership <- summarized_proteoforms_ioi$membership[match(peptide_coordinates_dt$peptide, summarized_proteoforms_ioi$peptide)]
  peptide_coordinates_dt$proteof_sig <- peptide_coordinates_dt$proteof_membership %in% summarized_proteoforms_ioi_sig$membership
  
  
  # Start, end of peptides
  peptide_coordinates_dt[, c('start', 'end')] <- t(sapply(1:nrow(peptide_coordinates_dt), function(i)  {
    pep_coords <- unlist(strsplit(peptide_coordinates_dt$chr_coords[i], ':'))
    pep_coords <- sapply(pep_coords, strsplit, split = '-')
    pep_start <- as.numeric(sapply(pep_coords, function(i) i[1]))
    pep_start <- min(pep_start)
    pep_end <- as.numeric(sapply(pep_coords, function(i) i[2]))
    pep_end <- max(pep_end)
    c(pep_start, pep_end)
  }))
  
  
  
  ## Order by genomic position
  if(all(peptide_coordinates_dt$strand == '-')) {
    
    idx_gen_order <- order(-peptide_coordinates_dt$end) 
  } else {
    idx_gen_order <- order(peptide_coordinates_dt$start)
  }
  peptide_coordinates_dt <- peptide_coordinates_dt[idx_gen_order, ]
  
  
  # Split by significance
  peptide_coordinates_sig_dt <- peptide_coordinates_dt[peptide_coordinates_dt$proteof_sig, ]
  
  # Color peptides
  peptide_col <- c()
  peptide_col <- rep('grey', length(unique(peptide_coordinates_dt$proteof_membership)))
  names(peptide_col) <- paste0('P',unique(peptide_coordinates_dt$proteof_membership))
  
  if(nrow(peptide_coordinates_sig_dt) > 0 ) {
    peptide_col[paste0('P', unique(peptide_coordinates_sig_dt$proteof_membership))]  <- proteoform_col[order(unique(peptide_coordinates_sig_dt$proteof_membership))]
  }
  
  if(all(grepl('_mod_', summarized_proteoforms_ioi$id))) {
    
    peptide_coordinates_rest_dt <- peptide_coordinates_dt[!peptide_coordinates_dt$proteof_sig, ]
    
  } else {
    peptide_coordinates_rest_dt <- peptide_coordinates_dt[grep('sim', peptide_coordinates_dt$proteof_membership), ]
    
  }
  
  # Create Granges object
  if(nrow(peptide_coordinates_sig_dt) == 0) {
    atrack_sig <- NULL
  } else {
    
    peptide_coordinates_sig <- GRangesList(lapply(1:nrow(peptide_coordinates_sig_dt), function(i)  {
      
      pep_chr_name <- as.character(peptide_coordinates_sig_dt$chr[i])
      pep_strand <- unique(peptide_coordinates_sig_dt$strand)
      pep_coords <- unlist(strsplit(peptide_coordinates_sig_dt$chr_coords[i], ':'))
      pep_coords <- sapply(pep_coords, strsplit, split = '-')
      pep_start <- as.numeric(sapply(pep_coords, function(i) i[[1]]))
      pep_end <- as.numeric(sapply(pep_coords, function(i) i[[2]]))
      
      
      GRanges(seqnames = pep_chr_name, 
              ranges = IRanges(start = pep_start,end =  pep_end),
              strand = pep_strand)
    }))
    
    names(peptide_coordinates_sig) <- peptide_coordinates_sig_dt$peptide
    
    # Unlist  GRangesList
    peptide_coordinates_sig <- unlist(peptide_coordinates_sig)
    
    # add peptide, id metadata
    peptide_coordinates_sig$peptide <- names(peptide_coordinates_sig)
    peptide_coordinates_sig$id <- paste0('P', peptide_coordinates_sig_dt$proteof_membership[match(names(peptide_coordinates_sig), peptide_coordinates_sig_dt$peptide)])
    
    
    atrack_sig <-  lapply(sort(unique(peptide_coordinates_sig$id)), function(j) {
      p <- peptide_coordinates_sig[peptide_coordinates_sig$id == j]
      AnnotationTrack(p, 
                      group =   p$peptide,
                      name = j,
                      background.title	= peptide_col[j],
                      col = 'black',
                      fill = 'black',
                      rotation.title = 90)
    })
  }
  
  
  if(nrow(peptide_coordinates_rest_dt) == 0) {
    atrack_rest <- NULL
  } else {
    peptide_coordinates_rest <- GRangesList(lapply(1:nrow(peptide_coordinates_rest_dt), function(i)  {
      # print(i)
      pep_chr_name <- as.character(peptide_coordinates_rest_dt$chr[i])
      pep_strand <- unique(peptide_coordinates_rest_dt$strand)
      pep_coords <- unlist(strsplit(peptide_coordinates_rest_dt$chr_coords[i], ':'))
      pep_coords <- sapply(pep_coords, strsplit, split = '-')
      pep_start <- as.numeric(sapply(pep_coords, function(i) i[[1]]))
      pep_end <- as.numeric(sapply(pep_coords, function(i) i[[2]]))
      
      
      
      GRanges(seqnames = pep_chr_name, 
              ranges = IRanges(start = pep_start,end =  pep_end),
              strand = pep_strand)
    }))
    
    names(peptide_coordinates_rest) <- peptide_coordinates_rest_dt$peptide
    
    # Unlist  GRangesList
    peptide_coordinates_rest <- unlist(peptide_coordinates_rest)
    
    # add peptide, id metadata
    peptide_coordinates_rest$peptide <- names(peptide_coordinates_rest)
    peptide_coordinates_rest$id <- paste0('P', peptide_coordinates_rest_dt$proteof_membership[match(names(peptide_coordinates_rest), peptide_coordinates_rest_dt$peptide)])
    
    
    # Remove duplicate peptides for vizualization
    atrack_rest <-  lapply(sort(unique(peptide_coordinates_rest$id)), function(j) {
      p <- peptide_coordinates_rest[peptide_coordinates_rest$id == j]
      AnnotationTrack(p, 
                      group =   p$peptide,
                      name = j,
                      background.title	= peptide_col[j],
                      col = 'black',
                      fill = 'black',
                      rotation.title = 90)
    })
  }
  
  ## Gene model track
  chr <- gsub('chr', '', as.character(unique(peptide_coordinates_dt$chr)))
  
  # itrack <- IdeogramTrack(genome = genome_version, chromosome = chr)
  gtrack <- GenomeAxisTrack()
  
  bm_list <-  list(unique(unlist(strsplit(peptide_coordinates_dt$db_match_ids, ','))))
  names(bm_list) <-  biomart_filters
  
  if(!grepl('^ENSG', ioi)) {
    
    biomart_ioi <- ioi_to_ensg$gene_id[match(ioi, ioi_to_ensg[[gtf_ioi_id]])]
  } else {
    biomart_ioi <- ioi
  }
  
  biomTrack <- BiomartGeneRegionTrack(genome = genome_version,
                                      gene = biomart_ioi,
                                      biomart = biomart_obj,
                                      filters = bm_list,
                                      # featureMap = fm,
                                      transcriptAnnotation = "transcript",
                                      name = ioi)
  chromosome(biomTrack) <- chr
  
  
  
  res <- c(
    # itrack,
           gtrack,
           biomTrack,
           atrack_sig,
           atrack_rest)
  
}

##############################################################################################################################
## uniprot tracks
createRangesUniprotKBtrack <- function(track_name) {
  
  
  if(!track_name %in% track_df$Name){
    stop(cat("Unknown track name."))
  }
  
  track <- dplyr::filter(track_df, Name %in% track_name)
  
  bed <- 1
  while (is.numeric(bed)) {
    bed <- tryCatch(
      data.table::fread(as.character(track$URL), header = FALSE, sep = "\t", 
                      quote = "", stringsAsFactors = FALSE,
                      data.table = FALSE, showProgress = FALSE), 
      error = function(e) {
        if (bed > 10) stop("cannot connect to url, too many failures", call.=FALSE)
        bed + 1
      })
  }
  
  if(is.null(bed)) {
    return(NA)
  }
  
  coords <- bed.to.coords(bed)
  
  ranges_obj <- IRanges(start = coords$start, end = coords$end)
  values(ranges_obj) <- DataFrame(uniprot =  coords$protein, id =  coords$id)
  
  return(ranges_obj)
  
}
##############################################################################################################################
bed.to.coords <- function(bed) {
  
  # Convert to ranges
  coords <- sapply(strsplit(bed$V14, ';'), function(i) i[[1]])
  track_ids <-  sapply(strsplit(bed$V14, ';'), function(i) paste(trimws(i[-1]), collapse = ';'))
  
  # Multiple sites
  if(any(grepl(',', coords))) {
    coords <- strsplit(coords, ',')
    times_to_mult <- sapply(coords, length)
    
    if(any(grepl('-', coords))) {
      coords <- as.data.frame(do.call(rbind, lapply(strsplit(unlist(coords), '-'), function(i) gsub('[^0-9]', '', i))))
      colnames(coords) <- c('start', 'end')
      coords$start <- as.numeric(coords$start)
      coords$end <- as.numeric(coords$end)
    } else {
      coords <- as.numeric(sapply(strsplit(unlist(coords), '-'), function(i) gsub('[^0-9]', '', i)))
      coords <- data.frame('start' = coords, 'end' = coords)
    }
    coords$protein <- rep(bed$V4, times_to_mult)
    coords$id <- rep(track_ids, times_to_mult)
    
  } else {
    # Single sites
    if(any(grepl('-', coords))) {
      coords <- as.data.frame(do.call(rbind, lapply(strsplit(coords, '-'), function(i) gsub('[^0-9]', '', i))))
      colnames(coords) <- c('start', 'end')
      coords$start <- as.numeric(coords$start)
      coords$end <- as.numeric(coords$end)
    } else {
      coords <- as.numeric(sapply(strsplit(coords, '-'), function(i) gsub('[^0-9]', '', i)))
      coords <- data.frame('start' = coords, 'end' = coords)
    }
    
    coords$protein <- bed$V4
    coords$id <- track_ids
  }
  
  return(coords)
}

##############################################################################################################################
## Annotate isoforms (unique group of ids in an outlying peptide)
annotate.isoforms <- function(summarized_proteoforms_sig, thres_qvalue, thres_fc, ioi_to_db_match, db_match_fasta) {
   
  unique_protef_sig_iois <- unique(summarized_proteoforms_sig$ioi)
  
  print(paste0('Unique iois: ', length(unique_protef_sig_iois)))
  pb  <- txtProgressBar(1, length(unique_protef_sig_iois), style=3)
  
  res <- do.call(rbind, lapply(1:length(unique_protef_sig_iois), function(i) {
    
    # print(i)
    
    setTxtProgressBar(pb, i)
    
    ioi = unique_protef_sig_iois[i]
    
    # Create peptide ranges
    summarized_proteoforms_sig_ioi <- summarized_proteoforms_sig[summarized_proteoforms_sig$ioi %in% ioi, ]
    
    summarized_proteoforms_sig_ioi$isoforms <- unlist(lapply(1:nrow(summarized_proteoforms_sig_ioi), function(i) {
    
    peptide <- gsub('[^A-Z]', '', summarized_proteoforms_sig_ioi$peptide[i])
    ioi <- summarized_proteoforms_sig_ioi$ioi[i]
    db_match_ioi <- ioi_to_db_match[[ioi]]
    
    db_match_seq <- db_match_fasta[db_match_ioi]
    
    matched_seq <- vmatchPattern(peptide, db_match_seq)
    db_match_ioi_included <- names(matched_seq)[sapply(matched_seq, length) == 1]
    
    # uniprot_id <- fasta_mapping_id$uniprot_id[fasta_mapping_id$ensembl_id %in% matched_ensembl_p]
    
    if(length(db_match_ioi_included) == 0) {
      return(NA)
    }
    
    
    # if(any(grepl('-', matched_ensembl_p))) {
    res <- paste(db_match_ioi_included, collapse = ';')
    return(res)
    # } else {
    #   return(NA)
    #   }
  }))
  
  # Subset to outlying and dominant proteoforms
  summarized_proteoforms_sig_ioi_outlying <- summarized_proteoforms_sig_ioi[c(which(summarized_proteoforms_sig_ioi$q_value_proteoform <= thres_qvalue & abs(summarized_proteoforms_sig_ioi$fc_proteoform) >= thres_fc),
                                                                              which(is.na(summarized_proteoforms_sig_ioi$q_value_proteoform)),
                                                                              grep('_mod_', summarized_proteoforms_sig_ioi$id)), ]
  
  
  # if(length(unique(summarized_proteoforms_sig_ioi_outlying$id)) == 0) {
  #   summarized_proteoforms_sig_ioi$isoforms <- NA
  #   return(summarized_proteoforms_sig_ioi)
  # }

  # Case of single proteoform
  if(length(unique(summarized_proteoforms_sig_ioi_outlying$id)) == 1) {
    summarized_proteoforms_sig_ioi$isoforms <- NA
    return(summarized_proteoforms_sig_ioi)
  }
  
  
  # Condition for calling putative spice variants if they appear exclusively in one proteoform
  isoforms <- unique(unlist(unique(summarized_proteoforms_sig_ioi_outlying$isoforms)))
  isoforms <- isoforms[!is.na(isoforms)]
  
  tokeep <- sapply(isoforms, function(j) {
    length(unique(summarized_proteoforms_sig_ioi_outlying$id[summarized_proteoforms_sig_ioi_outlying$isoforms %in% j])) == 1 &
      !any(unique(summarized_proteoforms_sig_ioi_outlying$membership[summarized_proteoforms_sig_ioi_outlying$isoforms %in% j]) %in% 'sim')
  })
  
  
  if(length(isoforms[tokeep]) > 0) {
    summarized_proteoforms_sig_ioi$isoforms[summarized_proteoforms_sig_ioi$q_value_proteoform > thres_qvalue | abs(summarized_proteoforms_sig_ioi$fc_proteoform) < thres_fc | is.na(summarized_proteoforms_sig_ioi$q_value_proteoform)] = NA 
    
    return(summarized_proteoforms_sig_ioi)
  } else {
    summarized_proteoforms_sig_ioi$isoforms <- NA
    return(summarized_proteoforms_sig_ioi)
  }
  }))
  
  return(res)
}

##############################################################################################################################
## Annotate with uniprot features
annotate.uniprot.features <- function(uniprot_features, 
                                      summarized_proteoforms_sig, 
                                      thres_qvalue,
                                      thres_fc,
                                      ioi_to_db_match, 
                                      fasta_mapping_id, 
                                      fasta_uniprot) {
  
  # Subset to significant proteoforms
  idx_proteof_sig <- which(summarized_proteoforms_sig$q_value_proteoform <= thres_qvalue & abs(summarized_proteoforms_sig$fc_proteoform) >= thres_fc)
  
  pb  <- txtProgressBar(1, length(uniprot_features), style=3)
  
  summarized_proteoforms_sig[idx_proteof_sig, uniprot_features]  <- do.call(cbind, lapply(1:length(uniprot_features), function(i) {
    
    uniprot_feature <- uniprot_features[i]
    
    print(uniprot_feature)
    setTxtProgressBar(pb, i)
    
    uniprot_ranges <- createRangesUniprotKBtrack(track_name = uniprot_feature)
  
    # Create peptide ranges
    res <- unlist(lapply(1:length(idx_proteof_sig), function(j) {
      # print(j)
      peptide <- gsub('[^A-Z]', '', summarized_proteoforms_sig$peptide[idx_proteof_sig[j]])
      
      ioi <- summarized_proteoforms_sig$ioi[idx_proteof_sig[j]]
      db_match_ioi <- ioi_to_db_match[[ioi]]
      
      
      uniprot_id <- fasta_mapping_id$uniprot_id[fasta_mapping_id$db_id %in% db_match_ioi]
      idx <- which(values(uniprot_ranges)[, 'uniprot'] %in% uniprot_id)
      
      if(length(idx) > 0) {
        id_ranges <- uniprot_ranges[idx]
        
        # Peptide ranges
        protein_seq <- unlist(fasta_uniprot[uniprot_id])
        
        peptide_range <- unlist(vmatchPattern(pattern = peptide, subject = protein_seq))
        
        overlapping_seq <- findOverlaps(peptide_range, id_ranges)
        if(length(overlapping_seq) > 0) {
          
          paste(unique(values(id_ranges)[subjectHits(overlapping_seq), 'id']), collapse = ';')
        } else {
          return(NA)
        }
      } else {
        return(NA)
      }
    }))
  }))
    
    return(summarized_proteoforms_sig)
  
  
}
##############################################################################################################################
## Generate uniprot feature tracks
generate.uniprot.feature.track <- function(ioi, 
                                           track_name, 
                                           summarized_proteoforms, 
                                           thres_qvalue, 
                                           ioi_to_db_match, 
                                           gtf_file,
                                           gtf_db_match_id,
                                           fasta_mapping_id,
                                           hg_chain_file) {
  
  if(!track_name %in% track_df$Name){
    print("Unknown track name")
    return(NULL)
    
  } else {
    track_df <- read.csv(paste0(path_data, 'FeaturesUniprotKB.csv'))
  }
  
  summarized_proteoforms_ioi <- summarized_proteoforms[summarized_proteoforms$ioi == ioi, ]
  
  idx_sig <- which(summarized_proteoforms_ioi$q_value_proteoform <= thres_qvalue)
  track_feature <- summarized_proteoforms_ioi[idx_sig, track_name]
  
  if(length(track_feature) == 0) {
    print(paste0('no overlapping feature for ', ioi, ' in ', track_name))
    return(NULL)
  }
  
  track <- dplyr::filter(track_df, Name %in% track_name)
  
  # Bed file
  bed <- data.table::fread(as.character(track$URL), header = FALSE, sep = "\t",
                           quote = "", stringsAsFactors = FALSE,
                           data.table = FALSE, showProgress = FALSE)
  
  
  # Find corresponding uniprot id
  db_id <- ioi_to_db_match[[gsub('\\..*', '', ioi)]]
  uniprot_id <- fasta_mapping_id$uniprot_id[fasta_mapping_id[['db_id']] %in% db_id]
  
  # Principal isoform index
  uniprot_id <- uniprot_id[!grepl('-', uniprot_id)]
  db_id <- fasta_mapping_id$db_id[fasta_mapping_id$uniprot_id %in% uniprot_id]
  
  
  idx <- which(bed$V4 == uniprot_id)
  
  # Chromosomal coordinates
  gene_ranges <- gtf_file[values(gtf_file)[, gtf_db_match_id] %in% db_id]
  gene_ranges <- gene_ranges[gene_ranges$type == 'CDS', ]
  seqlevelsStyle(gene_ranges) = "UCSC" 
  
  
  # Create iranges
  track_ioi_bed <- bed[idx, ]
  
  track_ioi_ranges <- GRanges(
    seqnames = track_ioi_bed$V1,
    strand = track_ioi_bed$V6,
    ranges = IRanges(start = track_ioi_bed$V2, end = track_ioi_bed$V3),
    name = track_ioi_bed$V14)
  
  if(!is.null(hg_chain_file)) {
    
    # Liftover
    ch = rtracklayer::import.chain(hg_chain_file)
    seqlevelsStyle(track_ioi_ranges) = "UCSC" 
    track_ioi_ranges <- unlist(liftOver(track_ioi_ranges, ch))
    
  }
  
  
  track_ioi_ranges_split <- split(track_ioi_ranges, track_ioi_ranges$name)
  res_ranges <- lapply(names(track_ioi_ranges_split), function(i) {
    # print(i)
    track_ioi_ranges_split_sub <- track_ioi_ranges_split[[i]]
    track_ioi_ranges_new <- disjoin(c(gene_ranges, track_ioi_ranges_split_sub))
    idx_start <- which(start(track_ioi_ranges_new) %in% start(track_ioi_ranges_split_sub))
    idx_end <- which(end(track_ioi_ranges_new) %in% end(track_ioi_ranges_split_sub))
    
    track_ioi_ranges_new <- unique(track_ioi_ranges_new[idx_start:idx_end])
    track_ioi_ranges_new <- intersect(track_ioi_ranges_new, gene_ranges)
    seqlevels(track_ioi_ranges_new) <- seqlevelsInUse(track_ioi_ranges_new)
    track_ioi_ranges_new$name <- i
    
    track_ioi_ranges_new
    
    
  })
  
  res_ranges <- unlist(GRangesList(res_ranges))
  
  res <- AnnotationTrack(res_ranges,  
                         group = res_ranges$name,
                         name = track_name,
                         background.title	= 'grey',
                         col = 'black',
                         fill = 'black',
                         rotation.title = 90)
  
  return(res)
  
  
}
##############################################################################################################################
generate.proteoform.quant <- function(peptide_quant, protein_quant, summarized_proteoforms, thres_qvalue, thres_fc) {
  
  idxtokeep <- c(which(summarized_proteoforms_sig$q_value_proteoform <= thres_qvalue &
                       abs(summarized_proteoforms_sig$fc_proteoform) >= thres_fc),
                 grep('sim', summarized_proteoforms_sig$id))
  

  summarized_proteoforms_sig_sub <- summarized_proteoforms_sig[idxtokeep, ]
  
  proteoforms_sig <- unique(summarized_proteoforms_sig_sub$ioi)
  
  res_all <- as.data.frame(do.call(rbind, lapply(proteoforms_sig, function(k) {
    
    idx <- summarized_proteoforms_sig_sub$ioi %in% k
    peps <- summarized_proteoforms_sig_sub$peptide[idx]
    
    peptide_quant_proteof <- peptide_quant[peps,]
    peptide_quant_proteof$id <- summarized_proteoforms_sig_sub$id[idx]
    
    peptide_quant_proteof <- peptide_quant_proteof %>% group_by(id) %>%
      summarise_all(median,na.rm = TRUE)
  })))
  rownames(res_all)  <- res_all$id
  res_all$id <- NULL
  
  toexclude <- unique(gsub('^.*_', '', rownames(res_all)))
  proteof_protein_quant <- rbind(protein_quant[setdiff(rownames(protein_quant), toexclude), ],
                                 res_all)
  
  return(proteof_protein_quant)
  
}
###