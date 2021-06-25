
# MAN_functions -----------------------------------------------------------

# v1.3 25/06/2021

# functions used for the analysis of Microbial association networks
# this script has to be in the sub-folder source, within the working directory

# the packages for this script are loaded by the main script

# load metadata
load_metadata <- function(filelist = input_file_list, physeqs = names_list){
  study_metadata_flag <- str_detect(input_file_list, "study_metadata") 
  if(sum(study_metadata_flag)!=1){
    # create a metadata file
    warning("No metadata file found, or multiple metadata files, I am creating one with ps names, but this may cause trouble later")
    studymetadata <- tibble(label = physeqs,
                            obj_type = rep(NA_character_, length(physeqs)),
                            studyId = rep(NA_character_, length(physeqs)),
                            FMBN = rep(NA_character_, length(physeqs)),
                            samples = rep(NA_real_, length(physeqs)),
                            short_descr = rep(NA_character_, length(physeqs)),
                            type = rep(NA_character_, length(physeqs)),
                            ref_short = rep(NA_character_, length(physeqs))
                            )
  } else {
    # will throw an error if the metadata is not either .tsv or .txt
    if("study_metadata.txt" %in% filelist){
      studymetadata <- read_tsv(file.path(data_folder,"study_metadata.txt"))
    } else {
      studymetadata <- read_tsv(file.path(data_folder,"study_metadata.tsv"))
    }
    studymetadata <-  arrange(studymetadata, label)
    # check if the names match otherwise throw a warning and create the metadata
    if(!(colnames(studymetadata)[1]=="label" && identical(pull(studymetadata,label),physeqs))){
      warning("Improperly formatted metadata file (unable to match names), I am creating one with ps names, but this might cause trouble later")
      studymetadata <- tibble(label = physeqs,
                              obj_type = rep(NA_character_, length(physeqs)),
                              studyId = rep(NA_character_, length(physeqs)),
                              FMBN = rep(NA_character_, length(physeqs)),
                              samples = rep(NA_real_, length(physeqs)),
                              short_descr = rep(NA_character_, length(physeqs)),
                              type = rep(NA_character_, length(physeqs)),
                              ref_short = rep(NA_character_, length(physeqs))
      )
    }
  }
  return(studymetadata)
}


# report step 0 -----------------------------------------------------------

report_step_0 <- function(my_physeq){
  step_0 <- c(stage = "original", 
              samples = nsamples(my_physeq), 
              sequences = sum(sample_sums(my_physeq)),
              taxa = ntaxa(my_physeq),
              prop_samples = 1,
              prop_seq = 1,
              prop_taxa = 1
  )
  return(step_0)
}


# prune samples by size ---------------------------------------------------

prune_samples_by_size <- function(physeq_obj, obj_name, plot_ecdf = F, minseqs = min_seqs){
  if(plot_ecdf){
    ssums <- data.frame(sampleSums= sample_sums(physeq_obj))
    g <- ggplot(ssums, mapping = aes(x=sampleSums))+stat_ecdf() +
      geom_hline(yintercept = 0.1, linetype = I(2)) +
      scale_y_continuous(breaks = seq(0,1,.2)) +
      labs(title = str_c("Cumulative distribution, sequences, ", obj_name, sep = ""),
           x = "sequences per sample", 
           y = "proportion of sequences") +
      theme(plot.title = element_text(hjust=0.5))
    print(g)
    }
  pruned_physeq <- prune_samples(sample_sums(physeq_obj)>minseqs, physeq_obj)
  return(pruned_physeq)
}


# remove objects with low sample number -----------------------------------

rem_low_sample_obj <- function(myphylist, ms = min_samples){
  physeq_to_keep <- sapply(myphylist, function(x) nsamples(x)>min_samples)
  if(length(myphylist[physeq_to_keep])==0){
    stop("no objects left to process, try reducing min_samples")
  } else {
    return(myphylist[physeq_to_keep])
  }
}


# report step n -----------------------------------------------------------

report_step_n <- function(my_physeq, my_physeq_o, stage_name){
  step_n <- c(stage = stage_name,
              samples = nsamples(my_physeq),
              sequences = sum(sample_sums(my_physeq)),
              taxa = ntaxa(my_physeq),
              prop_samples = nsamples(my_physeq)/nsamples(my_physeq_o),
              prop_seq = sum(sample_sums(my_physeq))/sum(sample_sums(my_physeq_o)),
              prop_taxa = ntaxa(my_physeq)/ntaxa(my_physeq_o)
  )
  return(step_n)
}


# set rank names ----------------------------------------------------------

# a function for getting taxonomic levels from a phyloseq object
# expects 7 levels or will throw an error
gset_rank_names <- function(myphyseq, tax_levels){
  rknames <- rank_names(myphyseq)
  if(length(rknames)!=7) stop("sorry, I can only handle physeq objects with 7 ranks, domain:species")
  if(identical(rknames, tax_levels)){
    cat("no change in rank names necessary\n")
  } else {
    colnames(tax_table(myphyseq))<-tax_levels
    cat("rank names changed\n")
  }
  return(myphyseq)
}


# remove chloroplasts and mitochondria ------------------------------------

remove_Chl_Mit <- function(myphyseq){
  myphyseq <- subset_taxa(myphyseq, class !="Chloroplast")
  myphyseq <- subset_taxa(myphyseq, order !="Chloroplast")
  myphyseq <- subset_taxa(myphyseq, family !="Mitochondria")
  return(myphyseq)
}


# tax_glom + change names -------------------------------------------------

# a function for performing taxonomic agglomeration and setting the taxa names 
# to the name of the taxon (this is necessary for both physeqs with ASVs and those
# taken from FMBN)

tax_glom_name_change <- function(myphyseq, taxa_glom = "genus"){
  if(!(taxa_glom %in% c("genus", "family"))) stop("only works for genus and family, sorry")
  ret_physeq <- tax_glom(myphyseq, taxrank = taxa_glom)
  taxa_table <- as(tax_table(ret_physeq), "matrix")
  tnames <- switch(taxa_glom,
                   genus = taxa_table[,6],
                   family = taxa_table[,5])
                   # order = taxa_table[,4]) etc.
  # fix duplicateee names if any
  if(anyDuplicated(tnames)){
    dupli <- which(duplicated(tnames))
    for(i in seq_along(dupli)){
      tnames[dupli[i]]<-paste(tnames[dupli[i]],i,sep="_")
    }
  }
  
  taxa_names(ret_physeq) <- tnames
  return(ret_physeq)
}

# prevalence and abundance filter -----------------------------------------

# a function for prevalence ad abundance filtering, takes a phyloseq object
# and other options and returns a list
# myphyseq is the phyloseq object to process
# name is the name (will be taken from the list during processing)
# filepref the prefix for filenames
# save_prev_ab_plot flag for saving the prev ab plot
# print_prev_ab_plot flag for printing the plot
# save_prev_table flag for saving the table with the top 50 most abundant and
# prevalent taxa
# the function works with defaults
filter_by_prev_ab <- function(myphyseq, 
                              name = "phyloseqobj",
                              prevfilter = prev_filter,
                              prevthreshold = prev_threshold,
                              passboth = pass_both,
                              abthreshold = ab_threshold,
                              filenm = out_filename_pref,
                              gres = dpi_option,
                              gtype = g_type,
                              saveplot = save_prev_ab_plot,
                              printplot = print_prev_ab_plot,
                              outfolder = output_folder,
                              savepat = save_prev_table){
  # creates list for the results, the object which will be returned by the function
  prev_ab_filter_results <- list(phyobj = NULL,
                                 phyobjname=NULL,
                                 prevabplot=NULL,
                                 prevabtable=NULL)
  # check if any taxon has 0 counts and removes it
  if(any(taxa_sums(myphyseq)==0)){
    myphyseq <- filter_taxa(myphyseq, function(x) sum(x)>0, prune = T)
  }
  # get the option on the OTU table
  taxa_are_rows_flag <- myphyseq@otu_table@taxa_are_rows
  # obtain a prevalence and abundance plot
  OTUmatrixf <- as(otu_table(myphyseq), "matrix")
  if (taxa_are_rows_flag) {
    OTUmatrixf <- t(OTUmatrixf)
  }
  OTUmatrixf_relab <- OTUmatrixf / rowSums(OTUmatrixf)
  # calculate the prevalence on all columns
  prevdf <- apply(
    X = OTUmatrixf,
    MARGIN = 2,
    FUN = function(x) {
      sum(x > 0)
    }
  )
  # calculate minimum relative abundance
  min_rel_ab <- apply(
    X = OTUmatrixf_relab,
    MARGIN = 2,
    FUN = function(x) {
      min(x)
    }
  )
  # calculate max. rel abundance
  max_rel_ab <- apply(
    X = OTUmatrixf_relab,
    MARGIN = 2,
    FUN = function(x) {
      max(x)
    }
  )
  prevdf <- data.frame(
    Prevalence = prevdf,
    TotalAbundance = colSums(OTUmatrixf),
    min_rel_ab = min_rel_ab,
    max_rel_ab = max_rel_ab
  )
  # merge taxonomic information
  taxa_metadata <-
    as.data.frame(as(tax_table(myphyseq),"matrix")) %>% rownames_to_column(var = "label")
  prevdf <- prevdf %>%
    rownames_to_column(var = "label") %>%
    left_join(., select(taxa_metadata, label, domain:species))
  # apply the filter
  if (prevfilter) {
    pass_prev_filter <-
      dplyr::filter(select(prevdf, label, Prevalence),
                    Prevalence > floor(nsamples(myphyseq) *
                                         prevthreshold)) %>%
      pull(label)
    if (passboth) {
      OTUtokeep <- intersect(names(which(max_rel_ab >= abthreshold)),
                             pass_prev_filter)
    } else {
      OTUtokeep <- union(names(which(max_rel_ab >= abthreshold)),
                         pass_prev_filter)
    }
  } else {
    OTUtokeep <- names(which(maxrelab >= ab_threshold))
  }
  prevdf <- prevdf  %>%
    mutate(
      relAbundance = TotalAbundance / sum(TotalAbundance),
      pass_filters = ifelse(label %in% OTUtokeep, "T", "F")
    )
  # apply the filter and put the result in the return list
  prev_ab_filter_results$phyobj <- prune_taxa(OTUtokeep, myphyseq)
  prev_ab_filter_results$phyobjname <- name
  
  # prevalence vs abundance plot
  OTUmatrixf <-
    OTUmatrixf[, which(colnames(OTUmatrixf) %in% OTUtokeep)]
  # fraction of remaining sequences
  
  f_seq_ret <-
    round(sum(sample_sums(prev_ab_filter_results$phyobj)) / sum(sample_sums(myphyseq)), 4)
  
  # make a plot
  title_text <- paste("Prevalence vs. abundance, by Phylum,", name, sep = " ")
  # original number of taxa
  ntaxa_prefilt <- ntaxa(myphyseq)
  subtitle_text <-
    paste(
      "using the filters you retain ",
      length(OTUtokeep),
      " taxa (triangles) out of ",
      ntaxa_prefilt,
      " (",
      f_seq_ret * 100,
      "% of init. seqs.)",
      sep = ""
    )
  # a prevalence and abundance plot change rank for color and facet as appropriate
  prev_ab_plot <-
    ggplot(prevdf,
           aes(
             x = TotalAbundance,
             y = Prevalence / nrow(OTUmatrixf),
             shape = as.factor(pass_filters),
             color = phylum
           )) +
    geom_point(size = 2, alpha = 0.7) +
    facet_wrap(~ phylum) +
    geom_hline(
      yintercept = ifelse(prev_filter, prev_threshold, 0),
      alpha = 0.5,
      linetype = 2
    ) +
    labs(
      x = "total abundance",
      y = "Prevalence [Frac. Samples]",
      shape = 'pass ab. threshold',
      title = title_text,
      subtitle = subtitle_text
    ) +
    scale_x_log10() +
    scale_y_continuous(minor_breaks = seq(0, 1, 0.05)) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90)
    )
  if(printplot) print(prev_ab_plot)
  # add the plot to the list
  prev_ab_filter_results$prevabplot <- prev_ab_plot
  
  if(saveplot) ggsave(file = paste(file.path(outfolder,out_filename_pref), name, "_prevabg.",gtype,sep=""),
                      prev_ab_plot,
                      device = gtype,
                      dpi = gres)
  # finish up with the table
  prevdf <- prevdf %>%
    mutate(relprev = Prevalence / nrow(OTUmatrixf)) %>%
    arrange(-relprev,-relAbundance)
  prev_ab_filter_results$prevabtable <- prevdf
  
  if(savepat) {
    write_tsv(prevdf, 
              file = paste(file.path(outfolder,out_filename_pref), name, "_prevabt.txt",sep=""))
  }
  return(prev_ab_filter_results)
}


# microbial association network inference ---------------------------------
# infers microbial association network with NetCoMi::netConstruct using 
# error catching (returning try-error object if it fails, the microNet object if
# successful); the method and the parameter list are provided with the function 
# (and are set as options in the main script)

infer_MAN <- function(myphyseq, inf_method, method_parameters){
  mnet_object <- try(netConstruct(myphyseq,
                                  measure = inf_method,
                                  normMethod = method_parameters$normmethodPar,
                                  zeroMethod = method_parameters$zeromethodPar,
                                  sparsMethod = method_parameters$sparMethod,
                                  measurePar = method_parameters$measureParList,
                                  alpha = method_parameters$alpha,
                                  dissFunc = method_parameters$dissFuncPar,
                                  verbose = method_parameters$verbosePar,
                                  cores = ncores,
                                  seed = 123456
  )
  )
  return(mnet_object)
}


# calculate network stats -------------------------------------------------
# calculates net stats if possible, otherwise returns a try-error object. 
# I am just using defaults for most parameters, 
# in the future I might want to set more parameters

calculate_net_stats <- function(microNet_obj, verbosePar = 1){
  netstats <- try(
    netAnalyze(
      microNet_obj,
      centrLCC = TRUE,
      avDissIgnoreInf = FALSE,
      sPathAlgo = "dijkstra",
      sPathNorm = TRUE,
      normNatConnect = TRUE,
      connectivity = TRUE,
      clustMethod = NULL,
      clustPar = NULL,
      clustPar2 = NULL,
      weightClustCoef = TRUE,
      hubPar = "eigenvector",
      hubQuant = 0.95,
      lnormFit = FALSE,
      weightDeg = FALSE,
      normDeg = TRUE,
      normBetw = TRUE,
      normClose = TRUE,
      normEigen = TRUE,
      verbose = verbosePar
    )
  )
  return(netstats)
}


# extract node properties -------------------------------------------------
# a function to extract, as a data frame, the node properties from a
# "microNetProps" object


# extracting the node stats 

extract_node_stats <- function(net_stat_list, nodestat){
    node_props <- data.frame(
      pos_degree = apply(net_stat_list$input$assoMat1, 2, function(x) sum(x>0)-1),
      neg_degree = apply(net_stat_list$input$assoMat1, 2, function(x) sum(x<0)),
      clust_memb = net_stat_list[["clustering"]][["clust1"]],
      degree = net_stat_list[["centralities"]][["degree1"]],
      between = net_stat_list[["centralities"]][["between1"]],
      close = net_stat_list[["centralities"]][["close1"]],
      eigenv = net_stat_list[["centralities"]][["eigenv1"]]
    )
    node_props <- node_props %>% 
      mutate(is_hub = rownames(node_props) %in% net_stat_list$hubs$hubs1) %>%
      rownames_to_column("label")
    node_props <- left_join(node_props, nodestat)
    return(node_props)
}


# converting microNet in tidygraph objects --------------------------------

# note that some of the conditions handled by if-else structures are redundant
# and unlikely to happen in practice. The function can probably be simplified

microNet_to_tidygraph <- function(micronet_obj,
                                  net_to_use = 1,
                                  add_names_to_edges = T,
                                  use_asso_matrix = T,
                                  fail_w_err = T
                                  ) {
  # assuming you have NetCoMi installed and loaded
  needed_packages <- c("dplyr" , "igraph", "tidygraph")
  if (!all(needed_packages %in% .packages())) {
    cat("installing/loading needed packages\n")
    .to_be_loaded <-
      needed_packages[!(needed_packages %in% .packages())]
    # check if all are in installed packages and install missing
    .inst2 <- .to_be_loaded %in% installed.packages()
    if (any(!.inst2))
      install.packages(.to_be_loaded[!.inst2])
    sapply(.to_be_loaded, require, character.only = TRUE)
  }
  # check the class of the object
  if ((class(micronet_obj) != "microNet")) {
    if (fail_w_err) {
      stop("this function only handles microNet objects\n")
    } else {
      warning("this function only handles microNet objects\n")
      tidygraph_from_micronet <- "not a microNet object"
      return(tidygraph_from_micronet)
    }
  } else {
    # get the adjacency matrix from microNet object
    if (net_to_use == 1) {
      adja_mat <- micronet_obj$adjaMat1
    } else {
      adja_mat <- micronet_obj$adjaMat2
    }
    # check if the adjacency matrix is NULL
    if (is.null(adja_mat) | !is.matrix(adja_mat)) {
      if (fail_w_err) {
        stop("cannot find the adjacency matrix\n")
      } else {
        warning("cannot find the adjacency matrix\n")
        tidygraph_from_micronet <- "not a microNet object"
        return(tidygraph_from_micronet)
      }
    } else {
      # create igraph object from the adjacency matrix
      igraph_from_micronet <- graph_from_adjacency_matrix(
        adja_mat,
        mode = "lower",
        weighted = T,
        diag = F,
        add.colnames = NULL
      )
      # create tidygraph object (from the adjacency matrix) note that this also has class igraph
      tidygraph_from_micronet <- as_tbl_graph(igraph_from_micronet)
      # check if there is at least one edge
      has_edges <- nrow((tidygraph_from_micronet %>% activate(edges) %>% as_tibble())) >= 1
      if(has_edges){
        # optionally add names to edges
        if (add_names_to_edges) {
          tidygraph_from_micronet <- tidygraph_from_micronet %>%
            activate(edges) %>%
            mutate(from_name = .N()$name[from],
                   to_name = .N()$name[to])
        }
        if (use_asso_matrix) {
          # create a network from the association matrix
          if (net_to_use == 1) {
            asso_mat <-
              micronet_obj$assoMat1
          } else {
            asso_mat <- micronet_obj$assoMat2
          }
          if (is.null(asso_mat) | !is.matrix(asso_mat)) {
            if (fail_w_err) {
              stop("cannot find the association matrix\n")
            } else {
              warning("cannot find the association matrix\n")
              return(tidygraph_from_micronet)
            }
          } else {
            asso_graph <- graph_from_adjacency_matrix(
              asso_mat,
              mode = "lower",
              weighted = T,
              diag = F,
              add.colnames = NULL
            )
            tidy_asso_graph <- as_tbl_graph(asso_graph)
            edge_tidy_asso_graph <- tidy_asso_graph %>%
              activate(edges) %>%
              as_tibble() %>%
              dplyr::rename(asso_est = weight) %>%
              mutate(asso_type = if_else(asso_est > 0, "copres", "mut_ex"))
            # join the association measure
            tidygraph_from_micronet <- tidygraph_from_micronet %>%
              activate(edges) %>%
              left_join(., edge_tidy_asso_graph)
          }
        }
        return(tidygraph_from_micronet)
      } else {
        warning("the network has 0 edges\n")
        tidygraph_from_micronet <- "no edges"
        return(tidygraph_from_micronet)
      }
    }
  }
}


# get node stats + calculate edge stats ----------------------------------------------------------
# A function for calculating edge stats (here it does only centrality_edge_betweennes)
# takes a tidygraph object and a further argument which is either a logical or a data frame
# with node stats to be merged with nodes. Returns a modified tidygraph object
merge_stats <- function(tg, node_stats = F, ebetw = calc_e_betw){
  if(!is.tbl_graph(tg)) stop("This is not a tidygraph")
  if(ebetw){
    tg <- tg %>% activate(edges) %>% mutate(edge_betw = centrality_edge_betweenness())
  }
  # check if node_stats is a data.frame to decide to merge nodes
  if("data.frame" %in% class(node_stats)){
    # renaming a bariable befor joining
    node_stats <- node_stats %>% rename(name = label)
    tg <- tg %>% activate(nodes) %>%
      left_join(.,node_stats)
  }
  return(tg)
}


# plot networks -----------------------------------------------------------
# c0l0r is either phylum or clust_memb
# this is patchy I should find a way to do it programmatically
plot_ggraph <- function(tidy_graph, name = "", method = "", 
                        c0l0r = "phylum", lp = "bottom", clp = "off"){
  g2plot <- tidy_graph %>% 
    activate(nodes) %>%
    mutate(t_deg = pos_degree + neg_degree) %>%
    dplyr::filter(pos_degree>0 | neg_degree>0) %>%
    mutate(clust_memb = as_factor(clust_memb))
  g2plot_title <- paste(name, method, sep = ", ")
  if(c0l0r == "phylum"){
  # note that using check_overlap = T may remove the names of some nodes
    ggraph_plot <- ggraph(g2plot, layout = 'fr', weights = weight) + 
    geom_edge_link(mapping = aes(edge_colour = asso_type, edge_width = weight),
                   alpha = I(0.5), show.legend = F) + 
    geom_node_point(mapping = aes(colour = phylum, size = t_deg)) +
    geom_node_text(mapping = aes(label = str_trunc(name, 15,"center", ellipsis = ".")), check_overlap = F) +
    labs(title = g2plot_title, size = "degree") +
    scale_edge_color_manual(values = (c("green","red"))) +
    scale_edge_width_continuous(range = c(1,4)) +
    coord_cartesian(clip = clp) +  
    theme_graph() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = lp)
  } else {
    ggraph_plot <- ggraph(g2plot, layout = 'fr', weights = weight) + 
      geom_edge_link(mapping = aes(edge_colour = asso_type, edge_width = weight),
                     alpha = I(0.5), show.legend = F) + 
      geom_node_point(mapping = aes(colour = clust_memb, size = t_deg)) +
      geom_node_text(mapping = aes(label = str_trunc(genus,15,"center", ellipsis = ".")), check_overlap = F) +
      labs(title = g2plot_title, size = "degree") +
      scale_edge_color_manual(values = (c("green","red"))) +
      scale_edge_width_continuous(range = c(1,4)) +
      theme_graph() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = lp) 
  }
  return(ggraph_plot)
}


# oddsratio_assotype ------------------------------------------------------

# A function to calculate odds ratio for copresence and mutual exclusion
# associations as a function of taxonomic relationships of the two nodes
# connected by an association edge. 
# Known issues: epi2by2() function has changed with epiR 2.0.26; 
# The odds ratio function

odds_ratio <- function(inputdf, taxo_level = "family"){
  epiR_version <- packageVersion("epiR")
  if(epiR_version < "2.0.26"){
    stop("\nYou must update epiR to version 2.0.26 or higher!")
  }
  mini_df <- switch(taxo_level,
                    family = select(inputdf, asso_type, sf),
                    order = select(inputdf, asso_type, so),
                    class = select(inputdf, asso_type, sc)
  )
  colnames(mini_df)[2] <- "same_taxon"
  # copresence
  epi_list_cop <- epiR::epi.2by2(xtabs(~same_taxon + asso_type, data = mini_df))
  # relevel factors for mutE
  mini_df_2 <- mini_df %>%
    mutate(
      asso_type = factor(asso_type, levels = c("mut_ex", "copres")),
      same_taxon = factor(as.character(same_taxon))
    )
  epi_list_mutE <- epiR::epi.2by2(xtabs(~same_taxon + asso_type, data = mini_df_2))
  # extract in a data frame the Wald incidence risk ratio
  which_test <- tibble(assort_test = c(str_c("cop_same",taxo_level),
                                        str_c("mutE_diff",taxo_level)))
  # extracts the appropriate chisq values
  chicop <- if(epi_list_cop$massoc.detail$chi2.correction){
    epi_list_cop$massoc.detail$chi2.strata.yates
  } else {
    epi_list_cop$massoc.detail$chi2.strata.uncor
  }
  chimutE <- if(epi_list_mutE$massoc.detail$chi2.correction){
    epi_list_mutE$massoc.detail$chi2.strata.yates
  } else {
    epi_list_mutE$massoc.detail$chi2.strata.uncor
  }
  OR <- rbind(
    cbind(epi_list_cop$massoc.detail$OR.strata.wald, chicop),
    cbind(epi_list_mutE$massoc.detail$OR.strata.wald, chimutE)
  )
  names(OR) <- str_c("OR", names(OR), sep = "_")
  RR <- rbind(
    epi_list_cop$massoc.detail$RR.strata.wald, 
    epi_list_mutE$massoc.detail$RR.strata.wald
  )
  names(RR) <- str_c("RR", names(RR), sep = "_")
  assort_results <- cbind(
    which_test,
    OR,
    RR
  )
  return(assort_results)
}

# df_return
# checks a list containing data frames or try-error objects, drops the try-error objects
# and binds the data frames by row
df_return <- function(input_list){
  classes <- map_chr(input_list, class)
  input_list<-input_list[which(classes == "data.frame")]
  df <- bind_rows(input_list, .id = "method")
  return(df)
}


# a report for the functions, a tribble
print_function_report <- function(){
ext_functions <- tribble(
  ~funct_name,              ~description,
  "load metadata",          "loads and checks the metadata",
  "report_step_0",          "create a report for step 0 (original phyloseqs)",
  "prune_samples_by_size",  "removes samples with less than min_seqs, optionally prints ECDF",
  "rem_low_sample_obj",     "removes phyloseq objects with less than min_samples",
  "report_step_n",          "create a report for step n (from sample pruning onward)",
  "gset_rank_names",        "get and set rank names",
  "remove_Chl_Mit",         "removes chloroplasts and mitochondria",
  "tax_glom_name_change",   "performs taxa agglomeratio and changes name",
  "filter_by_prev_ab",      "performs prevalence and abundance filtering, with optional graphs and tables",
  "infer_MAN",              "infers microbial association networks using netConstruct",
  "calculate_net_stats",    "calculate network and node statistics",
  "extract_node_stats",     "extracts node properties from a microNetProps object",
  "microNet_to_tidygraph",  "convert a microNet object to a tidygraph object",
  "merge_stats",            "merge node stats in the tidygraph object, optionally calculates edge betweenness",
  "plot_ggraph",            "use ggraph to plot the tidygraph, and returns the object",
  "odds_ratio",             "calculates odds ratios for taxonomic assortativity",
  "df_return",              "takes a list of data frames and try-error objects, removes the latter and binds the data frames"
)
return(ext_functions)
}


# Credits and copyright ---------------------------------------------------

# script created by Eugenio Parente, Università degli Studi della Basilicata
# 2021
# eugenio.parente@unibas.it 
# https://github.com/ep142
# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the \"Software\"), to deal 
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.

