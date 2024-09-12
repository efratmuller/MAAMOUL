#' MAAMOUL: A method for detecting microbiome-metabolome alterations in disease
#'   using metabolic networks
#'
#' MAAMOUL is a knowledge-based computational method that integrates metagenomic
#'   and metabolomic data to identify custom data-driven microbial metabolic
#'   modules associated with disease states. Unlike traditional statistical
#'   approaches, MAAMOUL leverages prior biological knowledge about bacterial
#'   metabolism to link genes to metabolites through a global, microbiome-wide
#'   metabolic network, and then projects genes' and metabolites' disease-
#'   association scores onto this network. The identified 'modules' are
#'   sub-networks in this graph that are significantly enriched with disease-
#'   associated features, both metagenomic and metabolomic.
#'
#' @param global_network_edges A path to a file holding the list of edges to be
#'   included in the global metabolic network. The file should be comma-
#'   delimited, with the first column listing EC's and the second column listing
#'   metabolites. Additional columns, if exits, will be ignored. Each row in the
#'   table indicates an edge between the EC and the metabolite.
#' @param ec_pvals A path to a file holding all metagenomic EC features and
#'   their corresponding p-values representing their association with disease.
#'   These do not have to be FDR corrected. The file should be tab-delimited,
#'   with a column named 'feature' holding EC codes in the same format as in
#'   the global network file, and a 'pval' column holding the p-values. Other
#'   columns will be ignored.
#' @param metabolite_pvals Similar to the `ec_pvals` file, but listing
#'   metabolite p-values. Metabolite codes/names should be in the same format
#'   as in the global network file.
#' @param out_dir A folder in which all output files will be saved.
#' @param SEED An integer to be used as a seed for result reproducibility.
#' @param NODE_FDR_THRESHOLD The FDR threshold to determine which nodes should
#'   be treated as 'anchors' (i.e. estimated to be disease-associated).
#'   Default: 0.1.
#' @param N_REPEATS The number of random coloring of nodes to perform.
#' @param MAX_DIST_BTWN_REDS A maximal distance between nodes for them to be
#'   considered as taking part in the same disease-associated module.
#'   Default: 4.
#' @param HCLUST_METHOD Either 'average', 'single' or 'complete'. Default:
#'   'average'. See `?hclust`.
#' @param CUTREE_H The height at which the hierarchical tree is cut to determine
#'   clusters.
#' @param MIN_MOD_SIZE The minimal size of a module to be outputted. Default: 3.
#' @param MIN_METS_IN_MOD Modules with less than this number of metabolite nodes
#'   will be discarded. Default: 0.
#' @param MIN_ECS_IN_MOD Modules with less than this number of EC nodes will be
#'   discarded. Default: 0.
#' @param N_VAL_PERM Number of node-weight permutations to perform for
#'   calculating the significance of each module.
#' @param MODULE_FDR_THRESHOLD The FDR threshold to determine which modules are
#'   significant.
#' @param N_THREADS Number of threads to use for parallel computing. Verify a
#'   sufficient number of cores with `parallel::detectCores()` first.
#'
#' @return The method outputs several tables and plots to the `out_dir` folder.
#'
#' @export
#'
#' @examples
#' write_test_files()
#' maamoul(
#'   global_network_edges = 'test_input/enzyme_compound_edges_kegg.csv',
#'   ec_pvals = 'test_input/ec_pvals.tsv',
#'   metabolite_pvals = 'test_input/mtb_pvals.tsv',
#'   out_dir = 'test_outputs',
#'   N_REPEATS = 100,
#'   N_VAL_PERM = 9,
#'   N_THREADS = 4
#'   )
maamoul <- function(
  global_network_edges,
  ec_pvals,
  metabolite_pvals,
  out_dir,
  # Misc. parameters
  SEED = 710,
  NODE_FDR_THRESHOLD = 0.1,
  N_REPEATS = 1000,
  MAX_DIST_BTWN_REDS = 4,
  HCLUST_METHOD = 'average',
  CUTREE_H = 0.8,
  MIN_MOD_SIZE = 3,
  MIN_ECS_IN_MOD = 0,
  MIN_METS_IN_MOD = 0,
  N_VAL_PERM = 99,
  MODULE_FDR_THRESHOLD = 0.2,
  N_THREADS = 1
) {

  # Check that all required packages are installed
  installed <- rownames(installed.packages())
  if (! "tidyr" %in% installed)        stop("Please install package 'tidyr'")
  if (! "dplyr" %in% installed)        stop("Please install package 'dplyr'")
  if (! "readr" %in% installed)        stop("Please install package 'readr'")
  if (! "igraph" %in% installed)       stop("Please install package 'igraph'")
  if (! "logger" %in% installed)       stop("Please install package 'logger'")
  if (! "BioNet" %in% installed)       stop("Please install package 'BioNet'")
  if (! "ggplot2" %in% installed)      stop("Please install package 'ggplot2'")
  if (! "RColorBrewer" %in% installed) stop("Please install package 'RColorBrewer'")
  if (! "dendextend" %in% installed) stop("Please install package 'dendextend'")
  if (! "cowplot" %in% installed) stop("Please install package 'cowplot'")
  if (! "foreach" %in% installed)      stop("Please install package 'conflicted'")
  if (! "doSNOW" %in% installed)       stop("Please install package 'stringr'")
  if (! "conflicted" %in% installed)   stop("Please install package 'conflicted'")

  # Required packages
  suppressMessages(require(tidyr))        # Tested with version:
  suppressMessages(require(dplyr))        # Tested with version:
  suppressMessages(require(readr))        # Tested with version:
  suppressMessages(require(igraph))       # Tested with version:
  suppressMessages(require(logger))       # Tested with version:
  suppressMessages(require(BioNet))       # Tested with version:
  suppressMessages(require(ggplot2))      # Tested with version:
  suppressMessages(require(RColorBrewer)) # Tested with version:
  suppressMessages(require(dendextend)) # Tested with version:
  suppressMessages(require(cowplot)) # Tested with version:
  suppressMessages(require(foreach))      # Tested with version:
  suppressMessages(require(doSNOW))       # Tested with version:
  suppressMessages(require(conflicted))   # Tested with version:
  conflict_prefer("select", "dplyr", quiet = T)
  conflict_prefer("filter", "dplyr", quiet = T)
  conflict_prefer("as_data_frame", "igraph", quiet = T)

  # Verify that all parameters are valid
  if (!file.exists(global_network_edges))      log_error('Invalid *global_network_edges* argument. File not found')
  if (!file.exists(ec_pvals))                  log_error('Invalid *ec_pvals* argument. File not found')
  if (!file.exists(metabolite_pvals))          log_error('Invalid *metabolite_pvals* argument. File not found')
  if (!is.numeric(NODE_FDR_THRESHOLD))         log_error('Invalid *NODE_FDR_THRESHOLD* argument. Should be a number between 0 and 1 (recommended: <= 0.1)')
  if (NODE_FDR_THRESHOLD <= 0 | NODE_FDR_THRESHOLD >= 1) log_error('Invalid *NODE_FDR_THRESHOLD* argument. Should be a number between 0 and 1 (recommended: <= 0.1)')
  if (!is.numeric(N_REPEATS))                  log_error('Invalid *N_REPEATS* argument. Should be an integer > 10')
  if (N_REPEATS < 10)                          log_error('Invalid *N_REPEATS* argument. Should be an integer > 10')
  if (!is.numeric(MAX_DIST_BTWN_REDS))         log_error('Invalid *MAX_DIST_BTWN_REDS* argument. Should be an integer (recommended values 2-5)')
  if (MAX_DIST_BTWN_REDS < 2 | MAX_DIST_BTWN_REDS > 6) log_error('Invalid *MAX_DIST_BTWN_REDS* argument. Should be an integer (recommended values 2-5)')
  if (!is.numeric(CUTREE_H))                   log_error('Invalid *CUTREE_H* argument. Should be a number between 0 and 1')
  if (CUTREE_H <= 0 | CUTREE_H >= 1)           log_error('Invalid *CUTREE_H* argument. Should be a number between 0 and 1')
  if (!is.numeric(MIN_MOD_SIZE))               log_error('Invalid *MIN_MOD_SIZE* argument. Should be an integer > 0')
  if (MIN_MOD_SIZE < 1)                        log_error('Invalid *MIN_MOD_SIZE* argument. Should be an integer > 0')
  if (!is.numeric(MIN_METS_IN_MOD))            log_error('Invalid *MIN_METS_IN_MOD* argument. Should be an integer >= 0')
  if (MIN_METS_IN_MOD < 0)                     log_error('Invalid *MIN_METS_IN_MOD* argument. Should be an integer >= 0')
  if (!is.numeric(MIN_ECS_IN_MOD))             log_error('Invalid *MIN_ECS_IN_MOD* argument. Should be an integer >= 0')
  if (MIN_ECS_IN_MOD < 0)                      log_error('Invalid *MIN_ECS_IN_MOD* argument. Should be an integer >= 0')
  if (!is.numeric(N_VAL_PERM))                 log_error('Invalid *N_VAL_PERM* argument. Should be an integer > 1')
  if (N_VAL_PERM <= 1)                         log_error('Invalid *N_VAL_PERM* argument. Should be an integer > 1')
  if (! HCLUST_METHOD %in% c('complete','average','single')) log_error('Invalid *HCLUST_METHOD* argument. Should be one of "average","single","complete"')
  if (N_THREADS >= parallel::detectCores() | N_THREADS < 1)  log_error('Invalid *N_THREADS* argument. Should be an integer between 1 and the number of available cores')
  if (!is.numeric(MODULE_FDR_THRESHOLD))              log_error('Invalid *MODULE_FDR_THRESHOLD* argument. Should be a number between 0 and 1 (recommended: <= 0.2)')
  if (MODULE_FDR_THRESHOLD <= 0 | MODULE_FDR_THRESHOLD >= 1) log_error('Invalid *MODULE_FDR_THRESHOLD* argument. Should be a number between 0 and 1 (recommended: <= 0.2)')

  # For rounding issues
  EPS = 0.000000001

  # Create required output folders
  if (dir.exists(out_dir)) log_info('Output directory "',out_dir,'" already exists. Files may be overriden.')
  dir.create(out_dir, showWarnings = F, recursive = T)

  # Start pipeline
  log_info('Working directory is: ', getwd(),'.')
  log_info('Starting module-identification pipeline.')
  set.seed(SEED)

  # ----------------------------------------------------------------------------
  # 1. Read input files
  # ----------------------------------------------------------------------------

  input <- read_inputs(global_network_edges, ec_pvals, metabolite_pvals)
  log_info('Loaded network information and feature p-values.')

  # Print some statistics about the overlap of observed features and network nodes
  log_info(sum(input$mtb_pvals$feature %in% input$all_net_nodes),
           ' of ', nrow(input$mtb_pvals),
           ' observed metabolite features are also in the network.')

  log_info(sum(input$ec_pvals$feature %in% input$all_net_nodes),
           ' of ', nrow(input$ec_pvals),
           ' observed EC features are also in the network.')

  log_info(sum(input$all_net_nodes %in%
                 c(input$ec_pvals$feature, input$mtb_pvals$feature)),
           ' of ', length(input$all_net_nodes),
           ' network nodes are observed in the data.')

  # ----------------------------------------------------------------------------
  # 2. Model p-values as beta-uniform mixture models and mark anchor nodes
  # ----------------------------------------------------------------------------

  # Fit a mixture model for the metabolite p-values
  #  (+ save plots describing the model's fit)
  bum_mtb <- fit_bum_model(
    input$mtb_pvals,
    node_type = 'mtb',
    plot_dir = out_dir
  )

  # Calculate a threshold for metabolite p-values corresponding to the desired FDR.
  # This threshold will be used for defining "anchor" nodes.
  # See explanation about the formula here:
  #  https://academic.oup.com/bioinformatics/article/19/10/1236/184434
  #  https://academic.oup.com/bioinformatics/article/24/13/i223/231653#394582999
  # Also, the relevant BioNet source code is here:
  #  https://rdrr.io/github/assaron/BioNet/src/R/Statistics.R
  mtb_thres <- fdrThreshold(NODE_FDR_THRESHOLD, bum_mtb) %>% round(4)
  if (mtb_thres < 0.001) mtb_thres <- 0.001 # We define a minimal threshold for cases where BUM-based thresholds are extremely low
  log_info('Metabolite p-value threshold based on BUM: ', mtb_thres, '.')

  # Mark metabolites as either anchor nodes (1) or not (2)
  input$mtb_pvals$anchor <- ifelse(input$mtb_pvals$pval < mtb_thres, 1, 2)
  n_mtb_anchors <- sum(input$mtb_pvals$anchor == 1)

  # Now apply the exact same steps - to the EC's
  bum_ec <- fit_bum_model(
    input$ec_pvals,
    node_type = 'ec',
    plot_dir = out_dir
  )
  ec_thres <- fdrThreshold(NODE_FDR_THRESHOLD, bum_ec) %>% round(4)
  if (ec_thres < 0.001) ec_thres <- 0.001
  log_info('EC p-value threshold based on BUM: ', ec_thres, '.')
  input$ec_pvals$anchor <- ifelse(input$ec_pvals$pval < ec_thres, 1, 2)
  n_ec_anchors <- sum(input$ec_pvals$anchor == 1)

  # Save mixture model parameters
  data.frame(
    node_type = c('EC','Metabolite'),
    BUM_param_a = c(bum_ec$a, bum_mtb$a),
    BUM_param_lambda = c(bum_ec$lambda, bum_mtb$lambda),
    FDR_THRESHOLD = NODE_FDR_THRESHOLD,
    anchor_pval_threshold = c(ec_thres, mtb_thres)
  ) %>%
    write_csv(file.path(out_dir, 'bum_parameters.csv'))

  if (n_mtb_anchors == 0)
    log_error('No anchor metabolites were found. Consider relaxing your FDR threshold.')
  if (n_ec_anchors == 0)
    log_error('No anchor ECs were found. Consider relaxing your FDR threshold.')
  log_info('Found ', n_ec_anchors, ' EC anchor nodes and ', n_mtb_anchors, ' metabolite anchor nodes.')

  # ----------------------------------------------------------------------------
  # 3. Initialize graph
  # ----------------------------------------------------------------------------

  # Initialize a graph with nodes marked as either anchors (1), non-anchors (2),
  #  or unknown (3)
  g_init <- init_graph(
    input$edges,
    input$mtb_pvals,
    input$ec_pvals
  )

  # Save node information for later
  anchor_nodes <- V(g_init)$name[V(g_init)$anchor == 1]
  g_nodes <- as_data_frame(g_init, 'vertices')

  # ----------------------------------------------------------------------------
  # 4. For each pair of anchor nodes, calculate the probability of them being
  #  connected in a disease-associated module.
  # ----------------------------------------------------------------------------

  # Takes a few minutes
  anchors_mat <- get_anchor_matrix(
    g = g_init,
    anchors = anchor_nodes,
    bum_mtb = bum_mtb,
    bum_ec = bum_ec,
    N_REPEATS = N_REPEATS,
    MAX_DIST_BETWEEN_RED_NODES = MAX_DIST_BTWN_REDS
  )

  # Organize probabilities in a table
  anchor_pairs_df <- as.data.frame(anchors_mat) %>%
    tibble::rownames_to_column('node1') %>%
    pivot_longer(cols = -node1, names_to = 'node2', values_to = 'prob')

  # ----------------------------------------------------------------------------
  # 5. Extract modules based on a hierarchical clustering
  # ----------------------------------------------------------------------------

  modules <- extract_modules_with_hclust(
    anchors_mat,
    g_nodes,
    CUTREE_H,
    MIN_MOD_SIZE,
    MIN_METS_IN_MOD,
    MIN_ECS_IN_MOD,
    HCLUST_METHOD,
    plot_outfile = file.path(out_dir, 'anchors_dendogram.svg')
  )
  module_assignments <- modules$module_assignments
  modules_overview <- modules$modules_overview

  log_info('Identified a total of ', nrow(modules_overview),' modules (before significance testing).')
  if (nrow(module_assignments) == 0) log_error('No modules identified')

  # ----------------------------------------------------------------------------
  # 6. Complete modules using Steiner tree approach
  #    (to assure each module forms a connected subgraph)
  # ----------------------------------------------------------------------------

  complete_modules <- complete_modules_with_steiner(
    g_init,
    modules_overview,
    module_assignments
    )

  # ----------------------------------------------------------------------------
  # 7. Also find modules in node-permuted graphs to later calculate true modules
  #    significance
  # ----------------------------------------------------------------------------

  # Setup local cluster for parallel computing
  cl <- makeCluster(N_THREADS)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = N_VAL_PERM, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)

  # Perform N_VAL_PERM permutations over node p-values, and then repeat all steps
  #  as before to get a "null distribution" of modules.
  log_info('Starting graph random coloring iterations - permuted graphs.')
  modules_perm <-
    foreach(j = 1:N_VAL_PERM,
            .combine='rbind',
            .options.snow = list(progress = progress),
            .packages = c('dplyr','igraph','logger')) %dopar% {

              # Permute p-values, within each node type (EC/metabolite)
              g_permuted <- permute_by_node_type(g_init)

              # Get updated list of anchor nodes
              anchor_nodes_permuted <- V(g_permuted)$name[V(g_permuted)$anchor == 1]

              # Now run module identification algorithm and return the resulting
              #  probability matrix
              anc_mat_shuffled <- get_anchor_matrix(
                g = g_permuted,
                anchors = anchor_nodes_permuted,
                bum_mtb = bum_mtb,
                bum_ec = bum_ec,
                N_REPEATS = N_REPEATS,
                MAX_DIST_BETWEEN_RED_NODES = MAX_DIST_BTWN_REDS
              )

              # Now use hierarchical clustering to get modules as before
              tmp <- extract_modules_with_hclust(
                anc_mat_shuffled,
                g_nodes = as_data_frame(g_permuted, 'vertices'),
                CUTREE_H,
                MIN_MOD_SIZE,
                MIN_METS_IN_MOD,
                MIN_ECS_IN_MOD,
                HCLUST_METHOD
              )

              # Note: For the 'null' modules we are only interested in the
              #  number of EC and metabolite anchors and their average p value
              #  per module, so we skip the steiner tree step.

              # # And lastly use steiner algorithm to complete each module into a
              # #  connected subgraph
              # complete_modules_perm <- complete_modules_with_steiner(
              #   g_permuted,
              #   tmp$modules_overview,
              #   tmp$module_assignments
              #   )

              # Return
              tmp$modules_overview %>% mutate(Permutation_ID = j)
  }
  log_info('Finished graph random coloring iterations - permuted graphs.')
  stopCluster(cl)

  modules_perm <- bind_rows(
    modules_overview %>% mutate(Permutation_ID = -1),
    modules_perm
  )

  # Finally, compute a p-value per module using the permuted modules.

  # Given a module of n anchor metabolite nodes, m anchor EC nodes, and an
  #  average p-value of p_hat, we ask:
  #  What are the chances of finding a module with at least n anchor metabolites
  #  and m anchor EC nodes that also has an average p value < p_hat?

  # Note: We slightly change the test described above to also take into account
  #  additional modules that were 'just as good' (>=size, <=p-value), in the same
  #  run. The motivation is that it is much more likely to get a single module
  #  meeting certain criteria than to get X modules with these properties...

  # For each module, compute how many modules *in the same run* were just as good
  tmp <- modules_perm %>%
    inner_join(modules_perm %>% rename(module_id.y = module_id),
               by = 'Permutation_ID',
               relationship = "many-to-many") %>%
    filter(n_anchors_ECs.x <= n_anchors_ECs.y) %>%
    filter(n_anchors_metabs.x <= n_anchors_metabs.y) %>%
    filter(mean_pval_anchors.x >= (mean_pval_anchors.y - EPS)) %>%
    group_by(Permutation_ID, module_id) %>%
    summarise(n_just_as_good_in_run = n(), .groups = 'drop')
  modules_perm <- modules_perm %>%
    left_join(tmp, by = c('module_id', 'Permutation_ID'))

  # For each module, compute how many modules *across all permutations* were
  #  just as good. Accordingly, compute the module's p-value.
  tmp <- cross_join(
    modules_perm,
    modules_perm %>% rename_with(function(s) paste0(s,'.y'))
    ) %>%
    filter(n_anchors_ECs <= n_anchors_ECs.y) %>%
    filter(n_anchors_metabs <= n_anchors_metabs.y) %>%
    filter(mean_pval_anchors >= (mean_pval_anchors.y - EPS)) %>%
    filter(n_just_as_good_in_run <= n_just_as_good_in_run.y) %>%
    # How many permutations had a module just as good?
    group_by(Permutation_ID, module_id) %>%
    summarise(nom = n_distinct(Permutation_ID.y), .groups = 'drop') %>% # Min: 1 (my run), max: N_VAL_PERM+1
    mutate(module_pval = nom / (N_VAL_PERM+1)) %>%
    # FDR correct
    group_by(Permutation_ID) %>%
    mutate(module_FDR = p.adjust(module_pval, 'fdr')) %>%
    ungroup() %>%
    select(-nom)

  modules_perm <- modules_perm %>%
    left_join(tmp, by = c('module_id', 'Permutation_ID'))

  log_info('Computed modules\' significance.')

  # Also add module-FDR info to main module overview
  modules_overview <- modules_overview %>%
    left_join(modules_perm %>%
                filter(Permutation_ID == -1) %>%
                select(module_id, module_pval, module_FDR),
              by = 'module_id')

  # Sanity check 1: Do permuted versions consistently result in fewer/smaller
  #  modules? (Note: these are not necessarily significant)
  p1 <- plot_true_vs_permuted_modules(
    modules_perm,
    N_VAL_PERM,
    MIN_MOD_SIZE,
    title = 'All modules'
  )

  # Sanity check 2: Do permuted versions consistently result in fewer/smaller
  #  *significant* modules?
  p2 <- plot_true_vs_permuted_modules(
    modules_perm %>%
      filter(module_FDR <= MODULE_FDR_THRESHOLD),
    N_VAL_PERM,
    MIN_MOD_SIZE,
    title = 'Significant modules only'
  )
  plot_outfile <- file.path(out_dir, 'true_vs_permuted_modules.png')
  ggsave(
    plot_outfile,
    plot_grid(p1 + theme(plot.margin = unit(c(0.5,0.1,0.2,0.2), "cm")),
              p2 + theme(plot.margin = unit(c(0.2,0.1,0.5,0.2), "cm")),
              ncol=1),
    width = 3.9,
    height = 5.5,
    dpi = 1200)

  # ----------------------------------------------------------------------------
  # 8. Save all results to files
  # ----------------------------------------------------------------------------

  # Save graph, module assignments, and module overview
  outfile1 <- file.path(out_dir, 'modules_overview.csv')
  outfile2 <- file.path(out_dir, 'complete_modules.csv')
  outfile3 <- file.path(out_dir, 'graph_and_data.rdata')
  write_csv(modules_overview, outfile1)
  write_csv(complete_modules, outfile2)
  save(g_init, complete_modules, modules_overview, file = outfile3)

  log_info('Done!')
  return()
}
