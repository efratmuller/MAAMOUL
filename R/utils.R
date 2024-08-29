#' Read all MAAMOUL inputs (internal function)
#'
#' @param global_network_edges A file listing all edges in the global network
#' @param ec_pvals A file with EC p-values
#' @param metabolite_pvals A file with metabolite p-values
#'
#' @noRd
read_inputs <- function(global_network_edges, ec_pvals, metabolite_pvals) {

  # Read clean network files
  edges <- read_csv(global_network_edges, show_col_types = FALSE)
  if (ncol(edges) > 2) log_warn('Only using the first 2 columns of the file: global_network_edges')
  edges <- edges %>% select(1:2) %>% rename(from = 1, to = 2)

  # Read node p-values
  mtb_pvals <- read_delim(
    metabolite_pvals,
    show_col_types = FALSE,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

  ec_pvals <- read_delim(
    ec_pvals,
    show_col_types = FALSE,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

  # Verify tables include needed columns
  if (! 'feature' %in% names(ec_pvals)) log_error('Missing a column named "feature" in the ec_pvals table')
  if (! 'feature' %in% names(mtb_pvals)) log_error('Missing a column named "feature" in the metabolite_pvals table')
  if (! 'pval' %in% names(ec_pvals)) log_error('Missing a column named "pval" in the ec_pvals table')
  if (! 'pval' %in% names(mtb_pvals)) log_error('Missing a column named "pval" in the metabolite_pvals table')

  # Verify that there's an overlap between network nodes and ec/metabolite features
  all_net_nodes <- unique(c(edges$from, edges$to))
  if (length(intersect(all_net_nodes, ec_pvals$feature)) == 0) log_error('No overlap identified between network nodes and EC features. Please make sure feature names and node names match')
  if (length(intersect(all_net_nodes, mtb_pvals$feature)) == 0) log_error('No overlap identified between network nodes and metabolite features. Please make sure feature names and node names match')

  # If a metabolite has more than one p-value, take the minimal one
  n_tmp <- nrow(mtb_pvals)
  mtb_pvals <- mtb_pvals %>%
    group_by(feature) %>%
    slice_min(order_by = pval, with_ties = FALSE) %>%
    ungroup()
  if (nrow(mtb_pvals) < n_tmp)
    log_info('Note that ', n_tmp-nrow(mtb_pvals),
             ' duplicated metabolites were identified and only the ones with minimal p-values are kept.')

  # Similarly, if an EC has more than one p-value, take the minimal one (can only happen in weird situations related to deprecated ECs)
  n_tmp <- nrow(ec_pvals)
  ec_pvals <- ec_pvals %>%
    group_by(feature) %>%
    slice_min(order_by = pval, with_ties = FALSE) %>%
    ungroup()
  if (nrow(ec_pvals) < n_tmp)
    log_info('Note that ', n_tmp-nrow(ec_pvals),
             ' duplicated ECs were identified and only the ones with minimal p-values are kept.')

  return(list(
    edges = edges,
    all_net_nodes = all_net_nodes,
    mtb_pvals = mtb_pvals,
    ec_pvals = ec_pvals
  ))
}


#' Initialize an igraph object with nodes marked as either anchors (1),
#'  non-anchors (2), or unknown (3).
#' @noRd
init_graph <- function(edges, mtb_df, ec_df) {
  require(igraph)

  # Node lists
  all_ecs <- unique(edges$from)
  all_metabs <- unique(edges$to)
  all_nodes <- c(all_ecs, all_metabs)

  # Organize a table with node assignments to categories 1-3 (as above)
  node_df <- bind_rows(
    # Actual metabolite scores
    mtb_df %>%
      mutate(type = 'Metabolite') %>%
      mutate(imputed = F),

    # Missing metabolite scores
    data.frame(
      feature = all_metabs[! all_metabs %in% mtb_df$feature],
      anchor = 3,
      type = 'Metabolite',
      imputed = T
    ),

    # Actual EC scores
    ec_df %>%
      mutate(type = 'EC') %>%
      mutate(imputed = F),

    # Missing EC scores
    data.frame(
      feature = all_ecs[! all_ecs %in% ec_df$feature],
      anchor = 3,
      type = 'EC',
      imputed = T
    )
  ) %>%
    # Remove nodes that aren't in the table
    filter(feature %in% all_nodes) %>%
    relocate(feature)

  # Add some node attributes for later plotting
  node_df <- node_df %>%
    mutate(shape = ifelse(type == 'Metabolite', 'circle', 'square'))

  # Build network with edge list + node weights
  g <- graph_from_data_frame(d = edges, directed = F, vertices = node_df)

  # A degree info per node
  g <- g %>% set_vertex_attr('degree', value = igraph::degree(g))

  log_info('Constructed a node-weighted network, with ', vcount(g),
           ' nodes and ', ecount(g),' edges.')
  return(g)
}

#' For eah pair of anchor nodes, compute the probability of them belonging to a
#'   single disease-associated metabolic module, using the random coloring of
#'   nodes by their p-values as described in the manuscript.
#'
#' @noRd
get_anchor_matrix <- function(
    g,
    anchors,
    bum_mtb,
    bum_ec,
    N_REPEATS,
    MAX_DIST_BETWEEN_RED_NODES
) {

  # Initialize a matrix of red (anchor) nodes X red (anchor) nodes
  anchors_mat <- matrix(0, nrow = length(anchors), ncol = length(anchors))
  rownames(anchors_mat) <- anchors
  colnames(anchors_mat) <- anchors

  # Flag grey metabolites of each type
  grey_mtbs <- (V(g)$anchor == 3) & (V(g)$type == 'Metabolite')
  grey_ecs <- (V(g)$anchor == 3) & (V(g)$type == 'EC')

  log_info('Starting graph random coloring iterations')
  for (i in 1:N_REPEATS) {
    if (i%%100 == 0) cat('.')
    if (i%%10000 == 0) cat('\n')

    # Create a temporary copy of the graph
    g_rand <- g

    # Randomly color nodes in the graph as either red or black using the
    #  following logic:
    #  1. For unobserved nodes - first randomly sample a p-value for them from
    #     the mixture model;
    V(g_rand)$pval[grey_mtbs] <- sample_from_bum(
      a = bum_mtb$a,
      lambda = bum_mtb$lambda,
      N = sum(grey_mtbs)
    )
    V(g_rand)$pval[grey_ecs] <- sample_from_bum(
      a = bum_ec$a,
      lambda = bum_ec$lambda,
      N = sum(grey_ecs)
    )

    #  2. Then, for each node: color it red with a probability (1-p) and black
    #     with probability p;
    #     (First correct p-values using FDR)
    V(g_rand)$fdr <- p.adjust(V(g_rand)$pval, method = 'fdr')

    V(g_rand)$color_rand <- sapply(
      1:vcount(g_rand),
      function(j, v_fdrs) {
        sample(1:2, prob = c(1-v_fdrs[j], v_fdrs[j]), size = 1)
      },
      v_fdrs = V(g_rand)$fdr
    )

    # Get the graph induced by red nodes only
    g_rand <- delete.vertices(g_rand, V(g_rand)[V(g_rand)$color_rand == 2])

    # Get shortest distances between every pair of anchor nodes
    anchors_in_g_rand <- anchors[anchors %in% V(g_rand)$name]
    dists <- distances(
      g_rand,
      v = anchors_in_g_rand,
      to = anchors_in_g_rand,
      algorithm = "unweighted"
    )

    # Mark pairs connected by short paths only with a 1, 0 otherwise.
    short_dists <- apply(
      dists,
      1:2,
      function(x, n) as.numeric(x <= n),
      n = MAX_DIST_BETWEEN_RED_NODES
    )

    # Record the resulting links between red nodes
    # (We use a trick using the reshape2 package so that we sum matrix cell by
    #  node names and not positions)
    # anchors_mat <- anchors_mat + short_dists
    anchors_mat <- reshape2::acast(rbind(
      reshape2::melt(anchors_mat),
      reshape2::melt(short_dists)
    ), Var1~Var2, sum)
  }
  log_info('End of graph random coloring iterations')

  # Convert to estimated probabilities
  anchors_mat <- anchors_mat / N_REPEATS

  return(anchors_mat)
}

#' Identify modules with hierarchical clustering
#'
#' @noRd
extract_modules_with_hclust <- function(
    anchors_mat,
    g_nodes,
    CUTREE_H,
    MIN_MODULE_SIZE,
    MIN_METABOLITES_IN_MODULE,
    MIN_ECS_IN_MODULE,
    HCLUST_METHOD = 'average',
    plot_outfile = NULL
) {
  require(dendextend)

  # Perform hierarchical clustering
  hc = hclust(as.dist(1-anchors_mat), method = HCLUST_METHOD)
  hc$height <- round(hc$height, 6) # See: https://stackoverflow.com/questions/46749829/the-height-component-of-tree-is-not-sorted-error-in-cutree

  # Cut the dendogram tree at a specified height
  modules <- cutree(hc, h = CUTREE_H)

  # Plot tree
  dhc <- hc %>% as.dendrogram(hang = 0.05)
  leaf_types <- 2-grepl('EC', hc$labels[hc$order]) # Marks EC leafs with 1 and metabolite leafs with 2

  # Plot
  if (!is.null(plot_outfile)) {
    p_width = max(10, round(ncol(anchors_mat) / 20))
    svg(file = plot_outfile, width = p_width, height = 5)
    p <- dhc %>%
      set("leaves_pch", c(17,19)[leaf_types]) %>%
      set("leaves_col", c("mediumorchid4", "darkgreen")[leaf_types]) %>%
      set("labels_cex", 0.3) %>%
      plot(); abline(h = CUTREE_H, col = adjustcolor("darkred", alpha = 0.5), lwd = 2)
    dev.off()
  } else { p <- NA }

  # Create a small data frame with module assignments
  modules_df <- data.frame(
    name = names(modules),
    tmp_module_id = unname(modules)
  ) %>%
    # Add node information
    left_join(g_nodes, by = 'name') %>%
    # Remove tiny modules, or those not meeting the MIN_METABOLITES_IN_MODULE /
    #  MIN_ECS_IN_MODULE criteria
    group_by(tmp_module_id) %>%
    mutate(n_anchors = n(),
           n_anchors_metabs = sum(type == 'Metabolite'),
           n_anchors_ECs = sum(type == 'EC'),
           mean_pval_anchors = mean(pval)) %>%
    ungroup() %>%
    filter(n_anchors >= MIN_MODULE_SIZE &
           n_anchors_metabs >= MIN_METABOLITES_IN_MODULE &
           n_anchors_ECs >= MIN_ECS_IN_MODULE) %>%
    # Rename modules from 1 to ...
    group_by(tmp_module_id) %>%
    mutate(module_id = cur_group_id()) %>%
    ungroup() %>%
    select(-tmp_module_id) %>%
    # Order nodes according to the hclust output
    mutate(name = factor(name, levels = hc$labels[hc$order]))

  modules_overview <- modules_df %>%
    select(module_id, n_anchors, n_anchors_metabs, n_anchors_ECs, mean_pval_anchors) %>%
    distinct()

  return(list(
    module_assignments = modules_df,
    modules_overview = modules_overview
    ))
}

#' Permutes node attributes within each node type
#'
#' Assumes that graph nodes have an attribute called 'type'
#'
#' @noRd
permute_by_node_type <- function(g) {
  # Extract node and edge lists of current graph
  edges_tmp <- as_data_frame(g, 'edges')
  nodes_tmp <- as_data_frame(g, 'vertices')

  # Permute node names, within each node type
  tmp <- nodes_tmp %>%
    select(name, type) %>%
    group_by(type) %>%
    mutate(new_name = sample(name)) %>%
    ungroup()
  node_name_map <- tmp$new_name
  names(node_name_map) <- tmp$name

  # Rename nodes
  nodes_tmp$name <- node_name_map[nodes_tmp$name]
  rownames(nodes_tmp) <- nodes_tmp$name

  # Re-construct the graph
  g_permuted <- graph_from_data_frame(d = edges_tmp, directed = F, vertices = nodes_tmp)
  return(g_permuted)
}

#' @noRd
plot_true_vs_permuted_modules <- function(
    modules_perm,
    N_VAL_PERM,
    MIN_MOD_SIZE,
    title = '') {
  tmp <- data.frame('Permutation_ID' = 1:N_VAL_PERM) %>% # Required so that permutations that resulted in no modules at all are also considered
    full_join(modules_perm, by = 'Permutation_ID') %>%
    group_by(Permutation_ID) %>%
    summarise(num_modules = sum(!is.na(n_anchors)),
              mean_module_size = mean(n_anchors),
              sd_module_size = sd(n_anchors),
              .groups = 'drop') %>%
    mutate(mean_module_size = ifelse(is.na(mean_module_size), MIN_MOD_SIZE, mean_module_size)) %>%
    mutate(Category = ifelse(Permutation_ID==-1, 'True', 'Permuted'))
  # log_info('Average number of modules in permuted graphs: ', round(mean(tmp$num_modules[-1]),2))

  ggplot(tmp, aes(x = num_modules, y = mean_module_size, fill = Category, size = Category)) +
    geom_jitter(color = 'black', shape = 21, alpha = 0.7, width = 0.2, height = 0.05) +
    theme_classic() +
    scale_fill_manual(values = c('grey70','darkred')) +
    scale_size_manual(values = c(2, 4)) +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    scale_x_continuous(expand = c(0.1, 0.1)) +
    ggtitle(title) +
    xlab('Total number of modules') +
    ylab('Average size of modules\n(anchors only)') +
    theme(legend.title = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
}
