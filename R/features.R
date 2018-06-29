#' Generate CBMS'15+ Attributes
#'
#' Generates CBMS'15+ attributes (evolution features) from a data frame and
#' the output of a clustering with \code{\link{clustering}}, according to
#' \insertCite{Niemann:CBMS15}{evoxploit}.
#'
#' @param df A data frame.
#' @param li_clustering The output of \code{\link{clustering}}.
#' @param label The class labels given as factor vector.
#'
#' @return A data frame with the following columns:
#' \itemize{
#' \item dist_s_\strong{x}_s\strong{y}: The distance of an instance to itself.
#' \item dev_from_pop_\strong{a}_s\strong{x}_\strong{y}: For numeric
#' attributes: The deviation from the population slope for \strong{a} between
#' \strong{x} and \strong{y}. For factors: whether the instance value AND at
#' least 50\% of the population have changed.
#' \item cluster_idx_s\strong{x}: The instance's cluster ID.
#' \item lof_s\strong{x}: The instance's \emph{Local Outlier Factor}
#' \insertCite{Breuning:LOF2000}{evoxploit}.
#' \item cluster_rep_s\strong{x}: Whether the instance is the cluster's
#' \emph{representative}. The cluster representative is the instance i, where
#' the sum of path distances between i and all cluster peers is minimal.
#' \item dist_to_rep_s\strong{x}: The instance's distance to its cluster
#' representative. \code{NA} if noise instance.
#' The distance matrix is inherited from \code{li_clustering}.
#' \item path_length_to_rep_s\strong{x}: The path length of direct-density
#' connected instances between an intance and its cluster representative.
#' \item dist_to_centroid_s\strong{x}: The (HEOM) distance between an instance
#' to its centroid.
#' \item dist_to_medoid_s\strong{x}: The (HEOM) distance between an instance to
#' its medoid.
#' \item frac_class_\strong{c}_in_neighborhood_s\strong{x}: fraction of
#' instances of class \strong{c} in the instance's \eqn{\epsilon}-neighborhood
#' \item silhouette_s\strong{x}: the instance's silhouette score. (noise
#' instances are treated as one separate cluster.)
#' \item real_diff_att_\strong{a}_s\strong{x}_\strong{y}:
#' the real difference of \strong{a} between \strong{y} and \strong{x} for an
#' instance. For numeric attributes only.
#' \item abs_diff_att_\strong{a}_s\strong{x}_\strong{y}:
#' the absolute difference of \strong{a} between \strong{y} and \strong{x} for
#' an instance. For numeric attributes only.
#' \item rel_diff_att_\strong{a}_s\strong{x}_\strong{y}:
#' the relative difference of \strong{a} between \strong{y} and \strong{x} for
#' an instance. For numeric attributes only.
#' \item has_changed_att_\strong{a}_s\strong{x}_\strong{y}:
#' whether the instance has changed in \strong{a} between \strong{x} and
#' \strong{y}. For factors only.
#' \item stays_outlier_s\strong{x}_\strong{y}:
#' whether the instance was noise in \strong{x} and remains noise in \strong{y}.
#' \item becomes_outlier_s\strong{x}_\strong{y}:
#' whether the instance was not noise in \strong{x} and becomes noise in
#' \strong{y}.
#' \item was_outlier_s\strong{x}_\strong{y}:
#' whether the instance was noise in \strong{x} and is not noise in \strong{y}.
#' \item never_outlier_s\strong{x}_\strong{y}:
#' whether the instance neither was noise in \strong{x} nor in \strong{y}.
#' \item frac_same_peers_s_\strong{x}_\strong{y}:
#' the fraction of cluster peers in both \strong{x} and \strong{y}.
#' \item same_minPts_nn_s_\strong{x}_\strong{y}:
#' the fraction of \eqn{minPts}-nearest neighbors in both \strong{x} and
#' \strong{y}.
#' \item real_diff_silhouette_s_\strong{x}_\strong{y}: Silhouette change.
#' \item real_diff_lof_s_\strong{x}_\strong{y}: LOF change.
#' \item diff_s_\strong{x}_\strong{y}: The instance's change w.r.t. the CFS
#' subset of the clustering.
#' }
#'
#' Legend:
#' \itemize{
#' \item \strong{x} - a study wave index
#' \item \strong{y} - a study wave index, \strong{y} > \strong{x}
#' \item \strong{a} - an attribute
#' \item \strong{c} - a class label
#' }
#'
#' @seealso \code{\link{create_IDA14_attributes}},
#' \code{\link{create_simple_attributes}}
#'
#' @references \insertAllCited{}
#'
#' @export
#'
create_CBMS15_attributes <- function(df = NULL, label = NULL,
                                     li_clustering = NULL) {

  if(is.null(df) | is.null(li_clustering)) {
    stop("Both df and li_clustering have to be provided.")
  }

  process_blocks_total <- 8
  process_block <- 1

  # data_norm_global: Apply z-score normalization to all numeric attributes
  data_norm_global <- df %>%
    map_if(is.numeric, ~ (. - min(., na.rm = TRUE))/
             (max(., na.rm = TRUE) - min(., na.rm = TRUE))) %>%
    map_if(is.matrix, as.numeric) %>%
    bind_cols()

  ti_new_atts <- NULL

  unique_waves <- get_unique_waves(names(df))


  # ~Regression ----
  # 1. Filter all attributes that occur in at least 2 waves
  # 2. For each 2-wave permutation, calculate the degree of how much different
  # an instance changes its values in comparison with the whole population.
  # 2.a. For numerical attributes, calculate the slope between consecutive
  # moments and then take the difference to the slope of the population (i.e.,
  # the slope of a regression curve).
  # 2.b. For categorical attributes, check whether the majority of the
  # population changes its values and then calculate the agreement of the
  # instance's change.
  attribute_stem <- get_global_attname_stem(names(df),
                                            min_wave_count = 2)

  # Distance between waves ----
  ti_new_atts %<>% bind_cols(
    evo_permutations(names(data_norm_global)) %>%
      map_dfc(function(x) {
        d <- 1:nrow(data_norm_global) %>%
          map_dbl(function(instance) {
            # print(instance)
            # if(instance == 77) browser()
            # browser()
            i1 <- data_norm_global[instance,
                                   li_clustering[[x[1]+1]]$subset_att_wave]
            i2 <- data_norm_global[instance,
                                   li_clustering[[x[2]+1]]$subset_att_wave]
            names(i1) <- names(i2)
            heom(bind_rows(i1,i2))[1,2]
          })
        tibble(value = d) %>%
          set_names(str_c("dist_s_", x[1], "_", x[2]))
      })
  )
  print(str_c("Create CBMS'15 Attributes: ",
              round(process_block/process_blocks_total*100), "% completed."))
  process_block <- process_block + 1

  ti_new_atts %<>% bind_cols(
    attribute_stem %>%
      map_dfc(function(att){
        data_sub <- df %>% select(matches(str_c("^", att, "_s\\d$")))
        att_order <- str_replace(names(data_sub), ".*_s(\\d)$", "\\1") %>%
          as.integer() %>% order()
        data_sub <- data_sub[,att_order]

        evo_permutations(att_order) %>%
          map2_dfc(evo_permutations(names(data_sub)), function(ao, att_perm){
            if(is.numeric(data_sub[[1]])) {
              data_sub[, ao] %>%
                set_names(c("att_1", "att_2")) %>%
                mutate(slope = att_2 - att_1) %>%
                transmute(dev_from_slope = slope - mean(slope, na.rm = TRUE)) %>%
                set_names(str_c("dev_from_pop_", att, "_s_",
                                att_perm[1], "_", att_perm[2]))
            } else {
              data_sub[, ao] %>%
                set_names(c("att_1", "att_2")) %>%
                mutate(has_changed = !att_2 == att_1) %>%
                mutate(pop_has_changed = ifelse(mean(has_changed, na.rm = TRUE)
                                                >= .5, TRUE, FALSE)) %>%
                transmute(dev_from_pop = ifelse(! has_changed == pop_has_changed,
                                                TRUE, FALSE)) %>%
                set_names(str_c("dev_from_pop_", att, "_s_",
                                att_perm[1], "_", att_perm[2]))
            }

          })
      })
  )

  print(str_c("Create CBMS'15 Attributes: ",
              round(process_block/process_blocks_total*100), "% completed."))
  process_block <- process_block + 1

  ### ~Features linked to one wave -----
  # Cluster idx ----
  ti_new_atts %<>% bind_cols(
    map_dfc(1:length(unique_waves), function(x) {
      li_clustering[[x]]$clustering_result$cluster %>% as.factor() %>%
        as_tibble() %>%
        set_names(str_c("cluster_idx_s", unique_waves[x]))
    }
    )
  )

  for(i in 1:length(li_clustering)) {
    # LOF ----
    lof_score <- dbscan::lof(li_clustering[[i]]$dist,
                             k = li_clustering[[i]]$clustering_result$minPts)

    ti_new_atts %<>% bind_cols(tibble(dummy = lof_score))
    names(ti_new_atts)[length(names(ti_new_atts))] <- str_c("lof_s",
                                                            unique_waves[i])

    cluster_assignment <- li_clustering[[i]]$clustering_result$cluster
    unique_clusters <- unique(cluster_assignment) %>% base::sort()
    # unique_clusters <- unique_clusters[unique_clusters > 0]

    eps_dist <- li_clustering[[i]]$dist
    eps_dist[eps_dist > li_clustering[[i]]$clustering_result$eps] <- +Inf
    eps_dist_matrix <- as.matrix(eps_dist)

    paths <- e1071::allShortestPaths(eps_dist)

    # Cluster representative ----
    # Graph distance to representative ----
    # Graph path length to representative ----
    ti_new_atts %<>% mutate(cluster_rep = FALSE,
                            dist_to_rep = NA,
                            path_length_to_rep = NA)

    for(ci in 1:length(unique_clusters)) {
      idx_cluster <- which(cluster_assignment == unique_clusters[ci])
      # idx_rep <- idx_cluster[which.min(lof_score[idx_cluster])]

      path_sums <- paths$length[idx_cluster, idx_cluster] %>% apply(1, sum)
      if(all(is.na(path_sums))) {
        set.seed(123)
        idx_rep <- sample(idx_cluster, 1)
      }
      else {idx_rep <- idx_cluster[which.min(path_sums)]}

      ti_new_atts$cluster_rep[idx_rep] <- TRUE

      pa <- map(idx_cluster, function(x) {
        if(is.na(paths$length[x,idx_rep])){return(NA)}
        else{e1071::extractPath(paths, idx_rep, x)}
      })

      ti_new_atts$dist_to_rep[idx_cluster] <- pa %>%
        map_dbl(function(x) {
          if(all(is.na(x))) {return(NA)}
          else{
            tibble(from = x[1:(length(x)-1)], to = x[2:length(x)]) %>%
              group_by(from, to) %>%
              mutate(distance = eps_dist_matrix[from[1], to[1]]) %>%
              pull(distance) %>%
              sum()
          }
        })
      ti_new_atts %<>% mutate(dist_to_rep = ifelse(is.infinite(dist_to_rep),
                                                   NA, dist_to_rep))

      ti_new_atts$path_length_to_rep[idx_cluster] <-
        map_int(pa, ~if(all(is.na(.))){NA}else{length(.)})
      ti_new_atts %<>% mutate(path_length_to_rep =
                                ifelse(is.infinite(dist_to_rep),
                                       NA, path_length_to_rep))
    }
    names(ti_new_atts)[which(names(ti_new_atts) == "cluster_rep")] <-
      str_c("cluster_rep_s", unique_waves[i])
    names(ti_new_atts)[which(names(ti_new_atts) == "dist_to_rep")] <-
      str_c("dist_to_rep_s", unique_waves[i])
    names(ti_new_atts)[which(names(ti_new_atts) == "path_length_to_rep")] <-
      str_c("path_length_to_rep_s", unique_waves[i])

    # Cluster centroid ----
    # Distance to centroid ----
    ti_new_atts %<>% mutate(dist_to_centroid = NA)


    data_norm <- data_norm_global %>%
      select(li_clustering[[i]]$subset_att_wave)
    for(ci in 1:length(unique_clusters)) {
      idx_cluster <- which(cluster_assignment == unique_clusters[ci])

      centroid <- data_norm %>%
        slice(idx_cluster) %>%
        map_if(is.numeric, mean, na.rm = TRUE) %>%
        map_if(is.factor, function(x) {
          tab <- table(x)
          most_occuring <- base::sort(tab, decreasing = TRUE) %>%
            names() %>% .[1]
          x[which(x == most_occuring)[1]]
        }) %>%
        bind_cols()

      ti_new_atts$dist_to_centroid[idx_cluster] <-
        heom(bind_rows(centroid, data_norm[idx_cluster, ]))[1,-1]
    }
    names(ti_new_atts)[which(names(ti_new_atts) == "dist_to_centroid")] <-
      str_c("dist_to_centroid_s", unique_waves[i])

    # Cluster medoid ----
    # Distance to medoid ----
    # browser()
    ti_new_atts %<>% mutate(dist_to_medoid = NA)
    for(ci in 1:length(unique_clusters)) {
      idx_cluster <- which(cluster_assignment == unique_clusters[ci])

      idx_medoid <- idx_cluster[
        li_clustering[[i]]$dist %>%
          as.matrix() %>%
          .[idx_cluster, idx_cluster] %>%
          apply(1,sum) %>%
          which.min()
        ]

      ti_new_atts$dist_to_medoid[idx_cluster] <-  li_clustering[[i]]$dist %>%
        as.matrix() %>%
        .[idx_medoid, idx_cluster] %>%
        as.double()
    }
    names(ti_new_atts)[which(names(ti_new_atts) == "dist_to_medoid")] <-
      str_c("dist_to_medoid_s", unique_waves[i])
    # browser()

    # fraction of instances of class x in eps neighborhood ----
    label_levels <- levels(label)
    minPts <- li_clustering[[i]]$clustering_result$minPts
    dist_matrix <- li_clustering[[i]]$dist %>% as.matrix()

    # browser()
    temp <-
      1:nrow(ti_new_atts) %>%
      map(function(x) {
        idx_minPtsNN <- dist_matrix[x, ] %>% order() %>% .[1+(1:minPts)]
        tab_labels <- label[idx_minPtsNN] %>% table()
        tab_labels <- tab_labels / minPts
        tab_labels %>% as.numeric()
        # tibble(var_names = names(tab_labels),
        #        values = tab_labels %>% as.integer()) %>%
        #   spread(var_names, values)
      })

    ti_new_atts %<>%
      bind_cols(
        matrix(unlist(temp), ncol = length(label_levels), byrow = TRUE) %>%
          as_tibble() %>%
          set_names(str_c("frac_class_", label_levels,
                          "_in_neighborhood_s", unique_waves[i]))

      )

    # Silhouette ----
    ti_new_atts$silhouette <-
      cluster::silhouette(x = li_clustering[[i]]$clustering_result$cluster,
                          dist = li_clustering[[i]]$dist)[ ,3]
    names(ti_new_atts)[which(names(ti_new_atts) == "silhouette")] <-
      str_c("silhouette_s", unique_waves[i])

  }

  print(str_c("Create CBMS'15 Attributes: ",
              round(process_block/process_blocks_total*100), "% completed."))
  process_block <- process_block + 1

  # ~Evolution Features ----
  # Attribute changes ----
  ti_new_atts %<>%
    bind_cols(
      evo_permutations(names(df)) %>%
        map(function(p) {
          get_global_attname_stem(names(df),
                                  min_wave_count = length(unique_waves)) %>%
            map(function(x) {
              df <- df[str_c(x, "_s", c(p[1],p[2]))]

              if(is.numeric(df[[1]])) {
                real_diff <- df[,2] - df[,1]
                abs_diff <- real_diff %>% abs()
                rel_diff <- real_diff / df[,1]
                return(
                  bind_cols(real_diff, abs_diff, rel_diff) %>%
                    set_names(str_c(c("real_diff", "abs_diff", "rel_diff"),
                                    "_att_", x, "_s_", p[1], "_", p[2]))
                )
              } else {
                has_changed <- (df[,1] != df[,2]) %>% as.logical()
                return(
                  tibble(has_changed) %>%
                    set_names(str_c("has_changed", "_att_", x,
                                    "_s_", p[1], "_", p[2]))
                )
              }
            }) %>%
            bind_cols()
        }) %>% bind_cols()
    )

  print(str_c("Create CBMS'15 Attributes: ",
              round(process_block/process_blocks_total*100), "% completed."))
  process_block <- process_block + 1

  # Outlierness change ----
  ti_new_atts %<>%
    bind_cols(
      evo_permutations(names(df)) %>%
        map(function(x) {
          bind_cols(
            c1 = li_clustering[[x[1]+1]]$clustering_result$cluster,
            c2 = li_clustering[[x[2]+1]]$clustering_result$cluster) %>%
            mutate(stays_outlier = ifelse(c1 == 0 & c2 == 0, TRUE, FALSE)) %>%
            mutate(becomes_outlier = ifelse(c1 == 0 & c2 == 1, TRUE, FALSE)) %>%
            mutate(was_outlier = ifelse(c1 == 1 & c2 == 0, TRUE, FALSE)) %>%
            mutate(no_outlier = ifelse(c1 == 1 & c2 == 1, TRUE, FALSE)) %>%
            select(-c(c1,c2)) %>%
            set_names(str_c(c("stays_outlier", "becomes_outlier",
                              "was_outlier", "never_outlier"), "_s_",
                            x[1], "_", x[2]))
        }) %>% bind_cols()
    )

  print(str_c("Create CBMS'15 Attributes: ",
              round(process_block/process_blocks_total*100), "% completed."))
  process_block <- process_block + 1

  # Cluster peer change ----
  ti_new_atts %<>%
    bind_cols(
      evo_permutations(names(df)) %>%
        map(function(x) {
          c1 = li_clustering[[x[1]+1]]$clustering_result$cluster
          c2 = li_clustering[[x[2]+1]]$clustering_result$cluster

          1:nrow(df) %>%
            map_dbl(function(o) {
              peers1 <- which(c1 == c1[o])
              peers2 <- which(c2 == c2[o])

              intersection_peers <- length(intersect(peers1, peers2))/
                max(length(peers1),length(peers2))
              # if(length(intersection_minPts_nn) == 0) intersection_minPts_nn <- NA
              return(intersection_peers)
            }) %>%
            as_tibble() %>%
            set_names(str_c("frac_same_peers_s_", x[1], "_", x[2]))
        }) %>% bind_cols()
    )

  print(str_c("Create CBMS'15 Attributes: ",
              round(process_block/process_blocks_total*100), "% completed."))
  process_block <- process_block + 1




  # minPts nearest neighbor change ----
  ti_new_atts %<>%
    bind_cols(
      evo_permutations(names(df)) %>%
        map(function(x) {
          c1 = li_clustering[[x[1]+1]]$clustering_result$cluster
          c2 = li_clustering[[x[2]+1]]$clustering_result$cluster

          minPts <- li_clustering[[x[1]+1]]$clustering_result$minPts

          d1 = li_clustering[[x[1]+1]]$dist %>% as.matrix()
          d2 = li_clustering[[x[2]+1]]$dist %>% as.matrix()

          1:nrow(df) %>%
            map_dbl(function(o) {
              idx1 <- d1[o, ] %>% order() %>% .[1+(1:minPts)]
              idx2 <- d2[o, ] %>% order() %>% .[1+(1:minPts)]

              intersection_minPts_nn <- length(intersect(idx1, idx2))/minPts
              # if(length(intersection_minPts_nn) == 0) intersection_minPts_nn <- NA
              return(intersection_minPts_nn)
            }) %>%
            as_tibble() %>%
            set_names(str_c("same_minPts_nn_s_", x[1], "_", x[2]))
        }) %>% bind_cols()
    )

  print(str_c("Create CBMS'15 Attributes: ",
              round(process_block/process_blocks_total*100), "% completed."))
  process_block <- process_block + 1

  # Silhouette and LOF change ----
  clm <- names(ti_new_atts)[str_detect(names(ti_new_atts), "(silhouette|lof)")]
  stem <- str_replace(clm, "^(.*)_s\\d$", "\\1") %>% unique()

  ti_new_atts %<>%
    bind_cols(
      stem %>%
        map_dfc(function(s){
          evo_permutations(names(df)) %>%
            map_dfc(function(x) {
              df <- ti_new_atts[, str_c(s, "_s", x)]
              tibble(value = (df[, 2] - df[, 1])[[1]] ) %>%
                set_names(str_c("real_diff_", s, "_s_", x[1], "_", x[2]))
            })
        })
    )

  print(str_c("Create CBMS'15 Attributes: ",
              round(process_block/process_blocks_total*100), "% completed."))
  process_block <- process_block + 1


  # Overall change of a participant ----
  ti_new_atts %<>%
    bind_cols(
      evo_permutations(names(df)) %>%
        map_dfc(function(x) {
          s1 <- li_clustering[[x[1]+1]]$subset_att_wave
          s2 <- li_clustering[[x[2]+1]]$subset_att_wave
          res <- 1:nrow(data_norm_global) %>%
            map_dbl(function(inst){
              heom(
                bind_rows(data_norm_global[inst, s1],
                          data_norm_global[inst, s2] %>%
                            set_names(s1))
              )[1,2]
            })

          tibble(value = res) %>% set_names(str_c("diff_s_", x[1], "_", x[2]))
        }
        ))

  # Clean
  ti_new_atts %<>% map_if(is.logical,
                          function(x) {
                            x %>% as.integer() %>% as.factor()
                          }) %>% bind_cols()

  ti_new_atts <- replace_nan_inf_with_na(ti_new_atts)

  return(ti_new_atts)
}

#' Generate IDA'14 Attributes
#'
#' Generate sequence features as described in
#' \insertCite{Hielscher:IDA14}{evoxploit}.
#'
#' @param df A data frame.
#' @param train_lgc A logical vector. \code{TRUE} elements represent training
#' instances of \code{df}. If not set, all instances of \code{df} are
#' treated as training instances.
#' @param label The class labels given as factor vector.
#' @param kdist_sample_size The number of randomly sampled \eqn{epsilon} values
#' drawn from \code{dist_to_line} used as candidates for clustering. Weighted
#' random sampling is used: the larger \code{dist_to_line}, the higher the
#' probability of being drawn.
#'
#'
#' @return A data frame. For each attribute \strong{a}, the data
#' frame contains an attribute \strong{a}_seq.
#'
#' @seealso \code{\link{create_CBMS15_attributes}},
#' \code{\link{create_simple_attributes}}
#'
#' @references \insertAllCited{}
#'
#' @export
#'
create_IDA14_attributes <- function(df, train_lgc = rep(TRUE, nrow(df)),
                                    label, kdist_sample_size = 50) {
  df <- df %>%
    map_if(is.numeric, ~ (. - min(., na.rm = TRUE))/
             (max(., na.rm = TRUE) - min(., na.rm = TRUE))) %>%
    map_if(is.matrix, as.numeric) %>%
    bind_cols()

  # Get name stems of attributes that occur in at least 2 waves
  attribute_stem <- get_global_attname_stem(names(df), min_wave_count = 2)

  # Create seq feature for each global attribute
  # attribute_stem
  # c("stea", "stea_alt75") %>%
  df_seq <- attribute_stem %>%
    map_dfc(function(att) {
      print(att)
      dist_obj <- df %>% select(matches(str_c("^", att, "_s\\d$"))) %>%
        heom() %>% as.dist()
      kdist <- kdist_info(dist_obj)
      # Randomly sample `kdist_sample_size` epsilons using dist_to_line as
      # probability weights. The higher dist_to_line, i.e. the closer to the
      # knee point, the higher the probability of being drawn.

      if(any(kdist$dist_to_line > 0)) {
        set.seed(123)
        s <- sample(1:length(kdist$knn_dist_sorted), size = kdist_sample_size,
                    prob = kdist$dist_to_line)

        # Perform various clusterings on the sample of epsilon values and
        # calculate gain ratio w.r.t label
        df_s <- s %>%
          map_df(function(x) {
            clustering_result <- dbscan::dbscan(dist_obj,
                                                eps = kdist$knn_dist_sorted[x],
                                                minPts = kdist$minPts)
            cluster_assignment <- clustering_result$cluster

            # Assess quality based on class distribution of train index only
            gain_ratio <- tibble(cluster_assignment_only_train =
                                   cluster_assignment[train_lgc]) %>%
              bind_cols(tibble(label = label) %>% slice(which(train_lgc)) %>%
                          set_names("label")) %>%
              FSelector::gain.ratio(label ~ cluster_assignment_only_train,
                                    data = .) %>% as.numeric()

            return(tibble(s = x, gain_ratio = gain_ratio))
          })

        # Select index with highest gain ratio and re-cluster with the
        # respective eps value
        if(max(df_s$gain_ratio) > 0) {
          s_opt <- arrange(df_s, desc(gain_ratio)) %>% pull(s) %>% .[1]
          print(max(df_s$gain_ratio))
          clustering_result <- dbscan::dbscan(dist_obj,
                                              eps = kdist$knn_dist_sorted[s_opt],
                                              minPts = kdist$minPts)
          cluster_assignment <- clustering_result$cluster %>% as.factor()

          tibble(value = cluster_assignment) %>%
            set_names(str_c(att, "_seq")) %>% return()
        } else {return(NULL)}
      } else {return(NULL)}
    })
  df_seq <- replace_nan_inf_with_na(df_seq)
  return(df_seq)
}

#' Generate Simple Temporal Attributes
#'
#' Generates simple temporal attributes, e.g. the mean and mode of an attribute
#' over all waves.
#'
#' @param df A data frame.
#'
#' @return A data frame with the following columns for each attribute \strong{a}
#' in \code{df}:
#' \itemize{
#' \item \strong{a}_desc_att_mean: The mean of \strong{a} over all waves.
#' For numeric attributes only.
#' \item \strong{a}_desc_att_sd: The st. deviation of \strong{a} over all waves.
#' For numeric attributes only.
#' #' \item \strong{a}_desc_att_median: The median of \strong{a} over all waves.
#' For numeric attributes only.
#' \item \strong{a}_desc_att_mode: The mode of \strong{a} over all waves.
#' \item \strong{a}_desc_att_min: The minimum of \strong{a} over all waves.
#' For numeric attributes only.
#' \item \strong{a}_desc_att_max: The maximum of \strong{a} over all waves.
#' For numeric attributes only.
#' }
#'
#' @seealso \code{\link{create_CBMS15_attributes}},
#' \code{\link{create_IDA14_attributes}}
#'
#' @export
#'
create_simple_attributes <- function(df) {
  attribute_stem <- get_global_attname_stem(names(df))

  df_new <- attribute_stem %>%
    map_dfc(function(att){
      data_sub <- df %>%
        select(matches(str_c("^", att, "_s\\d$"))) %>%
        mutate(r = row_number()) %>%
        gather(key, value, -r) %>%
        group_by(r)

      if(is.factor(df[[data_sub$key[1]]])) {
        data_sub$value <- factor(data_sub$value)
      }

      if(is.numeric(data_sub$value)) {
        suppressWarnings(
          df <- bind_cols(
            data_sub %>% summarize(att_mean = mean(value, na.rm = TRUE)) %>%
              select(-r),
            data_sub %>% summarize(att_sd = sd(value, na.rm = TRUE)) %>%
              select(-r),
            data_sub %>% summarize(att_median = median(value, na.rm = TRUE)) %>%
              select(-r),
            data_sub %>% summarize(att_mode = getmode(value)) %>%
              select(-r),
            data_sub %>% summarize(att_min = min(value, na.rm = TRUE)) %>%
              select(-r),
            data_sub %>% summarize(att_max = max(value, na.rm = TRUE)) %>%
              select(-r))
        )
      } else {
        df <- bind_cols(data_sub %>%
                          summarize(att_mode = getmode(value)) %>%
                          select(-r))
      }
      names(df) <- str_c(att, "_desc_",
                         str_replace(names(df), "^att_(.*)$", "\\1"))
      return(df)
    })

  df_new <- replace_nan_inf_with_na(df_new)

  return(df_new)
}
