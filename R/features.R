#' Generate CBMS'15+ Attributes
#'
#' Generates CBMS'15+ attributes (evolution features) from a data frame and
#' the output of a clustering with \code{\link{clustering}}, according to
#' \insertCite{Niemann:CBMS15}{evoxploit}.
#'
#' @param df A data frame.
#' @param li_clustering The output of \code{\link{clustering}}.
#' @param label The target variable, given as either factor or numeric variable.
#' @param suffix A string indicating the start of the wave index suffix.
#' @param verbose Whether or not to show some diagnostic messages. Defaults to FALSE.
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
#' \item \strong{c} - the target variable
#' }
#'
#' @seealso \code{\link{create_IDA14_attributes}},
#' \code{\link{create_simple_attributes}}
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @import dplyr
#' @import purrr
#' @import stringr
#'
create_CBMS15_attributes <- function(df, label, li_clustering,
                                     train_lgc = rep(TRUE, nrow(df)),
                                     suffix, verbose = FALSE, ...) {

  debug <- FALSE

  if(verbose) {
    process_blocks_total <- 10
    process_block <- 1
    message("START calculating features from or inspired by CBMS15")
  }
  if(debug) message("Debug mode ON")

  # data_norm_global: Apply range transformation to each numeric attributes
  data_norm_global <- range_trans(df)

  # store all evo features in this data frame
  ti_new_atts <- NULL

  # some helpers:
  # int vector of wave indeces
  unique_waves <- get_unique_waves(names(df), suffix = suffix)
  # list of int vectors of pairs of wave indices
  wave_pairs <- evo_permutations(unique_waves)
  # char vector of attribute name stems that occur in at least min_wave_count waves
  attribute_stem <- get_global_attname_stem(names(df), suffix = suffix,
                                            min_wave_count = 2)

  if(debug) browser()

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

  # Distance between waves ----
  ti_new_atts  <- ti_new_atts %>% bind_cols(
    wave_pairs %>%
      map_dfc(function(x) {
        d <- 1:nrow(data_norm_global) %>%
          map_dbl(function(instance) {
            i1 <- data_norm_global[instance,
                                   li_clustering[[x[1]+1]]$subset_att_wave]
            i2 <- data_norm_global[instance,
                                   li_clustering[[x[2]+1]]$subset_att_wave]
            names(i1) <- names(i2)
            calc_dist(bind_rows(i1,i2), trans = NULL)[1,2]
          })
        tibble(value = d) %>%
          set_names(paste0("dist", suffix, x[1], "_", x[2]))
      })
  )
  if(verbose) message(paste0(process_block, "/", process_blocks_total, " blocks completed."))
  if(verbose) process_block <- process_block + 1
  if(debug) browser()

  ti_new_atts  <- ti_new_atts %>% bind_cols(
    attribute_stem %>%
      map_dfc(function(att){
        data_sub <- df %>% select(matches(paste0("^", att, suffix, "\\d+$")))
        att_order <- str_replace(names(data_sub), paste0(".*", suffix, "(\\d+)$"), "\\1") %>%
          as.integer() %>% order()
        data_sub <- data_sub[,att_order]

        # print(att)
        # browser()
        map_dfc(evo_permutations(get_unique_waves(names(data_sub), suffix)),
            function(ao) {
              if(is.numeric(data_sub[[1]])) {
                data_sub[paste0(att, suffix, ao)] %>%
                  set_names(c("att_1", "att_2")) %>%
                  mutate(slope = att_2 - att_1) %>%
                  transmute(dev_from_slope = slope - mean(slope, na.rm = TRUE)) %>%
                  set_names(paste0("dev_from_pop_", att, suffix,
                                   ao[1], "_", ao[2]))
              } else {
                data_sub[paste0(att, suffix, ao)] %>%
                  set_names(c("att_1", "att_2")) %>%
                  mutate(has_changed = !att_2 == att_1) %>%
                  mutate(pop_has_changed = ifelse(mean(has_changed, na.rm = TRUE)
                                                  >= .5, TRUE, FALSE)) %>%
                  transmute(dev_from_pop = ifelse(! has_changed == pop_has_changed,
                                                  TRUE, FALSE)) %>%
                  set_names(paste0("dev_from_pop_", att, suffix,
                                   ao[1], "_", ao[2]))
              }
            })

        # evo_permutations(att_order) %>%
        #   map2_dfc(evo_permutations(att_order), function(ao, att_perm){
        #     if(is.numeric(data_sub[[1]])) {
        #       data_sub[, ao] %>%
        #         set_names(c("att_1", "att_2")) %>%
        #         mutate(slope = att_2 - att_1) %>%
        #         transmute(dev_from_slope = slope - mean(slope, na.rm = TRUE)) %>%
        #         set_names(paste0("dev_from_pop_", att, suffix,
        #                          att_perm[1], "_", att_perm[2]))
        #     } else {
        #       data_sub[, ao] %>%
        #         set_names(c("att_1", "att_2")) %>%
        #         mutate(has_changed = !att_2 == att_1) %>%
        #         mutate(pop_has_changed = ifelse(mean(has_changed, na.rm = TRUE)
        #                                         >= .5, TRUE, FALSE)) %>%
        #         transmute(dev_from_pop = ifelse(! has_changed == pop_has_changed,
        #                                         TRUE, FALSE)) %>%
        #         set_names(paste0("dev_from_pop_", att, suffix,
        #                          att_perm[1], "_", att_perm[2]))
        #     }
        #
        #   })
      })
  )

  if(verbose) message(paste0(process_block, "/", process_blocks_total, " blocks completed."))
  if(verbose) process_block <- process_block + 1
  if(debug) browser()

  ### ~Features linked to one wave -----
  # Cluster idx ----
  ti_new_atts  <- ti_new_atts %>% bind_cols(
    map_dfc(1:length(unique_waves), function(x) {
      li_clustering[[x]]$clustering_result$cluster %>% as.factor() %>%
        as_tibble() %>%
        set_names(str_c("cluster_idx_s", unique_waves[x]))
    }
    )
  )
  if(debug) browser()


  for(i in 1:length(li_clustering)) {
    # LOF ----
    lof_score <- dbscan::lof(li_clustering[[i]]$dist,
                             k = li_clustering[[i]]$clustering_result$minPts)

    ti_new_atts  <- ti_new_atts %>% bind_cols(tibble(dummy = lof_score))
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
    ti_new_atts  <- ti_new_atts %>% mutate(cluster_rep = FALSE,
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
      ti_new_atts  <- ti_new_atts %>% mutate(dist_to_rep = ifelse(is.infinite(dist_to_rep),
                                                   NA, dist_to_rep))

      ti_new_atts$path_length_to_rep[idx_cluster] <-
        map_int(pa, ~if(all(is.na(.))){NA}else{length(.)})
      ti_new_atts  <- ti_new_atts %>% mutate(path_length_to_rep =
                                ifelse(is.infinite(dist_to_rep),
                                       NA, path_length_to_rep))
    }
    names(ti_new_atts)[which(names(ti_new_atts) == "cluster_rep")] <-
      paste0("cluster_rep", suffix, unique_waves[i])
    names(ti_new_atts)[which(names(ti_new_atts) == "dist_to_rep")] <-
      paste0("dist_to_rep", suffix, unique_waves[i])
    names(ti_new_atts)[which(names(ti_new_atts) == "path_length_to_rep")] <-
      paste0("path_length_to_rep", suffix, unique_waves[i])
    if(debug) browser()

    # Cluster centroid ----
    # Distance to centroid ----
    ti_new_atts  <- ti_new_atts %>% mutate(dist_to_centroid = NA)

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
        calc_dist(bind_rows(centroid, data_norm[idx_cluster, ]), trans = NULL)[1,-1]
    }
    names(ti_new_atts)[which(names(ti_new_atts) == "dist_to_centroid")] <-
      paste0("dist_to_centroid", suffix, unique_waves[i])

    # Cluster medoid ----
    # Distance to medoid ----
    ti_new_atts  <- ti_new_atts %>% mutate(dist_to_medoid = NA)
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
      paste0("dist_to_medoid", suffix, unique_waves[i])
    if(debug) browser()

  # fraction of instances of class x in eps neighborhood ----
  if(is.factor(label)) {
    label_levels <- levels(label)
    minPts <- li_clustering[[i]]$clustering_result$minPts
    dist_matrix <- as.matrix(li_clustering[[i]]$dist)[, train_lgc]

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

    ti_new_atts  <- ti_new_atts %>%
      bind_cols(
        matrix(unlist(temp), ncol = length(label_levels), byrow = TRUE) %>%
          as_tibble() %>%
          set_names(str_c("frac_class_", label_levels,
                          "_in_neighborhood", suffix, unique_waves[i]))

      )
  }

    # Silhouette ----
    ti_new_atts$silhouette <-
      cluster::silhouette(x = li_clustering[[i]]$clustering_result$cluster,
                          dist = li_clustering[[i]]$dist)[ ,3]
    names(ti_new_atts)[which(names(ti_new_atts) == "silhouette")] <-
      str_c("silhouette", suffix, unique_waves[i])

    if(verbose) message(paste0(process_block, "/", process_blocks_total, " blocks completed."))
    if(verbose) process_block <- process_block + 1
    if(debug) browser()
  }

  # ~Evolution Features ----
  # Attribute changes ----
  ti_new_atts  <- ti_new_atts %>%
    bind_cols(
      wave_pairs %>%
        map(function(p) {
          get_global_attname_stem(names(df),
                                  min_wave_count = length(unique_waves)) %>%
            map(function(x) {
              # browser()
              df <- df[str_c(x, suffix, c(p[1], p[2]))]

              if(is.numeric(df[[1]])) {
                real_diff <- df[,2] - df[,1]
                abs_diff <- real_diff %>% abs()
                rel_diff <- real_diff / df[,1]
                return(
                  {quietly(bind_cols)}(real_diff, abs_diff, rel_diff) %>%
                    pluck("result") %>%
                    set_names(str_c(c("real_diff", "abs_diff", "rel_diff"),
                                    "_att_", x, suffix, p[1], "_", p[2]))
                )
              } else {
                has_changed <- (df[,1] != df[,2]) %>% as.logical()
                return(
                  # {quietly(tibble)}(has_changed) %>% pluck("result") %>%
                  tibble(has_changed) %>%
                    set_names(str_c("has_changed", "_att_", x,
                                    suffix, p[1], "_", p[2]))
                )
              }
            }) %>%
            bind_cols()
        }) %>% bind_cols()
    )

  if(verbose) message(paste0(process_block, "/", process_blocks_total, " blocks completed."))
  if(verbose) process_block <- process_block + 1
  if(debug) browser()

  # Outlierness change ----
  ti_new_atts  <- ti_new_atts %>%
    bind_cols(
      wave_pairs %>%
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
                              "was_outlier", "never_outlier"), suffix,
                            x[1], "_", x[2]))
        }) %>% bind_cols()
    )

  if(verbose) message(paste0(process_block, "/", process_blocks_total, " blocks completed."))
  if(verbose) process_block <- process_block + 1
  if(debug) browser()

  # Cluster peer change ----
  ti_new_atts  <- ti_new_atts %>%
    bind_cols(
      wave_pairs %>%
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
            set_names(str_c("frac_same_peers", suffix, x[1], "_", x[2]))
        }) %>% bind_cols()
    )

  if(verbose) message(paste0(process_block, "/", process_blocks_total, " blocks completed."))
  if(verbose) process_block <- process_block + 1
  if(debug) browser()



  # minPts nearest neighbor change ----
  ti_new_atts  <- ti_new_atts %>%
    bind_cols(
      wave_pairs %>%
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
            set_names(str_c("same_minPts_nn", suffix, x[1], "_", x[2]))
        }) %>% bind_cols()
    )

  if(verbose) message(paste0(process_block, "/", process_blocks_total, " blocks completed."))
  if(verbose) process_block <- process_block + 1

  # Silhouette and LOF change ----
  clm <- names(ti_new_atts)[str_detect(names(ti_new_atts), "(silhouette|lof)")]
  stem <- str_replace(clm, "^(.*)_s\\d$", "\\1") %>% unique()

  ti_new_atts  <- ti_new_atts %>%
    bind_cols(
      stem %>%
        map_dfc(function(s){
          wave_pairs %>%
            map_dfc(function(x) {
              df <- ti_new_atts[, str_c(s, suffix, x)]
              tibble(value = (df[, 2] - df[, 1])[[1]] ) %>%
                set_names(str_c("real_diff_", s, suffix, x[1], "_", x[2]))
            })
        })
    )

  if(verbose) message(paste0(process_block, "/", process_blocks_total, " blocks completed."))
  if(verbose) process_block <- process_block + 1
  if(debug) browser()

  # Overall change of a participant ----
  ti_new_atts  <- ti_new_atts %>%
    bind_cols(
      wave_pairs %>%
        map_dfc(function(x) {
          s1 <- li_clustering[[x[1]+1]]$subset_att_wave
          s2 <- li_clustering[[x[2]+1]]$subset_att_wave
          res <- 1:nrow(data_norm_global) %>%
            map_dbl(function(inst){
              calc_dist(
                bind_rows(data_norm_global[inst, s1],
                          data_norm_global[inst, s2] %>%
                            set_names(s1)), trans = NULL, ...
              )[1,2]
            })

          tibble(value = res) %>% set_names(str_c("diff", suffix, x[1], "_", x[2]))
        }
        ))
  if(debug) browser()

  # Clean
  # Convert logical features to dichotomous factors
  ti_new_atts  <- ti_new_atts %>% map_if(is.logical,
                          function(x) {
                            x %>% as.integer() %>% as.factor()
                          }) %>% bind_cols()

  ti_new_atts <- replace_nan_inf_with_na(ti_new_atts)
  if(debug) browser()

  if(verbose) message(paste0("FINISHED calculating CBMS15 features."))

  return(ti_new_atts)
}

#' Generate IDA'14 Attributes
#'
#' Generate sequence features as described in
#' \insertCite{Hielscher:IDA14}{evoxploit}.
#'
#' @param df A data frame.
#' @param label The class labels given as factor vector.
#' @param suffix A string indicating the start of the wave index suffix.
#' @param train_lgc A logical vector. \code{TRUE} elements represent training
#' instances of \code{df}. If not set, all instances of \code{df} are
#' treated as training instances.
#' @param kdist_sample_size The number of randomly sampled \eqn{epsilon} values
#' drawn from \code{dist_to_line} used as candidates for clustering. Weighted
#' random sampling is used: the larger \code{dist_to_line}, the higher the
#' probability of being drawn.
#' @param seed Random seed. Used to draw a random sample from the k-distances.
#' @param verbose Whether or not to show some diagnostic messages. Defaults to FALSE.
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
#' @import dplyr
#' @import purrr
#' @import stringr
#'
create_IDA14_attributes <- function(df, label, train_lgc = rep(TRUE, nrow(df)),
                                    suffix, kdist_sample_size = 50,
                                    seed = 123, verbose = FALSE) {

  if(verbose) message(paste0("START calculating IDA14 features."))

  df <- range_trans(df)

  # Get name stems of attributes that occur in at least 2 waves
  attribute_stem <- get_global_attname_stem(names(df), min_wave_count = 2,
                                            suffix = suffix)

  # Create seq feature for each global attribute
  # attribute_stem
  # c("stea", "stea_alt75") %>%
  df_seq <- attribute_stem %>%
    map_dfc(function(att) {
      # print(att)
      dist_obj <- df %>% select(matches(str_c("^", att, suffix, "\\d+$"))) %>%
        calc_dist(trans = NULL) %>% as.dist()
      kdist <- kdist_info(dist_obj)
      # Randomly sample `kdist_sample_size` epsilons using dist_to_line as
      # probability weights. The higher dist_to_line, i.e. the closer to the
      # knee point, the higher the probability of being drawn.

      if(any(kdist$dist_to_line > 0)) {
        set.seed(seed)
        s <- sample(1:length(kdist$knn_dist_sorted), size = kdist_sample_size,
                    prob = kdist$dist_to_line)

        # Perform various clusterings on the sample of epsilon values and
        # calculate gain ratio or linear correlation w.r.t label
        df_s <- s %>%
          map_df(function(x) {
            clustering_result <- dbscan::dbscan(dist_obj,
                                                eps = kdist$knn_dist_sorted[x],
                                                minPts = kdist$minPts)
            cluster_assignment <- clustering_result$cluster

            # Assess quality based on class distribution of train index only
              tmp <- tibble(cluster_assignment_only_train =
                                   cluster_assignment[train_lgc]) %>%
              bind_cols(tibble(label = label) %>% slice(which(train_lgc)) %>%
                          set_names("label"))

              if(is.factor(label)) {
                res <- FSelector::gain.ratio(label ~ cluster_assignment_only_train,
                                             data = tmp)
              } else {
                res <- FSelector::linear.correlation(label ~ cluster_assignment_only_train,
                                             data = tmp)
              }
            return(tibble(s = x, importance = as.numeric(res)))
          })

        # Select index with highest absolute gain ratio or correlation and
        # re-cluster with the respective eps value
        if(all(is.na(df_s$importance))) return(NULL)
        max_importance <- max(abs(df_s$importance), na.rm = TRUE)
        if(verbose) message(paste0("Max. importance of sequence feature derived from ", att, " = ", max_importance, "."))
        if(verbose & max_importance == 0) message(paste0("Drop ", att, " since importance == 0."))
        if(max_importance > 0) {
          s_opt <- arrange(df_s, desc(importance)) %>% pull(s) %>% .[1]
          clustering_result <- dbscan::dbscan(dist_obj,
                                              eps = kdist$knn_dist_sorted[s_opt],
                                              minPts = kdist$minPts)
          cluster_assignment <- clustering_result$cluster %>% as.factor()

          tibble(value = cluster_assignment) %>%
            set_names(paste0("seq_", att)) %>% return()
        } else {return(NULL)}
      } else {return(NULL)}
    })
  df_seq <- replace_nan_inf_with_na(df_seq)

  if(verbose) message(paste0("FINISHED calculating IDA14 features."))

  return(df_seq)
}

#' Generate Descriptive Temporal Attributes
#'
#' Generates descriptive temporal attributes, e.g. the mean and mode of an attribute
#' over all waves.
#'
#' @param df A data frame.
#' @param suffix A string indicating the start of the wave index suffix.
#' @param verbose Whether or not to show some diagnostic messages. Defaults to FALSE.
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
#' @import dplyr
#' @import purrr
#' @import stringr
#'
create_desc_attributes <- function(df, suffix, verbose = FALSE) {
  if(verbose) message(paste0("START calculating 'descriptive' evolution features."))

  attribute_stem <- get_global_attname_stem(names(df), suffix = suffix)

  df_new <- attribute_stem %>%
    map_dfc(function(att){
      data_sub <- df %>%
        select(matches(str_c("^", att, suffix, "\\d+$"))) %>%
        mutate(r = row_number()) %>%
        tidyr::gather(key, value, -r) %>%
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
      names(df) <- str_c("desc_", att, "_",
                         str_replace(names(df), "^att_(.*)$", "\\1"))
      return(df)
    })

  df_new <- replace_nan_inf_with_na(df_new)

  if(verbose) message(paste0("FINISHED calculating 'descriptive' evolution features."))

  return(df_new)
}
