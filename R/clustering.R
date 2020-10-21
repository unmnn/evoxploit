#' DBSCAN Clustering on Each Study Wave
#'
#' Performs feature selection and DBSCAN clustering on each study wave
#' independently.
#'
#' The function first employs CFS \insertCite{Hall:CFS2000}{evoxploit} to
#' obtain a small set of relevant,
#' non-redundant features for DBSCAN clustering. If for a variable only one wave
#' is selected (e.g. only som_bmi_s0 but not som_bmi_s1 and *._s2), it  expands
#' the feature space to include all possible realizations.
#'
#' Then, for each wave, it applies DBSCAN. If neccessary, it uses a heuristic
#' to find appropriate values for DBSCAN's parameters.
#'
#' Distance calculation is based on \code{\link{heom}}. If no values for
#' \eqn{minPts} and \eqn{\epsilon} are provided, the parameter are estimated
#' using a heuristic.
#'
#' @param df The data frame. TODO: Add requirements, e.g. w.r.t. suffix
#' @param label The class labels given as factor vector.
#' @param minPts \eqn{minPts}. If not set, it defaults to \eqn{log(n)},
#' where \eqn{n} is the number of objects, as suggested in
#' \insertCite{Kailing:RIS2003}{evoxploit}.
#' @param eps \eqn{\epsilon}.
#'
#' @seealso \code{\link{heom}}, \code{\link{best_att_subset_global}},
#' \code{\link{kdist_info}}
#'
#' @return A list containing the following elements for each wave:
#' \itemize{
#' \item \code{dist}: The \code{dist} matrix used for clustering.
#' \item \code{subset_att_wave}: A character vector of attribute names used
#' for clustering.
#' \item \code{kdist}: A list containing at least \eqn{minPts} and
#' \eqn{\epsilon}.
#' \item \code{clustering_result}: The DBSCAN clustering output.
#' }
#'
#' @references \insertAllCited{}
#'
#' @export
clustering <- function(df, label, minPts = NULL, eps = NULL, suffix = "_s", ...) {

  best_attributes_stem <- best_att_subset_global(df, label, suffix = suffix)
  unique_waves <- get_unique_waves(names(df), suffix = suffix)

  li_clustering_result <- vector("list", length(unique_waves))

  for(wave_idx in seq_along(li_clustering_result)) {
    # Subset of names of attributes to cluster on
    subset_att_wave <- paste0(best_attributes_stem, suffix, unique_waves[wave_idx])

    # Calculate Distance Matrix
    dist <- df %>%
      select(subset_att_wave) %>%
      calc_dist(...) %>%
      as.dist()

    # Get k-dist info to extract appropriate eps
    if(is.null(minPts) & is.null(eps)) li_kdist <- kdist_info(dist)
    if(!is.null(minPts) & is.null(eps)) li_kdist <- kdist_info(dist, k = minPts)
    if(is.null(minPts) & !is.null(eps)) li_kdist <- kdist_info(dist) #unrealistic case
    if(!is.null(minPts) & !is.null(eps)) li_kdist <- list(minPts = minPts, eps_opt = eps)

    clustering <- dbscan::dbscan(dist, eps = li_kdist$eps_opt,
                                           minPts = li_kdist$minPts)

    li_clustering_result[[wave_idx]] <-
      list(dist = dist, subset_att_wave = subset_att_wave, kdist = li_kdist,
           clustering_result = clustering)
  }

  return(li_clustering_result)
}

#' Attribute name stems from CFS attribute subset
#'
#' Performs CFS on all attributes that are available in all study waves and
#' returns a vector of (unique) attribute name stems.
#'
#' @param df A data frame (without class).
#' @param label The class labels given as factor vector.
#' @param suffix A string indicating the start of the wave index suffix.
#'
#' @return A character vector of unique attribute name stems.
#' @export
#'
best_att_subset_global <- function(df, label, suffix) {
  df_names <- names(df)

  # Get stems of global attribute names
  # should min_wave_count be an argument to this function?
  attribute_stem <- get_global_attname_stem(df_names, suffix = suffix)

  # Get waves
  num_suffix <- get_unique_waves(df_names, suffix = suffix)

  # expand attribute names with all possible suffixes
  attribute_subset <- str_c(rep(attribute_stem, length(num_suffix)),
                            suffix,
                            rep(num_suffix, each = length(attribute_stem)))

  # select these attributes + label
  df_cfs <- df %>%
    select(attribute_subset) %>%
    bind_cols(tibble(label = label))

  # get best attribute subset according to cfs
  cfs_result <- suppressWarnings(FSelector::cfs(label ~ ., df_cfs))

  # return unique attribute stems
  return(cfs_result %>% str_replace(paste0("^(.*)", suffix, "\\d+$"), "\\1") %>% unique())
}

#' k-dist graph calculations
#'
#' Estimates an appropriate value for the DBSCAN parameter
#' \eqn{\epsilon} based on the k-dist graph, as suggested in
#' \insertCite{Ester:DBSCAN1996}{evoxploit}.
#'
#' @param dist_obj A distance matrix of type \code{dist}.
#' @param k Used for k-dist graph. If not set, it defaults to \eqn{log(n)},
#' where \eqn{n} is the number of objects, as suggested in
#' \insertCite{Kailing:RIS2003}{evoxploit}.
#'
#' @return A list containing the following elements:
#' \itemize{
#' \item \code{knn_dist_sorted}: A double vector of sorted k-dist values.
#' Element names represent indices of \code{dist_obj}.
#' \item \code{dist_to_line}: A double vector of distances between
#' \code{knn_dist_sorted} and the line between its first and last element.
#' \item \code{idx_opt}: The index of the maximum value in \code{dist_to_line}.
#' \item \code{eps_opt}: The \code{idx_opt}th value in \code{knn_dist_sorted}.
#' \item \code{minPts}: k.
#' }
#'
#' @references \insertAllCited{}
#'
#' @export
#'
kdist_info <- function(dist_obj, k = round(log(attr(dist_obj, "Size")))) {
  k <- k - 1

  # Calculate matrix where j-th column depicts dist_objance to
  # j-nearest neighbor (knn$dist_obj)
  knn <- dbscan::kNN(dist_obj, k = k)

  num_knn_sorted <- base::sort(knn$dist[,k])

  # y = m*x + n
  intercept_point <- num_knn_sorted[1] # n
  slope <- (num_knn_sorted[length(num_knn_sorted)] - num_knn_sorted[1]) /
    length(num_knn_sorted) # m
  line <- intercept_point + slope * (1:length(num_knn_sorted))

  # ALTERNATIVE
  # p1 <- c(1,num_knn_sorted[1])
  # p2 <- c(length(num_knn_sorted),num_knn_sorted[length(num_knn_sorted)])
  #
  # dist_to_line <- abs((p2[2] - p1[2])*(1:length(num_knn_sorted)) -
  #                       (p2[1] - p1[1])*num_knn_sorted +
  #                       p2[1]*p1[2] - p2[2]*p1[1]) /
  #   sqrt((p2[2] - p1[2])^2 + (p2[1] - p1[1])^2)

  # Calculate distance to line
  dist_to_line <- abs(line - num_knn_sorted)

  idx_opt <- which.max(dist_to_line)
  eps_opt <- num_knn_sorted[idx_opt]

  return(list(knn_dist_sorted = num_knn_sorted,
              dist_to_line = dist_to_line,
              idx_opt = idx_opt,
              eps_opt = eps_opt,
              minPts = k + 1))
}
