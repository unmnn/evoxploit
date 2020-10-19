#' Object that contains original data and derived evolution features, in
#' addition to detailed information on the clustering.
#'
#' \code{Evoxploit} computes evolution features on longitudinal data.
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Usage:
#' \preformatted{
#' evo <- Evoxploit$new(data, label, wave_suffix = "_s", minPts = NULL,
#' eps = NULL, run = TRUE, verbose = FALSE)
#'
#' summary(evo)
#' evo$data
#' evo$label
#' evo$evo_features
#' evo$all_features
#' evo$clu
#' }
#'
#' @section Arguments:
#'
#' For Evoxploit$new():
#' \describe{
#' \item{data: }{(`data.frame`)\cr
#' The data set with all input features (predictors). The object (created with Predictor$new())
#' holding the machine learning model and the data.}
#' \item{target: }{(`factor` | `numeric`) \cr
#' The target variable (response).}
#' \item{wave_suffix: }{(`character(1)`)\cr
#' The wave suffix given as string.}
#' \item{minPts: }{(`integer(1)`)\cr
#' (optional) minPts parameter for DBSCAN clustering.}
#' \item{eps: }{(`double(1)`)\cr
#' (optional) eps parameter for DBSCAN clustering.}
#' \item{run: }{(`logical(1)`)\cr
#' Should the (evo) features extraction method be run?}
#' \item{verbose: }{(`logical(1)`)\cr
#' Whether or not to show some diagnostic messages. Defaults to FALSE. }
#' }
#'
#' @export
Evoxploit <-
  R6::R6Class(
    "Evoxploit",
    private = list(
      ..data = NULL,
      ..label = NULL,
      ..wave_suffix = NULL,
      ..wave_idx = NULL,
      ..clu = NULL,
      ..minPts = NULL,
      ..eps = NULL,
      ..train_lgc = NULL,
      ..data_evo = NULL
    ),
    public = list(
      initialize = function(data, label, wave_suffix = "_s", minPts = NULL, eps = NULL,
                            train_lgc = rep(TRUE, nrow(data)),
                            run = TRUE, verbose = FALSE){
        checkmate::assert_data_frame(data)
        checkmate::assert_true(all(purrr::map_lgl(data, ~ checkmate::test_numeric(.x) | checkmate::test_factor(.x))))

        checkmate::assert_true(is.factor(label) | is.numeric(label))
        checkmate::assert_true(nrow(data) == length(label))
        checkmate::assert_false(any(is.na(label)))
        checkmate::assert_true(length(train_lgc) == nrow(data))

        checkmate::assert_character(wave_suffix)
        if(!is.null(minPts)) checkmate::assert_integerish(minPts)
        if(!is.null(eps)) checkmate::assert_number(eps)
        checkmate::assert_logical(run)

        private$..data <- data
        private$..label <- label
        private$..wave_suffix <- wave_suffix
        private$..wave_idx <- get_unique_waves(names(data), suffix = private$..wave_suffix)
        private$..minPts <- minPts
        private$..eps <- eps
        private$..train_lgc <- train_lgc

        if(run) self$run(verbose)
      },
      run = function(verbose = FALSE) {
        private$..clu <- clustering(private$..data, private$..label,
                                    minPts = private$..minPts, eps = private$..eps)
        private$..data_evo <- vector("list", 3)
        private$..data_evo[[1]] <- create_CBMS15_attributes(private$..data,
                                                            private$..label,
                                                            private$..clu,
                                                            private$..train_lgc,
                                                            private$..wave_suffix,
                                                            verbose = verbose)
        private$..data_evo[[2]] <- create_IDA14_attributes(private$..data,
                                                           private$..label,
                                                           private$..train_lgc,
                                                           private$..wave_suffix,
                                                           verbose = verbose)
        private$..data_evo[[3]] <- create_desc_attributes(private$..data,
                                                            private$..wave_suffix,
                                                            verbose = verbose)
        names(private$..data_evo) <- c("cbms", "ida", "desc")

      },
      summary = function() {
        cat("=== Summary ===\n")
        cat("Data:\n")
        cat("  - observations: ", nrow(private$..data), "\n")
        cat("  - columns: ", ncol(private$..data), "\n")
        cat("    - #numeric variables: ", sum(purrr::map_int(private$..data, is.numeric)), "\n")
        cat("    - #factor variables: ", sum(purrr::map_int(private$..data, is.factor)), "\n")
        cat("  - waves: ", paste0(private$..wave_idx, collapse = ", "), "\n")
        cat("", "\n")
        cat("Target:", "\n")
        if(is.factor(private$..label)) {
          cat("  - #classes: ", length(unique(private$..label)), "\n")
          ta <- prop.table(table(private$..label))
          sub_max <- ta[ta == max(ta)]
          cat("  - predominant class: ", paste0(names(sub_max), collapse = ", "),
              ": ", round(100*sub_max[1]), "%", "\n")
        }else{
          cat("  - range:", min(private$..label), " - ", max(private$..label))
          cat("  - 25% quantile:", quantile(private$..label, 0.25))
          cat("  - mean:", mean(private$..label))
          cat("  - median:", median(private$..label))
          cat("  - 75% quantile:", quantile(private$..label, 0.75))
        }

        cat("", "\n")
        cat("Clustering:", "\n")
        if(is.null(private$..clu)) {
          cat("NULL", "\n")
        } else {
          for(i in seq_along(private$..clu)) {
            cat("Cluster ", i, "\n")
            cat("  - attributes: ", paste0(private$..clu[[i]]$subset_att_wave, collapse = ", "), "\n")
            cat("  - minPts: ", private$..clu[[i]]$kdist$minPts, "\n")
            cat("  - eps: ", private$..clu[[i]]$kdist$eps, "\n")
            cat("  - clusters: ", length(unique(private$..clu[[i]]$clustering_result)), "\n")
            clu_table <- prop.table(table(private$..clu[[i]]$clustering_result$cluster))
            cat("  - predominant cluster: ", paste0(names(clu_table[clu_table == max(clu_table)]), collapse = ", "),
                ": ", round(100*max(clu_table)), "%", "\n")
            cat("  - outlier: ", round(sum(private$..clu[[i]]$clustering_result$cluster == 0) /
                                         length(private$..clu[[i]]$clustering_result$cluster) * 100),
                "%", "\n")
          }
          cat("", "\n")
          cat("Evo features:", "\n")
          cat("  - sequence features: ", ncol(private$..data_evo$ida), "\n")
          cat("  - descriptive features: ", ncol(private$..data_evo$desc), "\n")
          cat("  - other evo features: ", ncol(private$..data_evo$cbms), "\n")

        }
        cat("======\n")
      }
    ),
    active = list(
      data = function(value) {
        if(missing(value)) {
          private$..data
        } else{
          stop("`$data` is read only", call. = FALSE)
        }
      },
      label = function(value) {
        if(missing(value)) {
          private$..label
        } else{
          stop("`$label` is read only", call. = FALSE)
        }
      },
      clu = function(value) {
        if(missing(value)) {
          private$..clu
        } else{
          stop("`$clu` is read only", call. = FALSE)
        }
      },
      evo_features = function() {
        if(is.null(private$..clu)) {
          warning("Did not create evo features yet. Use `$run()`.")
          return(NULL)
        }
        reduce(private$..data_evo, bind_cols)
      },
      all_features = function() {
        if(is.null(private$..clu)) {
          warning("Did not create evo features yet. Use `$run()`")
          return(NULL)
        }
        reduce(list(private$..data, private$..data_evo), bind_cols)
      }
    ))


#' Print summary of Evoxploit object
#'
#' Prints a comprehensive summary of the evoxploit object to the console.
#'
#' @param x The `Evoxploit` object
#'
#' @export
#'
summary.Evoxploit <- function(x) {x$summary()}
