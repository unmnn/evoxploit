#' Mode of a vector
#'
#' Calculate the mode of a vector.
#'
#' @param v A vector.
#'
#' @return The vector's mode.
#'
getmode <- function(v) {
  uniqv <- sort(unique(v))
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#' Replace NA/NaN/Inf values with a value
#'
#' Replace occurences of NA/NaN/Inf values with a given value.
#'
#' @param x A vector.
#'
#' @return A vector with replaced values.
#' !!! NA -> 0
replace_na_nan_inf_vec <- function(x, v = 0) {
  x[is.na(x)|(is.nan(x)|is.infinite(x))] <- v
  return(x)
}


#' Replace all NaN/Inf values with NA within a data frame
#'
#' Replaces all occurences of NaN and Inf within a data frame with NA.
#'
#' @param df A data frame.
#'
#' @return The modified data frame.
#'
#' @import purrr
#'
replace_nan_inf_with_na <- function(df) {
  map_dfc(df, function(x) {
    x[is.nan(x)|is.infinite(x)] <- NA
    return(x)
  })
}

#' Extract (unique) stems from attribute names
#'
#' Extracts unique stems from attribute names. Returns only the stems for
#' attributes which occur in at least \code{min_wave_count} waves.
#'
#' @param att_names A character vector of attribute names.
#' @param min_wave_count The minimum number of waves an attribute must be present
#' in for its name stem to be returned. If not set, the parameter defaults to
#' the total number of waves.
#' @param suffix A string indicating the start of the wave index suffix.
#'
#' @return A character vector of attribute name stems.
#'
#' @import stringr
#' @import dplyr
#'
#' @export
#'
get_global_attname_stem <- function(att_names,
                                    min_wave_count = length(get_unique_waves(att_names, suffix)),
                                    suffix = "_s") {

  unique_wave_idx <- get_unique_waves(att_names, suffix)

  attribute_stem <- str_replace(att_names, paste0("^(.*)", suffix, "\\d+"), "\\1") %>%
    as_tibble() %>%
    count(value) %>%
    filter(n >= min_wave_count) %>%
    pull(value)

  return(attribute_stem)
}

#' Extract unique wave numbers
#'
#' Extracts unique wave numbers (e.g. 0 for "som_bmi_s0") from a vector of
#' attribute names.
#'
#' @param att_names A character vector of attribute names.
#' @param suffix A string indicating the start of the wave index suffix.
#'
#' @return An integer vector of (sorted) wave numbers.
#'
#' @import stringr
#'
#' @export
#'
get_unique_waves <- function(att_names, suffix = "_s") {
  wave_idx <- as.integer(str_replace(att_names, paste0("^.*", suffix, "(\\d+)"), "\\1"))
  checkmate::assert_integer(wave_idx, any.missing = FALSE)

  unique_wave_idx <- sort(unique(wave_idx))
  checkmate::assert_true(all(diff(unique_wave_idx) > 0))

  return(unique_wave_idx)
}

#' Find pairwise permutations of attribute names
#'
#' Find all possible pairwise wave permuations of attribute names, e.g.
#' {(som_bmi_s0, som_bmi_s1), (som_bmi_s0, som_bmi_s2),
#' (som_bmi_s1, som_bmi_s2)}. Autmatically extracts the set of wave suffixes.
#'
#' @param int_waves An integer vector of wave indices.
#'
#' @return A list of pairs of wave indices
#'
#' @export
#'
# evo_permutations <- function(df_names, suffix = "_s") {
#   unique_waves <- get_unique_waves(df_names, suffix = suffix)
#
#   perm <- list()
#
#   for(i in unique_waves) {
#     for(j in unique_waves[-1]) {
#       if(i < j) perm[[length(perm) + 1]] <- c(i,j)
#     }
#   }
#
#   return(perm)
# }
evo_permutations <- function(int_waves) {
  perm <- combn(sort(int_waves), 2)
  return(split(perm, rep(1:ncol(perm), each = nrow(perm))))
}
