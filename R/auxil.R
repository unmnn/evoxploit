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
replace_nan_inf_with_na <- function(df) {
  map_dfc(df, function(x) {
    x[is.nan(x)|is.infinite(x)] <- NA
    return(x)
  })
}

#' An implementation of the HEOM metric.
#'
#' This creates distance matrix from given data.
#'
#' Please note that the function doesn't scale the variables.
#'
#' TODO: This breaks, if the data has only one column. Also, check if the sapply-
#' way to check numerics and factors is correct, since if there's only one column,
#' it returns logical vector for all cases, not for the column.
#'
#' For further reading, see \insertCite{Wilson:HEOM1997}{evoxploit}.
#'
#' Adapted from https://github.com/Tommytronic/Scatter-R/blob/master/R/heom.R.
#'
#' @param data The data frame to calculate the distance matrix from.
#' @return A distance matrix.
#'
#' @references \insertAllCited{}
#'
#' @export
heom <- function(data) {

  if(!is.data.frame(data)) data <- as.data.frame(data)

  # Index of factor columns
  factors <- which(map_lgl(data, is.factor))

  # Index of numerical columns
  numerics <- which(map_lgl(data, is.numeric))

  # Calculate the ranges of numeric variables; it's needed to calculate the distance
  # for numeric variables.
  # ranges <- map_dbl(data[numerics], max, na.rm = TRUE) - map_dbl(data[numerics], min, na.rm = TRUE)

  # Convert to numeric whole data
  # TODO: Why?
  data <- map_dfc(data, as.numeric) %>% as.matrix()

  #Distance matrix is an n by n matrix, where n is the number of rows in the data.
  dlen <- nrow(data)
  result <- matrix(nrow = nrow(data), ncol = nrow(data))
  row <- vector(mode = "numeric", length = ncol(data))

  # Two loops to calculate the distance between both
  for (i in 1:dlen) {
    for (j in 1:dlen) {

      # This is used to skip about half of the calculations, since the
      # distance between case x and y is same as between y and x.
      if(!is.na(result[i, j]))
        next

      # For factor attributes, the distance is 1, if they're equal, 0 if
      # not. What is inside the parentheses does this: if those factor-
      # columns are same, it results in 1 and if they're not same, the
      # result is zero (0). But since it must be just other way around,
      # we use exclamation mark to invert the result: the distance is
      # zero if they're same and one if they're not same.
      if(length(factors) > 0) {
        tryCatch(
          row[factors] <- as.double(!(data[i, factors] == data[j, factors])),
          error = function(e) browser()
        )
      }

      # For numeric type, distance is xi - yi / range. The code is just
      # elementwise calculations.
      if(length(numerics) > 0) {
        # row[numerics] <- ((data[i, numerics] - data[j, numerics]) / ranges) %>% as.double()
        row[numerics] <- as.double(data[i, numerics] - data[j, numerics])
      }

      # For NA, distance is 1
      row[is.na(row)] <- 1

      # Finally square the result and sum it
      ssum <- sum(row^2)
      result[i, j] <- ssum
      result[j, i] <- ssum
    }
  }

  # Since sqrt isn't really fast operation, this is done for all rows at once
  # just before returning the result.
  return(sqrt(result))
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
#'
#' @return A character vector of attribute name stems.
#' @export
#'
get_global_attname_stem <- function(att_names, min_wave_count = NULL) {
  if(is.null(min_wave_count)) {
    # Why is suppressWarnings() used here?
    suppressWarnings({
      min_wave_count <- str_replace(att_names, "^.*_s(\\d)$", "\\1") %>%
        as.integer() %>% unique() %>% length()
    })
  }

  num_suffix <- get_unique_waves(att_names)

  attribute_stem <- str_replace(att_names, "^(.*)_s.$", "\\1") %>%
    as_tibble() %>%
    count(value) %>%
    filter(n >= min_wave_count) %>%
    pull(value)

  return(attribute_stem)
}

#' Extract unique wave suffixes
#'
#' Extracts unique wave suffixes (e.g. 0 for "som_bmi_s0") from a vector of
#' attribute names.
#'
#' @param att_names A character vector of attribute names.
#'
#' @return An integer vector of (sorted) wave numbers.
#' @export
#'
get_unique_waves <- function(att_names) {
  num_suffix <-
    suppressWarnings(
      str_replace(att_names, "^.*_s(\\d)$", "\\1") %>%
        as.integer()
    ) %>%
    na.omit() %>%
    unique() %>%
    sort()

  return(num_suffix)
}
