#' Calculate a distance data frame
#'
#' Switchman function that forwards to a specific distance function or
#' uses a custom distance method.
#'
#' @param data Data frame.
#' @param fun The distance calculation method, either as a string or as custom function.
#' @param trans Optional transformation before distance calculation as string or custom function.
#'
#' @return nxn distance data frame
#'
#' @export
#'
calc_dist <- function(data, fun = "heom", trans = "range") {

  if(is.character(trans)) {
    # browser()
    checkmate::assert_choice(trans, c("range", "scale"))
    if(trans == "range") data <- range_trans(data)
    if(trans == "scale") data <- range_scale(data)
  }
  if(is.function(trans)) data <- trans(data)

  if(is.character(fun)) {
    checkmate::assert_choice(fun, c("heom", "euclidean"))
    if(fun == "heom") d <- heom(data)
    if(fun == "euclidean") stop("TODO: implement euclidean distance")
  }
  if(is.function(fun)) d <- fun(data)

  return(d)
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
#' @import purrr
#'
#' @export
#'
heom <- function(data) {

  if(is.matrix(data)) data <- as.data.frame(data)
  checkmate::assert_data_frame(data)

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


#' Range-transforms numeric columns of a data frame
#'
#' Range-transforms each numeric column of a data frame. Returns the complete
#' data frame.
#'
#' @param data Data frame.
#'
#' @return Data frame where each numeric column is range-transformed.
#'
#' @import purrr
#'
#' @export
#'
range_trans <- function(data) {
  data <- map_if(data, is.numeric,
                 ~ (. - min(., na.rm = TRUE))/
                   (max(., na.rm = TRUE) - min(., na.rm = TRUE))) %>%
    bind_cols()

  return(data)
}

#' Scales numeric columns of a data frame
#'
#' Scales each numeric column of a data frame. Returns the complete
#' data frame.
#'
#' @param data Data frame.
#'
#' @return Data frame where each numeric column is scaled.
#'
#' @import purrr
#'
#' @export
#'
range_scale <- function(data) {
  data <- map_if(data, is.numeric, scale) %>%
    map_if(is.matrix, as.numeric) %>%
    bind_rows()

  return(data)
}
