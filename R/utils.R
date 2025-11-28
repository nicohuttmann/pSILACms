#' Title
#'
#' @param data 
#' @param min_height 
#' @param topN 
#' @param plot 
#'
#' @returns
#' @export
#'
#' @examples
find_density_modes <- function(data, min_height = 0, topN = NULL, plot = F) {
  
  # Estimate the density
  dens <- density(data)
  
  # Find the x-values (modes) corresponding to the peaks
  # Local maxima occur where the derivative changes from positive to negative
  diff_sign <- diff(sign(diff(dens$y)))
  peaks_idx <- which(diff_sign == -2) + 1  # -2 indicates a peak (positive to negative)
  modes <- dens$x[peaks_idx]
  
  # Get the heights of the peaks
  peak_heights <- dens$y[peaks_idx]
  
  # Remove peaks under min height
  peaks_minh <- which(peak_heights / max(peak_heights) >= min_height)
  
  modes <- modes[peaks_minh]
  peak_heights <- peak_heights[peaks_minh]
  
  # Mode order by height 
  modes_sorted <- modes[order(peak_heights, 
                              decreasing = T)]
  
  peak_heights <- peak_heights[order(peak_heights, 
                                     decreasing = T)]
  
  # Keep only topN peaks 
  if (is.null(topN)) topN <- length(modes_sorted)
  
  modes_sorted <- modes_sorted[1:topN]
  #modes <- modes[1:topN]
  peak_heights <- peak_heights[1:topN]
  
  if (plot) {
    # Identify the highest peak
    highest_mode <- modes_sorted[which.max(peak_heights)]
    highest_peak_height <- max(peak_heights)
    
    # Plot the density with the modes marked
    plot(dens, main = "Bimodal Distribution with Modes")
    abline(v = modes_sorted, col = "red", lty = 2)
    text(modes_sorted, peak_heights, labels = round(modes_sorted, 2), pos = 3, col = "red")
    
    # Print results
    cat("All modes:", modes_sorted, "\n")
    cat("Highest mode:", highest_mode, "with height", highest_peak_height, "\n")
  }
  
  return(modes_sorted)
  
}

#' Test if element is a number between 0 and Inf
#'
#' @param x numeric vector 
#'
#' @returns
#' @export
#'
#' @examples
is_ratio <- function(x) {
  map_lgl(x, \(y) !is.null(y) && !is.na(y) && is.numeric(y) && y > 0 && y < Inf)
}


#' Test if element is a number between 0 and 1
#'
#' @param x numeric vector 
#'
#' @returns
#' @export
#'
#' @examples
is_SILAC_ratio <- function(x) {
  map_lgl(x, \(y) !is.null(y) && !is.na(y) && is.numeric(y) && y > 0 && y < 1)
}


#' Title
#'
#' @param x numeric vector 
#' @param p fraction of values to retain 
#'
#' @returns
#' @export
#'
#' @examples
in_inner_quantile <- function(x, p = 0.99, na.rm = T) {
  
  d <- (1 - p) / 2
  
  range <- quantile(x, seq(0, 1, d), na.rm = na.rm)
  
  keep <- x > range[2] & x < range[length(range) - 1]
  
  return(keep)
}

#' Title
#'
#' @param x numeric vector 
#' @param p fraction of values to retain 
#'
#' @returns
#' @export
#'
#' @examples
in_lower_quantile <- function(x, p = 0.99, na.rm = T) {
  
  d <- 1 - p
  
  range <- quantile(x, seq(0, 1, d), na.rm = na.rm)
  
  keep <- x < range[length(range) - 1]
  
  return(keep)
}

#' Title
#'
#' @param x numeric vector 
#' @param p fraction of values to retain 
#'
#' @returns
#' @export
#'
#' @examples
in_upper_quantile <- function(x, p = 0.99, na.rm = T) {
  
  d <- 1 - p
  
  range <- quantile(x, seq(0, 1, d), na.rm = na.rm)
  
  keep <- x > range[2]
  
  return(keep)
}

#' Title
#'
#' @param x numeric vector 
#' @param n number of standard deviations 
#' @param center method to determine center of distribution 
#'
#' @returns
#' @export
#'
#' @examples
in_n_sd <- function(x, n = 3, center = "mean", na.rm = T) {
  
  x_sd <- sd(x, na.rm = na.rm)
  
  if (center == "mean") x_center <- mean(x, na.rm = na.rm)
  else if (center == "median") x_center <- median(x, na.rm = na.rm)
  else stop('<center> must be "mean" or "median".')
  
  keep <- x > x_center - n * x_sd & x < x_center + n * x_sd
  
}


#' Title
#'
#' @param x numeric vector 
#' @param na.rm allow NAs
#'
#' @returns
#' @export
#'
#' @examples
fwhm <- function(x, na.rm = T) {
  
  if (na.rm) x <- na.omit(x)
  
  d <- density(x)
  
  r <- range(which(d$y - max(d$y) / 2 > 0))
  
  l <- d$x[r[2]] - d$x[r[1]]
  
  return(l)
  
}

#' Title
#'
#' @param sequence character vector 
#' @param mods named vector c("pattern" = "replacement")
#'
#' @returns
#' @export
#'
#' @examples
str_replace_all_m <- function(sequence, 
                              mods = c("SILAC-L" = "", 
                                       "SILAC-H" = "UniMod:188")) {
  for (i in names(mods)) 
    sequence <- sequence %>% 
      str_replace_all(i, mods[i]) %>% 
      str_remove_all("\\(\\)")
  
  return(sequence)
}
