#' Count neighboring timepoints for each element in a vector
#'
#' Helper function that counts how many neighboring values (differing by 1) exist
#' for each element in a vector. Used internally for consecutive timepoint filtering.
#'
#' @param x Numeric vector of timepoint indices
#' @param na.as Value to return for NA elements (default: 0)
#'
#' @returns Numeric vector of neighbor counts
#' @export
#'
#' @examples
.n_neighbors <- function(x, na.as = 0) {
  map_dbl(x, \(y) ifelse(is.na(y),
    na.as,
    sum(abs(na.omit(x) - y) == 1)
  ))
}


#' Filter pSILAC data by minimum number of replicates per timepoint
#'
#' Removes timepoints that have fewer than the specified minimum number of valid
#' replicate measurements. This ensures sufficient statistical power for each timepoint.
#'
#' @param data A data frame containing pSILAC data
#' @param min_values_per_timepoint Minimum number of valid replicates required per timepoint (default: 2)
#' @param return_all_rows If TRUE, returns all rows with a 'fit' column; if FALSE, returns only passing rows (default: TRUE)
#' @param value_col Name of the column containing the values to filter (default: "Fraction_Old")
#' @param variables_col Name of the column containing variable identifiers (default: "variables")
#' @param observations_col Name of the column containing observation identifiers (default: "observations")
#' @param pulse_time_col Name of the column containing pulse time information (default: "Pulse_time")
#' @param group_cols Character vector of additional grouping column names (default: c())
#'
#' @returns A data frame with the same structure as input, either filtered or with a 'fit' column
#' @export
#'
#' @examples
filter_minimum_datapoints <- function(data,
                                      min_values_per_timepoint = 2,
                                      return_all_rows = T,
                                      value_col = "Fraction_Old",
                                      variables_col = "variables",
                                      observations_col = "observations",
                                      pulse_time_col = "Pulse_time",
                                      group_cols = c()) {
  # Count data points per pulse time x group_cols
  data_reps <- data %>%
    dplyr::summarise(
      n_replicates = sum(is_SILAC_ratio(!!rlang::sym(value_col))),
      .by = all_of(c(variables_col, pulse_time_col, group_cols))
    ) %>%
    # Filter timepoints based on number of data points each
    dplyr::mutate(fit1 = ifelse(n_replicates >= min_values_per_timepoint,
      "fit",
      paste0("underrepresented_n", n_replicates)
    )) %>%
    dplyr::arrange(!!rlang::sym(pulse_time_col))


  # Combine filtered data
  if (return_all_rows) {
    data_return <- data %>%
      dplyr::full_join(data_reps,
        by = c(variables_col, pulse_time_col, group_cols)
      ) %>%
      dplyr::select(-n_replicates)

    if ("fit" %in% names(data_return)) {
      data_return <- data_return %>%
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ fit,
          fit1 == "fit" ~ fit,
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1),
          .default = fit1
        )) %>%
        dplyr::select(-fit1)
    } else {
      data_return <- data_return %>%
        dplyr::rename(fit = fit1)
    }

    # Filter out data points
  } else {
    data_return <- data %>%
      dplyr::full_join(data_reps,
        by = c(variables_col, pulse_time_col, group_cols)
      ) %>%
      dplyr::filter(!is.na(!!rlang::sym(value_col))) %>%
      dplyr::filter(fit1 == "fit") %>%
      dplyr::select(-c(n_replicates, fit1))
  }

  return(data_return)
}


#' Filter pSILAC data by consecutive timepoint coverage
#'
#' Ensures that data has sufficient consecutive timepoint coverage for reliable
#' kinetic analysis. Can filter based on whether any consecutive timepoints exist
#' or require all timepoints to be consecutive.
#'
#' @param data A data frame containing pSILAC data
#' @param min_consecutive_timepoints Minimum number of consecutive timepoints required (default: 2)
#' @param only_keep_consecutive_timepoints If TRUE, requires all timepoints to be consecutive;
#'                                         if FALSE, only requires the minimum number somewhere in the series (default: TRUE)
#' @param return_all_rows If TRUE, returns all rows with a 'fit' column; if FALSE, returns only passing rows (default: TRUE)
#' @param value_col Name of the column containing the values to filter (default: "Fraction_Old")
#' @param variables_col Name of the column containing variable identifiers (default: "variables")
#' @param observations_col Name of the column containing observation identifiers (default: "observations")
#' @param pulse_time_col Name of the column containing pulse time information (default: "Pulse_time")
#' @param group_cols Character vector of additional grouping column names (default: c())
#'
#' @returns A data frame with the same structure as input, either filtered or with a 'fit' column
#' @export
#'
#' @examples
filter_consecutive_timepoints <- function(data,
                                          min_consecutive_timepoints = 2,
                                          only_keep_consecutive_timepoints = T,
                                          return_all_rows = T,
                                          value_col = "Fraction_Old",
                                          variables_col = "variables",
                                          observations_col = "observations",
                                          pulse_time_col = "Pulse_time",
                                          group_cols = c()) {
  # Check input
  if (min_consecutive_timepoints > 2 && only_keep_consecutive_timepoints) {
    stop("<min_consecutive_timepoints> can only be used with max. 2 neighbors for the option <only_keep_consecutive_timepoints>. Sorry, I was lazy :)")
  }

  # Count data points per pulse time x group_cols
  if ("fit" %in% names(data)) {
    data_reps <- data %>%
      dplyr::filter(str_detect(fit, "(fit|plateau)")) %>%
      dplyr::summarise(.by = dplyr::all_of(c(
        variables_col,
        pulse_time_col,
        group_cols
      )))
  } else {
    data_reps <- data %>%
      dplyr::summarise(.by = dplyr::all_of(c(
        variables_col,
        pulse_time_col,
        group_cols
      )))
  }

  #
  data_neighbors <- data_reps %>%
    # Remove all pulse time points with less than min_values_per_timepoint
    dplyr::arrange(!!rlang::sym(pulse_time_col)) %>%
    dplyr::mutate(
      pulse_time_order =
        match(
          !!rlang::sym(pulse_time_col),
          sort(unique(!!rlang::sym(pulse_time_col)))
        )
    ) %>%
    dplyr::mutate(
      n_neighbors = .n_neighbors(pulse_time_order),
      .by = c(variables_col, group_cols)
    )


  # Remove time points
  if (only_keep_consecutive_timepoints) {
    data_final <- data_neighbors %>%
      dplyr::mutate(
        fit1 = any(n_neighbors + 1 >= min_consecutive_timepoints),
        .by = dplyr::all_of(c(
          variables_col,
          pulse_time_col,
          group_cols
        ))
      ) %>%
      dplyr::mutate(fit1 = ifelse(fit1,
        "fit",
        paste0("non_consecutive_n", n_neighbors + 1)
      ))
  } else {
    data_final <- data_neighbors %>%
      dplyr::mutate(
        fit1 = any(n_neighbors + 1 >= min_consecutive_timepoints),
        .by = dplyr::all_of(c(
          variables_col,
          group_cols
        ))
      ) %>%
      dplyr::mutate(fit1 = ifelse(fit1,
        "fit",
        paste0("non_consecutive_n", n_neighbors + 1)
      ))
  }


  # Combine filtered data
  if (return_all_rows) {
    data_return <- data %>%
      dplyr::full_join(
        data_final %>%
          dplyr::select(-c(
            pulse_time_order,
            n_neighbors
          )),
        by = c(variables_col, pulse_time_col, group_cols)
      )

    if ("fit" %in% names(data_return)) {
      data_return <- data_return %>%
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ fit,
          fit1 == "fit" ~ fit,
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1),
          .default = fit1
        )) %>%
        dplyr::select(-fit1)
    } else {
      data_return <- data_return %>%
        dplyr::rename(fit = fit1)
    }
  } else {
    data_return <- data %>%
      dplyr::inner_join(
        data_final %>%
          dplyr::select(-c(
            pulse_time_order,
            n_neighbors
          )),
        by = c(variables_col, pulse_time_col, group_cols)
      ) %>%
      dplyr::filter(fit1 == "fit") %>%
      dplyr::select(-fit1)
  }

  return(data_return)
}


#' Title
#'
#' @param data
#' @param S2N_threshold
#' @param return_all_rows
#' @param value_col
#' @param variables_col
#' @param observations_col
#' @param pulse_time_col
#' @param group_cols
#'
#' @returns
#' @export
#'
#' @examples
filter_SILAC_plateaus_0 <- function(data,
                                    S2N_threshold = 2,
                                    return_all_rows = T,
                                    value_col = "Fraction_Old",
                                    variables_col = "variables",
                                    observations_col = "observations",
                                    pulse_time_col = "Pulse_time",
                                    group_cols = c()) {
  data_plateaus <- data %>%
    dplyr::filter(sum(!!rlang::sym(value_col) == 0) == 2,
      .by = c(variables_col, pulse_time_col, group_cols)
    ) %>%
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) %>%
    dplyr::summarise(
      mode = find_density_modes(!!rlang::sym(value_col),
        topN = 1
      ),
      .by = c(pulse_time_col, group_cols)
    )


  data_final <- left_join(data, data_plateaus,
    by = c(pulse_time_col, group_cols)
  ) %>%
    dplyr::mutate(
      mean_Fraction_Old = mean(!!rlang::sym(value_col)),
      .by = c(variables_col, pulse_time_col, group_cols)
    ) %>%
    dplyr::mutate(fit1 = ifelse(mean_Fraction_Old > mode * S2N_threshold,
      "fit",
      "plateau_0"
    ))


  if (return_all_rows) {
    if ("fit" %in% names(data_final)) {
      data_return <- data_final %>%
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ fit,
          fit1 == "fit" ~ fit,
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1),
          .default = fit1
        )) %>%
        dplyr::select(-c(mean_Fraction_Old, mode, fit1))
    } else {
      data_return <- data_final %>%
        dplyr::rename(fit = fit1)
    }
  } else {
    data_return <- data_final %>%
      # Remove
      dplyr::filter(fit1 == "fit") %>%
      dplyr::select(-c(mean_Fraction_Old, mode, fit1))
  }

  return(data_return)
}


#' Title
#'
#' @param data
#' @param S2N_threshold
#' @param return_all_rows
#' @param value_col
#' @param variables_col
#' @param observations_col
#' @param pulse_time_col
#' @param group_cols
#'
#' @returns
#' @export
#'
#' @examples
filter_SILAC_plateaus_1 <- function(data,
                                    S2N_threshold = 2,
                                    return_all_rows = T,
                                    value_col = "Fraction_Old",
                                    variables_col = "variables",
                                    observations_col = "observations",
                                    pulse_time_col = "Pulse_time",
                                    group_cols = c()) {
  data_plateaus <- data %>%
    dplyr::filter(sum(!!rlang::sym(value_col) == 1) == 2,
      .by = c(variables_col, pulse_time_col, group_cols)
    ) %>%
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) %>%
    dplyr::summarise(
      mode = find_density_modes(!!rlang::sym(value_col),
        topN = 1
      ),
      .by = c(pulse_time_col, group_cols)
    )


  data_final <- left_join(data, data_plateaus,
    by = c(pulse_time_col, group_cols)
  ) %>%
    dplyr::mutate(
      mean_Fraction_Old = mean(!!rlang::sym(value_col)),
      .by = c(variables_col, pulse_time_col, group_cols)
    ) %>%
    dplyr::mutate(fit1 = ifelse(mean_Fraction_Old < 1 - ((1 - mode) * S2N_threshold),
      "fit",
      "plateau_1"
    ))


  if (return_all_rows) {
    if ("fit" %in% names(data_final)) {
      data_return <- data_final %>%
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ fit,
          fit1 == "fit" ~ fit,
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1),
          .default = fit1
        )) %>%
        dplyr::select(-c(mean_Fraction_Old, mode, fit1))
    } else {
      data_return <- data_final %>%
        dplyr::rename(fit = fit1)
    }
  } else {
    data_return <- data_final %>%
      # Remove
      dplyr::filter(fit1 == "fit") %>%
      dplyr::select(-c(mean_Fraction_Old, mode, fit1))
  }

  return(data_return)
}


#' Title
#'
#' @param data
#' @param S2N_threshold
#' @param return_all_rows
#' @param value_col
#' @param variables_col
#' @param observations_col
#' @param pulse_time_col
#' @param group_cols
#'
#' @returns
#' @export
#'
#' @examples
plot_SILAC_plateaus_0 <- function(data,
                                  S2N_threshold = 2,
                                  return_all_rows = F,
                                  value_col = "Fraction_Old",
                                  variables_col = "variables",
                                  observations_col = "observations",
                                  pulse_time_col = "Pulse_time",
                                  group_cols = c()) {
  data_plateaus <- data %>%
    dplyr::filter(sum(!!rlang::sym(value_col) == 0) == 2,
      .by = c(variables_col, pulse_time_col, group_cols)
    ) %>%
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col)))

  data_thresholds <- data_plateaus %>%
    dplyr::summarise(
      mode = find_density_modes(!!rlang::sym(value_col), topN = 1),
      .by = c(pulse_time_col, group_cols)
    )

  p <- data_plateaus %>%
    ggplot(aes(x = !!rlang::sym(value_col))) +
    geom_density(fill = "grey") +
    theme_classic() +
    coord_cartesian(expand = F, xlim = c(0, 1)) +
    facet_grid(
      rows = vars(!!rlang::sym(pulse_time_col)),
      cols = vars(!!rlang::sym(group_cols))
    )

  return(p)
}


#' Filter pSILAC data by coefficient of variation (CV)
#'
#' Filters data points based on the coefficient of variation (CV) within replicates
#' for each combination of variable, timepoint, and grouping columns. High CV indicates
#' high variability between replicates, suggesting unreliable measurements.
#'
#' @param data A data frame containing pSILAC data
#' @param max_cv Maximum allowed coefficient of variation (default: 0.3)
#' @param min_replicates Minimum number of replicates required to calculate CV (default: 2)
#' @param return_all_rows If TRUE, returns all rows with a 'fit' column indicating filter status.
#'                        If FALSE, returns only rows that pass the filter (default: TRUE)
#' @param value_col Name of the column containing the values to filter (default: "Fraction_Old")
#' @param variables_col Name of the column containing variable identifiers (default: "variables")
#' @param observations_col Name of the column containing observation identifiers (default: "observations")
#' @param pulse_time_col Name of the column containing pulse time information (default: "Pulse_time")
#' @param group_cols Character vector of additional grouping column names (default: c())
#'
#' @returns A data frame with the same structure as input, either filtered or with a 'fit' column
#' @export
#'
#' @examples
#' \dontrun{
#' # Filter with default 30% CV threshold
#' filtered_data <- filter_by_cv(data)
#'
#' # Use stricter 20% CV threshold
#' filtered_data <- filter_by_cv(data, max_cv = 0.2)
#'
#' # Return only passing data points
#' filtered_data <- filter_by_cv(data, return_all_rows = FALSE)
#' }
filter_by_cv <- function(data,
                         max_cv = 0.3,
                         min_replicates = 2,
                         return_all_rows = T,
                         value_col = "Fraction_Old",
                         variables_col = "variables",
                         observations_col = "observations",
                         pulse_time_col = "Pulse_time",
                         group_cols = c()) {
  # Calculate CV for each variable x timepoint x group combination
  data_cv <- data %>%
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) %>%
    dplyr::summarise(
      n_replicates = dplyr::n(),
      mean_value = mean(!!rlang::sym(value_col), na.rm = TRUE),
      sd_value = sd(!!rlang::sym(value_col), na.rm = TRUE),
      .by = dplyr::all_of(c(variables_col, pulse_time_col, group_cols))
    ) %>%
    dplyr::mutate(
      cv = sd_value / mean_value,
      fit1 = dplyr::case_when(
        n_replicates < min_replicates ~ paste0("insufficient_replicates_n", n_replicates),
        is.na(cv) | is.infinite(cv) ~ "cv_undefined",
        cv <= max_cv ~ "fit",
        TRUE ~ paste0("high_cv_", round(cv, 3))
      )
    )

  # Combine filtered data
  if (return_all_rows) {
    data_return <- data %>%
      dplyr::left_join(
        data_cv %>% dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), fit1),
        by = c(variables_col, pulse_time_col, group_cols)
      )

    if ("fit" %in% names(data_return)) {
      data_return <- data_return %>%
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ fit,
          fit1 == "fit" ~ fit,
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1),
          .default = fit1
        )) %>%
        dplyr::select(-fit1)
    } else {
      data_return <- data_return %>%
        dplyr::rename(fit = fit1)
    }
  } else {
    data_return <- data %>%
      dplyr::inner_join(
        data_cv %>% dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), fit1),
        by = c(variables_col, pulse_time_col, group_cols)
      ) %>%
      dplyr::filter(fit1 == "fit") %>%
      dplyr::select(-fit1)
  }

  return(data_return)
}

#' Plot density distribution for proteins plateauing at Fraction_Old = 1
#'
#' Creates density plots to visualize the distribution of Fraction_Old values for
#' proteins identified as potentially plateauing at 1. Useful for quality control
#' and threshold selection.
#'
#' @param data A data frame containing pSILAC data
#' @param S2N_threshold Signal-to-noise threshold for detecting plateaus (default: 2)
#' @param return_all_rows If TRUE, returns all rows with a 'fit' column; if FALSE, returns only passing rows (default: FALSE)
#' @param value_col Name of the column containing the values to filter (default: "Fraction_Old")
#' @param variables_col Name of the column containing variable identifiers (default: "variables")
#' @param observations_col Name of the column containing observation identifiers (default: "observations")
#' @param pulse_time_col Name of the column containing pulse time information (default: "Pulse_time")
#' @param group_cols Character vector of additional grouping column names (default: c())
#'
#' @returns A ggplot2 object showing density distributions
#' @export
#'
#' @examples
plot_SILAC_plateaus_1 <- function(data,
                                  S2N_threshold = 2,
                                  return_all_rows = F,
                                  value_col = "Fraction_Old",
                                  variables_col = "variables",
                                  observations_col = "observations",
                                  pulse_time_col = "Pulse_time",
                                  group_cols = c()) {
  data_plateaus <- data %>%
    dplyr::filter(sum(!!rlang::sym(value_col) == 1) == 2,
      .by = c(variables_col, pulse_time_col, group_cols)
    ) %>%
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col)))

  data_thresholds <- data_plateaus %>%
    dplyr::summarise(
      mode = find_density_modes(!!rlang::sym(value_col), topN = 1),
      .by = c(pulse_time_col, group_cols)
    )

  p <- data_plateaus %>%
    ggplot(aes(x = !!rlang::sym(value_col))) +
    geom_density(fill = "grey") +
    theme_classic() +
    coord_cartesian(expand = F, xlim = c(0, 1)) +
    facet_grid(
      rows = vars(!!rlang::sym(pulse_time_col)),
      cols = vars(!!rlang::sym(group_cols))
    )

  return(p)
}



#' Filter pSILAC data by removing specific timepoints with increasing values
#'
#' Identifies and filters out specific timepoints where the mean Fraction_Old value
#' increases compared to the previous timepoint. Only the problematic timepoint is removed,
#' not all timepoints. In pSILAC experiments, Fraction_Old should decrease over time as
#' old proteins are replaced by newly synthesized ones.
#'
#' @param data A data frame containing pSILAC data
#' @param max_increase Maximum allowed increase in mean value between consecutive timepoints
#'                     (default: 0, meaning no increase is allowed)
#' @param recursive If TRUE, recursively re-checks for increases after removing timepoints
#'                  (e.g., if timepoint 2 is removed, checks if timepoint 3 increases relative to timepoint 1)
#'                  (default: TRUE)
#' @param return_all_rows If TRUE, returns all rows with a 'fit' column indicating filter status.
#'                        If FALSE, returns only rows that pass the filter (default: TRUE)
#' @param value_col Name of the column containing the values to filter (default: "Fraction_Old")
#' @param variables_col Name of the column containing variable identifiers (default: "variables")
#' @param observations_col Name of the column containing observation identifiers (default: "observations")
#' @param pulse_time_col Name of the column containing pulse time information (default: "Pulse_time")
#' @param group_cols Character vector of additional grouping column names (default: c())
#'
#' @returns A data frame with the same structure as input, either filtered or with a 'fit' column
#' @export
#'
#' @examples
#' \dontrun{
#' # Filter timepoints with any increase, recursively (default)
#' filtered_data <- filter_increasing_timepoints(data)
#'
#' # Non-recursive: only check once
#' filtered_data <- filter_increasing_timepoints(data, recursive = FALSE)
#'
#' # Allow small increases up to 0.05
#' filtered_data <- filter_increasing_timepoints(data, max_increase = 0.05)
#' }
filter_increasing_trends <- function(data,
                                     max_increase = 0,
                                     recursive = TRUE,
                                     return_all_rows = T,
                                     value_col = "Fraction_Old",
                                     variables_col = "variables",
                                     observations_col = "observations",
                                     pulse_time_col = "Pulse_time",
                                     group_cols = c()) {
  data_return <- data %>%
    dplyr::arrange(
      !!rlang::sym(variables_col),
      !!rlang::sym(group_cols),
      !!rlang::sym(pulse_time_col)
    ) %>%
    dplyr::mutate(
      Fraction_mean = mean(!!rlang::sym(value_col)),
      .by = dplyr::all_of(c(
        variables_col,
        group_cols,
        pulse_time_col
      ))
    ) %>%
    dplyr::mutate(
      increase = any(
        diff(c(1, unique(Fraction_mean))) > max_increase
      ),
      .by = dplyr::all_of(c(variables_col, group_cols))
    ) %>%
    dplyr::mutate(fit = dplyr::case_when(
      fit == "fit" & increase ~ "increasing",
      fit != "fit" & increase ~ paste0(fit, "|increasing"),
      .default = fit
    )) %>%
    dplyr::select(-c(Fraction_mean, increase))

  if (!return_all_rows) data_return <- dplyr::filter(data_return, fit == "fit")

  return(data_return)
}

# filter_increasing_timepoints <- function(data,
#                                          max_increase = 0,
#                                          recursive = TRUE,
#                                          return_all_rows = T,
#                                          value_col = "Fraction_Old",
#                                          variables_col = "variables",
#                                          observations_col = "observations",
#                                          pulse_time_col = "Pulse_time",
#                                          group_cols = c()) {
#   # Helper function to identify increasing timepoints
#   identify_increases <- function(df) {
#     df %>%
#       dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) %>%
#       dplyr::arrange(!!rlang::sym(pulse_time_col)) %>%
#       dplyr::mutate(
#         mean_value = mean(!!rlang::sym(value_col), na.rm = TRUE),
#         .by = dplyr::all_of(c(variables_col, pulse_time_col, group_cols))
#       ) %>%
#       dplyr::arrange(!!rlang::sym(pulse_time_col)) %>%
#       dplyr::mutate(
#         # Get previous timepoint's mean value
#         prev_mean = {
#           unique_data <- data.frame(
#             tp = unique(!!rlang::sym(pulse_time_col)),
#             val = unique(mean_value)
#           ) %>% dplyr::arrange(tp)
#
#           current_tp <- unique(!!rlang::sym(pulse_time_col))[1]
#           current_idx <- match(current_tp, unique_data$tp)
#
#           if (current_idx > 1) {
#             unique_data$val[current_idx - 1]
#           } else {
#             1 # Baseline for first timepoint
#           }
#         },
#         .by = dplyr::all_of(c(variables_col, pulse_time_col, group_cols))
#       ) %>%
#       dplyr::mutate(
#         is_increasing = mean_value - prev_mean > max_increase
#       ) %>%
#       dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), is_increasing) %>%
#       dplyr::distinct()
#   }
#
#   # Iteratively identify and mark increasing timepoints
#   data_working <- data
#   all_flagged_timepoints <- data.frame()
#
#   repeat {
#     # Identify increases in current working dataset
#     increases <- identify_increases(data_working)
#
#     # Find timepoints to flag
#     new_flags <- increases %>%
#       dplyr::filter(is_increasing) %>%
#       dplyr::mutate(fit1 = "increasing_timepoint") %>%
#       dplyr::select(-is_increasing)
#
#     # If no new increases found, stop
#     if (nrow(new_flags) == 0) break
#
#     # Add to all flagged timepoints
#     all_flagged_timepoints <- dplyr::bind_rows(all_flagged_timepoints, new_flags)
#
#     # If not recursive, stop after first iteration
#     if (!recursive) break
#
#     # Remove flagged timepoints from working dataset for next iteration
#     data_working <- data_working %>%
#       dplyr::anti_join(new_flags %>% dplyr::select(-fit1),
#         by = c(variables_col, pulse_time_col, group_cols)
#       )
#   }
#
#   # Combine filtered data
#   if (return_all_rows) {
#     if (nrow(all_flagged_timepoints) > 0) {
#       data_return <- data %>%
#         dplyr::left_join(all_flagged_timepoints, by = c(variables_col, pulse_time_col, group_cols))
#
#       if ("fit" %in% names(data_return)) {
#         data_return <- data_return %>%
#           dplyr::mutate(fit = dplyr::case_when(
#             is.na(fit1) ~ fit,
#             fit1 == "increasing_timepoint" & fit == "fit" ~ fit1,
#             fit1 == "increasing_timepoint" & fit != "fit" ~ paste0(fit, "|", fit1),
#             .default = fit
#           )) %>%
#           dplyr::select(-fit1)
#       } else {
#         data_return <- data_return %>%
#           dplyr::mutate(fit1 = ifelse(is.na(fit1), "fit", fit1)) %>%
#           dplyr::rename(fit = fit1)
#       }
#     } else {
#       data_return <- data
#       if (!"fit" %in% names(data_return)) {
#         data_return <- data_return %>% dplyr::mutate(fit = "fit")
#       }
#     }
#   } else {
#     if (nrow(all_flagged_timepoints) > 0) {
#       data_return <- data %>%
#         dplyr::anti_join(all_flagged_timepoints %>% dplyr::select(-fit1),
#           by = c(variables_col, pulse_time_col, group_cols)
#         )
#     } else {
#       data_return <- data
#     }
#   }
#
#   return(data_return)
# }
#
#
# filter_increasing_timepoints <- function(data,
#                                          max_increase = 0,
#                                          recursive = TRUE,
#                                          return_all_rows = T,
#                                          value_col = "Fraction_Old",
#                                          variables_col = "variables",
#                                          observations_col = "observations",
#                                          pulse_time_col = "Pulse_time",
#                                          group_cols = c()) {
#
#   data_return <- data %>%
#     dplyr::arrange(!!rlang::sym(variables_col),
#                    !!rlang::sym(group_cols),
#                    !!rlang::sym(pulse_time_col)) %>%
#     dplyr::mutate(Fraction_mean = mean(!!rlang::sym(value_col)),
#                   .by = dplyr::all_of(c(variables_col,
#                                         group_cols,
#                                         pulse_time_col))) %>%
#     dplyr::mutate(increase = any(
#       diff(c(1, unique(Fraction_mean))) > max_increase),
#       .by = dplyr::all_of(c(variables_col, group_cols))) %>%
#     dplyr::mutate(fit = dplyr::case_when(
#       fit == "fit" & increase ~ "increasing",
#       fit != "fit" & increase ~ paste0(fit, "|increasing"),
#       .default = fit)) %>%
#     dplyr::select(-c(Fraction_mean, increase))
#
#   if (!return_all_rows) data_return <- dplyr::filter(data_return, fit == "fit")
#
#   return(data_return)
#
# }
#


#' Filter individual pSILAC data point outliers with increasing values
#'
#' Identifies and filters out individual replicate data points where the value
#' is higher than the previous timepoint's mean value. This removes outlier replicates
#' (e.g., 1 out of 3 replicates) while keeping other replicates at the same timepoint.
#' In pSILAC experiments, Fraction_Old should decrease over time.
#'
#' @param data A data frame containing pSILAC data
#' @param max_increase Maximum allowed increase relative to previous timepoint's mean
#'                     (default: 0, meaning no increase is allowed)
#' @param return_all_rows If TRUE, returns all rows with a 'fit' column indicating filter status.
#'                        If FALSE, returns only rows that pass the filter (default: TRUE)
#' @param value_col Name of the column containing the values to filter (default: "Fraction_Old")
#' @param variables_col Name of the column containing variable identifiers (default: "variables")
#' @param observations_col Name of the column containing observation identifiers (default: "observations")
#' @param pulse_time_col Name of the column containing pulse time information (default: "Pulse_time")
#' @param group_cols Character vector of additional grouping column names (default: c())
#'
#' @returns A data frame with the same structure as input, either filtered or with a 'fit' column
#' @export
#'
#' @examples
#' \dontrun{
#' # Filter individual outlier replicates with any increase (default)
#' filtered_data <- filter_increasing_datapoints(data)
#'
#' # Allow small increases up to 0.05
#' filtered_data <- filter_increasing_datapoints(data, max_increase = 0.05)
#'
#' # Return only passing data points
#' filtered_data <- filter_increasing_datapoints(data, return_all_rows = FALSE)
#' }
# filter_increasing_datapoints <- function(data,
#' #                                         max_increase = 0,
#' #                                         return_all_rows = T,
#' #                                         value_col = "Fraction_Old",
#' #                                         variables_col = "variables",
#' #                                         observations_col = "observations",
#' #                                         pulse_time_col = "Pulse_time",
#' #                                         group_cols = c()) {
#' #  # Calculate mean values per timepoint
#' #  data_with_means <- data %>%
#' #    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) %>%
#' #    dplyr::arrange(!!rlang::sym(pulse_time_col)) %>%
#' #    dplyr::mutate(
#' #      mean_value = mean(!!rlang::sym(value_col), na.rm = TRUE),
#' #      .by = dplyr::all_of(c(variables_col, pulse_time_col, group_cols))
#' #    )
#' #  # Get previous timepoint's mean for each data point
#' #  data_with_prev <- data_with_means %>%
#' #    dplyr::arrange(!!rlang::sym(pulse_time_col)) %>%
#' #    dplyr::mutate(
#' #      prev_mean = {
#' #        # Get unique timepoint-mean pairs
#' #        unique_data <- data.frame(
#' #          tp = unique(!!rlang::sym(pulse_time_col)),
#' #          val = unique(mean_value)
#' #        ) %>% dplyr::arrange(tp)
#' #        current_tp <- !!rlang::sym(pulse_time_col)
#' #        current_idx <- match(current_tp, unique_data$tp)
#' #        if (current_idx > 1) {
#' #          unique_data$val[current_idx - 1]
#' #        } else {
#' #          1 # Baseline for first timepoint
#' #        }
#' #      },
#' #      .by = dplyr::all_of(c(variables_col, group_cols))
#' #    )
#' #  # Identify individual data points that are higher than previous timepoint's mean
#' #  data_flagged <- data_with_prev %>%
#' #    dplyr::mutate(
#' #      is_increasing = !!rlang::sym(value_col) - prev_mean > max_increase,
#' #      fit1 = ifelse(is_increasing, "increasing_datapoint", "fit")
#' #    ) %>%
#' #    dplyr::select(dplyr::all_of(c(variables_col, observations_col, pulse_time_col, group_cols)), fit1)
#' #  # Combine filtered data
#' #  if (return_all_rows) {
#' #    data_return <- data %>%
#' #      dplyr::left_join(data_flagged, by = c(variables_col, observations_col, pulse_time_col, group_cols))
#' #    if ("fit" %in% names(data_return)) {
#' #      data_return <- data_return %>%
#' #        dplyr::mutate(fit = dplyr::case_when(
#' #          is.na(fit1) ~ fit,
#' #          fit1 == "fit" ~ fit,
#' #          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1),
#' #          .default = fit1
#' #        )) %>%
#' #        dplyr::select(-fit1)
#' #    } else {
#' #      data_return <- data_return %>%
#' #        dplyr::mutate(fit1 = ifelse(is.na(fit1), "fit", fit1)) %>%
#' #        dplyr::rename(fit = fit1)
#' #    }
#' #  } else {
#' #    data_return <- data %>%
#' #      dplyr::inner_join(data_flagged, by = c(variables_col, observations_col, pulse_time_col, group_cols)) %>%
#' #      dplyr::filter(fit1 == "fit") %>%
#' #      dplyr::select(-fit1)
#' #  }
#' #  return(data_return)
# }

#' Filter pSILAC data to keep only timepoints closest to Fraction_Old = 0.5
#'
#' Selects the 2 timepoints per variable and group whose mean Fraction_Old values
#' are closest to 0.5. This region often provides the most informative data for
#' kinetic analysis as it has the best signal-to-noise ratio for rate determination.
#'
#' @param data A data frame containing pSILAC data
#' @param n_timepoints Number of timepoints to keep (default: 2)
#' @param target_value Target Fraction_Old value to find closest timepoints to (default: 0.5)
#' @param return_all_rows If TRUE, returns all rows with a 'fit' column indicating filter status.
#'                        If FALSE, returns only rows that pass the filter (default: TRUE)
#' @param value_col Name of the column containing the values to filter (default: "Fraction_Old")
#' @param variables_col Name of the column containing variable identifiers (default: "variables")
#' @param observations_col Name of the column containing observation identifiers (default: "observations")
#' @param pulse_time_col Name of the column containing pulse time information (default: "Pulse_time")
#' @param group_cols Character vector of additional grouping column names (default: c())
#'
#' @returns A data frame with the same structure as input, either filtered or with a 'fit' column
#' @export
#'
#' @examples
#' \dontrun{
#' # Keep 2 timepoints closest to 0.5 (default)
#' filtered_data <- filter_closest_to_midpoint(data)
#'
#' # Keep 3 timepoints closest to 0.5
#' filtered_data <- filter_closest_to_midpoint(data, n_timepoints = 3)
#'
#' # Target a different value (e.g., 0.3)
#' filtered_data <- filter_closest_to_midpoint(data, target_value = 0.3)
#' }
filter_closest_to_midpoint <- function(data,
                                       n_timepoints = 2,
                                       target_value = 0.5,
                                       return_all_rows = T,
                                       value_col = "Fraction_Old",
                                       variables_col = "variables",
                                       observations_col = "observations",
                                       pulse_time_col = "Pulse_time",
                                       group_cols = c()) {
  # Calculate mean values per timepoint and distance from target
  data_distances <- data %>%
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) %>%
    dplyr::mutate(
      mean_value = mean(!!rlang::sym(value_col), na.rm = TRUE),
      .by = dplyr::all_of(c(variables_col, pulse_time_col, group_cols))
    ) %>%
    dplyr::mutate(
      distance_from_target = abs(mean_value - target_value),
      .by = dplyr::all_of(c(variables_col, pulse_time_col, group_cols))
    ) %>%
    dplyr::select(
      dplyr::all_of(c(variables_col, pulse_time_col, group_cols)),
      mean_value, distance_from_target
    ) %>%
    dplyr::distinct()

  # Select n closest timepoints per variable and group
  data_selected <- data_distances %>%
    dplyr::arrange(distance_from_target) %>%
    dplyr::mutate(
      rank = dplyr::row_number(),
      .by = dplyr::all_of(c(variables_col, group_cols))
    ) %>%
    dplyr::mutate(
      fit1 = ifelse(rank <= n_timepoints, "fit", "not_closest_to_midpoint")
    ) %>%
    dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), fit1)

  # Combine filtered data
  if (return_all_rows) {
    data_return <- data %>%
      dplyr::left_join(data_selected, by = c(variables_col, pulse_time_col, group_cols))

    if ("fit" %in% names(data_return)) {
      data_return <- data_return %>%
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ fit,
          fit1 == "fit" ~ fit,
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1),
          .default = fit1
        )) %>%
        dplyr::select(-fit1)
    } else {
      data_return <- data_return %>%
        dplyr::mutate(fit1 = ifelse(is.na(fit1), "fit", fit1)) %>%
        dplyr::rename(fit = fit1)
    }
  } else {
    data_return <- data %>%
      dplyr::inner_join(data_selected, by = c(variables_col, pulse_time_col, group_cols)) %>%
      dplyr::filter(fit1 == "fit") %>%
      dplyr::select(-fit1)
  }

  return(data_return)
}


#' Filter pSILAC data to keep one timepoint above and one below 0.5
#'
#' Selects exactly 2 timepoints per variable and group: the closest timepoint with
#' mean Fraction_Old above 0.5 and the closest below 0.5. This ensures balanced
#' sampling around the midpoint for kinetic analysis.
#'
#' @param data A data frame containing pSILAC data
#' @param target_value Target Fraction_Old value to bracket (default: 0.5)
#' @param allow_incomplete_brackets If TRUE, allows incomplete brackets (all timepoints above or below target)
#'                                  to still be marked as "fit" (default: FALSE)
#' @param return_all_rows If TRUE, returns all rows with a 'fit' column indicating filter status.
#'                        If FALSE, returns only rows that pass the filter (default: TRUE)
#' @param value_col Name of the column containing the values to filter (default: "Fraction_Old")
#' @param variables_col Name of the column containing variable identifiers (default: "variables")
#' @param observations_col Name of the column containing observation identifiers (default: "observations")
#' @param pulse_time_col Name of the column containing pulse time information (default: "Pulse_time")
#' @param group_cols Character vector of additional grouping column names (default: c())
#'
#' @returns A data frame with the same structure as input, either filtered or with a 'fit' column
#' @export
#'
#' @examples
#' \dontrun{
#' # Keep 1 timepoint above and 1 below 0.5 (default)
#' filtered_data <- filter_bracket_midpoint(data)
#'
#' # Bracket a different value (e.g., 0.3)
#' filtered_data <- filter_bracket_midpoint(data, target_value = 0.3)
#'
#' # Return only passing data points
#' filtered_data <- filter_bracket_midpoint(data, return_all_rows = FALSE)
#'
#' # Allow incomplete brackets (all above or all below target)
#' filtered_data <- filter_bracket_midpoint(data, allow_incomplete_brackets = TRUE)
#' }
filter_bracket_midpoint <- function(data,
                                    target_value = 0.5,
                                    allow_incomplete_brackets = FALSE,
                                    return_all_rows = T,
                                    value_col = "Fraction_Old",
                                    variables_col = "variables",
                                    observations_col = "observations",
                                    pulse_time_col = "Pulse_time",
                                    group_cols = c()) {
  # Calculate mean values per timepoint
  data_means <- data %>%
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) %>%
    dplyr::mutate(
      mean_value = mean(!!rlang::sym(value_col), na.rm = TRUE),
      .by = dplyr::all_of(c(variables_col, pulse_time_col, group_cols))
    ) %>%
    dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), mean_value) %>%
    dplyr::distinct()

  # Find closest timepoint above target
  data_above <- data_means %>%
    dplyr::filter(mean_value > target_value) %>%
    dplyr::mutate(
      distance_from_target = abs(mean_value - target_value)
    ) %>%
    dplyr::arrange(distance_from_target) %>%
    dplyr::slice(1, .by = dplyr::all_of(c(variables_col, group_cols))) %>%
    dplyr::mutate(fit1 = "fit") %>%
    dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), fit1)

  # Find closest timepoint below target
  data_below <- data_means %>%
    dplyr::filter(mean_value < target_value) %>%
    dplyr::mutate(
      distance_from_target = abs(mean_value - target_value)
    ) %>%
    dplyr::arrange(distance_from_target) %>%
    dplyr::slice(1, .by = dplyr::all_of(c(variables_col, group_cols))) %>%
    dplyr::mutate(fit1 = "fit") %>%
    dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), fit1)

  # Combine above and below selections
  data_selected <- dplyr::bind_rows(data_above, data_below)

  # Check if both above and below were found for each variable/group
  data_check <- data_selected %>%
    dplyr::summarise(
      n_selected = dplyr::n(),
      .by = dplyr::all_of(c(variables_col, group_cols))
    )

  # Mark variables/groups that don't have both above and below
  incomplete_groups <- data_check %>%
    dplyr::filter(n_selected < 2)

  if (nrow(incomplete_groups) > 0 && !allow_incomplete_brackets) {
    # Add flag for incomplete bracketing (only if not allowing incomplete brackets)
    data_selected <- data_selected %>%
      dplyr::left_join(
        incomplete_groups %>% dplyr::mutate(incomplete = TRUE),
        by = c(variables_col, group_cols)
      ) %>%
      dplyr::mutate(
        fit1 = ifelse(!is.na(incomplete), "incomplete_bracket", fit1)
      ) %>%
      dplyr::select(-incomplete, -n_selected)
  }

  # Combine filtered data
  if (return_all_rows) {
    data_return <- data %>%
      dplyr::left_join(data_selected, by = c(variables_col, pulse_time_col, group_cols))

    if ("fit" %in% names(data_return)) {
      data_return <- data_return %>%
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ ifelse(fit == "fit", "not_bracketing_midpoint", fit),
          fit1 == "fit" ~ fit,
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1),
          .default = fit1
        )) %>%
        dplyr::select(-fit1)
    } else {
      data_return <- data_return %>%
        dplyr::mutate(fit1 = ifelse(is.na(fit1), "not_bracketing_midpoint", fit1)) %>%
        dplyr::rename(fit = fit1)
    }
  } else {
    data_return <- data %>%
      dplyr::inner_join(data_selected, by = c(variables_col, pulse_time_col, group_cols)) %>%
      dplyr::filter(fit1 == "fit") %>%
      dplyr::select(-fit1)
  }

  return(data_return)
}


#' Filter pSILAC data by bracketing or closest timepoint to target value
#'
#' This function combines bracketing and closest-point selection strategies.
#' It first attempts to find one timepoint above and one below the target value
#' (bracketing). If bracketing is not possible (all timepoints are above or below
#' the target), it falls back to selecting the single closest timepoint to the
#' target value.
#'
#' @param data A data frame containing pSILAC data
#' @param target_value Target Fraction_Old value to bracket or approach (default: 0.5)
#' @param return_all_rows If TRUE, returns all rows with a 'fit' column indicating filter status.
#'                        If FALSE, returns only rows that pass the filter (default: TRUE)
#' @param value_col Name of the column containing the values to filter (default: "Fraction_Old")
#' @param variables_col Name of the column containing variable identifiers (default: "variables")
#' @param observations_col Name of the column containing observation identifiers (default: "observations")
#' @param pulse_time_col Name of the column containing pulse time information (default: "Pulse_time")
#' @param group_cols Character vector of additional grouping column names (default: c())
#'
#' @returns A data frame with the same structure as input, either filtered or with a 'fit' column
#' @export
#'
#' @examples
#' \dontrun{
#' # Bracket 0.5 or use closest timepoint if bracketing not possible
#' filtered_data <- filter_bracket_or_closest(data)
#'
#' # Target a different value (e.g., 0.3)
#' filtered_data <- filter_bracket_or_closest(data, target_value = 0.3)
#'
#' # Return only passing data points
#' filtered_data <- filter_bracket_or_closest(data, return_all_rows = FALSE)
#' }
filter_bracket_or_closest <- function(data,
                                      target_value = 0.5,
                                      return_all_rows = T,
                                      value_col = "Fraction_Old",
                                      variables_col = "variables",
                                      observations_col = "observations",
                                      pulse_time_col = "Pulse_time",
                                      group_cols = c()) {
  # Calculate mean values per timepoint
  data_means <- data %>%
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) %>%
    dplyr::mutate(
      mean_value = mean(!!rlang::sym(value_col), na.rm = TRUE),
      .by = dplyr::all_of(c(variables_col, pulse_time_col, group_cols))
    ) %>%
    dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), mean_value) %>%
    dplyr::distinct()

  # Find closest timepoint above target
  data_above <- data_means %>%
    dplyr::filter(mean_value > target_value) %>%
    dplyr::mutate(
      distance_from_target = abs(mean_value - target_value)
    ) %>%
    dplyr::arrange(distance_from_target) %>%
    dplyr::slice(1, .by = dplyr::all_of(c(variables_col, group_cols))) %>%
    dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), distance_from_target)

  # Find closest timepoint below target
  data_below <- data_means %>%
    dplyr::filter(mean_value < target_value) %>%
    dplyr::mutate(
      distance_from_target = abs(mean_value - target_value)
    ) %>%
    dplyr::arrange(distance_from_target) %>%
    dplyr::slice(1, .by = dplyr::all_of(c(variables_col, group_cols))) %>%
    dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), distance_from_target)

  # Combine above and below selections
  data_bracketing <- dplyr::bind_rows(data_above, data_below)

  # Check if both above and below were found for each variable/group
  data_check <- data_bracketing %>%
    dplyr::summarise(
      n_selected = dplyr::n(),
      .by = dplyr::all_of(c(variables_col, group_cols))
    )

  # Identify groups with complete brackets (both above and below)
  complete_groups <- data_check %>%
    dplyr::filter(n_selected == 2) %>%
    dplyr::select(dplyr::all_of(c(variables_col, group_cols)))

  # Identify groups with incomplete brackets (only above or only below)
  incomplete_groups <- data_check %>%
    dplyr::filter(n_selected < 2) %>%
    dplyr::select(dplyr::all_of(c(variables_col, group_cols)))

  # For complete brackets, use both timepoints
  data_complete <- data_bracketing %>%
    dplyr::inner_join(complete_groups, by = c(variables_col, group_cols)) %>%
    dplyr::mutate(fit1 = "fit") %>%
    dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), fit1)

  # For incomplete brackets, find the single closest timepoint
  if (nrow(incomplete_groups) > 0) {
    data_incomplete <- data_means %>%
      dplyr::inner_join(incomplete_groups, by = c(variables_col, group_cols)) %>%
      dplyr::mutate(
        distance_from_target = abs(mean_value - target_value)
      ) %>%
      dplyr::arrange(distance_from_target) %>%
      dplyr::slice(1, .by = dplyr::all_of(c(variables_col, group_cols))) %>%
      dplyr::mutate(fit1 = "fit") %>%
      dplyr::select(dplyr::all_of(c(variables_col, pulse_time_col, group_cols)), fit1)

    # Combine complete and incomplete selections
    data_selected <- dplyr::bind_rows(data_complete, data_incomplete)
  } else {
    data_selected <- data_complete
  }

  # Combine filtered data
  if (return_all_rows) {
    data_return <- data %>%
      dplyr::left_join(data_selected, by = c(variables_col, pulse_time_col, group_cols))

    if ("fit" %in% names(data_return)) {
      data_return <- data_return %>%
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ ifelse(fit == "fit", "not_bracket_or_closest", fit),
          fit1 == "fit" ~ fit,
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1),
          .default = fit1
        )) %>%
        dplyr::select(-fit1)
    } else {
      data_return <- data_return %>%
        dplyr::mutate(fit1 = ifelse(is.na(fit1), "not_bracket_or_closest", fit1)) %>%
        dplyr::rename(fit = fit1)
    }
  } else {
    data_return <- data %>%
      dplyr::inner_join(data_selected, by = c(variables_col, pulse_time_col, group_cols)) %>%
      dplyr::filter(fit1 == "fit") %>%
      dplyr::select(-fit1)
  }

  return(data_return)
}
