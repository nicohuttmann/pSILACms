#' Title
#'
#' @param x 
#' @param na.as 
#'
#' @returns
#' @export
#'
#' @examples
.n_neighbors <- function(x, na.as = 0) {
  map_dbl(x, \(y) ifelse(is.na(y), 
                         na.as,  
                         sum(abs(na.omit(x) - y) == 1)))
}


#' Title
#'
#' @param data 
#' @param min_values_per_timepoint 
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
    dplyr::summarise(n_replicates = sum(is_SILAC_ratio(!!rlang::sym(value_col))), 
                     .by = all_of(c(variables_col, pulse_time_col, group_cols))) %>% 
    # Filter timepoints based on number of data points each
    dplyr::mutate(fit1 = ifelse(n_replicates >= min_values_per_timepoint, 
                                "fit", 
                                paste0("underrepresented_n", n_replicates))) %>% 
    dplyr::arrange(!!rlang::sym(pulse_time_col)) 
  
  
  # Combine filtered data
  if (return_all_rows) {
    data_return <- data %>% 
      dplyr::full_join(data_reps, 
                       by = c(variables_col, pulse_time_col, group_cols)) %>% 
      dplyr::select(-n_replicates)
    
    if ("fit" %in% names(data_return)) {
      data_return <- data_return %>% 
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ fit, 
          fit1 == "fit" ~ fit, 
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1), 
          .default = fit1)) %>% 
        dplyr::select(-fit1)
    } else {
      data_return <- data_return %>% 
        dplyr::rename(fit = fit1)
    }
    
    # Filter out data points 
  } else {
    data_return <- data %>% 
      dplyr::full_join(data_reps, 
                       by = c(variables_col, pulse_time_col, group_cols)) %>% 
      dplyr::filter(!is.na(!!rlang::sym(value_col))) %>% 
      dplyr::filter(fit1 == "fit") %>% 
      dplyr::select(-c(n_replicates, fit1))
  }
  
  return(data_return)
  
}


#' Title
#'
#' @param data 
#' @param min_values_per_timepoint 
#' @param only_keep_complete_timepoints 
#' @param min_consecutive_timepoints 
#' @param only_keep_consecutive_timepoints 
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
  if (min_consecutive_timepoints > 2 && only_keep_consecutive_timepoints) 
    stop("<min_consecutive_timepoints> can only be used with max. 2 neighbors for the option <only_keep_consecutive_timepoints>. Sorry, I was lazy :)")
  
  # Count data points per pulse time x group_cols
  if ("fit" %in% names(data))
    data_reps <- data %>% 
      dplyr::filter(str_detect(fit, "(fit|plateau)")) %>% 
      dplyr::summarise(.by = dplyr::all_of(c(variables_col, 
                                             pulse_time_col, 
                                             group_cols)))
  else 
    data_reps <- data %>% 
      dplyr::summarise(.by = dplyr::all_of(c(variables_col, 
                                             pulse_time_col, 
                                             group_cols)))
  
  # 
  data_neighbors <- data_reps %>% 
    # Remove all pulse time points with less than min_values_per_timepoint
    dplyr::arrange(!!rlang::sym(pulse_time_col)) %>% 
    dplyr::mutate(pulse_time_order = 
                    match(!!rlang::sym(pulse_time_col), 
                          sort(unique(!!rlang::sym(pulse_time_col))))) %>% 
    dplyr::mutate(n_neighbors = .n_neighbors(pulse_time_order), 
                  .by = c(variables_col, group_cols))
  
  
  # Remove time points 
  if (only_keep_consecutive_timepoints) {
    data_final <- data_neighbors %>% 
      dplyr::mutate(fit1 = any(n_neighbors + 1 >= min_consecutive_timepoints), 
                    .by = dplyr::all_of(c(variables_col, 
                                          pulse_time_col, 
                                          group_cols))) %>% 
      dplyr::mutate(fit1 = ifelse(fit1, 
                                  "fit", 
                                  paste0("non_consecutive_n", n_neighbors + 1)))
  } else {
    data_final <- data_neighbors %>% 
      dplyr::mutate(fit1 = any(n_neighbors + 1 >= min_consecutive_timepoints), 
                    .by = dplyr::all_of(c(variables_col, 
                                          group_cols))) %>% 
      dplyr::mutate(fit1 = ifelse(fit1, 
                                  "fit", 
                                  paste0("non_consecutive_n", n_neighbors + 1)))
  }
  
  
  # Combine filtered data
  if (return_all_rows) {
    data_return <- data %>% 
      dplyr::full_join(data_final %>% 
                         dplyr::select(-c(pulse_time_order, 
                                          n_neighbors)), 
                       by = c(variables_col, pulse_time_col, group_cols)) 
    
    if ("fit" %in% names(data_return)) {
      data_return <- data_return %>% 
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ fit, 
          fit1 == "fit" ~ fit, 
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1), 
          .default = fit1)) %>% 
        dplyr::select(-fit1)
    } else {
      data_return <- data_return %>% 
        dplyr::rename(fit = fit1)
    }
  } else {
    data_return <- data %>% 
      dplyr::inner_join(data_final %>% 
                          dplyr::select(-c(pulse_time_order, 
                                           n_neighbors)), 
                        by = c(variables_col, pulse_time_col, group_cols)) %>% 
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
                  .by = c(variables_col, pulse_time_col, group_cols))  %>% 
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) %>% 
    dplyr::summarise(mode = find_density_modes(!!rlang::sym(value_col), 
                                               topN = 1), 
                     .by = c(pulse_time_col, group_cols)) 
  
  
  data_final <- left_join(data, data_plateaus,  
                          by = c(pulse_time_col, group_cols)) %>% 
    dplyr::mutate(mean_Fraction_Old = mean(!!rlang::sym(value_col)), 
                  .by = c(variables_col, pulse_time_col, group_cols)) %>% 
    dplyr::mutate(fit1 = ifelse(mean_Fraction_Old > mode * S2N_threshold, 
                                "fit", 
                                "plateau_0"))
  
  
  if (return_all_rows) {
    
    if ("fit" %in% names(data_final)) {
      data_return <- data_final %>% 
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ fit, 
          fit1 == "fit" ~ fit, 
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1), 
          .default = fit1)) %>% 
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
                  .by = c(variables_col, pulse_time_col, group_cols))  %>% 
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) %>% 
    dplyr::summarise(mode = find_density_modes(!!rlang::sym(value_col), 
                                               topN = 1), 
                     .by = c(pulse_time_col, group_cols)) 
  
  
  data_final <- left_join(data, data_plateaus,  
                          by = c(pulse_time_col, group_cols)) %>% 
    dplyr::mutate(mean_Fraction_Old = mean(!!rlang::sym(value_col)), 
                  .by = c(variables_col, pulse_time_col, group_cols)) %>% 
    dplyr::mutate(fit1 = ifelse(mean_Fraction_Old < 1 - ((1 - mode) * S2N_threshold), 
                                "fit", 
                                "plateau_1"))
  
  
  if (return_all_rows) {
    
    if ("fit" %in% names(data_final)) {
      data_return <- data_final %>% 
        dplyr::mutate(fit = dplyr::case_when(
          is.na(fit1) ~ fit, 
          fit1 == "fit" ~ fit, 
          fit1 != "fit" & fit != "fit" ~ paste0(fit, "|", fit1), 
          .default = fit1)) %>% 
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
                  .by = c(variables_col, pulse_time_col, group_cols)) %>% 
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) 
  
  data_thresholds <- data_plateaus %>% 
    dplyr::summarise(mode = find_density_modes(!!rlang::sym(value_col), topN = 1), 
                     .by = c(pulse_time_col, group_cols)) 
  
  p <- data_plateaus %>% 
    ggplot(aes(x = !!rlang::sym(value_col))) + 
    geom_density(fill = "grey") + 
    theme_classic() + 
    coord_cartesian(expand = F, xlim = c(0, 1)) + 
    facet_grid(rows = vars(!!rlang::sym(pulse_time_col)), 
               cols = vars(!!rlang::sym(group_cols)))
  
  return(p)
  
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
                  .by = c(variables_col, pulse_time_col, group_cols)) %>% 
    dplyr::filter(is_SILAC_ratio(!!rlang::sym(value_col))) 
  
  data_thresholds <- data_plateaus %>% 
    dplyr::summarise(mode = find_density_modes(!!rlang::sym(value_col), topN = 1), 
                     .by = c(pulse_time_col, group_cols)) 
  
  p <- data_plateaus %>% 
    ggplot(aes(x = !!rlang::sym(value_col))) + 
    geom_density(fill = "grey") + 
    theme_classic() + 
    coord_cartesian(expand = F, xlim = c(0, 1)) + 
    facet_grid(rows = vars(!!rlang::sym(pulse_time_col)), 
               cols = vars(!!rlang::sym(group_cols)))
  
  return(p)
  
}



