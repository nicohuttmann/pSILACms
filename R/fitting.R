#' Title
#'
#' @param data tibble containing fitted model 
#' @param x_col column containing fitted values 
#' @param y_col 
#' @param id_cols 
#' @param remove_cols columns to remove after mutation 
#' @returns
#' @export
#'
#' @examples
extract_fit_parameters_lm <- function(data, 
                                      x_col = "Fraction_Old", 
                                      y_col = "Pulse_time", 
                                      id_cols = "variables", 
                                      remove_cols = c("data", "model", "glance", "tidy")) {
  
  # Extract model 
  data_parameters <- data %>% 
    dplyr::filter(purrr::map_chr(model, class) == "lm") %>% 
    dplyr::mutate(glance = purrr::map(model, broom::glance), 
                  tidy = purrr::map(model, broom::tidy)) 
  
  
  # Calculate and extract parameters 
  if (!"(Intercept)" %in% data_parameters$tidy[[1]]$term) {
    
    # Fixed y-intercept 
    data_parameters <- data_parameters %>% 
      dplyr::mutate(rsq = glance %>% purrr::map_dbl("r.squared"), 
                    n = purrr::map_dbl(data, \(x) sum(!is.na(x[[x_col]]))), 
                    slope = purrr::map_dbl(tidy, \(x) x %>% 
                                             dplyr::filter(term == y_col) %>% 
                                             dplyr::pull("estimate")), 
                    t_half = log(0.5) / slope, 
                    t_mean = log( 1/exp(1) ) / slope)
    
  } else {
    
    # Fitted y-intercept 
    data_parameters <- data_parameters %>% 
      dplyr::mutate(rsq = glance %>% purrr::map_dbl("r.squared"), 
                    n = purrr::map_dbl(data, \(x) sum(!is.na(x[[x_col]]))), 
                    slope = purrr::map_dbl(tidy, \(x) x %>% 
                                             dplyr::filter(term == y_col) %>% 
                                             dplyr::pull("estimate")), 
                    y_intercept = purrr::map_dbl(tidy, \(x) x %>% 
                                                   dplyr::filter(term == "(Intercept)") %>% 
                                                   dplyr::pull("estimate")), 
                    P_0 = exp(y_intercept), 
                    t_half = ( log(0.5) - y_intercept ) / slope, 
                    t_mean = ( log( 1/exp(1) ) - y_intercept ) / slope)
    
    
  }
  
  # Rejoin data and remove columns 
  data_return <- data %>% 
    dplyr::select(-c(data, model)) %>% 
    dplyr::full_join(data_parameters, 
                     by = id_cols) %>% 
    dplyr::select(-dplyr::all_of(remove_cols))
  
  
  return(data_return)
  
}


#' Title
#'
#' @param data tibble of fitted models 
#' @param model_col column name containing fitted models 
#'
#' @returns
#' @export
#'
#' @examples
extract_fit_residuals <- function(data, model_col = "model") {
  
  data %>% 
    dplyr::mutate(augment = purrr::map(!!rlang::sym(model_col), broom::augment))
  
}

