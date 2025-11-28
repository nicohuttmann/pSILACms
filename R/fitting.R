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
extract_fit_parameters_sety <- function(data, 
                                        x_col = "Fraction_Old", 
                                        y_col = "Pulse_time", 
                                        id_cols = "variables", 
                                        remove_cols = c("data", "model", "glance", "tidy")) {
  
  data_parameters <- data %>% 
    dplyr::filter(purrr::map_chr(model, class) == "lm") %>% 
    dplyr::mutate(glance = purrr::map(model, broom::glance), 
                  tidy = purrr::map(model, broom::tidy)) %>% 
    dplyr::mutate(rsq = glance %>% purrr::map_dbl("r.squared"), 
                  n = purrr::map_dbl(data, \(x) sum(!is.na(x[[x_col]]))), 
                  slope = purrr::map_dbl(tidy, \(x) x %>% 
                                           dplyr::filter(term == y_col) %>% 
                                           dplyr::pull("estimate")), 
                  t_half = -log(2) / slope, 
                  t_mean = -1 / slope)
  
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
