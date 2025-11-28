

library(tidyverse)







# ---- Define analysis list ---- 
list_SILAC <- list(data_raw = tibble(), 
                   data_filtered = tibble(), 
                   data_filtered2 = tibble(), 
                   fit = list())


# ---- Add raw data ---- 

# MS1 
list_SILAC$data_raw_ms1 <- get_data_frame("Ms1.Area", 
                                          variables = !contam & SILAC_label != "R", 
                                          observations = SPP_fraction == "Total", 
                                          dataset = "SPP_Precursors") %>% 
  pivot_longer(-1, 
               names_to = "variables", 
               values_to = "Ms1.Area") %>% 
  filter(.is_ratio(Ms1.Area)) %>% 
  # Add precursor annotations 
  add_variables_data(c("SILAC.Precursor.Group", 
                       "SILAC_label", 
                       "Modified.Sequence", 
                       "PTM.Modified.Sequence"), 
                     dataset = "SPP_Precursors") %>% 
  # Calculate SILAC ratios 
  pivot_wider(id_cols = c("SILAC.Precursor.Group", 
                          "observations"), 
              names_from = "SILAC_label", 
              values_from = "Ms1.Area") %>% 
  mutate(L = replace_na(L, 0), 
         H = replace_na(H, 0)) %>% 
  mutate(Total = L + H, 
         Fraction_L = L / Total, 
         Fraction_H = H / Total) %>% 
  # Correct H fraction by estimate lysine pool H%
  add_observations_data("mode_Hp_ms1", dataset = "SPP_Precursors") %>% 
  mutate(Fraction_New = Fraction_H / mode_Hp_ms1, 
         Fraction_Old = 1 - Fraction_New) %>% 
  # Set values below 0 and above 1 to range 
  mutate(Fraction_Old = case_when(Fraction_Old < 0 ~ 0, 
                                  Fraction_Old > 1 ~ 1, 
                                  .default = Fraction_Old)) %>% 
  # Annotate samples 
  add_observations_data(c("pulse_time", 
                          "sex"), 
                        dataset = "SPP_Precursors") %>% 
  rename(variables = SILAC.Precursor.Group) %>% 
  select(observations, variables, pulse_time, sex, Fraction_L, Fraction_Old)

# MS2 
list_SILAC$data_raw_ms2 <- get_data_frame("Precursor.Quantity", 
                                          variables = !contam & SILAC_label != "R", 
                                          observations = SPP_fraction == "Total", 
                                          dataset = "SPP_Precursors") %>% 
  pivot_longer(-1, 
               names_to = "variables", 
               values_to = "Precursor.Quantity") %>% 
  filter(.is_ratio(Precursor.Quantity)) %>% 
  # Add precursor annotations 
  add_variables_data(c("SILAC.Precursor.Group", 
                       "SILAC_label", 
                       "Modified.Sequence", 
                       "PTM.Modified.Sequence"), 
                     dataset = "SPP_Precursors") %>% 
  # Calculate SILAC ratios 
  pivot_wider(id_cols = c("SILAC.Precursor.Group", 
                          "observations"), 
              names_from = "SILAC_label", 
              values_from = "Precursor.Quantity") %>% 
  mutate(L = replace_na(L, 0), 
         H = replace_na(H, 0)) %>% 
  mutate(Total = L + H, 
         Fraction_L = L / Total, 
         Fraction_H = H / Total) %>% 
  # Correct H fraction by estimate lysine pool H%
  add_observations_data("mode_Hp_ms1", dataset = "SPP_Precursors") %>% 
  mutate(Fraction_New = Fraction_H / mode_Hp_ms1, 
         Fraction_Old = 1 - Fraction_New) %>% 
  # Set values below 0 and above 1 to range 
  mutate(Fraction_Old = case_when(Fraction_Old < 0 ~ 0, 
                                  Fraction_Old > 1 ~ 1, 
                                  .default = Fraction_Old)) %>% 
  # Annotate samples 
  add_observations_data(c("pulse_time", 
                          "sex"), 
                        dataset = "SPP_Precursors") %>% 
  rename(variables = SILAC.Precursor.Group) %>% 
  select(observations, variables, pulse_time, sex, Fraction_L, Fraction_Old)






# ---- Filter data ---- 

## By consecutive data points ---- 
list_SILAC$data_filtered_ms1 <- list_SILAC$data_raw_ms1 %>% 
  filter_by_consecutive_datapoints()


## Remove datapoints that reach full turnover ---- 
list_SILAC$data_filtered2_ms1 <- list_SILAC$data_filtered_ms1 %>% 
  filter_by_SILAC_plateaus()



## By consecutive data points ---- 
list_SILAC$data_filtered_ms2 <- list_SILAC$data_raw_ms2 %>% 
  filter_by_consecutive_datapoints()


## Remove datapoints that reach full turnover ---- 
list_SILAC$data_filtered2_ms2 <- list_SILAC$data_filtered_ms2 %>% 
  filter_by_SILAC_plateaus()




# Visualize SILAC ratios for plateau identification ---- 
{
  list_SILAC$data_filtered %>% 
    ggplot(aes(x = Fraction_L)) + 
    geom_density(fill = "grey") + 
    theme_classic() + 
    coord_cartesian(expand = F, xlim = c(0, 1)) + 
    facet_grid(rows = vars(pulse_time), 
               cols = vars(sex))
  
  list_SILAC$data_filtered %>% 
    filter(is_SILAC_ratio(Fraction_L)) %>% 
    ggplot(aes(x = Fraction_L)) + 
    geom_density(fill = "grey") + 
    theme_classic() + 
    coord_cartesian(expand = F, xlim = c(0, 1)) + 
    facet_grid(rows = vars(pulse_time), 
               cols = vars(sex))
  
  
  list_SILAC$data_filtered %>% 
    ggplot(aes(x = Fraction_Old)) + 
    geom_density(fill = "grey") + 
    theme_classic() + 
    coord_cartesian(expand = F, xlim = c(0, 1)) + 
    facet_grid(rows = vars(pulse_time), 
               cols = vars(sex))
  
  list_SILAC$data_filtered %>% 
    filter(is_SILAC_ratio(Fraction_Old)) %>% 
    ggplot(aes(x = Fraction_Old)) + 
    geom_density(fill = "grey") + 
    theme_classic() + 
    coord_cartesian(expand = F, xlim = c(0, 1)) + 
    facet_grid(rows = vars(pulse_time), 
               cols = vars(sex))
  
  
  list_SILAC$data_filtered %>% 
    mutate(n_0 = sum(Fraction_Old == 0), 
           .by = c("variables", "pulse_time", "sex")) %>% 
    filter(is_SILAC_ratio(Fraction_Old)) %>% 
    ggplot(aes(x = Fraction_Old, 
               color = factor(n_0), 
               fill = factor(n_0))) + 
    geom_density(alpha = 0.3, linewidth = 1.2) + 
    theme_classic() + 
    scale_color_embl() + 
    scale_fill_embl() + 
    coord_cartesian(expand = F, xlim = c(0, 1)) + 
    facet_grid(rows = vars(pulse_time), 
               cols = vars(sex))
}





## Fit exponential model ---- 
list_SILAC$fit_ms1 <- list_SILAC$data_filtered2_ms1 %>% 
  filter(is_SILAC_ratio(Fraction_Old)) %>% 
  nest(.by = c("variables")) %>% 
  mutate(model = map(data, \(x) tryCatch({lm(log(Fraction_Old) ~ pulse_time + 0, x)}, 
                                         error = function(cond) NA, 
                                         warning = function(cond) NA)))

list_SILAC$parameters_ms1 <- list_SILAC$fit_ms1 %>% 
  extract_fit_parameters()


list_SILAC$fit_ms2 <- list_SILAC$data_filtered2_ms2 %>% 
  filter(is_SILAC_ratio(Fraction_Old)) %>% 
  nest(.by = c("variables")) %>% 
  mutate(model = map(data, \(x) tryCatch({lm(log(Fraction_Old) ~ pulse_time + 0, x)}, 
                                         error = function(cond) NA, 
                                         warning = function(cond) NA)))

list_SILAC$parameters_ms2 <- list_SILAC$fit_ms2 %>% 
  extract_fit_parameters()


# unfiltered plateaus 
list_SILAC$fit_np_ms1 <- list_SILAC$data_filtered_ms1 %>% 
  filter(is_SILAC_ratio(Fraction_Old)) %>% 
  nest(.by = c("variables")) %>% 
  mutate(model = map(data, \(x) tryCatch({lm(log(Fraction_Old) ~ pulse_time + 0, x)}, 
                                         error = function(cond) NA, 
                                         warning = function(cond) NA)))

list_SILAC$parameters_np_ms1 <- list_SILAC$fit_np_ms1 %>% 
  extract_fit_parameters()


list_SILAC$fit_np_ms2 <- list_SILAC$data_filtered_ms2 %>% 
  filter(is_SILAC_ratio(Fraction_Old)) %>% 
  nest(.by = c("variables")) %>% 
  mutate(model = map(data, \(x) tryCatch({lm(log(Fraction_Old) ~ pulse_time + 0, x)}, 
                                         error = function(cond) NA, 
                                         warning = function(cond) NA)))

list_SILAC$parameters_np_ms2 <- list_SILAC$fit_np_ms2 %>% 
  extract_fit_parameters()



# unfiltered 
list_SILAC$fit_n_ms1 <- list_SILAC$data_raw_ms1 %>% 
  filter(is_SILAC_ratio(Fraction_Old)) %>% 
  nest(.by = c("variables")) %>% 
  mutate(model = map(data, \(x) tryCatch({lm(log(Fraction_Old) ~ pulse_time + 0, x)}, 
                                         error = function(cond) NA, 
                                         warning = function(cond) NA)))

list_SILAC$parameters_n_ms1 <- list_SILAC$fit_n_ms1 %>% 
  extract_fit_parameters()


list_SILAC$fit_n_ms2 <- list_SILAC$data_raw_ms2 %>% 
  filter(is_SILAC_ratio(Fraction_Old)) %>% 
  nest(.by = c("variables")) %>% 
  mutate(model = map(data, \(x) tryCatch({lm(log(Fraction_Old) ~ pulse_time + 0, x)}, 
                                         error = function(cond) NA, 
                                         warning = function(cond) NA)))

list_SILAC$parameters_n_ms2 <- list_SILAC$fit_n_ms2 %>% 
  extract_fit_parameters()





list_SILAC$parameters_np_ms1 %>% 
  #filter(.in_lower_quantile(t_half)) %>% 
  ggplot(aes(x = t_half)) + 
  geom_density() + 
  scale_x_continuous(limits = c(0, 30))

patchwork::wrap_plots(
  imap(list_SILAC %>% 
        subset(., str_detect(names(.), "parameters")), 
      \(x, i) x %>% 
        filter(str_count(variables, "K") == 2) %>% 
        ggplot(aes(x = t_half)) + 
        geom_density() + 
        scale_x_continuous(limits = c(0, 30)) + 
        coord_cartesian(ylim = c(0, 0.15)) + 
        labs(title = i)), byrow = F
)


