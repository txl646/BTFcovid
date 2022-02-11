library(tidyverse)

covid_data <- readRDS(".//covid_data//covid_data.RData")

BTF <- function(target_country, path, 
                peakmatch=2, t_interval=75, t_gap=10,
                future_days=100)
  # peakmatch -  the order of the peak to match
  # t_interval - the length of matching time window
  # t_gap - washout interval
{
  ####### find which country match target country best (before surge) #######
  #### prepare candidate countries used by ARIMA ####
  # population information of different countries
  source(".//ARIMA//select_candidate.R")
  countries <- select_candidate(target_country, t_interval, peakmatch)
  countries_test <- countries$country
  start_date <- countries$start_date
  end_date <- countries$end_date
  N <- countries$pop
  
  #### fit ARIMA model for assistant countries ####
  source(".//ARIMA//covid_ARIMA_func.R")
  
  assist_ARIMA_fit <- lapply(1:length(countries_test), 
                             function(i)
                             {
                               fit_covid_arima(countries_test[i], N[i],
                                               start_date[i], end_date[i])
                             })
  
  #### fit ARIMA model for target country ####
  N_pre <- filter(readr::read_csv("covid_data/population.csv"),
                  country == target_country)$pop
  p_country_end <- find_valley(target_country)[peakmatch+1+1]
  projection_end <- find_valley(target_country)[peakmatch+2+1]  # end date of projection
  # p_country_end <- as.Date(end)
  target_fit <- fit_covid_arima(target_country, N_pre,
                                p_country_end-t_interval, p_country_end)
  
  
  ####  compare similarity between target and assistant countries (by MSE) ####
  MD <- rep(NA,length(countries_test))
  for(i in 1:length(countries_test))
  {
    # standardize by maximum
    target_standard <- target_fit$raw_data$I_fitted/max(target_fit$raw_data$I_fitted)
    test_standard <- assist_ARIMA_fit[[i]]$raw_data$I_fitted/max(assist_ARIMA_fit[[i]]$raw_data$I_fitted)
    
    MD[i] <- median((target_standard - test_standard)^2)
  }
  
  # rank MSE
  MD_rank <- cbind(countries_test, MD)[order(MD),]
  MD_rank <- data.frame(MD_rank) %>% 
    mutate(MD=signif(as.numeric(MD),4))  # 4 significant digits
  
  # select the country has minimum MSE (omit the country with all 0 data)
  index_sim <- which(MD==min(na.omit(MD[1:length(countries_test)])))
  
  # save rank table
  dir.create(paste0(path, target_country), showWarnings = F)
  path <- paste0(path, target_country, "//")
  library(gt)
  gt(MD_rank[1:10,]) %>% 
    gtsave(paste0(path, "rank.png"))

  ####### predict target country guided by other countries #######
  source(".//SIR//covid_SIR_func.R")
  assist_country_info <- countries[index_sim,]
  
  # future fit of SIR model for guiding countries
  f_assist_SIR_fit <- fit_covid_SIR(assist_country_info$country, N[index_sim], 
                                    as.Date(assist_country_info$end_date)+t_gap,
                                    as.Date(assist_country_info$end_date)+t_gap+future_days)
  
  # predict target by other countries
  predict_target <- 
    plot_covid_SIR_predict(beta = f_assist_SIR_fit$parameter$beta,
                           gamma = f_assist_SIR_fit$parameter$gamma,
                           model_country = assist_country_info$country,
                           predict_country = target_country,
                           pop = N_pre,
                           start_date = p_country_end,
                           time_gap = t_gap)
  
  ggsave(filename = paste0(path, "predict_", target_country,
                           "_by_", assist_country_info$country, ".png"),
         plot = predict_target$plot,
         width = 0.7*1200/100, height = 0.7*675/100, dpi = 100)
}