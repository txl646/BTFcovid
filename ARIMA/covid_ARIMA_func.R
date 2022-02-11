library(forecast)

fit_covid_arima <- function(country1, pop, start_date, end_date="2021-05-04",projection_days=10)
{
  covid_data <- readRDS(".//covid_data//covid_data.RData")
  dat0 <- covid_data %>% 
    filter(country==country1) %>% 
    select(-death,-recovered) %>% 
    na.exclude()
  
  dat <- covid_data %>% 
    filter(country==country1) %>% 
    select(-death,-recovered) %>% 
    slice(which(date == start_date):which(date == end_date)) %>% 
    na.exclude()
  
  # run arima model
  fit <- auto.arima(dat$confirmed, approximation = F)  # lambda = Lambda
  
  # calculate residuals
  # raw - fitted
  I_f <- fitted(fit)
  I_res <- dat$confirmed - I_f
  dat <- mutate(dat,
                I_fitted = I_f,
                I_residual = I_res)
  
  ####### projection data set #######
  projection <- as.data.frame(forecast(fit, h=projection_days))
  projection <- mutate(projection,
                       date=end_date+(1:nrow(projection)),
                       `Lo 95`=if_else(`Lo 95`>0,`Lo 95`,0))  # delete the negative error bar
  colnames(projection)[1] <- "confirmed"
  
  ####### draw plot #######
  p <- dat %>% 
    ggplot() +
    geom_line(data = dat0,
              color="black", inherit.aes = FALSE,
              aes(x=date, y=confirmed*pop),
              lwd=0.35, alpha=0.9) +
    geom_point(data = dat0,
               color = "grey30", inherit.aes = FALSE,
               aes(x=date, y=confirmed*pop),
               pch = 21, fill = "grey95", size = 1.25) +
    geom_line(aes(x=date, y=I_fitted*pop),
              color="red", alpha=0.5, size=1) +
    labs(title = paste0("ARIMA Model for ", country1)) +
    xlab("time") +
    ylab("cases") +
    theme_bw()
  
  result <- list(parameter = fit,
                 raw_data = dat, # within selected period
                 projection = projection,
                 plot = p)
}