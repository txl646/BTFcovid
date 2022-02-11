raw_data_by_country <- function(country)
{
  country1 <- country
  covid_data <- readRDS(".//covid_data//covid_data.RData")
  data <- covid_data %>% 
    filter(country==country1) %>% 
    mutate(death = if_else(death<0, as.integer(0), death)) %>% 
    na.exclude()
  
  result <- data
}

fit_covid_SIR <- function(country, pop, start_date, end_date="2021-05-04", 
                          prior_beta=0.5, prior_gamma=0.5)
{
  ####### prepare data #######
  country1 <- country
  covid_data <- readRDS(".//covid_data//covid_data.RData")
  dat0 <- covid_data %>% 
    filter(country==country1) %>% 
    mutate(death = if_else(death<0, as.integer(0), death)) %>% 
    na.exclude()
  
  dat <- covid_data %>% 
    filter(country==country1) %>% 
    mutate(death = if_else(death<0, as.integer(0), death)) %>% 
    slice(which(date == start_date):which(date == end_date)) %>% 
    na.exclude()
  
  N <- pop
  
  ####### find parameters estimation #######
  source(".//SIR//SIR_func.R")
  
  t <- 1:length(dat$confirmed)-1
  sse <- function(par)
  {
    out <- solve_sir(days = length(dat$confirmed),
                     beta = par[1],
                     gamma = par[2],
                     N = N,
                     S0 = 0.99*N-dat$confirmed[1]-
                       sum(dat0$death[1:which(dat0$date == start_date)])-
                       sum(dat0$recovered[1:which(dat0$date == start_date)]),
                     I0 = dat$confirmed[1],
                     R0 = dat$recovered[1]+dat$death[1],
                     time_step=0.1)
    
    I_esti <- out[which(out$time %in% t),]$I
    I_real <- dat$confirmed
    I_SSE <- sum((I_real-I_esti)^2)
    I_SSE
  }
  
  fit <- optim(fn=sse, par=c(prior_beta, prior_gamma))
  
  ####### calculate final result #######
  out <- solve_sir(days = length(dat$confirmed),
                   beta = fit$par[1],
                   gamma = fit$par[2],
                   N = N,
                   S0 = 0.99*N-dat$confirmed[1]-
                     sum(dat0$death[1:which(dat0$date == start_date)])-
                     sum(dat0$recovered[1:which(dat0$date == start_date)]),
                   I0 = dat$confirmed[1],
                   R0 = dat$recovered[1]+dat$death[1])
  
  # calculate residuals
  I_f <- out[which(out$time %in% t),]$I
  # raw - fitted
  I_res <- dat$confirmed-I_f
  dat <- mutate(dat,
                I_fitted = I_f,
                I_residual = I_res)
  
  ####### draw plot #######
  p <- out %>% 
    mutate(time = time + which(dat0$date == start_date)-1) %>% 
    mutate(time = dat0$date[1]+time) %>% 
    ggplot() +
    geom_line(data = dat0,
              color="black", inherit.aes = FALSE,
              aes(x=date, y=confirmed),
              lwd=0.35, alpha=0.9) +
    geom_point(data = dat0,
               color = "grey30", inherit.aes = FALSE,
               aes(x=date, y=confirmed),
               pch = 21, fill = "grey95", size = 1.25) +
    geom_line(aes(x=time, y=I),
              color="#0091F7", alpha=0.8, size=1.5) +
    labs(title = paste0("SIR Model for ", country1)) +
    xlab("Date") +
    ylab("Daily Cases") +
    theme_bw()
  
  result <- list(parameter = list(beta = fit$par[1],
                                  gamma = fit$par[2]),
                 raw_data = dat, # within selected period
                 output_data = out,
                 plot = p)
}

plot_covid_SIR_predict <- function(beta, gamma, 
                                     model_country, predict_country, 
                                     start_date, time_gap=0,
                                     pop, predict_days=100)
{
  start_date <- as.Date(start_date)+time_gap
  end_date <- "2021-05-04"
  ####### predict using other country's parameters #######
  covid_data <- readRDS(".//covid_data//covid_data.RData")
  dat0_ind <- covid_data %>% 
    filter(country==predict_country) %>% 
    mutate(death = if_else(death<0, as.integer(0), death)) %>% 
    na.exclude()
  
  dat_ind <- covid_data %>% 
    filter(country==predict_country) %>% 
    mutate(death = if_else(death<0, as.integer(0), death)) %>% 
    slice(which(date == start_date):which(date == end_date)) %>% 
    na.exclude()
  
  
  N_ind <- pop
  source(".//SIR//SIR_func.R")
  par <- c(beta, gamma, N_ind)
  out_ind <- solve_sir(days = length(dat_ind$confirmed)+predict_days,
                       beta = par[1],
                       gamma = par[2],
                       N = par[3],
                       S0 = 0.99*N_ind - dat_ind$confirmed[1]-
                         sum(dat0_ind$death[1:which(dat0_ind$date==start_date)])-
                         sum(dat0_ind$recovered[1:which(dat0_ind$date==start_date)]),
                       I0 = dat_ind$confirmed[1],
                       R0 = dat_ind$recovered[1]+dat_ind$death[1])
  
  out_ind <- out_ind %>%
    mutate(time=time+which(dat0_ind$date == start_date)-1) %>% 
    mutate(time = dat0_ind$date[1]+time)
  
  ####### specificity #######
  # jitter the parameters
  par <- cbind(beta+runif(10,-0.01,0.01), 
               gamma+runif(10,-0.01,0.01),
               N_ind+sample(1:100,10))
  
  SIR_err <- function(par)
  {
    out_err <- solve_sir(days = length(dat_ind$confirmed)+predict_days,
                         beta = par[1],
                         gamma = par[2],
                         N = par[3],
                         S0 = 0.99*N_ind - dat_ind$confirmed[1]-
                           sum(dat0_ind$death[1:which(dat0_ind$date==start_date)])-
                           sum(dat0_ind$recovered[1:which(dat0_ind$date==start_date)]),
                         I0 = dat_ind$confirmed[1],
                         R0 = dat_ind$recovered[1]+dat_ind$death[1]) %>%
      mutate(time=time+which(dat0_ind$date == start_date)-1) %>% 
      mutate(time = dat0_ind$date[1]+time) %>% 
      select(time, I)
  }
  
  out_er <- apply(par,1,SIR_err) %>% 
    reduce(full_join, by="time")
  out_err <- out_er %>% 
    mutate(I_max = apply(out_er[,-1],1,max),
           I_min = apply(out_er[,-1],1,min))
  
  ####### draw plot #######
  p <- out_ind %>%
    ggplot() +
    geom_line(data = dat0_ind,
              color="black", inherit.aes = FALSE,
              aes(x=date, y=confirmed),
              lwd=0.35, alpha=0.9) +
    geom_point(data = dat0_ind,
               color = "grey30", inherit.aes = FALSE,
               aes(x=date, y=confirmed),
               pch = 21, fill = "grey95", size = 1.25) +
    annotate("rect", xmin=start_date-time_gap, xmax=as.Date(Inf),
             ymin=0, ymax=Inf, alpha=0.1, fill="#0091F7") +
    # specificity shadow
    geom_ribbon(data = out_err, inherit.aes = FALSE,
                aes(ymin=I_min, ymax=I_max, x=time),
                fill="#0091F7", alpha=0.3) +
    # prediction line
    geom_line(aes(x=time, y=I),
              color="#0091F7", alpha=0.8, size=1.5) +
    
    labs(title = paste0("Projection of ", predict_country, " using ", model_country, " Model"),
         x = "Date",
         y = "Daily Cases") +
    theme_bw()
  
  result <- list(plot = p)
}
