find_peak <- function(country1)
{
  data <- covid_data %>% filter(country==country1) %>% 
    select(date, confirmed) %>% 
    na.exclude()
  if(nrow(data)!=0)
  {
    spline <- with(data, smooth.spline(date, confirmed, df=20))  # can change to cross-validation method here
    spline_p <- data %>% mutate(confirmed=fitted(spline))
    peak <- data %>% filter(ggpmisc:::find_peaks(spline_p$confirmed))
    
    result <- peak$date
    result
  }else{as.Date("0000-01-15")}
}

find_valley <- function(country1)
{
  data <- covid_data %>% filter(country==country1) %>% 
    select(date, confirmed) %>% 
    na.exclude()
  if(nrow(data)!=0)
  {
    spline <- with(data, smooth.spline(date, confirmed, df=20))
    spline_p <- data %>% mutate(confirmed=fitted(spline))
    peak <- data %>% filter(ggpmisc:::find_peaks(-spline_p$confirmed))
    
    result <- peak$date
    result
  }else{as.Date("0000-01-15")}
}

first_date <- function(country1)
{
  data <- covid_data %>% filter(country==country1) %>% 
    select(date, confirmed) %>% 
    na.exclude()
  data$date[1]
}

select_candidate <- function(target_country, t_interval, n_peak=2)
{
  population <- readr::read_csv("covid_data/population.csv") %>% 
    select(country, pop) %>% 
    filter(!is.na(pop)) %>% 
    filter(pop>1000000)
  
  # candidate countries
  peaks <- lapply(population$country, find_peak)
  peak_data <- plyr::ldply(peaks,rbind)
  peak_data[peak_data<0] <- NA_integer_
  peak_data <- cbind(population,peak_data)
  # change the name of target peak to more convenient name "peak"
  names(peak_data)[names(peak_data)==as.character(n_peak)] <- "peak"
  
  target_peak <- filter(peak_data, country==target_country)$peak  # 2nd peak of UK

  candidates <- peak_data %>% 
    filter(peak<target_peak) %>% 
    mutate(shift=target_peak-peak) %>% 
    mutate(peak = as.Date.numeric(peak, origin = "1970-01-01"))
  
  firstd <- sapply(candidates$country, first_date) %>% 
    as.Date.numeric(origin = "1970-01-01")
  
  candidates <- candidates %>% 
    mutate(end_date = peak+trunc(t_interval/2),
           start_date = peak-t_interval+trunc(t_interval/2)) %>% 
    filter(firstd<=start_date)  # delete those counties didn't suffice data
}