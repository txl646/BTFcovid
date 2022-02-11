sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)),
       {
         
         dS <- -beta * S * I / N
         dI <-  beta * S * I / N - gamma * I
         dR <-                     gamma * I
         
         return(list(c(dS, dI, dR)))
       }
       )
}

######## solve SIR by general ode function ########
solve_sir <- function(days, beta, gamma, N, S0, I0, R0=0, time_step=0.01)
{
  ## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
  init <- c(S = S0, 
            I = I0, 
            R = R0)
  
  ## beta: infection parameter; gamma: recovery parameter
  parameters <- c(beta = beta,
                  gamma = gamma,
                  N = N)
  
  ## Time frame
  times <- seq(0, days-1, by = time_step)
  
  ## default method : lsoda
  out <- deSolve::ode(y = init,
                      times = times, 
                      func = sir,
                      parms = parameters) %>% 
    as.data.frame()
}