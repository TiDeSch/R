library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

N <- 5000
gamma <- 1/7
R0 <- 2
beta = R0 * gamma

I0 <- 1
S0 <- N - I0

parameters <- c(N = N,
                gamma = gamma,
                R0 = R0,
                beta = beta)

state <- c(S = S0,
           I = I0,
           R = 0)
                

SIR_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    
    list(c(dS, dI, dR))
  })
}

times <- seq(0, 100, by = 1)

SIR_result <- ode(y = state, times = times, func = SIR_model, parms = parameters)
SIR_df <- as.data.frame(SIR_result)

print(
p <- ggplot(SIR_df, aes(x = time)) +
  geom_line(aes(y=S, colour="Susceptible")) + 
  geom_line(aes(y=R, colour="Recovered")) + 
  geom_line(aes(y=I, colour="Infected")) + 
  xlab(label="Time (days)") + 
  ylab(label="Proportion of Population") + 
  ggtitle("SIR Model") + 
  scale_colour_manual("Compartments",
  breaks=c("Susceptible","Infected","Recovered"),
  values=c("blue","red","darkgreen"))
)