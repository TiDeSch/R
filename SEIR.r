library(deSolve)
library(ggplot2)
library(dplyr)
library(gridExtra)

N <- 5000
gamma <- 0.1
beta <- 0.3
sigma <- 0.2
mu <- 0.01
R_0 <- beta / gamma

I0 <- 1 / N
R0 <- 0 / N
S0 <- (N - I0 - R0) / N
E0 <- 0 / N

parameters <- c(N = N,
                gamma = gamma,
                R_0 = R_0,
                beta = beta,
                sigma = sigma,
                mu = mu)

state <- c(S = S0,
           E = E0,
           I = I0,
           R = R0)

SEIR_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- mu - beta * S * I - mu * S
    dE <- beta * S * I - (mu + sigma) * E
    dI <- sigma * E - (mu + gamma) * I
    dR <- gamma * I - mu * R
    
    list(c(dS, dE, dI, dR))
  })
}

times <- seq(0, 200, by = 2)

SEIR_result <- ode(y = state, times = times, func = SEIR_model, parms = parameters)
SEIR_df <- as.data.frame(SEIR_result)

p1 <- ggplot(SEIR_df, aes(x = time)) +
  geom_line(aes(y=S, colour="Susceptible")) + 
  geom_line(aes(y=E, colour="Exposed")) + 
  geom_line(aes(y=I, colour="Infected")) + 
  geom_line(aes(y=R, colour="Recovered")) + 
  xlab(label="Time (days)") + 
  ylab(label="Proportion of Population") + 
  ggtitle("SIR Model") + 
  scale_colour_manual("Compartments",
  breaks=c("Susceptible", "Exposed", "Infected", "Recovered"),
  values=c("blue", "green", "red", "darkgreen"))

grid.arrange(p1, ncol = 1)