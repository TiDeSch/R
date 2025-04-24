library(deSolve)
library(ggplot2)
library(dplyr)
library(gridExtra)

N <- 5000
gamma <- 1/7
R_0_values <- c(0.5, 2, 5)
beta_values = R_0_values * gamma

I0 <- 1
R0 <- 0
S0 <- N - I0 - R0

parameters <- c(N = N,
               gamma = gamma)

state <- c(S = S0,
           I = I0,
           R = R0)
                
SIR_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    
    list(c(dS, dI, dR))
  })
}

times <- seq(0, 100, by = 1)

plots <- list()
z <- list()
for (i in 1:length(R_0_values)){
    R_0 <- R_0_values[i]
    beta <- beta_values[i]

    parameters$R_0 <- R_0
    parameters$beta <- beta

    SIR_result <- ode(y = state, times = times, func = SIR_model, parms = parameters)
    SIR_df <- as.data.frame(SIR_result)

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

    final_size <- function(R_inf, S0, R_0) {
        return(1 - R_inf - S0 * exp(-R_0 * R_inf))
    }

    final_size_parameters <- c(S0 = S0, R_0 = R_0)
    f_value <- uniroot(f=final_size, interval=c(0,1),
                       tol=0.0001,
                       S0=S0/N,
                       R_0=R_0)

    p <- p + geom_hline(yintercept = f_value$root, linetype = 'dashed', color = 'darkgreen')
    plots[[i]] <- p
    z[[i]] <- f_value$root
}

grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol = 1)
print(z)

overshoot <- numeric(length(z))
for (i in 1:length(z)){
    OS <- z[[i]] - 1 / R_0_values[i]
    overshoot[[i]] <- OS
}
print(overshoot)