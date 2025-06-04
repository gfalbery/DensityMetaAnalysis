
## Test out this for how I'm messing about with residuals

hist(rlnorm(n = 1000, meanlog = 2, sdlog = 0.5))

# Comparison of logged and not logged with homoscedasticity and linear effect ####

slo1 <- numeric()
slo2 <- numeric()

for (reps in 1:1000){

slope <- 2
y1 <- numeric()
y2 <- numeric()

for(i in 1:10){
  if(i == 1){
    y1 <- slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = 0.8))-exp(1))
    y2 <- slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1))
  }
  if(i>1){
    y1 <- c(y1, slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = 0.8))-exp(1)))
    y2 <- c(y2, slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1)))
  }
}

x <- rep(seq(1, 10, 1), each = 10)

m1 <- lm(y1~scale(x))
m2 <- lm(y2~scale(x))

slo1[reps] <- m1$coefficients[2]
slo2[reps] <- m2$coefficients[2]

}

hist(slo1, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2, add = TRUE, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1))

# Comparison of logged and not logged with heteroscedasticity (residuals more skewed with higher x) and linear effect ####

slo1i <- numeric()
slo2i <- numeric()

for (reps in 1:1000){
  
  slope <- 2
  y1 <- numeric()
  y2 <- numeric()
  
  for(i in 1:10){
    if(i == 1){
      y1 <- slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1*i))-exp(1))
      y2 <- slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1))
    }
    if(i>1){
      y1 <- c(y1, slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1*i))-exp(1)))
      y2 <- c(y2, slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1)))
    }
  }
  
  x <- rep(seq(1, 10, 1), each = 10)
  
  m1 <- lm(y1~scale(x))
  m2 <- lm(y2~scale(x))
  
  slo1i[reps] <- m1$coefficients[2]
  slo2i[reps] <- m2$coefficients[2]
  
}

hist(slo1i, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2i, add = TRUE, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1))

# Comparison of logged and not logged with heteroscedasticity (residuals more skewed with lower x) and linear effect ####

slo1d <- numeric()
slo2d <- numeric()

for (reps in 1:1000){
  
  slope <- 2
  y1 <- numeric()
  y2 <- numeric()
  
  for(i in 1:10){
    if(i == 1){
      y1 <- slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = (1.1-(0.1*i)))-exp(1)))
      y2 <- slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1))
    }
    if(i>1){
      y1 <- c(y1, slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = (1.1-(0.1*i)))-exp(1))))
      y2 <- c(y2, slope*i+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1)))
    }
  }
  
  x <- rep(seq(1, 10, 1), each = 10)
  
  m1 <- lm(y1~scale(x))
  m2 <- lm(y2~scale(x))
  
  slo1d[reps] <- m1$coefficients[2]
  slo2d[reps] <- m2$coefficients[2]
  
}

hist(slo1d, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2d, add = TRUE, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1))

# MOVE ONTO SATURATING VERSIONS ####

# Comparison of logged and not logged with homoscedasticity and saturating effect ####

slo1a <- numeric()
slo2a <- numeric()

for (reps in 1:1000){
  
  slope <- 2
  means <- c(2, 6, 9.5, 12.5, 15, 17, 18.5, 19.5, 20, 20)
  y1 <- numeric()
  y2 <- numeric()
  
  for(i in 1:10){
    if(i == 1){
      y1 <- means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = 0.8))-exp(1))
      y2 <- means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1))
    }
    if(i>1){
      y1 <- c(y1, means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = 0.8))-exp(1)))
      y2 <- c(y2, means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1)))
    }
  }
  
  x <- rep(seq(1, 10, 1), each = 10)
  
  m1 <- lm(y1~scale(x))
  m2 <- lm(y2~scale(x))
  
  slo1a[reps] <- m1$coefficients[2]
  slo2a[reps] <- m2$coefficients[2]
  
}

hist(slo1a, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2a, add = TRUE, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1))


hist(slo1a, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo1, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo2a, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

# Comparison of logged and not logged with heteroscedasticity (residuals more skewed with higher x) and saturating effect ####

slo1ia <- numeric()
slo2ia <- numeric()

for (reps in 1:1000){
  
  slope <- 2
  means <- c(2, 6, 9.5, 12.5, 15, 17, 18.5, 19.5, 20, 20)
  y1 <- numeric()
  y2 <- numeric()
  
  for(i in 1:10){
    if(i == 1){
      y1 <- means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1*i))-exp(1))
      y2 <- means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1))
    }
    if(i>1){
      y1 <- c(y1, means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1*i))-exp(1)))
      y2 <- c(y2, means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1)))
    }
  }
  
  x <- rep(seq(1, 10, 1), each = 10)
  
  m1 <- lm(y1~scale(x))
  m2 <- lm(y2~scale(x))
  
  slo1ia[reps] <- m1$coefficients[2]
  slo2ia[reps] <- m2$coefficients[2]
  
}

hist(slo1ia, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2ia, add = TRUE, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1))

hist(slo1ia, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo1i, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo2ia, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2i, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

# Comparison of logged and not logged with heteroscedasticity (residuals more skewed with lower x) and saturating effect ####

slo1da <- numeric()
slo2da <- numeric()

for (reps in 1:1000){
  
  slope <- 2
  y1 <- numeric()
  y2 <- numeric()
  
  for(i in 1:10){
    if(i == 1){
      y1 <- means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = (1.1-(0.1*i)))-exp(1)))
      y2 <- means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1))
    }
    if(i>1){
      y1 <- c(y1, means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = (1.1-(0.1*i)))-exp(1))))
      y2 <- c(y2, means[i]+((rlnorm(n = 10, meanlog = 1, sdlog = 0.1))-exp(1)))
    }
  }
  
  x <- rep(seq(1, 10, 1), each = 10)
  
  m1 <- lm(y1~scale(x))
  m2 <- lm(y2~scale(x))
  
  slo1da[reps] <- m1$coefficients[2]
  slo2da[reps] <- m2$coefficients[2]
  
}

hist(slo1da, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2da, add = TRUE, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1))

hist(slo1da, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo1d, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo2da, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2d, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

####################################################
####################################################

hist(slo2, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2i, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo2, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2d, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo2, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2a, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo2, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2ia, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo2, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo2da, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)


##########

hist(slo1, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo1i, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo1, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo1d, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo1, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo1a, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo1, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo1ia, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)

hist(slo1, col = adjustcolor("firebrick", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400))
hist(slo1da, col = adjustcolor("dodgerblue", 0.3), border = NA, breaks = seq(2, 10, 0.1), ylim = c(0, 400), add = TRUE)



