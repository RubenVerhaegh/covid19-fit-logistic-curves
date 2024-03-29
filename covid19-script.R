library(MASS)
library(viridis)
library(EstimationTools)
# Reading the input data
##########################################################
data <- read.csv("data/time_series_covid19_confirmed_global.csv", header=TRUE, stringsAsFactors = FALSE, sep=",")
attach(data)
data[nrow(data)+1,] <- NA
data[nrow(data),1] <- ""
data[nrow(data),2] <- "Global"
data[nrow(data),3] <- 0
data[nrow(data),4] <- 0
for (i in 5:ncol(data)) {
  data[nrow(data),i] <- 0
  data[nrow(data),i] <- sum(data[i])
}


# Returns the row index in the "data" dataframe
# of a country given by a string
country_to_index <- function(country) {
  country_seen <- FALSE
  for (i in 1:length(data$Lat)) {
    if (data$Country.Region[i] == country) {
      country_seen <- TRUE
      if (data$Province.State[i] == "") {
        return(i)
      }
    } 
  }
  
  if(country_seen) {
    stop("[!] Ambiguous location: ``", country, "''")
  } else {
    stop("[!] Location ``", country, "'' does not exist")
  }
}

# Returns a dataframe with one row containing the daily
# cumulative confirmed cases for a country as given by a string
daily_cases_for_country <- function(country, min_cases = 100) {
  index <- country_to_index(country)
  untruncated_row <- data[index, 5:ncol(data)]
  first_positive <- -1
  for (i in 1:ncol(untruncated_row)) {
    if (untruncated_row[1,i] > min_cases) {
      first_positive <- i
      break
    }
  }
  if (first_positive < 0) {
    stop("[!] Location ``", country, "'' exists but has 0 confirmed cases")
  }
  
  return(untruncated_row[, first_positive:ncol(untruncated_row)])
}

cumulative_to_increase <- function(cumulative) {
  cumulative_vector <- as.numeric(as.vector(cumulative))
  result <- numeric(length(cumulative_vector) - 1)
  for(i in 1 : length(result)) {
    result[i] <- cumulative_vector[i+1] - cumulative_vector[i]
  }
  return(result)
}

daily_increase_for_country <- function(country) {
  df <- daily_cases_for_country(country)
  print(ncol(df))
  for(i in 1 : (ncol(df)-1)) {
    print(i)
    df[1,i] <- df[1,i+1] - df[1,i]
  }
  return(df[,1:ncol(df)-1])
}

counts_to_observations <- function(counts) {
  observations <- numeric()
  for (i in 1:length(counts)) {
    if (counts[i] > 0) {
      observations <- c(observations, rep(i, counts[i]))
    }
  }
  return(observations)
}



statistical_distance <- function(df, vec) {
  df <- as.numeric(df)
  len <- length(vec)
  if (len != length(df)) {
    stop("Tried to compare to lists of different lengthts")
  }
  
  dist <- 0
  for (i in 1:length(vec)) {
    observed <- df[i]
    expected <- vec[i]
    dist <- dist + abs(observed - expected)
  }
  return(dist)
}

chisqr <- function(observed, estimated) {
  observed <- as.numeric(observed)
  observed_inc <- cumulative_to_increase(observed)
  estimated_inc <- cumulative_to_increase(estimated)
  
  dist <- 0
  for (i in 1:length(estimated_inc)) {
    if (is.na(estimated_inc[i])) {
      return(.Machine$integer.max)
    }
    if (observed_inc[i] != 0) {
      dist <- dist + observed_inc[i] * log(observed_inc[i] / estimated_inc[i])
    }
  }
  
  return(2 * dist)
}

power <- function(base, exponent) {
  (abs(base) ^ exponent) * sign(base)
}


# LOGISTIC model with 1 parameter:
#   r
logistic <- function(N0, parameters, days) {
  if (length(parameters) != 1) {
    stop(paste("Logistic model was called with", length(parameters), "parameters. (Must be 1)"))
  }
  r <- parameters[1]
  
  logistic <- integer(days)
  logistic[1] <- N0
  for (i in 2:days) {
    N  <- logistic[i-1]
    dN <- r * N
    logistic[i] <- N + dN
  }
  
  return(logistic)
}

# VERHULST model with 2 parameters:
#   r, K
verhulst <- function(N0, parameters, days) {
  if (length(parameters) != 2) {
    stop(paste("Logistic model was called with", length(parameters), "parameters. (Must be 2)"))
  }
  r <- parameters[1]
  K <- parameters[2]
  
  if(K == 0) {
    stop("[!] Verhulst can not be computed with K=0")
  }
  verhulst <- integer(days)
  verhulst[1] <- N0
  for (i in 2:days) {
    N  <- verhulst[i-1]
    dN <- r * N * (1 - N/K)
    verhulst[i] <- N + dN
  }
  return(verhulst)
}

# GENERIC model with 5 parameters:
#   r, K, alpha, beta, gamma
generic <- function(N0, parameters, days) {
  if (length(parameters) != 5) {
    stop(paste("Logistic model was called with", length(parameters), "parameters. (Must be 5)"))
  }
  r <- parameters[1]
  K <- parameters[2]
  alpha <- parameters[3]
  beta <- parameters[4]
  gamma <- parameters[5]
  
  if(K == 0) {
    stop("[!] Generic model can not be computed with K=0")
  }
  generic <- integer(days)
  generic[1] <- N0
  for (i in 2:days) {
    N  <- generic[i-1]
    dN <- r * power(N, alpha) * power((1 - power(N/K, beta)), gamma)
    # Line below is to fix discretisation errors
    dN <- max(dN, 0)
    generic[i] <- N + dN
  }
  return(generic)
}


# Draws a plot for a logistic model with different sets of input parameters
parameter_plots <- function(func, N0, parameters, days, var="", var_index=NA) {
  day <- c(1:days)
  cases <- data.frame()
  
  nr_inputs <- length(parameters[[1]])
  nr_parameters <- length(parameters)
  colors <- rainbow(nr_inputs)
  
  inputs <- data.frame()
  for (i in 1:nr_inputs) {
    input <- nr_parameters
    for (p in 1:nr_parameters) {
      input[p] <- parameters[[p]][i]
    }
    inputs <- rbind(inputs, input)
  }
  
  # Plot t against N
  xrange <- range(day)
  yrange <- range(func(N0, as.numeric(inputs[1,]), days))
  plot(xrange, yrange, type="n", xlab="t", ylab="N", cex.lab=1.1)
  for (i in 1:nr_inputs) {
    lines(day, func(N0, as.numeric(inputs[i,]), days), col = colors[i], lwd = 3)
  }
  if (!is.na(var_index)) {
    legend(0, yrange[2], legend=paste(var, "=", inputs[[var_index]]), lty = 1, lwd = 3, col = colors)
  }
  
  # Plot t against dN/dt
  day <- day[1:length(day) - 1]
  xrange <- range(day)
  yrange <- range(cumulative_to_increase(func(N0, as.numeric(inputs[1,]), days)))
  plot(xrange, yrange, type="n", xlab="t", ylab="dN/dt", cex.lab=1.1)
  for (i in 1:nr_inputs) {
    lines(day, cumulative_to_increase(func(N0, as.numeric(inputs[i,]), days)), col = colors[i], lwd = 3)
  }
  if (!is.na(var_index)) {
    legend(0, yrange[2], legend=paste(var, "=", inputs[[var_index]]), lty = 1, lwd = 3, col = colors)
  }
}



update_counters <- function(counter, goal) {
  nr_parameters <- length(counter)
  if (length(goal) != nr_parameters) {
    stop(paste("[!] `counter' has length", length(counter), "and `goal' has length", length(goal)))
  }
  
  for(i in 1:nr_parameters) {
    if (counter[i] < goal[i]) {
      counter[i] <- counter[i] + 1
      break
    } else {
      counter[i] <- 0
    }
  }
  
  return(counter)
}

update_search_bounds <- function(prev_min, prev_max, prev_steps, best_values, outer_min = prev_min, outer_max = prev_max){
  nr_parameters <- length(prev_min)
  if(length(prev_max) != nr_parameters || length(prev_steps) != nr_parameters || length(best_values) != nr_parameters) {
    stop(paste("[!] `min' has length", length(prev_min), 
               "`max' has length", length(prev_max), 
               "`steps' has length", length(prev_steps), 
               "`best' has length", length(best_values)))
  }
  
  new_min <- integer(nr_parameters)
  new_max <- integer(nr_parameters)
  
  for(i in 1:nr_parameters) {
    min <- prev_min[i]
    max <- prev_max[i]
    best <- best_values[i]
    steps <- prev_steps[i]
    step_size <- (max - min) / steps
    
    if (best == outer_min[i]) {
      new_min[i] <- best
      new_max[i] <- best + step_size
    } else if (best == outer_max[i]) {
      new_min[i] <- best - step_size
      new_max[i] <- best
    } else {
      new_min[i] <- best - step_size
      new_max[i] <- best + step_size
    }
  }
  
  output <- list()
  output[[1]] <- new_min
  output[[2]] <- new_max
  return(output)
}



test_fit_cumulative <- function(df, fit) {
  increase <- cumulative_to_increase(df[1,])
  increase_fit <- cumulative_to_increase(fit)
  probabilities <- increase_fit / sum(increase_fit)
  
  print(chisq.test(increase, p = probabilities))
}

# Plots fit (vector) to actual data (data frame)
plot_fit <- function(actual_data, fit, country = "??", func = NA, parameters = c("??")) {
  title_start <- "Cumulative daily confirmed cases in "
  # Plot cumulative cases
  day <- c(1:ncol(actual_data))
  cases <- actual_data[1,]
  
  start_date <- format(as.Date(colnames(actual_data)[1], format = "X%m.%d.%y"), "%b %d")
  end_date   <- format(as.Date(colnames(actual_data)[ncol(actual_data)], format = "X%m.%d.%y"), "%b %d")
  if (parameters != c("??")) parameters <- round(parameters, 3)
  estimates <- paste("(", parameters[1], sep="")
  for (i in 2:length(parameters)) {
    estimates <- paste(estimates, ", ", parameters[i], sep="")
  }
  estimates <- paste(estimates, ")", sep="")
  
  plot(day, cases, pch = 16, title(main = paste(title_start,
                                                country, " (", start_date, " - ", end_date, ")\n",
                                                "Estimates: ", estimates, sep = "")))
  lines(day, fit, pch = 1, col = "green", lwd = 2)
  legend('topleft', legend=c("Actual data", "Estimation"), pch = c(16, NA), lty = c(NA, 1), 
         col = c('black', 'green'), lwd = 2, inset = 0.03)
  
  
  # Plot daily increase
  title_start <- "Daily new confirmed cases in "
  actual_data <- as.data.frame(t(cumulative_to_increase(actual_data)))
  fit <- cumulative_to_increase(fit)
  day <- c(1:ncol(actual_data))
  cases <- actual_data[1,]
  
  if (parameters != c("??")) parameters <- round(parameters, 3)
  estimates <- paste("(", parameters[1], sep="")
  for (i in 2:length(parameters)) {
    estimates <- paste(estimates, ", ", parameters[i], sep="")
  }
  estimates <- paste(estimates, ")", sep="")
  
  plot(day, cases, pch = 16, title(main = paste(title_start,
                                                country, " (", start_date, " - ", end_date, ")\n",
                                                "Estimates: ", estimates, sep = "")))
  lines(day, fit, pch = 1, col = "green", lwd = 2)
  legend('topright', legend=c("Actual data", "Estimation"), pch = c(16, NA), lty = c(NA, 1), 
         col = c('black', 'green'), lwd = 2, inset = 0.03)
}

plot_fits <- function(actual_data, fits, legend_entries, colors = rainbow(length(fits))) {
  nr_fits <- length(fits)
  orig <- as.numeric(actual_data[1,])
  
  ndays <- length(orig)
  days <- c(1:ndays)
  
  xrange <- range(days)
  yrange <- range(orig)
  plot(xrange, yrange, xlab="Day", ylab="Cumulative confirmed cases", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  points(days, orig, pch = 16)
  for (i in 1:nr_fits) {
    lines(days, fits[[i]], col = colors[i], lwd = 2)
  }
  legend('topleft', legend = c("Actual data", legend_entries), 
         lty = c(NA, rep(1, nr_fits)), pch = c(16, rep(NA, nr_fits)), 
         col = c('black', colors), lwd = 2, cex = 1.1, inset = 0.03)
  
  
  incr <- cumulative_to_increase(orig)
  days <- c(1:(ndays-1))
  xrange <- range(days)
  yrange <- range(incr)
  plot(xrange, yrange, xlab="Day", ylab="Daily increase in cases", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  points(days, incr, pch = 16)
  for (i in 1:nr_fits) {
    lines(days, cumulative_to_increase(fits[[i]]), col = colors[i], lwd = 2)
  }
  legend('topright', legend = c("Actual data", legend_entries), 
         lty = c(NA, rep(1, nr_fits)), pch = c(16, rep(NA, nr_fits)), 
         col = c('black', colors), lwd = 2, cex = 1.1, inset = 0.03)
}

temp_plot <- function() {
  actual_data <- daily_cases_for_country("Netherlands")
  N0    <- actual_data[1,1]
  ndays <- ncol(actual_data)
  
  fit4  <- generic(N0, c(0.955, 69177, 0.782, 1.889, 5.173), ndays)
  fit5  <- generic(N0, c(1.973, 748584, 0.795, 0.719, 24.46), ndays)
  fit6  <- generic(N0, c(0.699, 303416, 0.863, 1.117, 26.78), ndays)
  fit7  <- generic(N0, c(1.069, 257882, 0.827, 0.951, 22.5), ndays)
  fit8  <- generic(N0, c(1.030, 144038, 0.782, 1.763, 22.37), ndays)
  fit9  <- generic(N0, c(1.232, 129666, 0.75, 2.183, 28.42), ndays)
  fit10 <- generic(N0, c(0.481, 150656, 0.878, 1.804, 28.32), ndays)
  fit11 <- generic(N0, c(0.626, 130507, 0.838, 2.027, 26.54), ndays)
  fit12 <- generic(N0, c(0.916, 127327, 0.788, 2.110, 26.18), ndays)
  
  fits <- list(fit4, fit5, fit9, fit12)
  legend_entries <- paste("Estimation (w=", c(4,5,9,12), ")", sep = "")
  
  plot_fits(actual_data, fits, legend_entries, c(viridis(4)[2], 'orange', viridis(4)[3], viridis(4)[4]))
}
temp_plot()



fit_function <- function(func, df, min, max, steps) {
  nr_parameters <- length(min)
  if(length(max) != nr_parameters || length(steps) != nr_parameters) {
    stop(paste("[!] min defined for", length(min), "parameters;",
               "max defined for", length(max), "parameters;",
               "steps defined for", length(steps), "parameters.",
               "All should be the same."))
  }
  
  days <- ncol(df)
  N0 <- df[1,1]
  
  min_dist <- .Machine$integer.max
  best_values <- rep(.Machine$integer.max, nr_parameters)
  
  nr_steps <- prod(steps + 1)
  step_sizes <- (max - min) / steps
  counter <- integer(nr_parameters)
  for(i in 1:nr_steps) {
    values <- min + counter * step_sizes
    dist <- statistical_distance(df, func(N0, values, days))
    if(is.na(dist)) {
      print(N0)
      print(values)
      print(days)
    }
    if (dist < min_dist) {
      min_dist = dist
      best_values = values
    }
    counter <- update_counters(counter, steps)
  }
  
  output <- list()
  output[[1]] <- best_values
  output[[2]] <- min_dist
  return(output)
}

fit_function_to_country <- function(func, country, min, max, search_width, search_depth, 
                                    print_progress=FALSE, make_plot=FALSE) {
  start_time <- Sys.time()
  nr_parameters <- length(min)
  if(length(max) != nr_parameters) {
    stop(paste("[!] `min' defined for", length(min), "parameters,",
               "`max' defined for", length(max), "parameters. (Should be the same)"))
  }
  outer_min <- min
  outer_max <- max
  
  df <- daily_cases_for_country(country)
  steps <- rep(search_width, nr_parameters)
  
  min_dist <- .Machine$integer.max
  best_values <- rep(.Machine$integer.max, nr_parameters)
  for (i in 1:search_depth) {
    if(print_progress) {
      print(i)
    }
    result <- fit_function(func, df, min, max, steps)
    best_values <- result[[1]]
    min_dist <- result[[2]]
    
    new_bounds <- update_search_bounds(min, max, steps, best_values, outer_min, outer_max)
    min <- new_bounds[[1]]
    max <- new_bounds[[2]]
  }
  
  end_time <- Sys.time()
  print(paste("Estimation took", end_time - start_time, "seconds"))
  
  if(make_plot) {
    fit <- func(df[1,1], best_values, ncol(df))
    plot_fit(df, fit, country, func, best_values)
    plot_fit(as.data.frame(t(cumulative_to_increase(df))), cumulative_to_increase(fit), 
             country, func, best_values, "Daily new confirmed cases in ")
  }
  
  output <- list()
  output[[1]] <- best_values
  output[[2]] <- min_dist
  return(output)
}



generic_poisson_regression <- function(cum_daily, K, beta) {
  if (beta == 0) return(NaN)
  inc_daily <- cumulative_to_increase(cum_daily)
  cum_daily <- cum_daily[1:length(cum_daily) - 1]
  logterm   <- log(cum_daily)
  logfrac   <- log(1 - power(cum_daily / K, beta))
  logfrac[which(is.infinite(logfrac))] <- NaN
  
  df <- data.frame(cum_daily, inc_daily, logterm, logfrac)
  names(df) <- c("cumulative", "increase", "log", "logfrac")
  
  fit <- glm(increase ~ logterm + logfrac, data = df, family = poisson())
  return(as.numeric(fit[[1]]))
}

generic_poisson_regression_grid <- function(cum_daily, min, max, nstep) {
  K_min   <- min[1]
  K_max   <- max[1]
  K_nstep <- nstep[1]
  if (K_nstep == 0) {
    K_step <- 0
  } else {
    K_step <- (K_max - K_min) / K_nstep
  }
  
  beta_min   <- min[2]
  beta_max   <- max[2]
  beta_nstep <- nstep[2]
  if (beta_nstep == 0) {
    beta_step <- 0
  } else {
    beta_step <- (beta_max - beta_min) / beta_nstep
  }
  
  
  min_dist <- .Machine$integer.max
  best_values <- rep(.Machine$integer.max, 5)
  for(i in 0:K_nstep) {
    K <- K_min + i * K_step
    for (j in 0:beta_nstep) {
      beta <- beta_min + j * beta_step
      result <- generic_poisson_regression(as.numeric(cum_daily), K, beta)
      
      if (!is.na(result)) {
        r     <- exp(result[1])
        alpha <- result[2]
        gamma <- result[3]
        
        curr_values <- c(r,K,alpha,beta,gamma)
        fit <- generic(cum_daily[1,1], curr_values, ncol(cum_daily))
        
        dist <- chisqr(cum_daily, fit)
        if (dist < min_dist) {
          min_dist <- dist
          best_values <- curr_values
        }
      }
    }
  }
  
  return(best_values)
}

generic_poisson_regression_country <- function(country, min, max, nstep, search_depth, print_progress = FALSE) {
  cum_daily <- daily_cases_for_country(country)[1,]
  
  outer_min <- min
  outer_max <- max
  best_values <- rep(.Machine$integer.max, 5)
  for (i in 1:search_depth) {
    if(print_progress) {
      print(i)
    }
    
    best_values <- generic_poisson_regression_grid(cum_daily, min, max, nstep)
    best_K <- best_values[2]
    best_beta <- best_values[4]
    new_bounds  <- update_search_bounds(min, max, nstep, c(best_K, best_beta), outer_min, outer_max)
    min <- new_bounds[[1]]
    max <- new_bounds[[2]]
  }
  
  print(best_values)
  fit <- generic(cum_daily[1,1], best_values, ncol(cum_daily))
  print(chisqr(cum_daily, fit))
  plot_fit(cum_daily, fit, country, generic, best_values)
}



# Example queries
##################################################################
fit_function_to_country(logistic, "Netherlands", c(0), c(10), 10, 10, TRUE, TRUE)
fit_function_to_country(verhulst, "Netherlands", c(0,1), c(1,1000000), 50, 10, TRUE, TRUE)
fit_function_to_country(verhulst, "US", c(0,1), c(3,1000000000),5, 8, TRUE, TRUE)
fit_function_to_country(generic,  "Global", c(0,1,0,0,0), c(0.5,1000000000,3,3,3),10, 8, TRUE, TRUE)
fit_function_to_country(generic,  "Netherlands", c(0,50000,0,0,0), c(2,250000,3,3,3), 3, 10, TRUE, TRUE)

parameter_plots(generic, 1, list(c(0.15, 0.15, 0.15), 
                                 c(1000, 1000, 1000), 
                                 c(1, 1, 1),
                                 c(1, 1, 1),
                                 c(1, 1.5, 2)), 100, "gamma", 5)
