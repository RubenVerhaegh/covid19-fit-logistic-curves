library(MASS)
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



statistical_distance <- function(vec, df) {
  len <- length(vec)
  if (len != ncol(df)) {
    stop("Tried to compare to lists of different lengthts")
  }
  
  dist <- 0
  for (i in 1:length(vec)) {
    observed <- df[1,i]
    expected <- vec[i]
    dist <- dist + abs(observed - expected)
  }
  return(dist)
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

update_search_bounds <- function(prev_min, prev_max, prev_steps, best_values){
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
    
    if (best == min) {
      new_min[i] <- best
      new_max[i] <- best + step_size
    } else if (best == max) {
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
plot_fit <- function(actual_data, fit, country = "??", func = NA, parameters = c("??"), title_start = "cumulative daily confirmed cases in ") {
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
  points(day, fit, pch = 1)
  legend(0, max(cases), legend=c("Actual data", "Estimation"), pch = c(16, 1))
}



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
    dist <- statistical_distance(func(N0, values, days), df)
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
                                    print_progress=FALSE, make_plot=FALSE, perform_test = FALSE) {
  start_time <- Sys.time()
  nr_parameters <- length(min)
  if(length(max) != nr_parameters) {
    stop(paste("[!] `min' defined for", length(min), "parameters,",
               "`max' defined for", length(max), "parameters. (Should be the same)"))
  }
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
    
    new_bounds <- update_search_bounds(min, max, steps, best_values)
    min <- new_bounds[[1]]
    max <- new_bounds[[2]]
#    print("min:")
#    print(min)
#    print("max:")
#   print(max)
#    print("---------------------------")
  }
  
  end_time <- Sys.time()
  print(paste("Estimation took", end_time - start_time, "seconds"))
  
  if(make_plot || perform_test) {
    fit <- func(df[1,1], best_values, ncol(df))
    if(make_plot) {
      plot_fit(df, fit, country, func, best_values)
      plot_fit(as.data.frame(t(cumulative_to_increase(df))), cumulative_to_increase(fit), country, func, best_values, "Daily new confirmed cases in ")
    }
    if(perform_test) test_fit_cumulative(df, fit)
  }
  
  output <- list()
  output[[1]] <- best_values
  output[[2]] <- min_dist
  return(output)
}


# Example queries
##################################################################
fit_function_to_country(logistic, "Netherlands", c(0), c(10), 10, 10, TRUE, TRUE, TRUE)
fit_function_to_country(verhulst, "Netherlands", c(0,1), c(1,1000000), 50, 10, TRUE, TRUE, TRUE)
fit_function_to_country(verhulst, "US", c(0,1), c(3,1000000000),5, 8, TRUE, TRUE, TRUE)
fit_function_to_country(generic,  "Global", c(0,1,0,0,0), c(0.5,1000000000,3,3,3),10, 8, TRUE, TRUE, TRUE)
fit_function_to_country(generic,  "Netherlands", c(0,50000,0,0,0), c(2,1000000,3,3,3), 10, 8, TRUE, TRUE, TRUE)

parameter_plots(generic, 1, list(c(0.15, 0.15, 0.15), 
                                 c(1000, 1000, 1000), 
                                 c(1, 1, 1),
                                 c(1, 1, 1),
                                 c(1, 1.5, 2)), 100, "gamma", 5)
