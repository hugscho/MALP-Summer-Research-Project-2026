library(MASS)
library(sn)
setwd('/Users/Hugo/Desktop/MALP-project-summer-2026/MALR project 2026/Writeup')

get_mean_vector <- function(y, x){ #takes two variables, y, and x, and creates a mean vector mew = [mew_y, mew_x]
  y_mean <- mean(y)
  x_mean <- colMeans(as.matrix(x))
  
  return(c(y_mean, x_mean))
}

get_cov_matrix <- function(y, x){ #takes a single variate y, and multivariatevariate x and generates covariance matrices s_yy, s_yx, s_xy, s_xx
  x <- as.matrix(x)
  n <- nrow(x)
  
  y_centered <- y - mean(y)
  x_centered <- sweep(x, 2, colMeans(x), '-')
  
  s_yy <- (1/n) * sum(y_centered^2)
  s_yx <- (1/n) * t(y_centered) %*% x_centered
  s_xy <- t(s_yx)
  s_xx <- (1/n) * t(x_centered) %*% x_centered
  
  return(list(s_yy = s_yy, s_yx = s_yx, s_xy = s_xy, s_xx = s_xx))
}

pcc <- function(y, x){ #takes single variate y and multivariate x and computes PCC between them
  cm <- get_cov_matrix(y, x)
  pearson <- sqrt(cm$s_yx %*% solve(cm$s_xx) %*% cm$s_xy / cm$s_yy)
  
  return(as.numeric(pearson))
}

ccc <- function(y, x){  #takes single variate y and single variate x and computes CCC between them
  if (ncol(as.matrix(x)) > 1){
    stop('ccc is only defined for univariate x')
  }
  
  cm <- get_cov_matrix(y, x)
  lin <- 2 * cm$s_yx / (cm$s_xx + cm$s_yy + (mean(x) - mean(y))**2)
  
  return(as.numeric(lin))
}

LSLP <- function(y, x){ #takes 
  means <- get_mean_vector(y, x)
  cm <- get_cov_matrix(y, x)
  
  intercept <- means[1] - as.numeric(cm$s_yx %*% solve(cm$s_xx) %*% means[-1])
  slope <- cm$s_yx %*% solve(cm$s_xx)
  
  return(list(intercept = intercept, slope = slope))
}

MALP <- function(y, x){
  means <- get_mean_vector(y, x)
  cm <- get_cov_matrix(y, x)
  
  lslp <- LSLP(y, x)
  gamma <- as.numeric(sqrt(cm$s_yx %*% solve(cm$s_xx) %*% cm$s_xy) / sqrt(cm$s_yy))
  
  intercept <- (lslp$intercept - means[1]) / gamma + means[1]
  slope <- lslp$slope / gamma
  
  return(list(intercept = intercept, slope = slope))
}

evaluator <- function(y, x, predictor){
  
  if(predictor == 'malp'){
    func <- MALP(y, x)
  } else if(predictor == 'lslp'){
    func <- LSLP(y, x)
  } else{
    stop('predictor must be malp or lslp')
  }
  
  new_x <- as.matrix(x) %*% matrix(func$slope, ncol = 1) + func$intercept
  
  PCC <- pcc(y, new_x)
  CCC <- ccc(y, new_x)
  MSE <- mean((new_x - y)**2)
  
  return(list(pcc = PCC, ccc = CCC, mse = MSE))
}

evaluator_five_fold <- function(y, x, predictor){
  
  x <- as.matrix(x)
  folds <- sample(rep(1:5, length.out = nrow(x)))
  results <- vector('list', 5)
  
  for (i in 1:5){
    train_y <- y[folds != i]
    test_y <- y[folds == i]
    train_x <- x[folds != i, drop = FALSE]
    test_x  <- x[folds == i, drop = FALSE]
    
    if(predictor == 'malp'){
      func <- MALP(train_y, train_x)
    } else if(predictor == 'lslp'){
      func <- LSLP(train_y, train_x)
    } else{
      stop('predictor must be malp or lslp')
    }
    
    new_x <- as.matrix(test_x) %*% matrix(func$slope, ncol = 1) + func$intercept
    
    results[[i]] <- list (
      PCC = pcc(test_y, new_x),
      CCC = ccc(test_y, new_x),
      MSE = mean((new_x - test_y)**2)
    )
  }
  
  mean_pcc <- mean(sapply(results, function(r) r$PCC))
  mean_ccc <- mean(sapply(results, function(r) r$CCC))
  mean_mse <- mean(sapply(results, function(r) r$MSE))
  
  return(list(pcc = mean_pcc, ccc = mean_ccc, mse = mean_mse))
}

evaluator_ten_fold <- function(y, x, predictor){
  
  x <- as.matrix(x)
  folds <- sample(rep(1:10, length.out = nrow(x)))
  results <- vector('list', 10)
  
  for (i in 1:10){
    train_y <- y[folds != i]
    test_y <- y[folds == i]
    train_x <- x[folds != i, , drop = FALSE]
    test_x  <- x[folds == i, , drop = FALSE]
    
    if(predictor == 'malp'){
      func <- MALP(train_y, train_x)
      error <- generate_error(test_y, test_x, 'malp')
    } else if(predictor == 'lslp'){
      func <- LSLP(train_y, train_x)
      error <- generate_error(test_y, test_x, 'lslp')
    } else{
      stop('predictor must be malp or lslp')
    }
    
    new_x <- as.matrix(test_x) %*% matrix(func$slope, ncol = 1) + func$intercept
    n_error <- length(error)
    iqr_error <- error[(ceiling(n_error/4) + 1):floor(3*n_error/4)]
    oqr_error <- c(error[1:ceiling(n_error/4)], error[(floor(3*n_error/4) + 1):n_error])
    
    
    results[[i]] <- list (
      PCC = pcc(test_y, new_x),
      CCC = ccc(test_y, new_x),
      MSE = mean((new_x - test_y)**2), 
      MSE_IQ = mean(iqr_error), 
      MSE_OQ = mean(oqr_error)
    )
  }
  
  mean_pcc <- mean(sapply(results, function(r) r$PCC))
  mean_ccc <- mean(sapply(results, function(r) r$CCC))
  mean_mse <- mean(sapply(results, function(r) r$MSE))
  mean_mse_iq <- mean(sapply(results, function(r) r$MSE_IQ))
  mean_mse_oq <- mean(sapply(results, function(r) r$MSE_OQ))
  
  return(list(pcc = mean_pcc, ccc = mean_ccc, mse = mean_mse, mse_iq = mean_mse_iq, mse_oq = mean_mse_oq))
}

gen_interval <- function(y, x, type, predictor, alpha = 0.05){
  mv <- get_mean_vector(y, x)
  cm <- get_cov_matrix(y, x)
  n <- length(y)
  p <- ncol(as.matrix(x))
  x0 <- seq(150, 350, length.out = 200)
  t <- qt(1 - alpha/2, n - p - 1)
  xbar <- colMeans(as.matrix(x))
  diffs <- sweep(as.matrix(x0), 2, xbar, '-')
  quad <- rowSums((diffs %*% solve(cm$s_xx)) * diffs)
  lslp <- LSLP(y, x)
  new_x <- as.matrix(x) %*% matrix(lslp$slope, ncol = 1) + lslp$intercept
  gamma <- pcc(y, new_x)
  s <- sqrt(cm$s_yy * (1 - gamma ** 2))
  s <- sqrt(cm$s_yy * (1 - gamma ** 2) * (n - 1) / (n - p - 1))

  if (predictor == 'LSLP'){
    if (type == 'CI'){
      fit <- as.matrix(x0) %*% matrix(lslp$slope, ncol = 1) + lslp$intercept
      margin <- t * s * sqrt(1/n * (1 + quad))
    }else if (type == 'PI') {
      fit <- as.matrix(x0) %*% matrix(lslp$slope, ncol = 1) + lslp$intercept
      margin <- t * s * sqrt(1 + 1/n * (1 + quad))
    }
  }else if (predictor == 'MALP') {
    malp <- MALP(y, x)
    D2ma <- 2/(1 + gamma) + 1/gamma**2 * quad - ((1 - gamma**2)/(cm$s_yy * gamma**4)) * as.numeric(diffs %*%(solve(cm$s_xx) %*% cm$s_xy)) ** 2
    if (type == 'CI'){
      fit <- as.matrix(x0) %*% matrix(malp$slope, ncol = 1) + malp$intercept
      margin <- t * s * sqrt((1/n) * (D2ma))
    }else if (type == 'PI') {
      fit <- as.matrix(x0) %*% matrix(lslp$slope, ncol = 1) + lslp$intercept
      margin <- t * s * sqrt(1 + (1/n) * (D2ma))
    }
  }
  return(list(x0 = x0, fit = fit, margin = margin))
}

generate_sample <- function(n, means, cm, rho = NULL){
  if (is.null(rho)){
    sigma <- rbind(cbind(cm$s_yy, cm$s_yx), cbind(cm$s_xy, cm$s_xx))
  
  }else{
    rho_dir <- as.numeric(cm$s_yx) / (sqrt(cm$s_yy) * sqrt(diag(cm$s_xx)))
    corm_x <- cov2cor(cm$s_xx)
    
    targetR2 <- 1 - (1 - rho ** 2) * (n - 1) / (n - (length(means) - 1) - 1)
    
    if (targetR2 <= 0 || targetR2 >= 1) {
      stop("Requested R2 is not achievable for this n and number of predictors.")
    }
    
    unscaledR2 <- as.numeric(t(rho_dir) %*% solve(corm_x) %*% rho_dir)
    rho_yx <- sqrt(targetR2 / unscaledR2) * rho_dir
    
    corm <- rbind(cbind(1, matrix(rho_yx, nrow = 1)), cbind(rho_yx, corm_x))
    sds <- c(sqrt(cm$s_yy), sqrt(diag(cm$s_xx))) 
    sigma <- diag(sds) %*% corm %*% diag(sds)
    
  }

  data <- mvrnorm(n, means, sigma)
  data_y <- data[, 1]
  data_x <- data[, 2:ncol(data), drop = FALSE]
  
  return(list(y = matrix(c(data_y), ncol = 1), x = matrix(c(data_x), ncol = ncol(data) - 1)))
}

generate_error <- function(y, x, predictor){
  if(predictor == 'malp'){
    func <- MALP(y, x)
  } else if(predictor == 'lslp'){
    func <- LSLP(y, x)
  } else{
    stop('predictor must be malp or lslp')
  }
  full_datapoints <- cbind(y, x)
  sorted_datapoints <- full_datapoints[order(full_datapoints[, 1]), ]
  y <- sorted_datapoints[, 1]
  x <- sorted_datapoints[, 2:ncol(sorted_datapoints), drop = FALSE]
  
  new_x <- as.matrix(x) %*% matrix(func$slope, ncol = 1) + func$intercept

  squared_error <- (y - new_x) ** 2

  return(squared_error)
}

generate_sample_skewed <- function(n, means, cm, rho = NULL, skew){
  if (is.null(rho)){
    sigma <- rbind(cbind(cm$s_yy, cm$s_yx), cbind(cm$s_xy, cm$s_xx))
  
  }else{
    rho_dir <- as.numeric(cm$s_yx) / (sqrt(cm$s_yy) * sqrt(diag(cm$s_xx)))
    corm_x <- cov2cor(cm$s_xx)
    
    targetR2 <- 1 - (1 - rho ** 2) * (n - 1) / (n - (length(means) - 1) - 1)
    
    if (targetR2 <= 0 || targetR2 >= 1) {
      stop("Requested R2 is not achievable for this n and number of predictors.")
    }
    
    unscaledR2 <- as.numeric(t(rho_dir) %*% solve(corm_x) %*% rho_dir)
    rho_yx <- sqrt(targetR2 / unscaledR2) * rho_dir
    
    corm <- rbind(cbind(1, matrix(rho_yx, nrow = 1)), cbind(rho_yx, corm_x))
    sds <- c(sqrt(cm$s_yy), sqrt(diag(cm$s_xx))) 
    sigma <- diag(sds) %*% corm %*% diag(sds)
    
  }

  sv <- c(0, rep(skew, length(means) - 1))
  sigma <- round(sigma, 10)

  data <- rmsn(n, means, sigma, sv)
  data_y <- data[, 1]
  data_x <- data[, 2:ncol(data), drop = FALSE]
  
  return(list(y = matrix(c(data_y), ncol = 1), x = matrix(c(data_x), ncol = ncol(data) - 1)))
}


##Scripts for studies  

Kim_etal_eye_replication <- function(){

  load('eye.rda')

  plot_LPs <- function(y, x, xlim, ylim, name){
    lslp <- LSLP(y, x)
    malp <- MALP(y, x)
    
    plot(x, y, xlim = xlim, ylim = ylim, xlab = 'cirrus', ylab = 'stratus', main = name)
    abline(a = lslp$intercept, b = lslp$slope, col = 'red') 
    abline(a = malp$intercept, b = malp$slope, col = 'blue')
  }

  par(mfrow = c(3, 2))

  x <- eye$Cirrus
  y <- eye$Stratus

  print('initial correlations between stratus and cirrus')
  print(pcc(y, x))
  print(ccc(y, x))

  #Corrections made to cirrus like in EJS 2026
  x_i <- x + (mean(y) - mean(x)) #equates the means of the variables 
  x_ii <- 0.76 * x - 0.51 #

  #correlations p, p_c between stratus and i_cirrus
  print('corrected correlations i) between stratus and cirrus')
  print(pcc(y, x_i))
  print(ccc(y, x_i))


  #correlations p, p_c between stratus and ii_cirrus
  print('corrected correlations ii) between stratus and cirrus')
  print(pcc(y, x_ii))
  print(ccc(y, x_ii))


  #predictor functions for OS
  lslp <- LSLP(y[eye$Eye == 'OS'], x[eye$Eye == 'OS'])
  malp <- MALP(y[eye$Eye == 'OS'], x[eye$Eye == 'OS'])
  print('predictor functions for OS')
  cat('MALP: ',malp[[1]],'x_0 + ',malp[[2]], '\n')
  cat('LSLP: ',lslp[[1]],'x_0 + ',lslp[[2]], '\n')
  plot_LPs(y[eye$Eye == 'OS'], x[eye$Eye == 'OS'], c(150, 350), c(100, 300), 'MALP and LSLP for OS values')


  #predictor functions for OD
  lslp <- LSLP(y[eye$Eye == 'OD'], x[eye$Eye == 'OD'])
  malp <- MALP(y[eye$Eye == 'OD'], x[eye$Eye == 'OD'])
  print('predictor functions for OD')
  cat('MALP: ',malp[[1]],'x_0 + ',malp[[2]], '\n')
  cat('LSLP: ',lslp[[1]],'x_0 + ',lslp[[2]], '\n')
  plot_LPs(y[eye$Eye == 'OD'], x[eye$Eye == 'OD'], c(150, 350), c(100, 300), 'MALP and LSLP for OD values')

  data <- data.frame()

  malp_os_ill <- evaluator(y[eye$Eye == 'OS'], x[eye$Eye == 'OS'], 'malp')
  lslp_os_ill <- evaluator(y[eye$Eye == 'OS'], x[eye$Eye == 'OS'], 'lslp')
  malp_od_ill <- evaluator(y[eye$Eye == 'OD'], x[eye$Eye == 'OD'], 'malp')
  lslp_od_ill <- evaluator(y[eye$Eye == 'OD'], x[eye$Eye == 'OD'], 'lslp')
  malp_os_5 <- evaluator_five_fold(y[eye$Eye == 'OS'], x[eye$Eye == 'OS'], 'malp')
  lslp_os_5 <- evaluator_five_fold(y[eye$Eye == 'OS'], x[eye$Eye == 'OS'], 'lslp')
  malp_od_5 <- evaluator_five_fold(y[eye$Eye == 'OD'], x[eye$Eye == 'OD'], 'malp')
  lslp_od_5 <- evaluator_five_fold(y[eye$Eye == 'OD'], x[eye$Eye == 'OD'], 'lslp')

  data <- rbind(data, data.frame(pcc = lslp_os_ill$pcc, ccc = lslp_os_ill$ccc, mse = lslp_os_ill$mse, predictor = 'lslp', eye = 'os', type = 'illustrative'))
  data <- rbind(data, data.frame(pcc = malp_os_ill$pcc, ccc = malp_os_ill$ccc, mse = malp_os_ill$mse, predictor = 'malp', eye = 'os', type = 'illustrative'))
  data <- rbind(data, data.frame(pcc = lslp_od_ill$pcc, ccc = lslp_od_ill$ccc, mse = lslp_od_ill$mse, predictor = 'lslp', eye = 'od', type = 'illustrative'))
  data <- rbind(data, data.frame(pcc = malp_od_ill$pcc, ccc = malp_od_ill$ccc, mse = malp_od_ill$mse, predictor = 'malp', eye = 'od', type = 'illustrative'))
  data <- rbind(data, data.frame(pcc = lslp_os_5$pcc, ccc = lslp_os_5$ccc, mse = lslp_os_5$mse, predictor = 'lslp', eye = 'os', type = 'five fold'))
  data <- rbind(data, data.frame(pcc = malp_os_5$pcc, ccc = malp_os_5$ccc, mse = malp_os_5$mse, predictor = 'malp', eye = 'os', type = 'five fold'))
  data <- rbind(data, data.frame(pcc = lslp_od_5$pcc, ccc = lslp_od_5$ccc, mse = lslp_od_5$mse, predictor = 'lslp', eye = 'od', type = 'five fold'))
  data <- rbind(data, data.frame(pcc = malp_od_5$pcc, ccc = malp_od_5$ccc, mse = malp_od_5$mse, predictor = 'malp', eye = 'od', type = 'five fold'))

  print(data)

  #generate the conversion formula and compare to Abedi et.al.

  malp_all <- MALP(y[eye$Eye == 'OS' | eye$Eye == 'OD'], x[eye$Eye == 'OS' | eye$Eye == 'OD'])
  print('predictor function from pooled OS and OD')
  cat('MALP: ',malp_all$slope,'x_0 + ',malp_all$intercept, '\n')

  print('predictor proposed by Abedi et.al')
  cat('MALP: ',0.760,'x_0 + ',0.510, '\n')

  plot(eye$Cirrus, eye$Stratus, xlim = c(150, 350), ylim = c(120, 260), main = 'MALPs with CIs and PIs')
  abline(a = 0.510, b = 0.760, col = 'red') 
  abline(a = malp_all$intercept, b = malp_all$slope, col = 'blue')
  LSLPpi <- gen_interval(y, x, 'PI', 'LSLP')
  lines(LSLPpi$x0, LSLPpi$fit + LSLPpi$margin, col = 'red', lty = 2)
  lines(LSLPpi$x0, LSLPpi$fit - LSLPpi$margin, col = 'red', lty = 2)
  MALPpi <- gen_interval(y, x, 'PI', 'MALP')
  lines(MALPpi$x0, MALPpi$fit + MALPpi$margin, col = 'blue', lty = 2)
  lines(MALPpi$x0, MALPpi$fit - MALPpi$margin, col = 'blue', lty = 2)
  MALPci <- gen_interval(y, x, 'CI', 'MALP')
  lines(MALPci$x0, MALPci$fit + MALPci$margin, col = 'green', lty = 2)
  lines(MALPci$x0, MALPci$fit - MALPci$margin, col = 'green', lty = 2)
  legend("topleft", legend = c("Abedi MALP", "Computed MALP","LSLP PI","MALP PI","MALP CI"), col = c("red", "blue", "red", "blue", "green"), lty = c(1, 1, 2, 2, 2), lwd = 2, bty = "n")

}

simstud1 <- function(dist){

  load('eye.rda')
  load('bodyFat.rda')
  
  data_raw <- vector('list', 2)
  
  y_1 <- eye$Stratus
  x_1 <- eye$Cirrus
  data_raw[[1]] <- list(y = y_1, x = x_1, name = 'eye data')
  
  y_2 <- bodyFat$PBF
  x_2 <- as.matrix(subset(bodyFat, select = c(Age)))
  data_raw[[2]] <- list(y = y_2, x = x_2, name = 'body fat data, p = 1')
  
  y_3 <- bodyFat$PBF
  x_3 <- as.matrix(subset(bodyFat, select = c(Age, WGT, HGT)))
  data_raw[[3]] <- list(y = y_3, x = x_3, name = 'body fat data, p = 3')
  
  y_4 <- bodyFat$PBF
  x_4 <- as.matrix(subset(bodyFat, select = c(Age, WGT, HGT, NCK, CST)))
  data_raw[[4]] <- list(y = y_4, x = x_4, name = 'body fat data, p = 5')
  #start data collection based on mv, cm, rho, n
  
  par(mfrow = c(4, 2))
  rho <- c(0.3, 0.5, 0.7, 0.9, 0.95)
  n <- c(100, 300, 1000)
  predictor <- c('lslp', 'malp')
  
  
  for (h in 1:length(data_raw)){
    y <- data_raw[[h]]$y
    x <- data_raw[[h]]$x
    data_name <- data_raw[[h]]$name
    mv <- get_mean_vector(y, x)
    cm <- get_cov_matrix(y, x)

    sim_stud_data <- data.frame()
    
    for (i in 1:length(rho)){
      
      for(j in 1:length(n)){
        
        for (k in 1:length(predictor)){
          n_reps <- 1000
          skew <- 20
          temp_data <- vector('list', n_reps)
          
          for (l in 1:n_reps){
            if (dist == 'normal'){
              sample <- generate_sample(n[j], mv, cm, rho[i])
            }else if (dist == 'skewed'){
               sample <- generate_sample_skewed(n[j], mv, cm, rho[i], skew)
            }else{
              stop('dist must be "normal" or "skewed"')
            }
            temp_data[[l]] <- evaluator(sample[[1]], sample[[2]], predictor[k])
            
          }
          avg_pcc <- mean(sapply(temp_data, function(x) x$pcc))
          avg_ccc <- mean(sapply(temp_data, function(x) x$ccc))
          avg_mse <- mean(sapply(temp_data, function(x) x$mse))
          
          sim_stud_data <- rbind(sim_stud_data, data.frame(
            predictor = predictor[k], 
            rho = rho[i], 
            n = n[j], 
            pcc = avg_pcc, 
            ccc = avg_ccc, 
            mse = avg_mse))
        }
      }
    }
    
    print(sim_stud_data)
    
    ###plotting simulated data
    
    plot(x = NULL, type = 'o', main = paste('Simulation using ', data_name, ' metrics - CCC'), xlim = c(0.2, 1), ylim = c(0, 1), xlab = 'pcc', ylab = 'ccc')
    points(sim_stud_data$rho[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[1]], sim_stud_data$ccc[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[1]], type = 'b',lty = 3, col = 'red')
    points(sim_stud_data$rho[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[2]], sim_stud_data$ccc[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[2]], type = 'b',lty = 2, col = 'red')
    points(sim_stud_data$rho[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[3]], sim_stud_data$ccc[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[3]], type = 'b',lty = 1, col = 'red')
    
    points(sim_stud_data$rho[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[1]], sim_stud_data$ccc[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[1]], type = 'b',lty = 3, col = 'blue')
    points(sim_stud_data$rho[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[2]], sim_stud_data$ccc[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[2]], type = 'b',lty = 2, col = 'blue')  
    points(sim_stud_data$rho[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[3]], sim_stud_data$ccc[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[3]], type = 'b',lty = 1, col = 'blue')
    
    legend(0.2, 1, legend = c(paste('lslp', n[[1]]), paste('lslp', n[[2]]), paste('lslp', n[[3]]), paste('malp', n[[1]]), paste('malp', n[[2]]), paste('malp', n[[3]])), col = c('red', 'red', 'red', 'blue', 'blue', 'blue'),lty = c(3, 2, 1, 3, 2, 1), cex = 0.8)
    
    plot(x = NULL, type = 'o', main = paste('Simulation using ', data_name, ' metrics - MSE'), xlim = c(0.2, 1), ylim = c(0, max(sim_stud_data$mse)), xlab = 'pcc', ylab = 'mse')
    points(sim_stud_data$rho[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[1]], sim_stud_data$mse[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[1]], type = 'b',lty = 3, col = 'red')
    points(sim_stud_data$rho[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[2]], sim_stud_data$mse[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[2]], type = 'b',lty = 2, col = 'red')
    points(sim_stud_data$rho[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[3]], sim_stud_data$mse[sim_stud_data$predictor == 'lslp' & sim_stud_data$n == n[3]], type = 'b',lty = 1, col = 'red')
    
    points(sim_stud_data$rho[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[1]], sim_stud_data$mse[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[1]], type = 'b',lty = 3, col = 'blue')
    points(sim_stud_data$rho[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[2]], sim_stud_data$mse[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[2]], type = 'b',lty = 2, col = 'blue')  
    points(sim_stud_data$rho[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[3]], sim_stud_data$mse[sim_stud_data$predictor == 'malp' & sim_stud_data$n == n[3]], type = 'b',lty = 1, col = 'blue')
    
    legend(0.8, max(sim_stud_data$mse), legend = c(paste('lslp', n[[1]]), paste('lslp', n[[2]]), paste('lslp', n[[3]]), paste('malp', n[[1]]), paste('malp', n[[2]]), paste('malp', n[[3]])), col = c('red', 'red', 'red', 'blue', 'blue', 'blue'),lty = c(3, 2, 1, 3, 2, 1), cex = 0.8)
  }
  
}

simstud2 <- function(dist){

  load('eye.rda')
  load('bodyFat.rda')
  
  data_raw <- vector('list', 2)
  
  y_1 <- eye$Stratus
  x_1 <- eye$Cirrus
  data_raw[[1]] <- list(y = y_1, x = x_1, name = 'eye data')
  
  y_2 <- bodyFat$PBF
  x_2 <- as.matrix(subset(bodyFat, select = c(Age)))
  data_raw[[2]] <- list(y = y_2, x = x_2, name = 'body fat data p = 1')
  
  y_3 <- bodyFat$PBF
  x_3 <- as.matrix(subset(bodyFat, select = c(Age, WGT, HGT)))
  data_raw[[3]] <- list(y = y_3, x = x_3, name = 'body fat data p = 3')
  
  y_4 <- bodyFat$PBF
  x_4 <- as.matrix(subset(bodyFat, select = c(Age, WGT, HGT, NCK, CST)))
  data_raw[[4]] <- list(y = y_4, x = x_4, name = 'body fat data p = 5')
  #start data collection based on mv, cm, rho, n
  
  par(mfrow = c(3, 3))
  rho <- c(0.6, 0.75, 0.9)
  n <- c(100, 300, 1000)
  predictor <- c('lslp', 'malp')
  
  
  #for (h in 1:length(data_raw)){
  for (h in 4:4){
    y <- data_raw[[h]]$y
    x <- data_raw[[h]]$x
    data_name <- data_raw[[h]]$name
    mv <- get_mean_vector(y, x)
    cm <- get_cov_matrix(y, x)
    print(data_name)
    sim_stud_data <- data.frame()
    
    for (i in 1:length(rho)){
      
      for(j in 1:length(n)){
        results <- list()

        for (k in 1:length(predictor)){
          n_reps <- 10000
          skew <- 20
          error <- vector('list', n_reps)
          mse <- vector('list', n_reps)
          ten_fold <- vector('list', n_reps)
          for (l in 1:n_reps){
            if (dist == 'normal'){
              sample <- generate_sample(n[j], mv, cm, rho[i])
            }else if (dist == 'skewed'){
               sample <- generate_sample_skewed(n[j], mv, cm, rho[i], skew)
            }else{
              stop('dist must be "normal" or "skewed"')
            }
            error[[l]] <- t(as.vector(generate_error(sample$y, sample$x, predictor[k])))
            mse[[l]] <- evaluator(sample[[1]], sample[[2]], predictor[k])$mse
            ten_fold[[l]] <- evaluator_ten_fold(sample[[1]], sample[[2]], predictor[k])
          }
          error <- do.call(rbind, error)
          avg_error <- colMeans(error)
          n_error <- length(avg_error)
          iqr_error <- avg_error[(ceiling(n_error/4) + 1):floor(3*n_error/4)]
          oqr_error <- c(avg_error[1:ceiling(n_error/4)], avg_error[(floor(3*n_error/4) + 1):n_error])

          avg_mse <- mean(sapply(mse, function(x) x))
          avg_iqr_error <- mean(sapply(iqr_error, function(x) x))
          avg_oqr_error <- mean(sapply(oqr_error, function(x) x))

          ten_fold_mse <- mean(sapply(ten_fold, function(x) x$mse))
          ten_fold_iqr <- mean(sapply(ten_fold, function(x) x$mse_iq))
          ten_fold_oqr <- mean(sapply(ten_fold, function(x) x$mse_oq))


          results[[predictor[k]]] <- avg_error

          sim_stud_data <- rbind(sim_stud_data, data.frame(
            predictor = predictor[k], 
            rho = rho[i], 
            n = n[j],
            mse = avg_mse,
            iqr_error = avg_iqr_error,
            oqr_error = avg_oqr_error,     
            ten_mse = ten_fold_mse, 
            ten_iqr = ten_fold_iqr, 
            ten_oqr = ten_fold_oqr
            ))
        }

      vec <- 1:length(results[['lslp']])
      plot(x = NULL, type = 'o', main = paste(data_name, ', rho =', rho[i], ', n =', n[j], ', ', dist), xlim = c(1, length(vec)), ylim = range(c(results[['lslp']], results[['malp']])), xlab = 'index of predictand, in order', ylab = 'error')
      points(vec, results[['lslp']], type = 'b', lty = 1, col = 'red')
      points(vec, results[['malp']], type = 'b', lty = 1, col = 'blue')
      }
    }
    
    options(width = 200)
    print(sim_stud_data)
    
  }
  
}

real_data_analysis <- function(){

  real <- read.csv('real_data.csv')
  
  par(mfrow = c(3, 2))

  qqnorm(real$mGFR, main = 'Normal Q-Q Plot for mGFR')
  qqline(real$mGFR, col = "blue")

  qqnorm(real$creatinine, main = 'Normal Q-Q Plot for creatinine')
  qqline(real$creatinine, col = "blue")

  qqnorm(real$cystatinC, main = 'Normal Q-Q Plot for cystatinC')
  qqline(real$cystatinC, col = "blue")

  qqnorm(real$age, main = 'Normal Q-Q Plot for age')
  qqline(real$age, col = "blue")


  data_raw <- vector('list', 2)
  
  set.seed(1)

  y_1 <- real$mGFR
  x_1 <- real$cystatinC
  data_raw[[1]] <- list(y = y_1, x = x_1, name = 'CystatinC')

  y_2 <- real$mGFR
  x_2 <- real$creatinine
  data_raw[[2]] <- list(y = y_2, x = x_2, name = 'creatinine')

  y_3 <- real$mGFR
  x_3 <- as.matrix(subset(real, select = c('creatinine', 'cystatinC')))
  data_raw[[3]] <- list(y = y_3, x = x_3, name = 'Creatinine, CystatinC')
  
  y_4 <- real$mGFR
  x_4 <- as.matrix(subset(real, select = c('creatinine', 'cystatinC', 'age')))
  data_raw[[4]] <- list(y = y_4, x = x_4, name = 'Creatinine, CystatinC, Age')

  #start data collection based on mv, cm, rho, n
  
  predictor <- c('lslp', 'malp')
  
  analysed_data <- data.frame()

  for (i in 1:length(data_raw)){
    y <- data_raw[[i]]$y
    x <- data_raw[[i]]$x
    data_name <- data_raw[[i]]$name
    
    for (j in 1:length(predictor)){
      ev_data <- evaluator(y, x, predictor[[j]])
      error <- generate_error(y, x, predictor[[j]])
      ten_fold_reps <- replicate(50, evaluator_ten_fold(y, x, predictor[[j]]), simplify = FALSE)
      ten_mse <- mean(sapply(ten_fold_reps, function(r) r$mse))
      ten_iqr <- mean(sapply(ten_fold_reps, function(r) r$mse_iq))
      ten_oqr <- mean(sapply(ten_fold_reps, function(r) r$mse_oq))

      n_error <- length(error)
      iqr_errors <- error[(ceiling(n_error/4) + 1):floor(3*n_error/4)]
      oqr_errors <- c(error[1:ceiling(n_error/4)], error[(floor(3*n_error/4) + 1):n_error])

      mse_iq <- mean(sapply(iqr_errors, function(x) x))
      mse_oq <- mean(sapply(oqr_errors, function(x) x))
          
      analysed_data <- rbind(analysed_data, data.frame(
        predictor = predictor[j], 
        pcc = ev_data$pcc, 
        ccc = ev_data$ccc, 
        mse = ev_data$mse,
        mse_iq = mse_iq,
        mse_oq = mse_oq,
        ten_mse = ten_mse, 
        ten_iqr = ten_iqr, 
        ten_oqr = ten_oqr, 
        prediction_data = data_name 
        ))
        }
      }

    options(width = 200)
    print(analysed_data)

}

#Kim_etal_eye_replication()
#simstud1('skewed') #either 'normal' or 'skewed'
#simstud2('skewed') #either 'normal' or 'skewed'
real_data_analysis()