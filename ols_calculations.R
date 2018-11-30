library(RConics)
library(tidyverse)
library(broom)

# http://faculty.cas.usf.edu/mbrannick/regression/regma.htm
# http://www.statmethods.net/advstats/matrix.html
# https://onlinecourses.science.psu.edu/stat501/node/382

# create data
y <- c(1, 2, 3, 3, 4)
x0 <- c(1, 1, 1, 1, 1)
x1 <- c(6, 7, 8, 9, 7)
x2 <- c(5, 7, 10, 8, 9)
df <- data.frame(y, x1, x2)

# alternate data from duke spreadsheet tutorial on multiple regression
y <- c(45, 58, 50, 54, 62, 53)
x0 <- c(1, 1, 1, 1, 1, 1)
x1 <- c(18, 25, 15, 22, 24, 20)
df <- data.frame(y, x1)

# create y matrix
y_mat <- as.matrix(data.frame(y = y))

# create x matrix
x_mat <- as.matrix(data.frame(x = x0, x1, x2))
x_mat <- as.matrix(data.frame(x = x0, x1))

# population variance of y
pop_var_y <- sum((y - mean(y))^2) / length(y)
pop_var_y

# sampel variance of y
samp_var_y <- sum((y - mean(y))^2) / (length(y) - 1)
samp_var_y

# degrees of freedom = # observations - # coefficients
degrees_of_freedom <- length(y) - ncol(x_mat)
degrees_of_freedom 

# schaums outline of econmetrics has derivation of equation for B
# goal is to find values of B that minimize SSE, in matrix form to accomodate multiple B: (Y - BX)(Y - BX)
# take derivative w respect to B and set equal to 0
# d/dB = Y^2 - YBX - YBX + (BX)^2; = Y^2 - 2YBX + B^2X^2
# d/dB = -2YX + 2BX^2
# 0 = -2YX + 2BX^2
# 2YX = 2BX; YX = BX^2; YX / X^2 = B; or in matrix notation YX * (X'X)^-1 = B 


# find X'Y
xy <- t(x_mat) %*% y_mat

# find X'X
# x_mat %*% t(x_mat)
xx <- t(x_mat) %*% x_mat

# find determinant of x'x
xx_det <- det(xx)


# find cofactor of x'x
cofactor_mat <- matrix(NA, 3, 3)
for(row in 1:nrow(xx)){
        for(col in 1:nrow(xx)){
                cofactor <- cofactor(xx, row, col)
                print(cofactor)
                cofactor_mat[row, col] <- cofactor
        }
}

# find adjoint of cofactor matrix of x'x
# note, RConics adjoint function provides transpose of cofactors for a square matrix, 
# so basically skips step of having to find cofactors with for loop above
# adjoint(xx)
# note that transpose of this matrix just returns the same matrix
adjoint <- t(cofactor_mat)

# inverse of A equals adjoint of A times 1/determinant of A: A^-1 = AdjA * 1/|A| 
xx_inv <- adjoint * 1 / xx_det

# find inverse of xx
xx_inv <- solve(xx)

# check if xx_inv * xx = I
xx_inv %*% xx

# calculate b = xx_inv * xy
b <- xx_inv %*% xy

# manually calculate b
b0 <- sum(xx_inv[1, ] * xy)
b1 <- sum(xx_inv[2, ] * xy)
b2 <- sum(xx_inv[3, ] * xy)

# check b using lm
m1 <- lm(y ~ x1 + x2, data = df)
m1 <- lm(y ~ x1, data = df)
m1 %>% tidy()

# get y-hat predictions with X * B
x_mat
b
y_hat <- x_mat %*% b
y_hat

# check predicitons with predict()
predict(object = m1, newdata = df)

# get errors (actual - predicted)
errors <- y - y_hat
errors

# get squared errors
sq_errors <- errors^2
sq_errors

# get simple average sq_error
simp_avg_sq_error <- mean(sq_errors)
simp_avg_sq_error

# get r-squared, which = 1 - (avg_sq_error / pop_var_y)
r_squared <- 1 - (simp_avg_sq_error / pop_var_y)
r_squared

# get mean squared error = sum(sq_errors) / degrees_of_freedom
mean_sq_error <- sum(sq_errors) / degrees_of_freedom
mean_sq_error

# adjusted r-squared = 1 - (mean_sq_error / samp_var_y)
adj_r_squared <- 1 - (mean_sq_error / samp_var_y)
adj_r_squared

# "standard error of regression" = RMSE (root mean_sq_error)
root_mean_sq_error <- sqrt(mean_sq_error)
root_mean_sq_error

# "covariance matrix of the coefficent estimates" = X'X^-1 * mean_sq_error
covariance_matrix <- xx_inv * mean_sq_error
covariance_matrix

# standard errors of coefficient estiamtes are sq.roots of diagonal elements of covariance matrix
standard_errors <- c()
for(i in 1:nrow(covariance_matrix)) {
        standard_errors <- c(standard_errors, sqrt(covariance_matrix[i, i]))
}
standard_errors

# t_stats for coefficient estimates = coefficients / standard_errors
t_stats <- b / standard_errors
t_stats

# get critical_value 
critical_value <- qt(p = .975, df = degrees_of_freedom)
critical_value

# get p_values
p_value <- 2 * pt(q = -abs(t_stats), df = degrees_of_freedom)
p_value

# get conf_int for coefficients (useful if using ols for t_test on difference of means w/ single dummy variable)
conf_int_upper <- (standard_errors * critical_value) + b
conf_int_upper
conf_int_lower <- (standard_errors * critical_value) - b
conf_int_lower

# get new data to predict 
z0 <- c(1) 
z1 <- c(30)
Z <- as.matrix(data.frame(x = z0, z1))
Z

# prediction is new data * coefficients
prediction <- Z %*% b
prediction

# "variance of (prediction) the mean" (? not sure what this refers to) = (Z*covariance_matrix) * transpose(Z)
variance_of_prediction_mean <- (Z %*% covariance_matrix) %*% t(Z)
variance_of_prediction_mean

# "variance of the prediction" = variance_of_prediction_mean + mean_sq_error
variance_of_prediction <- variance_of_prediction_mean + mean_sq_error
variance_of_prediction

# "standard error of the (prediction) mean" = sqrt(variance_of_prediction_mean)
standard_error_of_prediction_mean <- sqrt(variance_of_prediction_mean)
standard_error_of_prediction_mean

# "standard error of the prediction" = sqrt(variance_of_prediction)
standard_error_of_prediction <- sqrt(variance_of_prediction)
standard_error_of_prediction

# conf_int for prediction
conf_int_upper <- prediction + (critical_value * standard_error_of_prediction)
conf_int_upper

conf_int_lower <- prediction - (critical_value * standard_error_of_prediction)
conf_int_lower

# check predicitons with predict()
# ?predict.lm
predict(object = m1, newdata = tibble(x1 = 30), interval = "predict")





