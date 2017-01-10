library(RConics)

# http://faculty.cas.usf.edu/mbrannick/regression/regma.htm
# http://www.statmethods.net/advstats/matrix.html
# https://onlinecourses.science.psu.edu/stat501/node/382

# create data
y <- c(1, 2, 3, 3, 4)
x0 <- c(1, 1, 1, 1, 1)
x1 <- c(6, 7, 8, 9, 7)
x2 <- c(5, 7, 10, 8, 9)
df <- data.frame(y, x1, x2)

# create y matrix
y_mat <- as.matrix(data.frame(y = y))

# create x matrix
x_mat <- as.matrix(data.frame(x = x0, x1, x2))

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
# note, RConics adjoint function doesn't seem to return transpose for some reason?? use t()
# adjoint(cofactor_mat)
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

df2 <- data.frame(y, x1)
m2 <- lm(y ~ x1, data = df2)
m2

