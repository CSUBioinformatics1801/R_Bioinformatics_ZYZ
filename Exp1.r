# ------------encode utf-8-------------------
#datetime: 20210322
#author: Sean Peldom Zhang
#R.version()=4.0.3
# ------------Homework------------------------
str_test <- "hello world"
print(str_test)

A <- c(1:3)
B <- seq(from = 1, to = 5, by = 1)
C <- rep(0, 5)
print(A + B)

a <- factor(c(1, 1, 2, 2, 3, 3, 3, 1, 2))
b <- a
levels(b) <- c("Male", "Female", "unknown")
print(a)
print(b)

levels(b) <- c("a", "b", "c", "d")
print(b)

mat1 <- matrix(1:12, 3, 4)
# mat2 <- t(mat1)
mat2= matrix(1:12, 4, 3)
print(2*mat1)
print(mat1*mat1)
print(mat1 %*% mat2)

b <- factor(c(1, 1, 2, 2, 3, 3, 3, 1, 2));
levels(b) <- c("Male", "Female", "unknown");
dat <- data.frame(Sex = b,
                  GLU = c(8.8, 6.5, 5.4, 5.6, 6.7, 7.8, 4.5, 9.7, 5.0), 
                  GHb = c(5.5, 3.4, 5.5, 4.3, 3.5, 3.4, 5.5, 7.0, 4.3));
print(dat)

x <- 7
if (x > 6.1) {
    cat("GLU<U+8FC7><U+9AD8>")
} else if (x < 3.9) {
    cat("GLU<U+8FC7><U+4F4E>")
}

x <- 0
for (i in 1:9) {
    x <- x + dat$GLU[i]
}
print(x)

i=1
y=0
while(i<=nrow(dat)){
    y=y+dat$GHb[i]
    i=i+1
}
print(y)

Male_num <- 0
for (i in c(1:nrow(dat))) {
    if (dat[i, ]["Sex"] == "Male") {
        Male_num <- Male_num + 1
    }
}
print(Male_num)

Female_num <- 0
for (i in c(1:nrow(dat))) {
    rotator <- FALSE
    while (dat[i, ]["Sex"] == "Female" && dat[i, ]["GLU"] < 6.1
    && dat[i, ]["GLU"] > 3.9 && rotator == FALSE) {
        Female_num <- Female_num + 1
        rotator <- TRUE
    }
}
print(Female_num)
