x <- rnorm(100);
y <- rexp(100)
result <- data.frame(x = x, y = y);
write.csv(result, "test.csv")
