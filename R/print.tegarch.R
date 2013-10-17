print.tegarch <-
function(x, ...)
{
vcovmat <- vcov(x)
out1 <- rbind(x$par, sqrt(diag(vcovmat)))
rownames(out1) <- c("", "s.e.")
out2 <- rbind(x$objective, x$sic)
rownames(out2) <- c("Log-likelihood:", "SIC:")
colnames(out2) <- ""

cat("Date:", x$date, "\n")
cat("Message:", x$message, "\n")
if(!is.null(x$NOTE)){
  cat("NOTE:", x$NOTE, "\n")
}
cat("\n")
cat("Coefficients:\n")
print(out1)
print(out2)
}
