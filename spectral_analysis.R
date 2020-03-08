setwd("/home/iamu/Github/s-p500-volatility")

allData <- read.csv("SP500_Weekly_Preprocessed.csv", header=TRUE)
offset <- 600
trainingSize <- 1241
validationSize <- 12
training <- allData[(offset + 1):(offset + trainingSize),]
validation <- allData[(offset + trainingSize + 1):(offset + trainingSize + validationSize),]
fullSet <- rbind(training, validation)
m = 12

SP500_Growth_Adjusted_Weekly_Close <- training$D_Close

raw.spec <- spec.pgram(SP500_Growth_Adjusted_Weekly_Close, taper = 0)
raw.spec <- spec.pgram(SP500_Growth_Adjusted_Weekly_Close, kernel("daniell", c(m,m)), taper = 0)
specDataFrame<-data.frame("freq"=raw.spec$freq,"spec"=raw.spec$spec)

require(data.table)
topSpecs<-specDataFrame[order(-specDataFrame$spec),]
require(pracma)
specPeaks<-findpeaks(raw.spec$spec)

smoothed_close <- kernapply(training$D_Close, kernel("daniell", c(m,m)))
interp_smoothed_close <- approx(x=approx(range(1:trainingSize),n=length(smoothed_close))$y, y=smoothed_close, n=trainingSize)

X <- data.frame(Week=training$Week,
                Close = interp_smoothed_close$y
)

for (peak in specPeaks[,2]) {
  X <- cbind(X, sin(2*pi*X$Week/peak))
  colnames(X)[length(colnames(X))] <- paste("sin", toString(peak), sep="", collapse="")
  X <- cbind(X, cos(2*pi*X$Week/peak))
  colnames(X)[length(colnames(X))] <- paste("cos", toString(peak), sep="", collapse="")
}


mod <- lm(Close ~ . - Week, data = X)  # Regress Close on everything (but Week)
summary(mod)

library(ggplot2)
ggplot(X, aes(x=Week, y=Close)) + geom_point()

X$resid <- residuals(mod)
X$pred <- predict(mod)
ggplot(data = X) + labs(title="Weekly S&P 500 close (Training)", subtitle="Normalized by variance and 6% APR growth, 24 week moving average smoothed") + geom_line(aes(x = Week, y = Close, color="Observed")) + geom_line(aes(x = Week, y = pred, color="Predicted"))

smoothed_close <- kernapply(training$D_Close, kernel("daniell", c(m,m)))
interp_smoothed_close <- approx(x=approx(range(1:trainingSize),n=length(smoothed_close))$y, y=smoothed_close, n=trainingSize)
XFuture <- data.frame(Week=training$Week, Close=interp_smoothed_close$y)

for (peak in specPeaks[,2]) {
  XFuture <- cbind(XFuture, sin(2*pi*XFuture$Week/peak))
  colnames(XFuture)[length(colnames(XFuture))] <- paste("sin", toString(peak), sep="", collapse="")
  XFuture <- cbind(XFuture, cos(2*pi*XFuture$Week/peak))
  colnames(XFuture)[length(colnames(XFuture))] <- paste("cos", toString(peak), sep="", collapse="")
}

XFuture$pred <- predict(mod, newdata=XFuture)
ggplot(data = XFuture) + labs(title="Weekly S&P 500 close, (Full Set)", subtitle="Normalized by variance and 6% APR growth, 24 week moving average smoothed") + geom_line(aes(x = Week, y = Close, color="Observed")) + geom_line(aes(x = Week, y = pred, color="Predicted"))

