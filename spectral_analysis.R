######################################################################################################################
# First harmonic of historical data
#
# We assume wave-like behavior in S&P 500 price volatility, as investors "buy low" and "sell high." The lowest
# frequency sustainable DAMPED wave decay wate should have an inflection point at the at rate of exponential growth
# of the market volume, "flat-lining." However, DRIVEN waves can sustain above this inflection point. The below script
# isolates a driven wave, with amplification rate on order of the first harmonic of the entire history of the market.
# Its driving APR is approximately 7.6% in the training set, by which we calibrate the approximate historical growth
# rate of the market and the highest end of a DISPERSION RELATION for market volatility, as waves in the price of the
# S&P 500. (In general, the rate of exponential growth of the market is not exactly constant in time, though close,
# hence we see some compression on the future wave as exponential growth rate ostensibly increases.)

setwd("/home/iamu/Github/s-p500-volatility")

offset <- 0
trainingSize <- 3206
validationSize <- 1603
m <- 280 #Avoid changing - smooths apparent volatility and couples to decay rate 
pCutoff <- 1

wavePeriod <- 2500
waveOffset <- 2*pi*0.25

allData <- read.csv("SP500_Weekly_Preprocessed.csv", header=TRUE)
training <- allData[(offset + 1):(offset + trainingSize),]
validation <- allData[(offset + trainingSize + 1):(offset + trainingSize + validationSize),]
fullSet <- rbind(training, validation)

SP500_Growth_Adjusted_Weekly_Close <- training$D_Close

raw.spec <- spec.pgram(SP500_Growth_Adjusted_Weekly_Close, taper = 0)
raw.spec <- spec.pgram(SP500_Growth_Adjusted_Weekly_Close, kernel("daniell", c(m,m)), taper = 0)
specDataFrame<-data.frame("freq"=raw.spec$freq,"spec"=raw.spec$spec)

require(data.table)
topSpecs<-specDataFrame[order(-specDataFrame$spec),]
require(pracma)
specPeaks<-findpeaks(raw.spec$spec)

smoothed_close <- kernapply(training$D_Close, kernel("daniell", c(m,m)))
interp_smoothed_close <- approx(x=approx(seq(from=(offset+1), to=(offset+trainingSize)),n=length(smoothed_close))$y, y=smoothed_close, n=trainingSize)

X <- data.frame(Week=training$Week,
                Close = interp_smoothed_close$y
)

library(ggplot2)
ggplot(X, aes(x=Week, y=Close)) +  geom_point()

#for (peak in specPeaks[,2]) {
#  X <- cbind(X, sin(2*pi*X$Week/peak))
#  colnames(X)[length(colnames(X))] <- paste("sin", toString(peak), sep="", collapse="")
#  X <- cbind(X, cos(2*pi*X$Week/peak))
#  colnames(X)[length(colnames(X))] <- paste("cos", toString(peak), sep="", collapse="")
#}

X <- cbind(X, sin(2*pi*X$Week/wavePeriod + waveOffset))
colnames(X)[length(colnames(X))] <- paste("sin", toString(wavePeriod), sep="", collapse="")

mod <- lm(Close ~ . - Week, data = X)  # Regress Close on everything (but Week)
summary(mod)
coeffs <- summary(mod)$coefficients

omitPredictors <- c()
for (c in 1:nrow(coeffs)) {
  pvalue = coeffs[c,ncol(coeffs)]
  if (is.na(pvalue) || pvalue > pCutoff) {
    omitPredictors <- append(omitPredictors, rownames(coeffs)[c])
  }
}

newOmitCount <- length(omitPredictors)
while (newOmitCount > 0) {

  X <- data.frame(Week=training$Week,
                  Close = interp_smoothed_close$y
  )
  
  #for (peak in specPeaks[,2]) {
  #  predictorName <- paste("sin", toString(peak), sep="", collapse="")
  #  
  #  if (!(predictorName %in% omitPredictors)) {
  #    X <- cbind(X, sin(2*pi*X$Week/peak))
  #    colnames(X)[length(colnames(X))] <- predictorName
  #  }
  #  
  #  predictorName <- paste("cos", toString(peak), sep="", collapse="")
  #  
  #  if (!(predictorName %in% omitPredictors)) {
  #    X <- cbind(X, cos(2*pi*X$Week/peak))
  #    colnames(X)[length(colnames(X))] <- predictorName
  #  }
  #}
  
  X <- cbind(X, sin(2*pi*X$Week/wavePeriod + waveOffset))
  colnames(X)[length(colnames(X))] <- paste("sin", toString(wavePeriod), sep="", collapse="")
  
  mod <- lm(Close ~ . - Week, data = X)  # Regress Close on everything (but Week)
  coeffs <- summary(mod)$coefficients
  
  newOmitCount <- 0
  for (c in 1:nrow(coeffs)) {
    pvalue = coeffs[c,ncol(coeffs)]
    if (is.na(pvalue) || pvalue > pCutoff) {
      omitPredictors <- append(omitPredictors, rownames(coeffs)[c])
      newOmitCount <- newOmitCount + 1
    }
  }

}
  
summary(mod)

library(ggplot2)
ggplot(X, aes(x=Week, y=Close)) +  geom_point()

X$resid <- residuals(mod)
X$pred <- predict(mod)
ggplot(data = X) + labs(title=paste("Weekly S&P 500 close (Training, R^2=", toString(summary(mod)$r.squared), ")", sep="", collapse=""), subtitle="Normalized by variance and APR volatility amplification, 280 week moving average smoothed") + geom_line(aes(x = Week, y = Close, color="Observed")) + geom_line(aes(x = Week, y = pred, color="Predicted"))

smoothed_close <- kernapply(fullSet$D_Close, kernel("daniell", c(m,m)))
interp_smoothed_close <- approx(x=approx(seq(from=offset+1, to=(offset+trainingSize+validationSize)),n=length(smoothed_close))$y, y=smoothed_close, n=(trainingSize+validationSize))
XAll <- data.frame(Week=fullSet$Week, Close=interp_smoothed_close$y)

#for (peak in specPeaks[,2]) {
#  XAll <- cbind(XAll, sin(2*pi*XAll$Week/peak))
#  colnames(XAll)[length(colnames(XAll))] <- paste("sin", toString(peak), sep="", collapse="")
#  XAll <- cbind(XAll, cos(2*pi*XAll$Week/peak))
#  colnames(XAll)[length(colnames(XAll))] <- paste("cos", toString(peak), sep="", collapse="")
#}

XAll <- cbind(XAll, sin(2*pi*XAll$Week/wavePeriod + waveOffset))
colnames(XAll)[length(colnames(XAll))] <- paste("sin", toString(wavePeriod), sep="", collapse="")

XFuture <- XAll[(trainingSize + 1):(trainingSize + validationSize),]
XFuture$pred <- predict(mod, newdata=XFuture)
SS.total <- sum((XFuture$Close - mean(XFuture$Close))^2)
SS.residual <- sum((XFuture$pred- XFuture$Close)^2)
ggplot(data = XFuture) + labs(title=paste("Weekly S&P 500 close, (Validation, R^2=", toString(1-SS.residual/SS.total), ")", sep="", collapse=""), subtitle="Normalized by variance and ~30% APR volatility decay, 6 week moving average smoothed") + geom_line(aes(x = Week, y = Close, color="Observed")) + geom_line(aes(x = Week, y = pred, color="Predicted"))

XAll$pred <- predict(mod, newdata=XAll)
ggplot(data = XAll) + labs(title="Weekly S&P 500 close, (Full Set)", subtitle="Normalized by variance and APR volatility amplification, 280 week moving average smoothed") + geom_line(aes(x = Week, y = Close, color="Observed")) + geom_line(aes(x = Week, y = pred, color="Predicted"))

