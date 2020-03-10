######################################################################################################################
# First and second harmonics of historical data
#
# We assume wave-like behavior in S&P 500 price volatility, as investors "buy low" and "sell high." The lowest
# frequency sustainable DAMPED wave decay wate should have an inflection point at the at rate of exponential growth
# of the market volume, "flat-lining." However, DRIVEN waves can sustain above this inflection point. The below script
# isolates a driven wave, with amplification rate on order of the first harmonic of the entire history of the market.
# Its driving APR is approximately 7.6% in the training set, by which we calibrate the approximate historical growth
# rate of the market and the highest end of a DISPERSION RELATION for market volatility, as waves in the price of the
# S&P 500. (In general, the rate of exponential growth of the market is not exactly constant in time, though close,
# hence we see some compression on the future wave as exponential growth rate ostensibly increases.)
#
# This script attempts to first the first and second harmonics of driven waves, in overall historical S&P 500 data.
# Fit-tuning is on-going. The model must be validated.

setwd("/home/iamu/Github/s-p500-volatility")

offset <- 0
trainingSize <- 4809
validationSize <- 0
m <- 180 #Avoid changing - smooths apparent volatility and couples to decay rate 
pCutoff <- 1

wavePeriod <- 2500
waveOffset <- 2*pi*0.25
waveDrive <- 1.00140482691777 #Weekly, compounding to 7.6% APR
wavePeriod2 <- wavePeriod / 2
waveOffset2 <- waveOffset + pi / 4
waveDrive2 <- 1.00140482691777 #Weekly, compounding to 15% APR

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

X <- cbind(X, waveDrive^X$Week * sin(2*pi*X$Week/wavePeriod + waveOffset))
colnames(X)[length(colnames(X))] <- paste("sin", toString(wavePeriod), sep="", collapse="")
X <- cbind(X, waveDrive2^X$Week * sin(2*pi*X$Week/wavePeriod2 + waveOffset2))
colnames(X)[length(colnames(X))] <- paste("sin", toString(wavePeriod2), sep="", collapse="")

mod <- lm(Close ~ . - Week, data = X)  # Regress Close on everything (but Week)
summary(mod)
coeffs <- summary(mod)$coefficients

#omitPredictors <- c()
#for (c in 1:nrow(coeffs)) {
#  pvalue = coeffs[c,ncol(coeffs)]
#  if (is.na(pvalue) || pvalue > pCutoff) {
#    omitPredictors <- append(omitPredictors, rownames(coeffs)[c])
#  }
#}

#newOmitCount <- length(omitPredictors)
#while (newOmitCount > 0) {

  #X <- data.frame(Week=training$Week,
  #                Close = interp_smoothed_close$y
  #)
  
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
  
  #X <- cbind(X, sin(2*pi*X$Week/wavePeriod + waveOffset))
  #colnames(X)[length(colnames(X))] <- paste("sin", toString(wavePeriod), sep="", collapse="")
  #X <- cbind(X, sin(2*pi*X$Week/wavePeriod2 + waveOffset2))
  #colnames(X)[length(colnames(X))] <- paste("sin", toString(wavePeriod2), sep="", collapse="")
  
  #mod <- lm(Close ~ . - Week, data = X)  # Regress Close on everything (but Week)
  #coeffs <- summary(mod)$coefficients
  
  #newOmitCount <- 0
  #for (c in 1:nrow(coeffs)) {
  #  pvalue = coeffs[c,ncol(coeffs)]
  #  if (is.na(pvalue) || pvalue > pCutoff) {
  #    omitPredictors <- append(omitPredictors, rownames(coeffs)[c])
  #    newOmitCount <- newOmitCount + 1
  #  }
  #}

#}
  
summary(mod)

library(ggplot2)
ggplot(X, aes(x=Week, y=Close)) +  geom_point()

X$resid <- residuals(mod)
X$pred <- predict(mod)
ggplot(data = X) + labs(title=paste("Weekly S&P 500 close (Training, R^2=", toString(summary(mod)$r.squared), ")", sep="", collapse=""), subtitle="Normalized by variance and APR volatility amplification, 280 week moving average smoothed") + geom_line(aes(x = Week, y = Close, color="Observed")) + geom_line(aes(x = Week, y = pred, color="Predicted"))