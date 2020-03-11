setwd("/home/iamu/Github/s-p500-volatility")

p0 <- 2500
c0 <- 1/2

weeklyPeriodDamp <- function (wavePeriod) { ( 1 + ( atan( c0 * log( wavePeriod / p0 ) ) / pi ) ) ^ ( 1 / 52.1429 ) }

offset <- 1300
trainingSize <- 413
validationSize <- 138
m <- 6 #Avoid changing - smooths apparent volatility and couples to decay rate 
pCutoff <- 0.1

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
interp_smoothed_close <- approx(x=approx(range((offset+1):(offset+trainingSize)),n=length(smoothed_close))$y, y=smoothed_close, n=trainingSize)

X <- data.frame(Week=training$Week,
                Close = interp_smoothed_close$y
)

for (peak in specPeaks[,2]) {
  damping <- weeklyPeriodDamp(peak)
  X <- cbind( X, damping ^ X$Week * sin( 2 * pi * X$Week/peak ) )
  colnames(X)[length(colnames(X))] <- paste("sin", toString(peak), sep="", collapse="")
  X <- cbind( X, damping ^ X$Week * cos( 2 * pi * X$Week/peak ) )
  colnames(X)[length(colnames(X))] <- paste("cos", toString(peak), sep="", collapse="")
}


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
  
  for (peak in specPeaks[,2]) {
    damping <- weeklyPeriodDamp(peak)

    predictorName <- paste("sin", toString(peak), sep="", collapse="")
    
    if (!(predictorName %in% omitPredictors)) {
      X <- cbind( X, damping ^ X$Week * sin( 2 * pi * X$Week/peak ) )
      colnames(X)[length(colnames(X))] <- predictorName
    }
    
    predictorName <- paste("cos", toString(peak), sep="", collapse="")
    
    if (!(predictorName %in% omitPredictors)) {
      X <- cbind( X, damping ^ X$Week * cos( 2 * pi * X$Week/peak ) )
      colnames(X)[length(colnames(X))] <- predictorName
    }
  }
  
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
ggplot(data = X) + labs(title=paste("Weekly S&P 500 close (Training, R^2=", toString(summary(mod)$r.squared), ")", sep="", collapse=""), subtitle="Normalized by variance and ~30% APR volatility decay, 6 week moving average smoothed") + geom_line(aes(x = Week, y = Close, color="Observed")) + geom_line(aes(x = Week, y = pred, color="Predicted"))

smoothed_close <- kernapply(fullSet$D_Close, kernel("daniell", c(m,m)))
interp_smoothed_close <- approx(x=approx(range((offset+1):(offset+trainingSize+validationSize)),n=length(smoothed_close))$y, y=smoothed_close, n=(trainingSize+validationSize))
XAll <- data.frame(Week=fullSet$Week, Close=interp_smoothed_close$y)

for (peak in specPeaks[,2]) {
  damping <- weeklyPeriodDamp(peak)
  XAll <- cbind( XAll, damping ^ XAll$Week * sin( 2 * pi * XAll$Week/peak ) )
  colnames(XAll)[length(colnames(XAll))] <- paste("sin", toString(peak), sep="", collapse="")
  XAll <- cbind( XAll, damping ^ XAll$Week * cos( 2 * pi * XAll$Week/peak ) )
  colnames(XAll)[length(colnames(XAll))] <- paste("cos", toString(peak), sep="", collapse="")
}

XFuture <- XAll[(trainingSize + 1):(trainingSize + validationSize),]
XFuture$pred <- predict(mod, newdata=XFuture)
SS.total <- sum((XFuture$Close - mean(XFuture$Close))^2)
SS.regression <- sum((XFuture$pred- mean(XFuture$Close))^2)
ggplot(data = XFuture) + labs(title=paste("Weekly S&P 500 close, (Validation, R^2=", toString(1-summary(mod)$r.squared), ")", sep="", collapse=""), subtitle="Normalized by variance and ~30% APR volatility decay, 6 week moving average smoothed") + geom_line(aes(x = Week, y = Close, color="Observed")) + geom_line(aes(x = Week, y = pred, color="Predicted"))

XAll$pred <- predict(mod, newdata=XAll)
ggplot(data = XAll) + labs(title="Weekly S&P 500 close, (Full Set)", subtitle="Normalized by variance and ~30% volatility decay, 6 week moving average smoothed") + geom_line(aes(x = Week, y = Close, color="Observed")) + geom_line(aes(x = Week, y = pred, color="Predicted"))

