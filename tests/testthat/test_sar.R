context('SAR works')
## values of S0, N0, A0 and predicted S for SAR from Newman et al. 


test_that('predicted SAR values are correct', {
  bell.sar <- meteSAR(S0=32, N0=920, Amin=1, A0=16)$pred$S
  newmanS <- c(12.2560715, 15.8991283, 20.2158684, 25.7856639, 32)
  expect_true(all(round(bell.sar, 3) == round(newmanS, 3)))
})

test_that('SAR works for absolute and scaled area in the same way', {
  expect_equal(meteSAR(S0=30, N0=1000, Amin=1, A0=16)$pred$S,
               meteSAR(S0=30, N0=1000, Amin=1/16, A0=1)$pred$S)
})

test_that('EAR works for absolute and scaled area in the same way', {
  expect_equal(meteSAR(S0=30, N0=1000, Amin=1, A0=16, EAR=TRUE)$pred$S,
               meteSAR(S0=30, N0=1000, Amin=1/16, A0=1, EAR=TRUE)$pred$S)
})


test_that('predicted EAR values are correct', {
  areas <- 2^(-1*(0:8)) * 256 # Serpentine areas
  S0 <- 24
  N0 <- 37182
  A0 <- 256
  
  rawAxis <- c(48, 192, 334, 476, 614, 760, 900, 1040, 1182, 1322, 1466)
  actualAxis <- c(4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6)
  rawMETE <- c(166, 492, 622, 746, 864, 970, 1070, 1170, 1268)
  actualMETE <- predict(lm(actualAxis ~ rawAxis), newdata=data.frame(rawAxis=rawMETE))
  
  serpEAR <- downscaleSAR(meteESF(S0=S0, N0=N0), areas, A0=A0, EAR=TRUE)
  log(serpEAR$S) - actualMETE
  
  this.esf <- meteESF(S0=S0, N0=N0, E0=10^10)
  this.esf$La
  
  
  
  b <- 5.92e-05 #sum(this.esf$La)
  n <- 1:N0
  
  eq7.74 <- sapply(areas, function(A) S0/log(1/b) * sum(exp(-b*n)/n * A0/(n*A + A0) * (n*A/(n*A + A0))^n))
  log(eq7.74)
  
  earIdentity <- downscaleSAR(this.esf, A=A0, A0=A0)$S - downscaleSAR(this.esf, A0-areas[-1], A0=A0)$S
  eq7.74[-1]
  serpEAR$S[-1]
  
  par(mfrow=c(1, 4))
  plot(actualMETE, log(eq7.74), xlim=c(-4, 2), ylim=c(-4, 2)); abline(0, 1)
  plot(actualMETE, log(eq7.75), xlim=c(-4, 2), ylim=c(-4, 2)); abline(0, 1)
  plot(actualMETE, log(eq7.76), xlim=c(-4, 2), ylim=c(-4, 2)); abline(0, 1)
  plot(actualMETE, log(eq7.77), xlim=c(-4, 2), ylim=c(-4, 2)); abline(0, 1)
  
  plot(areas, actualMETE, log='x')
  lines(areas, log(serpEAR$S))
})
