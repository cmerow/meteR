context('test of power distribution predictions')

S0 <- 31
N0 <- 877
E0 <- 917872

this.esf <- meteESF(S0=S0, N0=N0, E=E0/1.134)


test_that('predicted rank ipd is correct', {
  yRaw <- c(0.072, 0.498, 0.929, 1.356, 1.787, 2.208, 2.639)
  yAct <- 6:0
  newmanPsi <- predict(lm(yAct~yRaw), 
                       newdata=data.frame(yRaw=c(0.555, 0.637, 0.801, 1.171, 1.186, 
                                                 1.202, 1.217, 1.73, 1.766, 1.802, 
                                                 1.843, 2.567, 2.603, 2.634)))
  newmanRank <- c(1, 2, 7, 65, 65, 65, 65, 605, 605, 605, 605, 877, 877, 877)
  
  this.ipd <- ipd(this.esf)
  this.ipdRank <- meteDist2Rank(this.ipd)
  
  expect_true(all(abs(log(this.ipdRank[newmanRank[1:3]], 10) - newmanPsi[1:3]) < 0.06))
  
  minmax <- do.call(rbind, tapply(newmanPsi, newmanRank, range))
  
  expect_true(all(log(this.ipdRank[unique(newmanRank[-(1:3)])], 10) <= minmax[-(1:3), 2] & 
                    log(this.ipdRank[unique(newmanRank[-(1:3)])], 10) >= minmax[-(1:3), 1]))
})


test_that('predicted rank sipd is correct', {
  
})