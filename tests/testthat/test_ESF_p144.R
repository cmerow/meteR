context('ESF_p144')
# tests the analytically computed results from Harte 2011 on page 144

d=data.frame(sp=c('a','b','c','c','d','d','d'),
             abun=c(1,1,1,1,1,3,8),
             pow=c(23,9,6,9,3,2,1))
esf1 <- meteESF(spp=d$sp, abund=d$abun, power=d$pow)
sad1=sad(esf1)
plot(sad1)
#sad1$d(c(1,2,12))
#meteDist2Rank(sad1)
#predictESF(esf1,c(12,12,12,2,2,1,1),c(1,2,3,6,9,9,23))


test_that('ESF computed correctly',{
  expect_is(esf1,"meteESF")

})


test_that('SAD computed correctly',{
  expect_is(sad1,"sad")
})


test_that('IPD computed correctly',{
  expect_equal(1,1)
  
  
})

