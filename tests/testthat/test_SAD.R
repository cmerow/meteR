context('test SAD')
## values from Newman et al.

valid <- data.frame(species = c('TAROFF', 'UNKSP', 'AQUCOE', 'OSMOCC', 'VIGMUL', 'CASOCC',
                                'GERRIC', 'SENAMP', 'ERISPE', 'EUCENG', 'DRASPE', 'AGOGLA', 'ERIELA',
                                'CASSUL', 'SOLMUL', 'LATLEU', 'ERIFOR', 'DELBAR', 'POLDOU', 'SENCRA',
                                'NOCMON', 'IPOAGG', 'VICAME', 'LIGPOR', 'POTGRA', 'FRAVES', 'HYMHOO',
                                'ANDSEP', 'LUPARG', 'HELQUI', 'VIONUT', 'BOEDRU'), 
                    abund = c(1, 1, 2, 2, 2, 6, 7, 7, 8, 8, 10, 11, 11, 12, 12, 13, 15, 17, 17, 19, 
                              20, 23, 33, 40, 47, 55, 55, 80, 86, 95, 101, 104), 
                    newmanPred = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5, 6, 6, 8, 9, 11, 13, 15, 
                                   18, 21, 25, 31, 37, 46, 58, 74, 97, 138, 240))

ourSAD <- sad(meteESF(valid$species, valid$abund))

test_that('SAD from empirical data works', {
  expect_true(all((sort(meteDist2Rank(ourSAD)) - valid$newmanPred)/valid$newmanPred < 0.05))
})

