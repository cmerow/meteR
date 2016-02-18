context('SAR works')

data(bellview)
head(bellview)
bell.sar <- meteSAR(bellview$species, abund=rep(1, nrow(bellview)), x=bellview$x, y=bellview$y, Amin=1, A0=16)
