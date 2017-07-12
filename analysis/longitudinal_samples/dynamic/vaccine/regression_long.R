library(data.table)
library(stargazer)
library(plm)
library(lmtest)

multiply.1e2 = function(x){
	return(x*1e2)
}

make.fit = function(data){
	FITS = list()
	i=1

	fit = plm(I ~ lag(V, 0:4) + VR01 + VR05 + VR10 + VR20, data = dat, index = c('hostId', 'time'), model='random')

	print(summary(fit))
	fit2 = coeftest(fit, vcov.=function(x) vcovHC(x,cluster='group'))

	FITS[[i]] = fit
	FITS[[2]] = fit2
	stargazer(FITS, type = 'text', apply.coef = multiply.1e2, apply.se = multiply.1e2, dep.var.caption = 'times $10^{-2}$')
	FITS = NA
}


dat = data.frame(fread('dynamic_subsampled.csv'))

runIds = seq(0,99)
rates = rep(c(.05,.1,.2,0.0,0.01), each=20)

uniqueIds = unique(dat$hostId)
#uniqueIds = uniqueIds[sample(length(uniqueIds), 0.01*length(uniqueIds))]

allRunIds = runIds

nhosts = length(uniqueIds)

dat = dat[dat$hostId %in% uniqueIds,]
dat = data.table(dat)

dat$VR00 = ifelse(dat$VR==0, 1, 0)
dat$VR01 = ifelse(dat$VR==0.01, 1, 0)
dat$VR05 = ifelse(dat$VR==0.05, 1, 0)
dat$VR10 = ifelse(dat$VR==0.10, 1, 0)
dat$VR20 = ifelse(dat$VR==0.20, 1, 0)
dat$VRnot0 = ifelse(dat$VR.init==0,1 , 0)

make.fit(data = dat[dat$time>4])
