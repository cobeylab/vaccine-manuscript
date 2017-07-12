library(ggplot2)
library(cowplot)

textSize = 12
pointSize = 1.0
lineSize = 1
plotDirectory='./plots'
plot_themes  = 	theme_classic() +
  theme(axis.line = element_line(size=1)) +
  theme(axis.ticks = element_line(size=0.5)) +
  theme(axis.ticks.length = unit(-0.1,'cm')) +
  theme(axis.title.x=element_text(size=textSize)) + 
  theme(axis.text.x=element_text(size= textSize, margin=margin(5,5,5,5,'pt'))) + 
  theme(axis.title.y=element_text(size= textSize)) +
  theme(axis.text.y=element_text(size= textSize, margin=margin(5,5,5,5,'pt'))) +
  theme(plot.title=element_text(size=textSize+2)) +
  theme(plot.margin=unit(c(5,5,5,5),'mm')) +
  theme(legend.title=element_text(size=textSize)) +
  theme(legend.text=element_text(size=textSize)) +
  theme(legend.position ='bottom') +
  theme(legend.direction='horizontal') +
  theme(legend.margin = unit(0,'cm')) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.line = element_blank())

getOR = function(a,b){
  b=-b
  OR = 1-((a-b)/(1-(a-b)))/(a / (1-a))
  return(OR)
}

getve = function(a,b){
  b=-b
  RR = 1-((a-b)/a)
  return(RR)
}

fix.data = function(groups, means, ses, tstat){
  dat = data.frame(groups, means, ses)
  dat$lower = means-dat$se*tstat
  dat$upper = means+dat$se*tstat
  dat$VE = getve(means[1],means)
  dat$VE.lower = getve(dat$lower[1], dat$lower)
  dat$VE.upper = getve(dat$upper[1], dat$upper)
  dat$order = seq(nrow(dat),1)
  return(dat)
}

n=987500
groups = c('Constant','Vaccinated this season', 'last season', '2 seasons ago', '3 seasons ago', '4 seasons ago',
            'Vaccination rate = 1%', 'rate = 5%', 'rate = 10%')
means = c(.0994, -.0334, -.0278,-.0205,-.0174,-.0108, -.00926, -.0575,NA)
ses = c(.00230,.00324,.00330,.00236,.00217,.00157,.00467,.00645,NA)

tstat = qt(1-.025,n-2)

dyndat = fix.data(groups, means, ses, tstat)
dyndat$order = dyndat$order*2-1
dyndat$type = 'dynamic'

means = c(.0991, -.0465, -.0362, -.0265, -.0147, -.0128, -.00954, -.0275, -.0242)
ses = c(.00346, .00200, .00182, .00127, .00193, .00204, .00463, .00503, .00351)
n=1627600
tstat = qt(1-.025,n-2)

statdat = fix.data(groups,means,ses,tstat)
statdat$order = statdat$order*2
statdat$type = 'static'

pooldat = rbind(statdat,dyndat)
pooldat = pooldat[order(pooldat$order),]
pooldat$groups = ifelse(pooldat$type=='dynamic', paste(pooldat$groups,""),pooldat$groups)
pooldat$groups = factor(pooldat$groups, levels = pooldat$groups[order(pooldat$order)])

plot = ggplot(pooldat[1:(nrow(pooldat)-2),], aes(x=VE, y=groups)) + 
  geom_point(aes(color = type)) + 
  geom_errorbarh(aes(xmax = VE.upper, xmin=VE.lower, color = type), height=0) +
  scale_y_discrete(breaks = c(6,7,9,3,2,1,5,8), labels= c('rate = 10%','rate = 5%','Vaccination rate = 1%',
                                                             '4 seasons ago', '3 seasons ago', '2 seasons ago', 'last season','Vaccinated this season')) +
  xlab('1 - Risk ratio') +
  ylab('') + plot_themes + xlim(c(0,1.0)) +
  theme(legend.title=element_blank())


save_plot('VE.pdf', plot, ncol=2, nrow = 1, base_aspect_ratio = .8)
