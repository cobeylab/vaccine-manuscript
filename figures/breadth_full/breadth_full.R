library(RSQLite)
library(ggplot2)
library(sensitivity)
library(gridExtra)
library(grid)
library(scales)
library(viridis)
library(plyr)
library(cowplot)


textSize = 14
pointSize = 1.0
lineSize = 1
plot_themes  = 	theme_classic() +
  theme(axis.line = element_line(size=1)) +
  theme(axis.ticks = element_line(size=0.5)) +
  theme(axis.ticks.length = unit(-0.1,'cm')) +
  theme(axis.title.x=element_text(size=textSize)) + 
  theme(axis.text.x=element_text(size= textSize, margin=margin(5,5,5,5,'pt'))) + 
  theme(axis.title.y=element_text(size= textSize)) +
  theme(axis.text.y=element_text(size= textSize, margin=margin(5,5,5,5,'pt'))) +
  theme(plot.title=element_text(size=textSize+2)) +
  theme(plot.margin=unit(c(4,4,0,4),'mm')) +
  theme(legend.title=element_text(size=textSize)) +
  theme(legend.text=element_text(size=textSize)) +
  theme(legend.position ='bottom') +
  theme(legend.direction='horizontal') +
  theme(legend.box.margin = margin(0,0,0,0,'mm')) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.line = element_blank())

percentile = function(percentile){
  function(x) quantile(x, percentile, na.rm=T)
}

makeDriftDensity = function(vaccineDF,nbreaks=41){
  rates = sort(unique(vaccineDF$vaccinationRate))
  dens = data.frame(matrix(nrow = length(rates) * (nbreaks-1), ncol=4))
  names(dens) = c('vaccinationRate','cumulativeDrift','density','mean')
  for(i in 1:length(rates)){
    rate = rates[i]
    breaks = seq(0,35,length= nbreaks)
    subDF = vaccineDF[vaccineDF$vaccinationRate == rate,]
    drift.hist = hist(subDF$cumulativeDrift, breaks=breaks,xlim=c(0,c(max(vaccineDF$cumulativeDrift))), plot=FALSE)
    drift.dens = drift.hist$counts/(sum(drift.hist$counts))
    tmp = cbind(rep(rate,length(breaks)-1),drift.hist$breaks[1:length(breaks)-1],drift.dens, mean(subDF$cumulativeDrift[subDF$fluLike==1]))
    start.index = ((i-1)*(length(breaks)-1))+1
    dens[start.index:(start.index+length(breaks)-2),] = tmp
  }
  return(dens)
}

makeIncDensity = function(vaccineDF, nbreaks=41){
  rates = sort(unique(vaccineDF$vaccinationRate))
  dens = data.frame(matrix(nrow = length(rates) * (nbreaks-1), ncol=4))
  names(dens) = c('vaccinationRate','cumulativeIncidence','density','mean')
  for(i in 1:length(rates)){
    rate = rates[i]
    breaks = seq(0,3.2,length=nbreaks)
    subDF = vaccineDF[vaccineDF$vaccinationRate == rate,]
    drift.hist = hist(subDF$cumulativeIncidence, breaks=breaks,xlim=c(0,c(max(vaccineDF$cumulativeIncidence))), plot=FALSE)
    drift.dens = drift.hist$counts/(sum(drift.hist$counts))
    tmp = cbind(rep(rate,length(breaks)-1),drift.hist$breaks[1:length(breaks)-1],drift.dens,mean(subDF$cumulativeIncidence[subDF$fluLike==1]))
    start.index = ((i-1)*(length(breaks)-1))+1
    dens[start.index:(start.index+length(breaks)-2),] = tmp
  }
  return(dens)
}

makeDensityPlot = function(driftDens, incDens, plotName, meanDF, vaccineDF){
  
  plot1 = ggplot(data = driftDens, aes(x=vaccinationRate,y=cumulativeDrift)) + 
    xlab('Vaccination coverage') +
    ylab('Cumulative antigenic evolution') +
    geom_tile(aes(fill=density)) + 
    geom_point(data = vaccineDF, stat = 'summary',
               fun.y = 'mean',
               color = 'white', size = 1) + 
    # stat_summary(data = vaccineDF, aes(y = cumulativeDrift),
    #              fun.data = "mean_cl_boot",
    #              geom = 'errorbar',
    #              size = .5,
    #              color = 'white') + 
    stat_summary(data = vaccineDF, geom='errorbar',
                  fun.min = percentile(.05),
                  fun.max = percentile(.95),
                  width = 0,
                  color = 'white', size = .3) +
    scale_fill_viridis(option='plasma') +#, limits = c(0,1)) +
    guides(fill=guide_colorbar(barwidth=11,barheight=0.5, title.position = 'top')) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    theme(axis.ticks = element_line(size=.1))  +
    #geom_smooth(data=vaccineDF, aes(x=vaccinationRate,y=cumulativeDrift), color='white', size=0.6, se=TRUE, alpha=0.7) +
    plot_themes +
    ggtitle(paste("Breadth = ",unique(vaccineDF$vaccineImmuneBreadth),sep='')) +
    theme(plot.title = element_text(hjust=0))
  
  plot2 = ggplot(data = incDens, aes(x=vaccinationRate,y=cumulativeIncidence)) + 
    xlab('Vaccination coverage') +
    ylab('Cumulative incidence') +
    geom_tile(aes(fill=density)) + 
    geom_point(data = vaccineDF, stat = 'summary',
               fun.y = 'mean',
               color = 'white', size = 1) + 
    # stat_summary(data = vaccineDF, aes(y = cumulativeIncidence),
    #              fun.data = "mean_cl_boot",
    #              geom = 'errorbar',
    #              size = .5,
    #              color = 'white') + 
    stat_summary(data = vaccineDF, geom='errorbar',
                 fun.min = percentile(.05),
                 fun.max = percentile(.95),
                 width = 0,
                 color = 'white', size = .3) +
    scale_fill_viridis(option='plasma')+ #, limits = c(0,1)) +
    guides(fill=guide_colorbar(barwidth=11,barheight=0.5, title.position = 'top')) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    theme(axis.ticks = element_line(size=.1))  +
    #geom_smooth(data=vaccineDF, aes(x=vaccinationRate,y=cumulativeIncidence), color='white', size=0.6, se=TRUE, alpha=0.7) +
    #geom_smooth(data=vaccineDF[vaccineDF$fluLike==1,], aes(x=vaccinationRate,y=cumulativeIncidence), color='white', size=0.6, se=TRUE, alpha=0.7, linetype='dotted') +
    #geom_smooth(data=vaccineDF[vaccineDF$extinct==1,], aes(x=vaccinationRate,y=cumulativeIncidence), color='white', size=0.6, se=TRUE, alpha=0.7, linetype='dashed') +
    plot_themes
  
    return(list(plot1=plot1, plot2=plot2))
  
}

summarize.df = function(df){
  zeroDF = df[df$vaccinationRate==0,]
  out.df = ddply(df, .(vaccinationRate, vaccineImmuneBreadth), here(summarise),
                 meanDrift = mean(cumulativeDrift),
                 meanInc = mean(cumulativeIncidence),
                 count = sum(cumulativeDrift)/mean(cumulativeDrift),
                 drift.p.greater = wilcox.test(cumulativeDrift, zeroDF$cumulativeDrift, alternative='greater')$p.value,
                 drift.p.less = 1-drift.p.greater,
                 inc.p.greater = wilcox.test(cumulativeIncidence, zeroDF$cumulativeIncidence, alternative='greater')$p.value,
                 inc.p.less = 1-inc.p.greater,
                 drift.p = as.numeric(min(drift.p.greater, drift.p.less) < 0.05 & count >= 5),
                 inc.p = as.numeric(min(drift.p.greater, drift.p.less) < 0.05 & count >= 5)
  )
  return(out.df)
}

format.data = function(vaccineDF){
  vaccineDF = vaccineDF[vaccineDF$vaccineImmuneBreadth %in% c(1, 0.3,0.2, 0.1, 0.05,.5,.7),]
  
  vaccineDF$cumulativeIncidence = vaccineDF$cumulativeIncidence/50000000
  vaccineDF$vaccinationRate = as.numeric(as.character(vaccineDF$vaccinationRate))
  return(vaccineDF)
}

resultsDirs = c('../../analysis/breadth_1_density/vaccine/','../../analysis/breadth_low_density/vaccine/')

for(dir in resultsDirs){
  resultsDb = paste(dir,'results.sqlite',sep='')
  comboDb = dbConnect(SQLite(), dbname = resultsDb)
  initExtension(comboDb)
  if(dir == resultsDirs[1]){
    vaccineDF = dbGetQuery(comboDb, 'SELECT * FROM pooled_results')
  }
  else{
    df2 = dbGetQuery(comboDb, 'SELECT * FROM pooled_results')
  }
}

completeDF = rbind(vaccineDF, df2)
completeDF = format.data(completeDF)
vaccineDF = completeDF[completeDF$tmrcaLimit == 0,]

summaryDF = summarize.df(vaccineDF)
fluDF = summarize.df(vaccineDF[vaccineDF$fluLike==1,])

BREADTHS = sort(unique(vaccineDF$vaccineImmuneBreadth))
full.breadth.plots = list()
i=0
for(breadth in BREADTHS){
  vaccineSubDF = vaccineDF[vaccineDF$vaccineImmuneBreadth == breadth,]
  driftDens = makeDriftDensity(vaccineSubDF[vaccineSubDF$tmrcaLimit==0,])
  incDens = makeIncDensity(vaccineSubDF[vaccineSubDF$tmrcaLimit==0,])
  summary = ddply(vaccineDF[vaccineDF$tmrcaLimit==0,], .(vaccinationRate), summarise, 
                  mean.drift = mean(cumulativeDrift), 
                  mean.incidence = mean(cumulativeIncidence),
                  var.drift = var(cumulativeDrift),
                  var.incidence = var(cumulativeIncidence),
                  extinct = sum(extinct),
                  tmrcaLimit = sum(tmrcaLimit),
                  fluLike = sum(fluLike))
  
  plots = makeDensityPlot(driftDens,incDens,paste('density_breadth=',breadth,sep=''),summaryDF, vaccineSubDF)
  full.breadth.plots[[2*i+1]] = plots$plot1
  full.breadth.plots[[2*i+2]] = plots$plot2

  i=i+1
}

plot.labels = rep(NA, length(BREADTHS))
plot.labels[0:(length(BREADTHS)-1)*2+1] = LETTERS[c(1:length(BREADTHS))]
plot = plot_grid(plotlist=full.breadth.plots, labels = plot.labels, align = 'h', ncol = 4)
plotName = 'breadth_full'
save_plot(paste(plotName,'.pdf',sep=''), plot, ncol=4, nrow = 5, base_aspect_ratio = 1.1)
