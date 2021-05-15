library(RSQLite)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(viridis)
library(cowplot)
library(tidyverse)
library(Hmisc)

textSize = 12
pointSize = 1.0
lineSize = 1
dynamicDir = '../../analysis/breadth_1_density/vaccine/'
staticDir = '../../analysis/static_density/vaccine/'
resultsDirs = c(dynamicDir,staticDir)

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

makePlot = function(summaryDF, vaccineDF, plotName){
  pd = position_dodge(width=.005)
  
  plot1 = ggplot(data = summaryDF, aes(x=vaccinationRate,y=meanDrift)) + 
    stat_summary(data = vaccineDF, aes(color = sim, y = cumulativeDrift),
                 fun.data = "mean_cl_boot",
                 geom = 'errorbar',
                 size = .5) + 
    # geom_errorbar(data = vaccineDF, aes(color = sim, y = cumulativeDrift),
    #               stat = 'summary',
    #               fun.ymin = percentile(.05),
    #               fun.ymax = percentile(.95),
    #               width = 0,
    #               size = .4,
    #               position = pd
    # ) +
    geom_point(aes(color = sim), size=.7) + 
    # geom_ribbon(data = vaccineDF, aes(color = sim, fill = sim, y = cumulativeDrift),
    #               stat = 'summary',
    #               fun.ymin = percentile(.05),
    #               fun.ymax = percentile(.95),
    #               #width = 0,
    #               size = .3,
    #               #position = pd,
    #               alpha= .05
    #               ) +
    xlab('Vaccination coverage') +
    ylab('Cumulative antigenic evolution') +
    #geom_smooth(data = vaccineDF, aes(colour = sim, y=cumulativeDrift), size=0.5) +
    scale_color_brewer(palette='Set1') +
    guides(colour= guide_legend('',title.position='top')) +
    #geom_smooth(data = fluDF, aes(colour=factor(sim), y=fluDrift), size=0.5, linetype = 'dashed', show_guide=FALSE) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    theme(axis.ticks = element_line(size=.1)) +
    ylim(c(0,30)) +
    plot_themes
  
  plot2 = ggplot(data = summaryDF, aes(x=vaccinationRate,y=meanInc)) + 
    stat_summary(data = vaccineDF, aes(color = sim, y = cumulativeIncidence),
                 fun.data = "mean_cl_boot",
                 geom = 'errorbar',
                 size = .5) + 
    # geom_errorbar(data = vaccineDF, aes(color = sim, y = cumulativeIncidence),
    #               stat = 'summary',
    #               fun.ymin = percentile(.05),
    #               fun.ymax = percentile(.95),
    #               width = 0,
    #               size = .4,
    #               position = pd
    # ) +
    geom_point(aes(color = sim), size=.7) + 
    xlab('Vaccination coverage') +
    ylab('Cumulative incidence') +
    #geom_smooth(data = vaccineDF, aes(colour = sim, y=cumulativeIncidence), size=0.5) +
    scale_color_brewer(palette='Set1') +
    #scale_color_viridis(discrete=TRUE, labels = c('All simulations', 'Surviving only')) +
    guides(color = guide_legend('',title.position='top', override.aes = list(colour='black', linetype=c('solid','dashed')))) +
    #geom_smooth(data = fluDF, aes(colour=factor(sim), y=fluInc), size=0.5, linetype = 'dashed',  show_guide=FALSE) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    theme(axis.ticks = element_line(size=.1))  +
    ylim(c(0,3)) +
    plot_themes 
	legend = get_legend(plot1)
	plot = plot_grid(plot1 + theme(legend.position = 'none'), plot2 + theme(legend.position = 'none'), labels = c('A','B'), align = 'h', ncol = 2)
	plot = plot_grid(plot, legend, nrow =2, rel_heights = c(1, .2))
	print(plot)
	save_plot(paste(plotName,'.pdf',sep=''), plot, ncol=2, nrow = 1, base_aspect_ratio = .9)
}

for(dir in resultsDirs){
  resultsDb = paste(dir,'results.sqlite',sep='')
  comboDb = dbConnect(SQLite(), dbname = resultsDb)
  initExtension(comboDb)
  if(dir == resultsDirs[1]){
    vaccineDF = dbGetQuery(comboDb, 'SELECT * FROM pooled_results WHERE tmrcaLimit==0')
    vaccineDF$sim = 'Dynamic'
  }
  else{
    df2 = dbGetQuery(comboDb, 'SELECT * FROM pooled_results WHERE tmrcaLimit == 0')
    df2$sim = 'Static'
  }
}

vaccineDF = rbind(vaccineDF, df2)

baseLineFlux = vaccineDF$meanFluxRate[vaccineDF$vaccinationRate==0 & vaccineDF$fluLike==1]
vaccineDF$cumulativeIncidence = vaccineDF$cumulativeIncidence/50000000

vaccineDF$acceleration = vaccineDF$meanFluxRate - mean(baseLineFlux)

summaryDF = vaccineDF %>%
  group_by(vaccinationRate, vaccineImmuneBreadth, sim) %>%
  summarise(meanDrift = mean(cumulativeDrift),
            meanInc = mean(cumulativeIncidence))


fluDF = vaccineDF %>%
  filter(fluLike == 1) %>%
  group_by(vaccinationRate, vaccineImmuneBreadth,sim) %>%
  summarise(fluDrift = mean(cumulativeDrift),
              fluInc = mean(cumulativeIncidence))

zerovacrow = summaryDF[summaryDF$vaccinationRate == 0,]
zerovacrow$sim = 'Static'
summaryDF = rbind(summaryDF, zerovacrow)

zerovac_append = vaccineDF %>% filter(vaccinationRate == 0) %>% mutate(sim = 'Static')
vaccineDF = rbind(vaccineDF, zerovac_append)

makePlot(summaryDF, vaccineDF, 'stat_dyn')

stsumm = summaryDF %>% filter(sim == 'Static') %>% rename(sdrift = meanDrift, stinc = meanInc)
dynsumm = summaryDF %>% filter(sim=='Dynamic') %>% rename(ddrift = meanDrift, dinc = meanInc)

outsumm = merge(stsumm,dynsumm, by='vaccinationRate') %>% 
  select(vaccinationRate,ddrift,sdrift,stinc,dinc) %>%
  mutate(dbig = ddrift > sdrift, ddiff = ddrift - sdrift, idiff=dinc-stinc)
outsumm

showdiff = function(rate){
  dat = vaccineDF %>% filter(vaccinationRate == toString(rate))
  a = dat %>% filter(sim=='Dynamic') %>% pull(cumulativeDrift) %>% smean.cl.boot()
  b = dat %>% filter(sim=='Dynamic') %>% pull(cumulativeIncidence) %>% smean.cl.boot()
  
  c = dat %>% filter(sim=='Static') %>% pull(cumulativeDrift) %>% smean.cl.boot()
  d = dat %>% filter(sim=='Static') %>% pull(cumulativeIncidence) %>% smean.cl.boot()
  
  out = do.call(rbind,list(a,b,c,d)) %>% data.frame()
  out$sim = c('Dynamic','Dynamic', 'Static','Static')
  out$stat = c('evo','inc','evo','inc')
  return(out)
}

showdiff(0.085)

cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeDrift, method = 'spearman')
cor.test(vaccineDF$vaccinationRate, vaccineDF$cumulativeIncidence, method = 'spearman')
wilcox.test(vaccineDF$cumulativeDrift[vaccineDF$sim=='Dynamic'], vaccineDF$cumulativeDrift[vaccineDF$sim=='Static'], alternative = 'less')
wilcox.test(vaccineDF$cumulativeIncidence[vaccineDF$sim=='Dynamic'], vaccineDF$cumulativeIncidence[vaccineDF$sim=='Static'], alternative = 'less')
