
# test for skewness of distributions
library(tidyverse)
library(RSQLite)
library(cowplot)


get.data = function(filename){
  comboDb = dbConnect(SQLite(), dbname = filename)
  initExtension(comboDb)
  vaccineDF = dbGetQuery(comboDb, 'SELECT * FROM pooled_results')
  return(vaccineDF)
}

dir = '../../analysis/breadth_low_density/vaccine/results.sqlite'
vaccineDF = get.data(dir)

dir2 = '../../analysis/breadth_1_density/vaccine/results.sqlite'
df2 = get.data(dir2) %>% filter(vaccinationRate >0)

vaccineDF = rbind(vaccineDF,df2)

start = 21; end=21
drift.ranges = seq(start, end, length=(end-start)+1)
tests = vector('list', length(drift.ranges))
i=1
for(drift.cutoff in drift.ranges){
  y = vaccineDF %>% filter(cumulativeDrift >= drift.cutoff | tmrcaLimit == 1) %>% 
    group_by(vaccineImmuneBreadth, vaccinationRate) %>%
    dplyr::summarize(counts = n()) %>%
    arrange(vaccinationRate, vaccineImmuneBreadth) %>% data.frame()
  
  novac.counts = y$counts[y$vaccinationRate==0]

  z = y %>% group_by(vaccineImmuneBreadth, vaccinationRate) %>%
    mutate(pval = chisq.test(matrix(c(novac.counts,500- novac.counts, counts, 500-counts), nrow=2))$p.value) %>%
    data.frame()
  z$cutoff = drift.cutoff
  z$pval[z$counts < novac.counts] = 1
  z$novac = novac.counts
  
  tests[[i]] = z
  i=i+1
}

out = do.call(rbind, tests) %>% 
  select(vaccineImmuneBreadth, vaccinationRate, counts, novac, pval, cutoff) %>%
  mutate(significant = ifelse(pval<.05, 1, 0))
         
CUTOFFS = unique(out$cutoff)
for(this.cutoff in CUTOFFS){
  out.plot = ggplot(out %>%filter(cutoff == this.cutoff, vaccinationRate>0) %>% mutate(vaccineImmuneBreadth = paste('b =', vaccineImmuneBreadth)), aes(x=vaccinationRate, y = counts)) + 
    geom_bar(stat='identity', aes(fill = factor(significant))) + 
    scale_fill_manual(values = c('grey50', 'red'), guide=FALSE) +
    facet_wrap(~vaccineImmuneBreadth) + 
    geom_hline(aes(yintercept = novac)) + 
    #ggtitle(paste('>', this.cutoff,' a.u. cumulative evolution, including TMRCA > 10 years\nLine shows counts from unvaccinated sims\nRed bars show significantly more frequent excess evolution compared to no vaccination', sep='')) + 
    ylab('Number of simulations with excessive evolution\n 
         (>21 a.u. cumulative evolution or TMRCA > 10 years)')+
    xlab('Annual vaccination rate') 
  
  save_plot(paste('cutoffs_withTMRCA/',this.cutoff,'_cutoff.pdf', sep=''),out.plot, base_aspect_ratio = 1.2, base_height =7.5)
  
}


out %>% filter(cutoff == 21) %>% filter(vaccineImmuneBreadth==1, vaccinationRate==0)
vaccineDF %>% filter(vaccineImmuneBreadth==0.3) %>%
  group_by(vaccinationRate) %>%
  summarise(cumulativeIncidence = mean(cumulativeIncidence)) %>%
  filter(vaccinationRate > 0.01)


head(vaccineDF)
tmrcatest = vaccineDF %>% group_by(vaccinationRate, vaccineImmuneBreadth) %>%
  summarize(tmrcaLimit =log( sum(fluLike)/sum(tmrcaLimit))) %>% data.frame() %>% mutate(vaccinationRate = vaccinationRate*365)
fit = lm(tmrcaLimit ~ vaccinationRate*vaccineImmuneBreadth, data= tmrcatest %>% filter(is.finite(tmrcaLimit)))
summary(fit)
breadthplot = ggplot(tmrcatest, aes(x=vaccineImmuneBreadth, y=tmrcaLimit)) +
  geom_point() + #facet_wrap(~vaccineImmuneBreadth, scale = 'free') +
  geom_smooth(method= 'lm') +
  ylab('Log ratio Flu-like:TMRCA fraction')
rateplot = ggplot(tmrcatest, aes(x=vaccinationRate, y=tmrcaLimit)) +
  geom_point() + #facet_wrap(~vaccineImmuneBreadth, scale = 'free') 
  geom_smooth(method= 'lm') +
  ylab('Log ratio Flu-like:TMRCA fraction')
tab = ggtexttable(round(summary(fit)$coefficients,2))

plot_grid(rateplot, breadthplot, tab, nrow=2, rel_widths = c(1,1,2))
unique(vaccineDF$vaccinationRate)[1]

novac  = vaccineDF %>% filter(vaccinationRate==0, tmrcaLimit ==0, cumulativeDrift > 21)
vac2 = vaccineDF%>% filter(vaccinationRate == unique(vaccineDF$vaccinationRate)[10], 
                           vaccineImmuneBreadth==0.05, cumulativeDrift > 21)
qqplot(novac$cumulativeDrift, vac2$cumulativeDrift, xlim=c(20,30), ylim=c(20,30))
abline(a=0,b=1)
