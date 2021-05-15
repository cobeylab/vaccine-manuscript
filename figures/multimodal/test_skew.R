# test for skewness of distributions
library(tidyverse)
library(RSQLite)
library(cowplot)


get_data = function(filename){
  comboDb = dbConnect(SQLite(), dbname = filename)
  initExtension(comboDb)
  vaccineDF = dbGetQuery(comboDb, 'SELECT * FROM pooled_results')
  dbDisconnect(comboDb)
  return(vaccineDF)
}

resultsDirs = c('../../analysis/breadth_low_density/vaccine/','../../analysis/breadth_1_density/vaccine/')

vaccineDF = resultsDirs %>%
  paste0(.,'results.sqlite') %>%
  map_dfr(get_data)
n_replicates = vaccineDF %>%
  group_by(vaccinationRate, vaccineImmuneBreadth) %>%
  summarise(counts = n()) %>%
  pull(counts) #%>%
  unique()

start = 21; end=21
drift.ranges = seq(start, end, length=(end-start)+1)
tests = vector('list', length(drift.ranges))
i=1

novac_data = vaccineDF %>% 
  filter(vaccinationRate==0)

mean_novac_drift = novac_data %>% filter(tmrcaLimit == 0) %>% pull(cumulativeDrift) %>% mean()
mean_novac_inc = novac_data %>% filter(tmrcaLimit == 0) %>% pull(cumulativeIncidence) %>% mean()

for(drift.cutoff in drift.ranges){
  test_data = vaccineDF %>% 
    group_by(vaccineImmuneBreadth, vaccinationRate) %>%
    summarize(replicates = n(),
              counts = sum(cumulativeDrift >= drift.cutoff | tmrcaLimit == 1)) %>%
    arrange(vaccinationRate, vaccineImmuneBreadth) %>%
    ungroup()
  
  
  novac.counts = novac_data %>%
    filter(cumulativeDrift >= drift.cutoff | tmrcaLimit == 1) %>% 
    nrow()
  novac.replicates = nrow(novac_data)
  

  test_data$pval = NA
  for(j in 1:nrow(test_data)){
    test_data$pval[j] = fisher.test(matrix(c(novac.counts, novac.replicates - novac.counts, 
                                     test_data$counts[j], test_data$replicates[j]-test_data$counts[j]), nrow=2),
                            alternative = 'less')$p.value 

  }

  test_data$cutoff = drift.cutoff
  #z$pval[z$counts < novac.counts] = 1
  test_data$novac = novac.counts
  
  tests[[i]] = test_data
  i=i+1
}

out = do.call(rbind, tests) %>% 
  select(vaccineImmuneBreadth, vaccinationRate, counts, replicates, novac, pval, cutoff) %>%
  mutate(significant = ifelse(pval<.05, 1, 0))

CUTOFFS = unique(out$cutoff)
for(this.cutoff in CUTOFFS){
  out.plot = ggplot(out %>% 
                      filter(cutoff == this.cutoff, vaccinationRate>0) %>% 
                      mutate(vaccineImmuneBreadth = paste('b =', vaccineImmuneBreadth)), 
                    aes(x=vaccinationRate, y = counts/replicates)) + 
    geom_bar(stat='identity', aes(fill = factor(significant))) + 
    scale_fill_manual(values = c('grey50', 'red'), guide=FALSE) +
    facet_wrap(~vaccineImmuneBreadth) + 
    geom_hline(aes(yintercept = novac/novac.replicates)) + 
    #ggtitle(paste('>', this.cutoff,' a.u. cumulative evolution, including TMRCA > 10 years\nLine shows counts from unvaccinated sims\nRed bars show significantly more frequent excess evolution compared to no vaccination', sep='')) + 
    ylab('Number of simulations with excessive evolution\n 
         (>21 a.u. cumulative evolution or TMRCA > 10 years, n=500)')+
    xlab('Annual vaccination coverage') 
  
  save_plot(paste('cutoffs_withTMRCA/',this.cutoff,'_cutoff.pdf', sep=''),out.plot, base_aspect_ratio = 1.2, base_height =7.5)
  
}


# 
# head(vaccineDF)
# tmrcatest = vaccineDF %>% group_by(vaccinationRate, vaccineImmuneBreadth) %>%
#   summarize(tmrcaLimit =log( sum(fluLike)/sum(tmrcaLimit))) %>% data.frame() %>% mutate(vaccinationRate = vaccinationRate*365)
# fit = lm(tmrcaLimit ~ vaccinationRate*vaccineImmuneBreadth, data= tmrcatest %>% filter(is.finite(tmrcaLimit)))
# summary(fit)
# breadthplot = ggplot(tmrcatest, aes(x=vaccineImmuneBreadth, y=tmrcaLimit)) +
#   geom_point() + #facet_wrap(~vaccineImmuneBreadth, scale = 'free') +
#   geom_smooth(method= 'lm') +
#   ylab('Log ratio Flu-like:TMRCA fraction')
# rateplot = ggplot(tmrcatest, aes(x=vaccinationRate, y=tmrcaLimit)) +
#   geom_point() + #facet_wrap(~vaccineImmuneBreadth, scale = 'free') 
#   geom_smooth(method= 'lm') +
#   ylab('Log ratio Flu-like:TMRCA fraction')
# tab = ggtexttable(round(summary(fit)$coefficients,2))
# 
# plot_grid(rateplot, breadthplot, tab, nrow=2, rel_widths = c(1,1,2))
# unique(vaccineDF$vaccinationRate)[1]
# 
# novac  = vaccineDF %>% filter(vaccinationRate==0, tmrcaLimit ==0, cumulativeDrift > 21)
# vac2 = vaccineDF%>% filter(vaccinationRate == unique(vaccineDF$vaccinationRate)[10], 
#                            vaccineImmuneBreadth==0.05, cumulativeDrift > 21)
# qqplot(novac$cumulativeDrift, vac2$cumulativeDrift, xlim=c(20,30), ylim=c(20,30))
# abline(a=0,b=1)
