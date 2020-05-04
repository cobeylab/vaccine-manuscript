library(RSQLite)
library(viridis)
library(cowplot)
library(tidyverse)
library(scales)

textSize = 12
pointSize = 1.0
lineSize = 1
plotDirectory='./plots'
plot_themes  = 	theme_classic() +
  theme(axis.line = element_line(size=1)) +
  theme(axis.ticks = element_line(size=0.5)) +
  theme(axis.ticks.length = unit(-0.1,'cm')) +
  #			theme(axis.ticks.margin = unit(0.5,'cm')) +
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
# 


theory.plotter = function(R0){
  gamma=.2
  mu = .000091
  beta = R0*(gamma+mu)
  left = -(gamma+2*mu)/2
  right = (gamma^2 + 4*mu*beta)^.5/2
  
  p = left+right
  return(p*365)
}

geteqI = function(beta,gamma,mu,p){
  R0 = beta/(gamma+mu+p)
  I_eq = mu/beta*(R0 - 1) - p/beta
  return(I_eq)
}

theoretical.threshold = function(beta,gamma, mu){
  left = -(gamma+2*mu)/2
  right = (gamma^2 + 4*mu*beta)^.5/2
  
  p = left+right
  return(p)
  
}

makePlot = function(vaccineDF, thresholds){
  
  plot.out = ggplot(data = vaccineDF, aes(x=vaccinationRate,y=extinct)) + 
    xlab('Vaccination rate') +
    ylab('Fraction extinct') +
    geom_vline(data = thresholds, aes(xintercept = rate, color = threshold))+
    geom_point(size=1, alpha=0.7) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
    theme(axis.ticks = element_line(size=.1))  +
    scale_x_continuous(breaks= trans_breaks(identity, identity, n = 3)) + 
    ggtitle(bquote(paste(italic('R'[0]),' = ', .(round(unique(vaccineDF$R0),2)))))+
    plot_themes
  
  return(plot.out)
}

get.vaccineDF = function(fileName){
  resultsDb = fileName
  comboDb = dbConnect(SQLite(), dbname = resultsDb)
  initExtension(comboDb)
  
  statusDF = dbGetQuery(comboDb, 'SELECT * FROM status')
  paramsDF = dbGetQuery(comboDb, 'SELECT beta, vaccinationRate, runId FROM parameters')
  vaccineDF = merge(statusDF, paramsDF)
  
  
  vaccineDF$vaccinationRate = vaccineDF$vaccinationRate*365
  dbDisconnect(comboDb)
  return(vaccineDF)
}

vaccineDF = get.vaccineDF('results.sqlite') %>% mutate(totalN = 5e7)
mu = 1/(30*365)
gamma = 0.2

vaccineDF  = vaccineDF %>%
  mutate(R0 = round(beta/(gamma + 1/(30*365)),2))
for(i in 1:nrow(vaccineDF)){
  vaccineDF$expected.threshold[i] = theoretical.threshold(vaccineDF$beta[i], gamma, mu) *365
}

summaryDF = vaccineDF %>% group_by(R0, vaccinationRate, expected.threshold, totalN) %>%
  summarise(extinct = sum(extinct)/length(runId)) %>% data.frame()

simplots = vector('list', length(nrow(out.df)))

out.df = summaryDF %>% 
  group_by(R0, totalN) %>% 
  filter(extinct==1) %>%  #note that some R0s don't satisfy this condition
  filter(vaccinationRate == min(vaccinationRate)) %>%
  data.frame()%>%
  rename('simulated' = 'vaccinationRate', 'theory' = 'expected.threshold') 
  
summaryDF = merge(summaryDF, out.df %>%select(R0, simulated))

x = vaccineDF[vaccineDF$R0==1.8,]
resolution = mean(diff(unique(x$vaccinationRate)))

out.df = out.df %>%
  mutate(resolution = resolution) %>%
  gather(key='Threshold', value = 'rate', 2:3) 
out.df$resolution[out.df$Threshold=='theory'] = NA
  

plotter.df = out.df %>% filter(Threshold == 'simulated')

summary.plot = ggplot(out.df, aes(x=R0, y= rate, color = Threshold, shape = Threshold)) + 
  geom_point(size=1) +
  xlab(expression(italic('R'[0]))) +
  #scale_shape_manual(values=c(1,3,4)) +
  ylab('Threshold annual vaccination rate') +
  scale_color_discrete(name = "Threshold") +
  stat_function(data = out.df %>% filter(Threshold == 'theory'), fun=theory.plotter) +
  geom_errorbar(aes(ymin = rate-resolution, ymax=rate+resolution), size=.5, width = .05) +
  scale_shape_manual(values= c(19,46)) +
  scale_color_manual(values = c('red','black'))
summary.plot

head(out.df)
summaryDF$labels = paste('R[0] == ',(summaryDF$R0))
out.df$labels = paste('R[0] == ',(out.df$R0))

breakfun = function(x) {
  len = 3
  out = pretty(x, n=len)
  return(out)       
}

sim.grid = ggplot(summaryDF, aes(x=vaccinationRate, y = extinct)) +
  geom_point(size= .5) +
  geom_vline(data = out.df, aes(xintercept = rate, color = Threshold))+
  facet_wrap(~labels, scale = 'free_x', labeller = label_parsed, nrow = 2) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black')) +
  scale_x_continuous(breaks = breakfun) +
  xlab('Annual vaccination rate') + ylab('Fraction extinct')


save_plot('./plots/threshold_summary.pdf', summary.plot, base_aspect_ratio=1.4)
save_plot('./plots/threshold_profile.pdf', sim.grid, base_height = 3.5, base_aspect_ratio = 2.3)

makeplot(vaccineDF, )