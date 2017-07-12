library(ggplot2)
library(viridis)
library(cowplot)
textSize = 12
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
  theme(plot.margin=unit(c(5,5,5,5),'mm')) +
  theme(legend.title=element_text(size=textSize)) +
  theme(legend.text=element_text(size=textSize)) +
  theme(legend.position ='bottom') +
  theme(legend.direction='horizontal') +
  theme(legend.margin = unit(0,'cm')) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.line = element_blank())

beta = 0.36 #transmission rate (day^-1)
gamma = 0.2 #recovery rate (day^-1)
nu = 1/(30*365) #birth and death rate (day^-1)
R0 = beta/(gamma+nu)
c = .07

get.gamma = function(d){
    delta.mean = .6
    delta.sd = .3
    alpha = (delta.mean*delta.mean)/(delta.sd*delta.sd)
    beta = (delta.sd*delta.sd)/(delta.mean)
    
    p.d = dgamma(d, shape = alpha, rate=1/beta)
    return(p.d*.08)
}

get.threshold = function(P){
  threshold = uniroot(f = function(d) beta*(1-get.R.mut(d,P,R0,beta,gamma,c))-(gamma+nu), c(-5,10))
  return(threshold$root)
}


get.R.mut = function(d,P,R0,beta,gamma,c){
  cd=c*d
  cd[cd>1] = 1
  
  S.eq = 1/R0
  I.eq = nu*(R0*(1-P)-1)/beta
  R.eq = 1-S.eq -I.eq
  
  S.eq[I.eq < 0] = 1-P[I.eq<0]
  R.eq[I.eq < 0] = P[I.eq<0]
  I.eq[I.eq < 0] = 0
   
  R.mut = (1 - cd)*R.eq 
  R.mut[R.mut>1] = 1
  return(R.mut)
}

s.w = function(d,P,gamma,nu,R0,beta,c,delta){

  R.mut = get.R.mut(d,P,R0,beta,gamma,c)
  S.mut = 1-R.mut
  s =  beta*S.mut - (gamma+nu)
  s[s<0] = NA
  mu = get.gamma(d)
  return(s*mu)
}

P = seq(0,1, length=100)
D = seq(0,10, length=100)

df = expand.grid(P=P, D=D)
df$R0 = beta/(gamma+nu)
#eqs = get.R.mut(df$D,df$P,gamma,df$R0,beta,c)
df$R.mut =  get.R.mut(df$D,df$P,df$R0,beta,gamma,c)
df$S.mut = 1 - df$R.mut
df$s = beta*df$S.mut - (gamma+nu)
df$s[df$s<0] = 0
df$threshold = NA
for(i in 1:nrow(df)){
  df$threshold[i] = get.threshold(df$P[i])
}

sampleDistanceIndex = 3

plot1 = ggplot(data = df, aes(x=P,y=D)) + 
  geom_tile(aes(fill=s)) + scale_fill_viridis() +
  xlab(expression(paste('Fraction of new hosts vaccinated ', italic('p')))) +
  ylab('Antigenic distance \nbetween mutant and resident') +
  guides(fill=guide_colorbar(expression(paste('Invasion fitness ', italic('s'))),barwidth=11,barheight=0.5, title.position = 'top', trans='log')) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_line(aes(y=threshold), color = 'white', size=1) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
  theme(axis.ticks = element_line(size=.1))  +
  geom_hline(yintercept = unique(df$D)[sampleDistanceIndex], color='red') +
  plot_themes

plot2 = ggplot(data = df, aes(x=D)) + 
  stat_function(aes(D),fun=get.gamma) + 
  xlab('Antigenic distance') +
  ylab('Probability density') +
  plot_themes

df$s.weighted = df$s*get.gamma(df$D)
df$s.weighted = ifelse(df$s.weighted<0, 0, df$s.weighted)

sampleRates = unique(df$P)[c(1,11,12,20,82)]
sampledf = df[df$P>0.42 & df$P<0.5,]

sub = df[df$D==unique(df$D)[sampleDistanceIndex],]
plot3 = ggplot(sub,aes(x=P,y=s)) + 
          geom_line() +
          xlab(expression(paste('Fraction of new hosts vaccinated ', italic('p')))) +
          ylab(expression(paste('Invasion fitness ', italic('s')))) +
          plot_themes

plot = plot_grid(plot1,plot3,plot2, labels = c('A','B','C'), nrow=1)
plotName = 'invasion'
save_plot(paste('',plotName,'.pdf',sep=''), plot, ncol=3, nrow = 1, base_aspect_ratio = 2)


