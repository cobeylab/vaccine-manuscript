#!/usr/bin/env Rscript
library(data.table)


pretty= function(text){
	years = seq(1:20)
	a = paste(text, years, sep='.')
	return(a)
}

args = commandArgs(TRUE)
print(getwd())
runId = args[1]

file = paste('results/',runId,'/out.long',sep='')

#file = 'out.long'

dat = data.frame(fread(file, header=F, sep = ':', na.strings='n', colClasses=list(character=1:4)))

dat = dat[sample(nrow(dat),0.01*nrow(dat)),]
names(dat) = c('hostId', 'infectionDates', 'vaccinationDates', 'birthDate')

years = seq(1,20)

out = data.frame(matrix(nrow = nrow(dat), ncol = 41))
names(out) = c('hostId', paste('I', years, sep='.'), paste('V', years, sep='.'))

hostIds = paste(runId,dat$hostId,sep='.')
dat$hostId = hostIds

for(i in 1:nrow(dat)){
	print(i/nrow(dat))
	row = dat[i,]
	
	hostId = row$hostId
	birthDate = row$birthDate
	
	infectionDates = row$infectionDates
	vaccinationDates = row$vaccinationDates

	if(!is.na(infectionDates)){
		I = ifelse(years %in% (floor(as.numeric(strsplit(infectionDates,';')[[1]]) - 420/365) + 1), 1, 0)
	}else{
		I = rep(0,20)
	}
	
	if(!is.na(vaccinationDates)){
		V = ifelse(years %in% (floor(as.numeric(strsplit(vaccinationDates,';')[[1]]) - 300/365) + 1), 1, 0)
	}else{
		V = rep(0,20)
	}
		
	# V.5 = rep(min(1, sum(V[1:5])),20)
	
	# I.5 = rep(min(1, sum(I[1:5])),20)

	x = (c(hostId,I,V))
	
	out[i, ] = x
}

out = reshape(out, varying = c( pretty('I'), pretty('V')), idvar = 'hostId', direction='long')

write.csv(out, paste('results/',runId,'/long.csv', sep=''),row.names=F)

#write.csv(out, 'long.csv',row.names=F)