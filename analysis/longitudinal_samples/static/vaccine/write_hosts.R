library(ggplot2)
library(plyr)
library(reshape)
library(scales)
library(grid)
library(gridBase)
library(gridExtra)
library(RSQLite)
library(doParallel)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

#setwd('/Volumes/frankwen/vaccination/dist/vaccine')

fill.extinct = FALSE
subsample = as.numeric(args[1])
out.file = args[2]
# if(subsample){
	# out.file = "dynamic_subsampled.csv"
# }else{
	# out.file = 	"dynamic.csv"
# }
parallel= FALSE

get.hosts = function(runId, vaccineDF){
	print(runId)

	runDir = paste('./results/',runId,sep='')
	
	hostsFile = paste(runDir,'/long.csv',sep='')
	
	hosts = fread(hostsFile, sep=',')
	hosts = data.frame(hosts)
	hosts$hostId = rep(seq(runId * 25000, (1+runId)*25000-1), 20 )
	
	hostIds = unique(hosts$hostId)
	
	
	
	if(subsample){
		hostIds = hostIds[sample(length(hostIds), 0.10*length(hostIds))]
		hosts = hosts[hosts$hostId %in% hostIds,]
	}
	
	extinct = vaccineDF$extinct[vaccineDF$runId == runId]
	
	hosts$extinct = extinct
	hosts$VR.init = vaccineDF$vaccinationRate[vaccineDF$runId == runId]
	hosts$VR = vaccineDF$vaccinationRate[vaccineDF$runId == runId]
	
	
	if(extinct & fill.extinct){
		lastDate = vaccineDF$lastDate[vaccineDF$runId == runId]
		timeLeft = 20.0 - lastDate
		P.survive = exp(-1/30 * timeLeft)
		
		reborn.hostIds = hostIds[sample(length(hostIds), (1-P.survive)*(length(hostIds)))]
		
		hosts[,c('I','V')][hosts$hostId %in% reborn.hostIds,] = 0
		#hosts$birthDate[hosts$hostId %in% reborn.hostIds] = lastDate
		hosts$VR[hosts$time > lastDate] = 0
		
	}
	
	if(extinct & !fill.extinct){
		lastDate = vaccineDF$lastDate[vaccineDF$runId == runId]

		hosts[,c('I','V','VR')][hosts$time -1 + 300/365 > lastDate,] = NA
		
	}
		
	hosts$runId = runId
	#hosts = data.frame(data.matrix(hosts))
	#dbDisconnect(inDb)
	hosts$extinctNow = ifelse(is.na(hosts$VR), 1, 0)
	
	hosts$I = as.numeric(hosts$I)
	hosts$V = as.numeric(hosts$V)
	hosts$VR = as.numeric(hosts$VR)
		
	return(hosts[,c('hostId','time','I','V','VR','VR.init','extinctNow','runId')])
}


get.vaccineDF = function(resultsDb){
	comboDb = dbConnect(SQLite(), dbname = resultsDb)
	initExtension(comboDb)

	vaccineDF = dbGetQuery(comboDb, 'SELECT * FROM pooled_results')
	vaccineDF$cumulativeIncidence = vaccineDF$cumulativeIncidence/50000000
	vaccineDF$vaccinationRate = round(vaccineDF$vaccinationRate*365,2)
	vaccineDF = vaccineDF[vaccineDF$tmrcaLimit==0,]
	dbDisconnect(comboDb)
	return(vaccineDF)
}

resultsDb = './results.sqlite'
vaccineDF = get.vaccineDF(resultsDb)
RUNIDS = vaccineDF$runId

if(parallel){
	registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
	out = foreach(i = RUNIDS,
				.packages = c('data.table'),
				.combine = rbind,
				.export=c('get.hosts','make.fit'))%dopar%{
		print(i)
		get.hosts(i,vaccineDF)
	}
	print('writing...')
	write.csv(out, file = out.file, row.names=FALSE)
}else{
	for(i in RUNIDS){
		if(i == RUNIDS[1]){
			out = get.hosts(i,vaccineDF)
			write.table(out, file = out.file, sep = ',', na = '',row.names=FALSE)
		}else{
			out = get.hosts(i,vaccineDF)
			write.table(out, file = out.file, sep = ',', na='', row.names=FALSE, col.names=FALSE, append=TRUE)
		}
	}	
}
# 
# print('Fitting...')
# fit = glm(data = out, infected ~ vaccinated + vaccinationRate + vaccinated:vaccinationRate)
# coeffs = data.frame(fit$coefficients)
# print(summary(fit))
# write.csv(coeffs , file = 'fit.csv')
