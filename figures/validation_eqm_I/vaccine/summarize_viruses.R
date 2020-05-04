library(RSQLite)
library(data.table)
library(plyr)

fn = 'samples.sqlite'

df = dbConnect(SQLite(), dbname = fn)
initExtension(df)

samples = dbGetQuery(df, 'SELECT * FROM hosts')

samples = count(samples,c('deme','date','ag1','ag2'))
names(samples) = c('deme','date','ag1','ag2','count')
samples = data.table(samples)

outdf = dbConnect(SQLite(), dbname = 'abundances.sqlite')
initExtension(outdf)

for(demeName in unique(samples$deme)){
	tmp = samples[samples$deme == demeName,]
	tmp[,abundance:=count/sum(count),by=date]
	tmp = data.frame(tmp)
	tmp = tmp[,c('deme','date','ag1','ag2','abundance')]
	dbBegin(outdf)
	dbWriteTable(outdf, 'refStrains', tmp,append=TRUE)
	dbCommit(outdf)
}

