#
# Limma analysis of CHD1 isogenic screens
# jamesc@icr.ac.uk, 1st Dec 2015
#

library(limma)

setwd("~/Dropbox/182_CHD1_Isogenic_siScreen_limma_analysis")


kinome_etc_data <- read.table(
	file="CHD1_isogenics_KS_TS_CGC_TNKS_WNT_libs_summary_with_rep_zscores.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

epi_lib_data <- read.table(
	file="CHD1_isogenics_EPI_lib_summary_with_rep_zscores.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)
mypastefun <- function(x){
	paste0(x, collapse="_")
}

kinome_etc_data_22RV1 <- as.matrix(kinome_etc_data[,4:9])
rownames(kinome_etc_data_22RV1) <- apply(kinome_etc_data[,1:3],1,mypastefun)

kinome_etc_data_RWPE1 <- as.matrix(kinome_etc_data[,10:15])
rownames(kinome_etc_data_RWPE1) <- apply(kinome_etc_data[,1:3],1,mypastefun)

epi_lib_data_22RV1 <- as.matrix(epi_lib_data[,4:9])
rownames(epi_lib_data_22RV1) <- gsub(" ", "", apply(epi_lib_data[,1:3],1,mypastefun))

epi_lib_data_RWPE1 <- as.matrix(epi_lib_data[,10:15])
rownames(epi_lib_data_RWPE1) <- gsub(" ", "", apply(epi_lib_data[,1:3],1,mypastefun))





# scaterplot kinome etc data for 22RV1 models
pdf("kinome_etc_data_22RV1_reps_scatter_plot.pdf",5,5)
plot(
	kinome_etc_data_22RV1[,1],
	kinome_etc_data_22RV1[,4],
	pch=19,
	col=rgb(1,0,0,0.25),
	xlim=c(-10,3),
	ylim=c(-10,3),
	xlab="Parental Z-scores",
	ylab="CHD1 null Z-scores",
	main="22RV1 models in\nkinome and other libraries"
	)
points(
	kinome_etc_data_22RV1[,2],
	kinome_etc_data_22RV1[,5],
	pch=19,
	col=rgb(0,1,0,0.25)
	)
points(
	kinome_etc_data_22RV1[,3],
	kinome_etc_data_22RV1[,6],
	pch=19,
	col=rgb(0,0,1,0.25)
	)
lines(c(-100,100),c(0,0), lty=2, col=rgb(0,0,0,0.5))
lines(c(0,0),c(-100,100), lty=2, col=rgb(0,0,0,0.5))
abline(0,1)
dev.off()

# scaterplot epi_lib data for 22RV1 models
pdf("epi_lib_data_22RV1_reps_scatter_plot.pdf",5,5)
plot(
	epi_lib_data_22RV1[,1],
	epi_lib_data_22RV1[,4],
	pch=19,
	col=rgb(1,0,0,0.25),
	xlim=c(-10,3),
	ylim=c(-10,3),
	xlab="Parental Z-scores",
	ylab="CHD1 null Z-scores",
	main="22RV1 models in\nepigentics library"
	)
points(
	epi_lib_data_22RV1[,2],
	epi_lib_data_22RV1[,5],
	pch=19,
	col=rgb(0,1,0,0.25)
	)
points(
	epi_lib_data_22RV1[,3],
	epi_lib_data_22RV1[,6],
	pch=19,
	col=rgb(0,0,1,0.25)
	)
lines(c(-100,100),c(0,0), lty=2, col=rgb(0,0,0,0.5))
lines(c(0,0),c(-100,100), lty=2, col=rgb(0,0,0,0.5))
abline(0,1)
dev.off()



# scaterplot kinome etc data for RWPE1 models
pdf("kinome_etc_data_RWPE1_reps_scatter_plot.pdf",5,5)
plot(
	kinome_etc_data_RWPE1[,1],
	kinome_etc_data_RWPE1[,4],
	pch=19,
	col=rgb(1,0,0,0.25),
	xlim=c(-10,3),
	ylim=c(-10,3),
	xlab="Parental Z-scores",
	ylab="CHD1 null Z-scores",
	main="RWPE1 models in\nkinome and other libraries"
	)
points(
	kinome_etc_data_RWPE1[,2],
	kinome_etc_data_RWPE1[,5],
	pch=19,
	col=rgb(0,1,0,0.25)
	)
points(
	kinome_etc_data_RWPE1[,3],
	kinome_etc_data_RWPE1[,6],
	pch=19,
	col=rgb(0,0,1,0.25)
	)
lines(c(-100,100),c(0,0), lty=2, col=rgb(0,0,0,0.5))
lines(c(0,0),c(-100,100), lty=2, col=rgb(0,0,0,0.5))
abline(0,1)
dev.off()

# scaterplot epi_lib data for RWPE1 models
pdf("epi_lib_data_RWPE1_reps_scatter_plot.pdf",5,5)
plot(
	epi_lib_data_RWPE1[,1],
	epi_lib_data_RWPE1[,4],
	pch=19,
	col=rgb(1,0,0,0.25),
	xlim=c(-10,3),
	ylim=c(-10,3),
	xlab="Parental Z-scores",
	ylab="CHD1 null Z-scores",
	main="RWPE1 models in\nepigentics library"
	)
points(
	epi_lib_data_RWPE1[,2],
	epi_lib_data_RWPE1[,5],
	pch=19,
	col=rgb(0,1,0,0.25)
	)
points(
	epi_lib_data_RWPE1[,3],
	epi_lib_data_RWPE1[,6],
	pch=19,
	col=rgb(0,0,1,0.25)
	)
lines(c(-100,100),c(0,0), lty=2, col=rgb(0,0,0,0.5))
lines(c(0,0),c(-100,100), lty=2, col=rgb(0,0,0,0.5))
abline(0,1)
dev.off()


#
# limma analysis
#

# kinome_etc_data_22RV1
# kinome_etc_data_RWPE1
# epi_lib_data_22RV1
# epi_lib_data_RWPE1

treatment<-factor(
	rep(
		c(
			"parental",
			"null"
			),
		c(3,3)),
		levels=c("parental","null")
		)

design<-model.matrix(~treatment)

run_limma <- function(x, design, prefix="limma_analysis"){
	fit<-lmFit(x, design)
	fit<-eBayes(fit)
	result<-topTable(
		fit,
		coef="treatmentnull",
		adjust="BH",
		number=nrow(x)
		)
	
	scores.withMeans <- cbind(
		x,
		apply(x[,1:3],1,mean),
		apply(x[,4:6],1,mean)
		)
	colnames(scores.withMeans) <- c(
		colnames(scores.withMeans)[1:6],
		"scores.par.mean",
		"scores.null.mean"
		)

	scores.withIDs <- data.frame(
		ID=rownames(scores.withMeans),
		scores.withMeans,
		row.names=NULL,
		stringsAsFactors=FALSE
		)
	colnames(scores.withIDs)[1] <- "ID"
	
	result <- data.frame(
		ID=rownames(result),
		result
		)
	
	scores.result <- merge(
		scores.withIDs,
		result,
		by.x="ID",
		by.y="ID"
		)

	pdf(file=paste(prefix, "significant_results_scatterplot.pdf", sep="_"), 5,5)
	plot(
		scores.result$scores.par.mean,
		scores.result$scores.null.mean,
		xlab="CHD1 parental mean z-score",
		ylab="CHD1 null mean z-score",
		pch=19,
		xlim=c(-10,3),
		ylim=c(-10,3),
		col=rgb(0,0,0,0.25)
		)
	points(
		scores.result$scores.par.mean[which(
			scores.result$adj.P.Val <= 0.05
			)],
		scores.result$scores.null.mean[which(
			scores.result$adj.P.Val <= 0.05
			)],
		pch=19,
		col=rgb(1,0,0,0.5)
		)
	lines(
		c(-2,-2),c(-100,100),
		lty=2
		)
	lines(
		c(-100,100),c(-2,-2),
		lty=2
		)
	i <- NULL
	for(i in 1:nrow(scores.result)){
		if(
			scores.result$scores.null.mean[i] <= -2 &
			scores.result$scores.par.mean[i] > -2 &
			scores.result$adj.P.Val[i] <= 0.05
			){
			text(
				scores.result$scores.null.mean[i],
				scores.result$scores.par.mean[i],
				strsplit(as.character(scores.result$ID[i]), split="_")[[1]][1],
				pos=2,
				cex=0.8
				)
		}
	}
	dev.off()
	write.table(
		scores.result,
		file=paste(prefix, "Limma_toptable.txt", sep="_"),
		col.names=TRUE,
		row.names=FALSE,
		sep="\t",
		quote=FALSE
		)

}





# kinome_etc_data_22RV1
run_limma(
	kinome_etc_data_22RV1,
	design=design,
	prefix="kinome_etc_data_22RV1"
	)

# kinome_etc_data_RWPE1
run_limma(
	kinome_etc_data_RWPE1,
	design=design,
	prefix="kinome_etc_data_RWPE1"
	)

# epi_lib_data_22RV1
run_limma(
	x=epi_lib_data_22RV1,
	design=design,
	prefix="epi_lib_data_22RV1"
	)

# epi_lib_data_RWPE1
run_limma(
	x=epi_lib_data_RWPE1,
	design=design,
	prefix="epi_lib_data_RWPE1"
	)

