#
# Limma analysis of drug screen data with Helen
# 9th May 2016
#

setwd("~/Dropbox/DrJCampbell/drug-screen-limma/drug-data")

require(limma)

p1112_22RV1_mut_file <- "150921_p11_12_IS_DS_PROSTATE_22RV1_1112_CHD1_PE6_sA_zscores_summary_with_rep_zscores.txt"
p1112_22RV1_par_file <- "150921_p11_12_IS_DS_PROSTATE_22RV1_1112_CHD1_par_sA_zscores_summary_with_rep_zscores.txt"
p1112_RWPE1_mut_file <- "150921_p11_12_IS_DS_PROSTATE_RWPE1_1112_74F_sA_zscores_summary_with_rep_zscores.txt"
p1112_RWPE1_par_file <- "150921_p11_12_IS_DS_PROSTATE_RWPE1_1112_par_sA_zscores_summary_with_rep_zscores.txt"
p13_22RV1_par_file <- "IS_PROSTATE_22RV1_CHD1_PARENTAL_2015-09-23_complib_p13_sA_zscores_summary_with_rep_zscores.txt"
p13_22RV1_mut_file <- "IS_PROSTATE_22RV1_CHD1_PE6-LOF_2015-09-23_complib_p13_sA_zscores_summary_with_rep_zscores.txt"
p13_RWPE1_mut_file <- "IS_PROSTATE_RWPE1_CHD1_CHD1LOF7F4_2015-09-23_complib_p13_sA_zscores_summary_with_rep_zscores.txt"
p13_RWPE1_par_file <- "IS_PROSTATE_RWPE1_CHD1_PARENTAL_2015-09-23_complib_p13_sA_zscores_summary_with_rep_zscores.txt"



#
# limma analysis
#


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
	

	write.table(
		result,
		file=paste(prefix, "Limma_toptable.txt", sep="_"),
		col.names=TRUE,
		row.names=FALSE,
		sep="\t",
		quote=FALSE
		)

}


p1112_22RV1_par <- read.table(
	p1112_22RV1_par_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)
	
p1112_22RV1_mut <- read.table(
	p1112_22RV1_mut_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

p1112_RWPE1_par <- read.table(
	p1112_RWPE1_par_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

p1112_RWPE1_mut <- read.table(
	p1112_RWPE1_mut_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

p13_22RV1_par <- read.table(
	p13_22RV1_par_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

p13_22RV1_mut <- read.table(
	p13_22RV1_mut_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

p13_RWPE1_par <- read.table(
	p13_RWPE1_par_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

p13_RWPE1_mut <- read.table(
	p13_RWPE1_mut_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)
	

# Join up the columns we need from the par and
# mut data frames


# join_the compound (eg. FLUOROURACIL_50000) and
# function (e.g. anti-metabolite) to make row ids

p1112_rownames <- paste(
	p1112_22RV1_par[,"compound"],
	p1112_22RV1_par[,"Function"]
	)

p13_rownames <- paste(
	p13_22RV1_par[,"compound"],
	p13_22RV1_par[,"Function"]
	)


#
#
#

p1112_22RV1 <- as.matrix(
	cbind(
		p1112_22RV1_par[,c(24,25,26)],
		p1112_22RV1_mut[,c(24,25,26)]
		)
	)
colnames(p1112_22RV1) <- c(
	"par1",
	"par2",
	"par3",
	"mut1",
	"mut2",
	"mut3"
	)
rownames(p1112_22RV1) <- p1112_rownames


p1112_RWPE1 <- as.matrix(
	cbind(
		p1112_RWPE1_par[,c(24,25,26)],
		p1112_RWPE1_mut[,c(24,25,26)]
		)
	)
colnames(p1112_RWPE1) <- c(
	"par1",
	"par2",
	"par3",
	"mut1",
	"mut2",
	"mut3"
	)
rownames(p1112_RWPE1) <- p1112_rownames	

p13_22RV1 <- as.matrix(
	cbind(
		p13_22RV1_par[,c(24,25,26)],
		p13_22RV1_mut[,c(24,25,26)]
		)
	)
colnames(p13_22RV1) <- c(
	"par1",
	"par2",
	"par3",
	"mut1",
	"mut2",
	"mut3"
	)
rownames(p13_22RV1) <- p13_rownames

p13_RWPE1 <- as.matrix(
	cbind(
		p13_RWPE1_par[,c(24,25,26)],
		p13_RWPE1_mut[,c(24,25,26)]
		)
	)
colnames(p13_RWPE1) <- c(
	"par1",
	"par2",
	"par3",
	"mut1",
	"mut2",
	"mut3"
	)
rownames(p13_RWPE1) <- p13_rownames




# 
run_limma(
	na.omit(p1112_22RV1), # table of data
	design=design,
	prefix="p1112_22RV1"
	)

run_limma(
	na.omit(p1112_RWPE1), # table of data
	design=design,
	prefix="p1112_RWPE1"
	)

run_limma(
	na.omit(p13_22RV1), # table of data
	design=design,
	prefix="p13_22RV1"
	)

run_limma(
	na.omit(p13_RWPE1), # table of data
	design=design,
	prefix="p13_RWPE1"
	)



pdf(
	"scatter_plots_of_reps_22RV1_and_RWPE1_in_plate_11_12_13_screen_zscores_160509.pdf",
	width=9,
	height=9
	)
plot(
	as.data.frame(p1112_22RV1),
	col=rgb(0,0,0,0.5),
	pch=19,
	cex=0.5,
	main="22RV1 (plate 11/12)"
	)
plot(
	as.data.frame(p1112_RWPE1),
	col=rgb(0,0,0,0.5),
	pch=19,
	cex=0.5,
	main="RWPE1 (plate 11/12)"
	)
plot(
	as.data.frame(p13_22RV1),
	col=rgb(0,0,0,0.5),
	pch=19,
	cex=0.5,
	main="22RV1 (plate 13)"
	)
plot(
	as.data.frame(p13_RWPE1),
	col=rgb(0,0,0,0.5),
	pch=19,
	cex=0.5,
	main="RWPE1 (plate 13)"
	)

dev.off()