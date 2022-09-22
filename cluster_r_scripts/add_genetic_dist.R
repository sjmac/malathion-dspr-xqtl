args = commandArgs(trailingOnly=TRUE)
filename = args[1]

#blah.txt release 6 recombination from flybase
#cat blah.txt | awk '{split($1,a,"-\\["); split($2,b,":"); cM=a[2]; gsub("\\]","",cM);\
#   chr=b[1]; myrange=b[2]; print chr, myrange, cM}' | sed 's:\.\.: :' |\
#   awk '{printf("chr%s\t%s\t%s\n",$1,($2+$3)/2,$4)}' > flymap.r6.txt

library(tidyverse)
add_genetic = function(df){
	df$cM = rep(NA,nrow(df))
	fm=read.table("helperfiles/flymap.r6.txt",header=FALSE)
	colnames(fm)=c("CHROM","POS","cM")
	library(splines)
	for(chrs in c("chrX","chr2L","chr2R","chr3L","chr3R")){
		fmX = fm %>% filter(CHROM==chrs)
		out = ksmooth(fmX$POS,fmX$cM,kernel="normal",bandwidth=3e6)
		f_of_x = splinefun(out$x,out$y)
		temp = f_of_x(df$POS[df$CHROM==chrs])
		df$cM[df$CHROM==chrs] = temp
		}
	df
	}

xx=read.table(filename)
xx2=add_genetic(xx)
filenameout = gsub(".txt",".cM.txt",filename)
write.table(xx2,filenameout)


