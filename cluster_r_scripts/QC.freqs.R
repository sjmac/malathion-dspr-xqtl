args = commandArgs(trailingOnly=TRUE)
filename = args[1]
Freqout = paste0("FREQ_",filename)
Totout = paste0("TOT_",filename)
NumberOfFounders = as.numeric(args[2])
t=read.table(filename,header=TRUE)
# biallelic
tt = sapply(t$REF,function(x) nchar(as.character(x)) == 1) & sapply(t$ALT,function(x) nchar(as.character(x)) == 1)
t=t[tt,]
# just genotypes
t2 = t[,-c(1:4,ncol(t))]
# add chromosome and position back later
chrpos = t[,c(1:2)]
# Total coverage
T = t2[,seq(1,ncol(t2),2)] + t2[,seq(2,ncol(t2),2)]
# observed frequency of reference allele
F = t2[,seq(1,ncol(t2),2)] / T
temp = colnames(T)
temp2 = gsub('_R', '', temp)
colnames(T) = temp2
colnames(F) = temp2
# bad coverage
badCoverage = apply(T[,1:NumberOfFounders],1,function(x) sum(x < 10))
T = data.frame(T[!badCoverage,])
F = data.frame(F[!badCoverage,])
chrpos = data.frame(chrpos[!badCoverage,])

####  check het for A2, it doesn't look bad...
#temp1 = F$A2[chrpos$CHROM=="chr3R"]
#temp2 = F$A2ST[chrpos$CHROM=="chr3R"]
#plot(temp1*(1-temp1),temp2*(1-temp2),pch=".")
#plot(temp1,temp2,pch=".")
#hist(temp2,breaks=100)
#### end checks

# fix A2
F$A2[chrpos$CHROM=="chr3R"] = F$A2ST[chrpos$CHROM=="chr3R"]

# fix B5
temp = F$B5Xiso1[chrpos$CHROM=="chr2L"]
# temp2 = F$B5[chrpos$CHROM=="chr2L"]
# plot(temp,temp2,pch=".")
# this is sort of brutal, but I am forcing calls to zero or 1
# if "strongly" heterozygous I call it alt, if strongly homo for 1 I call it a 1
# since this is B5 against iso1, any ALT alleles should suggest that is the B5 chromosome
temp3 = temp
temp3[temp3>0.95] = 1
temp3[temp3<=0.95] = 0
# table(temp3)
F$B5[chrpos$CHROM=="chr2L"] = temp3

# drop columns...not needed
F = F[,-which(names(F) %in% c("A2ST","B5Xiso1"))]
T = T[,-which(names(T) %in% c("A2ST","B5Xiso1"))]
NumberOfFounders = NumberOfFounders - 2

# drop SNPs that have high heterozygosity
H = t(apply(F[,1:NumberOfFounders],1,function(x) 2*x*(1-x)))
# summary(H)
# barplot(apply(H,2,function(x) sum(x>0.1)))
badHeterozygosity = apply(H,1,function(x) sum(x>0.10))
# table(badHeterozygosity)
T = T[badHeterozygosity==0,]
F = F[badHeterozygosity==0,]
chrpos = chrpos[badHeterozygosity==0,]

# drop SNPs that are not informative, that is near zero or 1 over all founders
badInformativeness = apply(F[,1:NumberOfFounders],1,function(x) (sum(x)<0.98 | sum(x)>(NumberOfFounders-0.98)))
# table(badInformativeness)
T = T[!badInformativeness,]
F = F[!badInformativeness,]
chrpos = chrpos[!badInformativeness,]

# output as standard R table (easier to read in)
write.table(cbind(chrpos,F),Freqout)
write.table(cbind(chrpos,T),Totout)


