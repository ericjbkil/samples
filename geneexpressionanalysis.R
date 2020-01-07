##ALL OF THIS STUFF IS FROM HW6. CHANGE ACCORDINGLY FOR THE PRACTICAL.
#step 3 
##import in the file with proper header and row names
##read in annotations as delimitated file with row.names also set
help("read.table")
data <- read.table("GSE9727_series_matrix-clean.txt", sep = "\t", row.names = 1, as.is = T, header = T)
annotations <- read.delim("GPL339-14792-clean.txt", row.names = 1)
nrow(data)
ncol(data)
View(data)
##take a log base 2 of your data to start making comparisons with more normalized data
Logdata <- log2(data)

#step 4 
##load in the proper package to perform the statistical analysis
##set the correct categorical values for the given groupings with "factor"
library("genefilter")
View(Logdata)
fac <- factor(c(0,0,0,0,1,1,1,2,2,2))
##need to set unequal variances and as a matrix
Anovaall <- rowFtests(as.matrix(Logdata), fac, var.equal = F)

help("factor")
View(Anovaall)

Logdata0and2 <- Logdata[,1:7]
View(Logdata0and2)
fac0and2 <- factor(c(0,0,0,0,1,1,1))
Anova0and2 <- rowFtests(as.matrix(Logdata0and2), fac0and2, var.equal = F)

Logdata0 <- Logdata[,1:4]
Logdata6 <- Logdata[,8:10]
Logdata0and6 <- cbind(Logdata0, Logdata6)
View(Logdata0and6)
fac0and6 <- factor(c(0,0,0,0,1,1,1))
Anova0and6 <- rowFtests(as.matrix(Logdata0and6), fac0and6, var.equal = F)

Allwithpvalues <- cbind(Logdata, Anovaall)
View(Allwithpvalues)
signpvalues <- Allwithpvalues[Allwithpvalues$p.value<0.05,]
nrow(signpvalues)
View(signpvalues)

twohrspvalues <- cbind(Logdata0and2, Anova0and2)
signpvalues2 <- twohrspvalues[twohrspvalues$p.value<0.05,]
nrow(signpvalues2)
View(signpvalues2)

sixhrspvalues <- cbind(Logdata0and6, Anova0and6)
signpvalues6 <- sixhrspvalues[sixhrspvalues$p.value<0.05,]
nrow(signpvalues6)
View(signpvalues6)


#step 5
##calculate averages for your logged data and set the variable with a given name
Logdata$"AveUntreated" <- rowMeans(Logdata[,1:4])
Logdata$"Ave2hrsTreated" <- rowMeans(Logdata[,5:7])
Logdata$"Ave6hrsTreated" <- rowMeans(Logdata[,8:10])
colnames(Logdata) 
View(Logdata)
str(Logdata)
##now calculate the differences between the log scores between experimental and control with a specific name
Logdata$"Logdiff2and0" <- as.matrix(Logdata[,12] - Logdata[,11])
Logdata$"Logdiff6and0" <- as.matrix(Logdata[,13] - Logdata[,11])
colnames(Logdata)

#step 6
##calculating fold changes by using the log diff scores
Logdata$"FoldChange2and0" <- ifelse(Logdata$'Logdiff2and0'<0,
                               1/2^Logdata$'Logdiff2and0'*-1,
                               2^Logdata$'Logdiff2and0')
Logdata$"FoldChange6and0" <- ifelse(Logdata$'Logdiff6and0'<0,
                                    1/2^Logdata$'Logdiff6and0'*-1,
                                    2^Logdata$'Logdiff6and0')
View(Logdata)

#step 7
##combine all of your logged data, with the annotations and the statistical test values

Combined2and0 <- cbind(Logdata, Anova0and2, annotations)
View(Combined2and0)
##ordering the data
#Orderedannotations <- annotations[order(row.names(annotations)),]

#step 8
##re-set the upregulated values as a new variable
upreg2and0 <- Combined2and0[Combined2and0$p.value< 0.05 & Combined2and0$FoldChange2and0>2,]
View(upreg2and0)
nrow(upreg2and0)
##see how many of these upregulated calls are unique or non-duplicated
#nrow(upreg[(unique(upreg$Gene.Symbol)),])
nrow(upreg2and0[(!duplicated(upreg2and0$Gene.Symbol)),])
##re-set the downregulated values as a new variable (be careful for </> arrows and spacing)
downreg2and0 <- Combined2and0[Combined2and0$p.value< 0.05 & Combined2and0$FoldChange2and0< -2,]
nrow(downreg2and0)
nrow(downreg2and0[(!duplicated(downreg2and0$Gene.Symbol)),])
View(downreg2and0)


Combined6and0 <- cbind(Logdata, Anova0and6, annotations)
upreg6and0 <- Combined6and0[Combined6and0$p.value< 0.05 & Combined6and0$FoldChange6and0>2,]
nrow(upreg6and0)
View(upreg6and0)
nrow(upreg6and0[(!duplicated(upreg6and0$Gene.Symbol)),])
downreg6and0 <- Combined6and0[Combined6and0$p.value< 0.05 & Combined6and0$FoldChange6and0< -2,]
nrow(downreg6and0)
nrow(downreg6and0[(!duplicated(downreg6and0$Gene.Symbol)),])



