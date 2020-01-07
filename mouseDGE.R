#install and load required packages
install.package("dplyr")
library(dplyr)
library(genefilter)
library(multtest)
library(tibble)
library(reshape2)

#pull .tsv files from working directory and set as a table with 
#headers and tab separators
KO_10 <- read.table('abundance10_KO.tsv', sep = '\t', header = TRUE)
CTRL_2 <- read.table('abundance2_CTRL.tsv', sep = '\t', header = TRUE)
CTRL_5 <- read.table('abundance5_CTRL.tsv', sep = '\t', header = TRUE)
KO_5 <- read.table('abundance5_KO.tsv', sep = '\t', header = TRUE)
KO_8 <- read.table('abundance8_KO.tsv', sep = '\t', header = TRUE)
CTRL_9 <- read.table('abundance9_CTRL.tsv', sep = '\t', header = TRUE)

#make sure the data looks alright at first-glance
View(KO_10)

#change tpm columns names for individual files
names(KO_10)[names(KO_10) == "tpm"] <- "tpm_KO_10"
names(KO_5)[names(KO_5) == "tpm"] <- "tpm_KO_5"
names(KO_8)[names(KO_8) == "tpm"] <- "tpm_KO_8"
names(CTRL_2)[names(CTRL_2) == "tpm"] <- "tpm_CTRL_2"
names(CTRL_5)[names(CTRL_5) == "tpm"] <- "tpm_CTRL_5"
names(CTRL_9)[names(CTRL_9) == "tpm"] <- "tpm_CTRL_9"

#create new datasets with just tpm columns
tpm10kow <- KO_10[c(1,5)]
tpm10ko <- KO_10[c(5)]
tpm5ko <- KO_5[c(5)]
tpm8ko <- KO_8[c(5)]
tpm2ctrl <- CTRL_2[c(5)]
tpm5ctrl <- CTRL_5[c(5)]
tpm9ctrl <- CTRL_9[c(5)]

#new dataset with gene ids and tpms
combined <- cbind(tpm10kow, tpm5ko, tpm8ko, tpm2ctrl, tpm5ctrl, tpm9ctrl)
combinedwoid <- cbind(tpm10ko, tpm5ko, tpm8ko, tpm2ctrl, tpm5ctrl, tpm9ctrl)
View(combined)



#run the filter and keep rows with tpm greater than 1
filtered <- filter_at(combined, vars(starts_with("tpm")), all_vars(. >= 1))
filtered2 <- filter_at(combinedwoid, vars(starts_with("tpm")), all_vars(. >= 1))
geneid <- filtered[c(1)]
View(geneid)
View(filtered2)

#perform log2 transformations
log2data <- log2(filtered2)
View(log2data)

#compute unequal variance F test - wrong (didn't log2 transform)
fac <- factor(c(0,0,0,1,1,1))
Anova2 <- rowFtests(as.matrix(log2data), fac, var.equal = FALSE)
View(Anova2)

#combine data with p values
wpvalues2 <- cbind(log2data, Anova2)

#perform BH correction and convert to dataframe
newANOVA2 <- mt.rawp2adjp(Anova2$p.value, proc = c("BH"), alpha = 0.05, na.rm = FALSE)
adjustedpvalue2<-as.data.frame(newANOVA2[c(1,2)])

#merge p values with rest of dataset
total2 <- cbind(wpvalues2, adjustedpvalue2)

#identify number of sign p-values
signftest2 <- total2[total2$p.value<= 0.05,]
signbh2 <- total2[total2$adjp.BH <= 0.05,]

View(total2)
#change exponent format to decimal
#options(scipen = 999)

#split total into without geneids
totalwoid2 <- total2[c(2,3,4,5,6,7)]
View(totalwoid2)

#put gene ids back with log2data
newlog <- cbind(geneid, log2data)
View(newlog)

#compute averages
newlog$"AveKO" <- rowMeans(newlog[,2:4])
newlog$"AveCTRL" <- rowMeans(newlog[,5:7])
colnames(newlog) 
View(newlog)
str(newlog)

##now calculate the differences between the log scores between experimental and control with a specific name
newlog$"Logdiff" <- as.matrix(newlog[,8] - newlog[,9])

#calculate fold change
newlog$"FoldChange" <- ifelse(newlog$'Logdiff'<0,
                               1/2^newlog$'Logdiff'*-1,
                               2^newlog$'Logdiff')

#determine up- and down-regulated
final <- cbind(newlog, Anova2)
upreg <- final[final$p.value< 0.05 & final$FoldChange>2,]
nrow(upreg)
View(upreg)
downreg <- final[final$p.value< 0.05 & final$FoldChange< -2,]
nrow(downreg)
View(downreg)

withBH <- cbind(newlog, adjustedpvalue2)
withBHanno <- merge(withBH, annotations, by = "target_id")

#read in annotations file
annotations <- read.delim("mart_export.txt")
View(annotations)
View(final)

#and merge annotations to the final, upreg, downreg file
names(annotations)[names(annotations) == "Transcript.stable.ID.version"] <- "target_id"
finalwanno <- merge(final, annotations, by = "target_id")
View(finalwanno)

upregwanno <- merge(upreg, annotations, by = "target_id")
downregwanno <- merge(downreg, annotations, by = "target_id")

View(upregwanno)
View(downregwanno)

#determine number of uniquely mapped
nrow(upregwanno[(!duplicated(upregwanno$Gene.name)),])
nrow(downregwanno[(!duplicated(downregwanno$Gene.name)),])

#combine upreg and downreg
upreganddownreg <- rbind(upregwanno, downregwanno)
View(upreganddownreg)

#write the table
write.table(upreganddownreg, file = "MouseDE.txt", sep = "\t")


#create a heat map
install.packages("gplots")
library("gplots")
install.packages("RColorBrewer")
library(RColorBrewer)

#all DE
#create sorted differential g.e. list
upreganddownregsort <- upreganddownreg[order(upreganddownreg$FoldChange),]
#check if any duplicates in Gene.name column and the index of these duplicates
duplicated(upreganddownregsort$Gene.name)
#remove indeces with a duplicate value
keepfrombiglist <- upreganddownregsort[-c(4,8,10,12,17,34,60,79,85,109,130),]
#check to see if all values are FALSE
View(keepfrombiglist)
duplicated(keepfrombiglist$Gene.name)
#rename column names to something more concise/readable
names(keepfrombiglist)[names(keepfrombiglist) == "tpm_CTRL_2"] <- "2_CTRL"
names(keepfrombiglist)[names(keepfrombiglist) == "tpm_CTRL_5"] <- "5_CTRL"
names(keepfrombiglist)[names(keepfrombiglist) == "tpm_CTRL_9"] <- "9_CTRL"
names(keepfrombiglist)[names(keepfrombiglist) == "tpm_KO_8"] <- "8_KO"
names(keepfrombiglist)[names(keepfrombiglist) == "tpm_KO_5"] <- "5_KO"
names(keepfrombiglist)[names(keepfrombiglist) == "tpm_KO_10"] <- "10_KO"
#keep data corresponding to tpm values and gene name
bigdf <- keepfrombiglist[c(2,3,4,5,6,7,14)]
bigdf2 <- keepfrombiglist[c(2,3,4,5,6,7)]
#create a vector of string elements with Gene names
bigvector <- bigdf[,7]
#make this vector of strings the rownames for dataframe without Gene.name column (df2)
rownames(bigdf2) <- bigvector
View(bigdf2)
#make data into a matrix
bigmat_data <- data.matrix(bigdf2[,1:ncol(bigdf2)])
#scalar to help form heatmap
bigdf2 <- scale(bigmat_data)
#heatmap command
heatmap.2(bigmat_data, scale = "none", xlab = "Condition", cexRow = .3, cexCol = 1, 
          ylab = "Gene", keysize = 0.5, lhei=c(1,5), lwid=c(1,3), col = bluered(100), trace = "none", 
          offsetRow = 0, offsetCol = 0, density.info = "none") 

#top 10
top5 <- final[final$p.value< 0.05 & final$FoldChange>4,]
bottom5 <- final[final$p.value< 0.05 & final$FoldChange< -3.92,]
select <- rbind(top5, bottom5)
selectwanno <- merge(select, annotations, by = "target_id")
View(selectwanno)
keep <- selectwanno[-c(7),]

names(keep)[names(keep) == "tpm_CTRL_2"] <- "2_CTRL"
names(keep)[names(keep) == "tpm_CTRL_5"] <- "5_CTRL"
names(keep)[names(keep) == "tpm_CTRL_9"] <- "9_CTRL"
names(keep)[names(keep) == "tpm_KO_8"] <- "8_KO"
names(keep)[names(keep) == "tpm_KO_5"] <- "5_KO"
names(keep)[names(keep) == "tpm_KO_10"] <- "10_KO"

View(keep)
df <- keep[c(2,3,4,5,6,7,14)]
df2 <- keep[c(2, 3, 4,5,6,7)]
View(df2)
avector <- df[,7]

rownames(df2) <- avector
View(df2)

mat_data <- data.matrix(df2[,1:ncol(df2)])
df2 <- scale(mat_data)
heatmap.2(mat_data, scale = "none", xlab = "Condition", cexRow = .8, cexCol = .8, keysize = 1,
          ylab = "Gene", srtCol = 50, col = bluered(100), trace = "none", 
          density.info = "none" )

heatmap.2(mat_data, scale = "none", xlab = "Condition", cexRow = 1, cexCol = 1, 
          ylab = "Gene", keysize = 1, lhei=c(1,5), lwid=c(1,3), col = bluered(100), trace = "none", 
          offsetRow = 0, offsetCol = 0, density.info = "none") 

