help("rowttests")
help(library)
help(factor)
help(as.matrix)
help(merge)
help(hist)
help(rnorm)
help(boxplot)
help(qqnorm)

library(genefilter)

label <- factor(c(0,0,0,0,1,1,1,1))
ttests <- rowttests(as.matrix(Pmlogexdata),label)
ttests
Pmlogexdata.stats <- merge(ttests, Pmlogexdata, by = 0)
head(Pmlogexdata.stats)
Pmlogexdata.stats
hist(Pmlogexdata.stats[,3])

randomnormaldist <- rnorm(12625, mean = 7.868)
hist(randomnormaldist)
boxplot(randomnormaldist)

Exdata
hist(Exdata$Patient.3.Exercise.8.pm)
hist(Logexdata$Patient.3.Exercise.8.pm)

qqnorm(randomnormaldist); qqline(randomnormaldist)
qqnorm(Exdata$Patient.3.Exercise.8.pm); qqline(Exdata$Patient.3.Exercise.8.pm)
qqnorm(Logexdata$Patient.3.Exercise.8.pm); qqline(Logexdata$Patient.3.Exercise.8.pm)
