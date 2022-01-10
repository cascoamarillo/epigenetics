
file=read.table("~/Documents/ArkLab/Methylome/collinearity//As11.collinearity.PBmet.rnaseq.csv", stringsAsFactors=FALSE, header=TRUE, sep=",")

head(file)


######
log2.rpkm<- log2(file$chr_rpkm   + 1)
log2.rpkm.1<- log2(file$chr1_rpkm   + 1)
log2.rpkm.2<- log2(file$chr2_rpkm   + 1)
log2.rpkm.diff<- log2(abs(file$chr1_rpkm - file$chr2_rpkm))
chr_PB<- file$chr_4mC + file$chr_6mA
chr1_PB<- file$chr1_4mC + file$chr1_6mA
chr2_PB<- file$chr2_4mC + file$chr2_6mA
sqrt.chr_PB<- sqrt(abs(chr1_PB - chr2_PB))

#####ggplot2
library(ggplot2)

#basic boxplot
ggplot(file, aes(x=chr1_PBmet, y=chr2_PBmet,shape=collinearity, color=collinearity)) +
  geom_point(size=1) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(file, aes(x=log2.rpkm.1, y=log2.rpkm.2,shape=collinearity, color=collinearity)) +
  geom_point(size=0.5)

ggplot(file, aes(x=log2.rpkm.1, y=log2.rpkm.2)) + geom_point(size=1, shape=23)

# Add the regression line
geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95)

ggplot(file, aes(x=log2.rpkm.1, y=log2.rpkm.2)) + 
  geom_point()+
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95)

ggplot(file, aes(x=log2.rpkm.1, y=log2.rpkm.2,shape=collinearity, color=collinearity))+
  geom_point(size=0.1) +
  geom_smooth(method=lm)


ggplot(file, aes(x = sqrt.chr_PB, y = log2.rpkm.diff, color = collinearity)) + 
  geom_point(size=0.1)

# Create a Multiple Boxplot Importing
ggplot(file, aes(x = chr_PBmet, y = log2.rpkm.diff, fill = collinearity)) + 
  geom_boxplot()

ggplot(file, aes(x = TE, y = log2.rpkm.diff, fill = collinearity)) + 
  geom_boxplot()

ggplot(file, aes(x = chr_TE, y = log2.rpkm, fill = chr_PBmet)) + 
  geom_boxplot()

p10 <- ggplot(file, aes(x = log2.rpkm.1, y = log2.rpkm.2)) + geom_boxplot()
p10

##labels
p10 <- p10 + scale_x_discrete(name = "H3K methylation") +  scale_y_continuous(name = "Mean log2 RPKM")
p10

##add title
p10 <- p10 + ggtitle("Boxplot of Mean log2 RPKM")
p10
png(filename = "lg2RPKM_H3K-met.png")
dev.off()

#  Add p-value
#install.packages("ggpubr")
library(ggpubr)
pH3K + stat_summary()
pH3K + stat_compare_means()
pH3K + stat_compare_means(method = "t.test")

library(plyr)
ddply(file, ~ IP4mC, wilcox.test, log2.rpkm.avg = log2.rpkm.avg)
ddply(file,"IP4mC",
      function(x) {
        w <- wilcox.test(log2.rpkm.avg~IP4mC,data=x)
        with(w,data.frame(statistic,p.value))
      })

# Create a Multiple Boxplot Importing
pH3K <- ggplot(file, aes(x = file$H3K4.9.27, y = log2.rpkm.avg, fill = AI)) + 
  geom_boxplot()
pH3K  <- pH3K  + scale_x_discrete(name = "H3K methylation") +  scale_y_continuous(name = "Mean log2 RPKM")
pH3K
pIP <- ggplot(file, aes(x = file$IP, y = log2.rpkm.avg, fill = H3K4.9.27)) + 
  geom_boxplot()
IP  <- pIP  + scale_x_discrete(name = "IP methylation") +  scale_y_continuous(name = "Mean log2 RPKM")

png(filename = "lg2RPKM_H3K_AI.png")
IP
dev.off()

# Create a Boxplot Importing
ggplot(file, aes(x = file$IP, y = log2.rpkm.avg, fill = H3K4.9.27)) + 
  geom_boxplot() +
  facet_wrap(~ H3K4.9.27, scale = "free")


####compare
head(diamonds)
head(diam)
diam <- diamonds[diamonds$cut=="Fair"|diamonds$cut=="Ideal",]
boxplots <- ggplot(diam, aes(x=cut, price)) + geom_boxplot(aes(fill=cut)) + facet_wrap(~ color)
print(boxplots)
head(diam)
library(plyr)
ddply(diam,"color",
      function(x) {
        w <- wilcox.test(price~cut,data=x)
        with(w,data.frame(statistic,p.value))
      })



#######BOOTSTRAP
####
library(boot)

# function to obtain R-Squared from the data
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
}
# bootstrapping with 1000 replications
head(file)
results <- boot(data=file, statistic=rsq,
                R=100, formula= AvRNA.rpkmAVG ~ IP)
# view results
results
plot(results)
# get 95% confidence interval
boot.ci(results, type="bca")

p11 <- ggplot(results, aes(x = file$IP, y = AvRNA.rpkmAVG)) + geom_boxplot()
p11
####
# Bootstrap 95% CI for R-Squared
library(boot)
# function to obtain R-Squared from the data
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
}
# bootstrapping with 1000 replications
head(mtcars)
results <- boot(data=mtcars, statistic=rsq,
                R=1000, formula=mpg~wt+disp)
head(results)

# view results
results
plot(results)

# get 95% confidence interval
boot.ci(results, type="bca")