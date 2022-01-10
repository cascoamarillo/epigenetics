
file=read.table("~/Documents/ArkLab/Methylome/motifs/base-modification/Av11.original_PB.euk_15cells_m4C-m6A.csv", stringsAsFactors=FALSE, header=TRUE, sep=",")

file=read.table("~/Documents/ArkLab/Methylome/motifs/base-modification/Av11.original_PB.euk_15cells_m6A-10x.csv", stringsAsFactors=FALSE, header=TRUE, sep=",")

file=read.table("~/Documents/ArkLab/Methylome/motifs/base-modification/Av11.original_PB.euk_15cells_m4C-m6A-10x_frac0.1.csv", stringsAsFactors=FALSE, header=TRUE, sep=",")


head(file)

#####ggplot2
library(ggplot2)

ggplot(file, aes(x=frac, color=plus.one)) +
  geom_density()+
  xlim(0, 50)

p <- ggplot(file, aes(x=coverage, color=modification)) +
  geom_density()+
  xlim(0, 50)
p
p + scale_colour_manual( values = c("red","blue"))
p + scale_colour_manual( values = c("blue","red")) + geom_vline(xintercept = 10, linetype="dotted", 
                 color = "black", size=1)

aggregate(file$V10, by=list(Category=file$V3), FUN=sum)

library(dplyr)
count(file, plus.one, plus.two)

ggplot(data=file, aes(x = plus.one, y = frac, fill = plus.two)) +
  geom_bar(stat = "identity")

p <-ggplot(data=file, aes(x = fraction, fill = modification)) +
  geom_histogram(stat = "bin", binwidth = 0.2, position=position_dodge())+
  theme_minimal()
p + scale_fill_manual(values=c('blue','red')) + scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1))

ggplot(data=file, aes(x = frac, y = coverage, fill = modification)) +
  geom_boxplot()

ggplot(file, aes(x = plus.one, y = frac, fill = plus.two )) + 
  geom_boxplot()

ggplot(data=file, aes(x = frac, fill = minus.two)) +
  geom_bar(stat = "bin", binwidth = 0.2, position=position_dodge())
