#This script produces a set of graphs comparing GC3 and stop codon usage and a table of Spearman's rank coefficients
#Location of file is set
#The folder contains 1 CSV file with the raw data

setwd("File_location")

if (file.exists("GCvAAusage.xls")) {l
file.remove("GCvAAusage.xls")
}


getwd()

list.files()


filelist <- list.files(path = ".", pattern = "\\.csv$")

i <- filelist[1]


data1 <- read.csv(i, header = TRUE, fill = TRUE)

GC3.inR <- unlist(data1$GlobalGC3)

#The columns containing the the stop codon usage are known  

for (j in c(7: 9)) {

aa.perc <- unlist(data1[j])

aa <- colnames(data1)[j]

aa.1 <- gsub("Perc", "", aa)

fig.name <- gsub("Perc", ".pdf", aa)

plot.title <- sprintf("%s versus GC3", aa.1)


#now do the plot

pdf(fig.name)

plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim =c(0,100))

abline(lm(aa.perc~GC3.inR), col = "black")
Rsquare <- lm(aa.perc~GC3.inR)


mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.8)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 3, adj = 0.5, line = -1.3, cex = 0.8)


dev.off()

}



#For all the four length stop codons

for (j in c(22: 33)) {

aa.perc <- unlist(data1[j])

aa <- colnames(data1)[j]

aa.1 <- gsub("Perc", "", aa)

fig.name <- gsub("Perc", ".pdf", aa)

plot.title <- sprintf("%s versus GC3", aa.1)


#now do the plot

pdf(fig.name)

plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim = c(0,100))


abline(lm(aa.perc~GC3.inR), col = "black")
Rsquare <- lm(aa.perc~GC3.inR)


mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.8)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 3, adj = 0.5, line = -1.3, cex = 0.8)


dev.off()


}



#One figure with subplots for the four length codons

pdf("stop_matrix.pdf", width = 21/2.54, height = 27.7/2.54)

par(mfrow=c(4,3), oma = c(4,2,2,2) + 0.1, mar = c(3,2,3,2) + 0.1)


for (j in c(22: 33)) {

aa.perc <- unlist(data1[j])

aa <- colnames(data1)[j]

aa.1 <- gsub("Perc", "", aa)

fig.name <- gsub("Perc", ".pdf", aa)

plot.title <- sprintf("%s versus GC3", aa.1)

#now do the plot


par(pty="s")

plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", mgp = c(2, 1, 0), col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim= c(0,100))

abline(lm(aa.perc~GC3.inR), col = "black")
Rsquare <- lm(aa.perc~GC3.inR)


mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.6)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 4, adj = 0.6, line = -1.3, cex = 0.6)


}

dev.off()


#One figure with all of the 3 and four length stop codons

pdf("all_stop_matrix.pdf", width = 21/2.54, height = 27.7/2.54)

par(mfrow=c(4,3), oma = c(4,2,2,2) + 0.1, mar = c(3,2,3,2) + 0.1)

z= 7
for (j in c(22: 33)) {

while(z<10) {
for (k in c(7: 9)) {

aa.perc <- unlist(data1[k])

aa <- colnames(data1)[k]

aa.1 <- gsub("Perc", "", aa)


plot.title <- sprintf("%s versus GC3", aa.1)


par(pty="s")

plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", mgp = c(2, 1, 0), col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim= c(0,100))



abline(lm(aa.perc~GC3.inR), col = "black")
Rsquare <- lm(aa.perc~GC3.inR)


mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.6)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 4, adj = 0.6, line = -1.3, cex = 0.6)
z = z+1
}

for (j in c(22:33)){
aa.perc <- unlist(data1[j])

aa <- colnames(data1)[j]

aa.1 <- gsub("Perc", "", aa)


plot.title <- sprintf("%s versus GC3", aa.1)



par(pty="s")
plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", mgp = c(2, 1, 0), col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim= c(0,100))


abline(lm(aa.perc~GC3.inR), col = "black")
Rsquare <- lm(aa.perc~GC3.inR)


mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.6)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 4, adj = 0.6, line = -1.3, cex = 0.6)

}
}
dev.off()
}



#Speaman's rank coefficient table

topline <- c("aa", "Rho","slope","P")
out1 <- print(topline, quote=FALSE)
write(out1, file = "GCvStopusage.xls", ncolumns=4, append =TRUE, sep="\t")
z = 7

for (j in c(22: 33)) {

while (z<10) {
for (k in c(7:9)) {

aa.perc <- unlist(data1[k])

aa <- colnames(data1)[k]

aa.1 <- gsub("Perc", "", aa)

cor.gc3.aa <- cor.test(GC3.inR, aa.perc, method = "spearman", exact =FALSE)

rho.gc3.aa <- round(cor.gc3.aa$estimate, digits = 4)
Pval.gc3.aa  <- cor.gc3.aa$p.value
slope<- round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)

out <-print(c(aa.1, rho.gc3.aa, slope, Pval.gc3.aa), quote=FALSE)

write(out, file = "GCvStopusage.xls", ncolumns=4, append = TRUE, sep="\t")
z= z+1
}
}



aa.perc <- unlist(data1[j])

aa <- colnames(data1)[j]

aa.1 <- gsub("Perc", "", aa)

cor.gc3.aa <- cor.test(GC3.inR, aa.perc, method = "spearman", exact = FALSE)

rho.gc3.aa <- round(cor.gc3.aa$estimate, digits = 4)
Pval.gc3.aa  <- cor.gc3.aa$p.value
slope<- round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)

out <-print(c(aa.1, rho.gc3.aa, slope, Pval.gc3.aa), quote=FALSE)

write(out, file = "GCvStopusage.xls", ncolumns=4, append = TRUE, sep="\t")

}








