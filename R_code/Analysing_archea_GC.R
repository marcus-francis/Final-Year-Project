# a hash means that a line is a comment and R will not read it - but you can
#we will make a file call GCVAAusage.xls - so as to keep things clean we will start by removing any such file if it already exists


#the first thing to do is to tell R where your data file is.  Best is to have any script (like this) in the same folder.  This is easiest done by push and click and dragging files.  You will be shown how. You can always check where the working directory currently is:
setwd("C:/Users/Mole/Documents/Uni work/UNI WORK/Year 4/The project/learningR_2021_1styrs/Archea_CSV")

if (file.exists("GCvAAusage.xls")) {l
file.remove("GCvAAusage.xls")
}


getwd()

#you can also ask what files are in there just to check:

list.files()

# so now you are in the right folder/directory.  Now you need to read in the data so that R knows about it

# if you have multiple data files (lets say they all end .csv) then:

filelist <- list.files(path = ".", pattern = "\\.csv$")

#path= "." says "in this folder", 

#The pattern instruction says don't just look for all files, just the ones that look like this.  Notice how the csv is in quotes.  This is very important. Things in quotes mean literal, things not in quotes are things R knows about
#the $ after csv tells R that the filename ends there. The \\ means a literal dot.  a dot on it own means any charatcer.  

#you now have a variable "filelist" which is a vector containing the name of all your files that end ".csv"

#you can now run through each file and open each in turn and process each in turn.  This is one thing that makes R just amazing - what you can do to one file you can do to many so long as they are the same structure:
#can you work out where the oppening curly brack finishes?  Clue - it is right at the bottom.  Why?

#for (i in filelist) {
i <- filelist[1]

#data2 <- read.table(i, header = TRUE, fill = TRUE)

#if you don't have a list of files just one you can still do this - the list will be just one long.  i takes the value of the name of the first file, then the second file etc
#read.table is a very generic way of opening a data file and is prone to failure as it tries to work out the structure of your data - and often gets its wrong.  Better is to have the file in csv format - each column seperated by a comma, and use:

data1 <- read.csv(i, header = TRUE, fill = TRUE)

#header=TRUE says that the file has a topline which gives the names of the columns - not all input files are like this.  I suspect you can guess what you write if it isn't true.

#we can use either the column name or the column number to read in data.  We have a column called GC3Perc for example:

GC3.inR <- unlist(data1$Global_GC3)

#alternatively if you know the column number, use that. If GCgenomic is in column 3:

GCgen.inR <- unlist(data1$Global_GC)

#in both cases we call the data that has been read in by a name we invent GC3.inR and GCgen.inR in this case.  It is really good idea to name your variables after what they represent.  R really likes "." to break up bits of names
#now we have two variables we can start to do some tests and plots

#if you want to know about the mean of a vector:

#mean.gc3 <- mean(GC3.inR)

#sd.gc3 <- sd(GC3.inR)

#var.gc3 <- sd.gc3^2

#sem.gc3 <- sd.gc3/sqrt(length(GC3.inR))

#here we are making se of R's built in functions mean and sd.  There is no variance function or standard error of the mean so we just calculate these ourselves

#perhaps you want to correlate two variables:

cor.gc3.gcgen <- cor.test(GC3.inR, GCgen.inR, method = "spearman")

#here we are saving the result of R's built in cor.test as a new variable, cor.gc3.gcgen

rho.gc3.gcgen <- round(cor.gc3.gcgen$estimate, digits = 4)
Pval.gc3.gcgen  <- cor.gc3.gcgen$p.value

#the above pulls out the rho value and P value

#and why not make a pdf of a plot

pdf("Figure1.pdf")

plot(GC3.inR, GCgen.inR)

dev.off()

#this last line is like the close file instruction
#this will work but is very basic.  You can do also sorts with plots in R

pdf("Figure1_fancy.pdf")

plot(GC3.inR, GCgen.inR, xlab="GC3%", ylab = "Genomic GC content", col = "blue", cex=0.5, pch=18, main = "GC3% v GC genomic")

#cex gives the size of the points, pch is an option to decide the shape of the points

#this plots a black regression line
abline(lm(GCgen.inR~GC3.inR), col = "black")
Rsquare <- lm(GCgen.inR~GC3.inR)


#and this adds text saying what the slope of the line is (you could add the correlation and P value if you prefer)
mtext(paste0("Gradient = ", round(coef(lm(GCgen.inR~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.8)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 3, adj = 0.5, line = -1.3, cex = 0.8)
#mtext(paste0("P-value =", round(Pval.gc3.gcgen, digits = 4)), side =3, adj = 0.75, line = -1.3, cex = 0.8)


dev.off()



#this is the simple case where there are just two columns you want to look at.  This file contains % usage of all amino acids as well.  How could we look at the correlation of GC3 and each % amino acid usage?  we could do them one by one - which would be fine. 
# an alternative is to use knowledge of header titles and run through each automatically.  

for (j in c(5: 7)) {

aa.perc <- unlist(data1[j])

aa <- colnames(data1)[j]

aa.1 <- gsub("Perc", "", aa)

fig.name <- gsub("Perc", ".pdf", aa)

plot.title <- sprintf("%s versus GC3",aa.1 )

if(j == 5) {
stop.inR <- unlist(data1$TAAPerc)
} else if(j ==6){
stop.inR <- unlist(data1$TGAPerc)
} else {
stop.inR <- unlist(data1$TAGPerc)
}

#now do the plot

pdf(fig.name)

plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim =c(0,100))


abline(lm(aa.perc~GC3.inR), col = "black")
Rsquare <- lm(aa.perc~GC3.inR)

cor.gc3.stop <- cor.test(GC3.inR, stop.inR, method = "spearman", exact = FALSE)
rho.gc3.stop <- round(cor.gc3.stop$estimate, digits = 4)
Pval.gc3.stop  <- cor.gc3.stop$p.value


#and this adds text saying what the slope of the line is (you could add the correlation and P value if you prefer)
mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.8)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 3, adj = 0.5, line = -1.3, cex = 0.8)
#mtext(paste0("P-value =", round(Pval.gc3.stop, digits = 4)), side =3, adj = 0.75, line = -1.3, cex = 0.8)



dev.off()


}

#For all the four length stop codons

for (j in c(20: 31)) {

aa.perc <- unlist(data1[j])

aa <- colnames(data1)[j]

aa.1 <- gsub("Perc", "", aa)

fig.name <- gsub("Perc", ".pdf", aa)

plot.title <- sprintf("%s versus GC3", aa.1)

#if(j == 5) {
#stop.inR <- unlist(data1$TAAPerc)
#} else if(j ==6){
#stop.inR <- unlist(data1$TGAPerc)
#} else {
#stop.inR <- unlist(data1$TAGPerc)
#}

#now do the plot

pdf(fig.name)

plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim = c(0,100))


abline(lm(aa.perc~GC3.inR), col = "black")
Rsquare <- lm(aa.perc~GC3.inR)

#cor.gc3.stop <- cor.test(GC3.inR, stop.inR, method = "spearman", exact = FALSE)
#rho.gc3.stop <- round(cor.gc3.stop$estimate, digits = 4)
#Pval.gc3.stop  <- cor.gc3.stop$p.value


#and this adds text saying what the slope of the line is (you could add the correlation and P value if you prefer)
mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.8)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 3, adj = 0.5, line = -1.3, cex = 0.8)
#mtext(paste0("P-value =", round(Pval.gc3.stop, digits = 4)), side =3, adj = 0.75, line = -1.3, cex = 0.8)



dev.off()


}



#but you might prefer one figure with all the plots as subfigures:

pdf("stop_matrix.pdf", width = 21/2.54, height = 27.7/2.54)

par(mfrow=c(4,3), oma = c(4,2,2,2) + 0.1, mar = c(3,2,3,2) + 0.1)

#to be pretty I'll do the first 12 on one page...
for (j in c(20: 31)) {

aa.perc <- unlist(data1[j])

aa <- colnames(data1)[j]

aa.1 <- gsub("Perc", "", aa)

#fig.name <- gsub("Perc", ".pdf", aa)

plot.title <- sprintf("%s versus GC3", aa.1)

#now do the plot


par(pty="s")
#the above forces the plot to be square
plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", mgp = c(2, 1, 0), col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim= c(0,100))


abline(lm(aa.perc~GC3.inR), col = "black")
Rsquare <- lm(aa.perc~GC3.inR)


mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.6)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 4, adj = 0.6, line = -1.3, cex = 0.6)





}

dev.off()


#but you might prefer one figure with all the plots as subfigures:

pdf("all_stop_matrix.pdf", width = 21/2.54, height = 27.7/2.54)

par(mfrow=c(4,3), oma = c(4,2,2,2) + 0.1, mar = c(3,2,3,2) + 0.1)

#to be pretty I'll do the first 12 on one page...
z= 5
#for (j in c(20: 31)) {

while(z<8) {
for (k in c(5: 7)) {

aa.perc <- unlist(data1[k])

aa <- colnames(data1)[k]

aa.1 <- gsub("Perc", "", aa)

#fig.name <- gsub("Perc", ".pdf", aa)

plot.title <- sprintf("%s versus GC3", aa.1)

#now do the plot


par(pty="s")
#the above forces the plot to be square
plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", mgp = c(2, 1, 0), col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim= c(0,100))


abline(lm(aa.perc~GC3.inR), col = "black")
Rsquare <- lm(aa.perc~GC3.inR)


mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.6)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 4, adj = 0.6, line = -1.3, cex = 0.6)
z = z+1
}

for (j in c(20:31)){
aa.perc <- unlist(data1[j])

aa <- colnames(data1)[j]

aa.1 <- gsub("Perc", "", aa)

#fig.name <- gsub("Perc", ".pdf", aa)

plot.title <- sprintf("%s versus GC3", aa.1)

#now do the plot


par(pty="s")
#the above forces the plot to be square
plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", mgp = c(2, 1, 0), col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim= c(0,100))


abline(lm(aa.perc~GC3.inR), col = "black")
Rsquare <- lm(aa.perc~GC3.inR)


mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.6)
mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 4, adj = 0.6, line = -1.3, cex = 0.6)

}
}
dev.off()



#and what about making a table of the slopes, rho and P values for each amino acid?

topline <- c("aa", "Rho","slope","P")
out1 <- print(topline, quote=FALSE)
write(out1, file = "GCvStopusage.xls", ncolumns=4, append =TRUE, sep="\t")
z = 5

for (j in c(20: 31)) {

while (z<8) {
for (k in c(5:7)) {

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





