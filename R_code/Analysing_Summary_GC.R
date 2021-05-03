# a hash means that a line is a comment and R will not read it - but you can
#we will make a file call GCVAAusage.xls - so as to keep things clean we will start by removing any such file if it already exists


#the first thing to do is to tell R where your data file is.  Best is to have any script (like this) in the same folder.  This is easiest done by push and click and dragging files.  You will be shown how. You can always check where the working directory currently is:
setwd("C:/Users/Mole/Documents/Uni work/UNI WORK/Year 4/The project/learningR_2021_1styrs/Summary_CSV")

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

Names.inR <- unlist(data1$Names)

Stops <- matrix(c(data1$TAAPerc, data1$TGAPerc, data1$TAGPerc), nrow = 3, ncol = 3, byrow = TRUE)
colours <- c("green","blue","red")
Codons <- c("TAA","TGA","TAG")

pdf(file = "Stop_Codons.pdf")
barplot(Stops,names.arg=Names.inR, main="Stop Codon Usage", xlab="Domain", ylab = 'Percentage Stop Codon Usage', col = colours)
#legend('topright',  Codons, cex = 1.3, fill = colours,xpd=TRUE, inset=c(0,-0.15))
dev.off()

Stops <- matrix(c(data1$TAAAPerc, data1$TAATPerc, data1$TAAGPerc, data1$TAACPerc), nrow = 4, ncol = 3, byrow = TRUE)
colours <- c("green","blue","red", "cyan")
Codons <- c("TAAA","TAAT","TAAG", "TAAC")


pdf(file = "TAA_Codons.pdf")
barplot(Stops,names.arg=Names.inR, main="Stop Codon Usage", xlab="Domain", ylab = 'Percentage Stop Codon Usage', col = colours)
#legend('topright',  Codons, cex = 1.3, fill = colours,xpd=TRUE, inset=c(0,-0.15))
dev.off()


Stops <- matrix(c(data1$TGAAPerc, data1$TGATPerc, data1$TGAGPerc, data1$TGACPerc), nrow = 4, ncol = 3, byrow = TRUE)
colours <- c("green","blue","red", "cyan")
Codons <- c("TGAA","TGAT","TGAG", "TGAC")


pdf(file = "TGA_Codons.pdf")
barplot(Stops,names.arg=Names.inR, main="Stop Codon Usage", xlab="Domain", ylab = 'Percentage Stop Codon Usage', col = colours)
#legend('topright',  Codons, cex = 1.3, fill = colours,xpd=TRUE, inset=c(0,-0.15))
dev.off()




Stops <- matrix(c(data1$TAGAPerc, data1$TAGTPerc, data1$TAGGPerc, data1$TAGCPerc), nrow = 4, ncol = 3, byrow = TRUE)
colours <- c("green","blue","red", "cyan")
Codons <- c("TAGA","TAGT","TAGG", "TAGC")


pdf(file = "TAG_Codons.pdf")
barplot(Stops,names.arg=Names.inR, main="Stop Codon Usage", xlab="Domain", ylab = 'Percentage Stop Codon Usage', col = colours)
#legend('topright',  Codons, cex = 1.3, fill = colours,xpd=TRUE, inset=c(0,-0.15))
dev.off()


#but you might prefer one figure with all the plots as subfigures:

pdf("stop_matrix.pdf", width = 21/2.54, height = 27.7/2.54)

par(mfrow=c(4,3), oma = c(4,2,2,2) + 0.1, mar = c(3,2,3,2) + 0.1)

#to be pretty I'll do the first 12 on one page...
j =1
while (j <5)) {

aa.perc <- unlist(data1[j])

aa <- colnames(data1)[j]

aa.1 <- gsub("Perc", "", aa)

fig.name <- gsub("Perc", ".pdf", aa)

plot.title <- sprintf("GC3 versus %s", aa.1)

#now do the plot


par(pty="s")
#the above forces the plot to be square
plot(GC3.inR, aa.perc, xlab="GC3%", ylab = "% use of stop codon", mgp = c(2, 1, 0), col = "blue", cex=0.5, pch=18, main = plot.title, xlim = c(0,100), ylim= c(0,100))
lines(GC3.inR, aa.perc, col = "black")


#abline(lm(aa.perc~GC3.inR), col = "black")
#Rsquare <- lm(aa.perc~GC3.inR)


#mtext(paste0("Gradient = ", round(coef(lm(aa.perc~GC3.inR))[2], digits = 4)), side = 3, adj = 0.05, line = -1.3, cex = 0.6)
#mtext(paste0("R-Squared = ", round(summary(Rsquare)$r.squared, digits = 4)), side = 4, adj = 0.6, line = -1.3, cex = 0.6)





}

dev.off()


