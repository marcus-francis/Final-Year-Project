#This script produces a set of graphs comparing the overall usage of each stop codon
#Location of file is set
#The folder contains 1 CSV file with the raw data

setwd("File Location")

if (file.exists("GCvAAusage.xls")) {l
file.remove("GCvAAusage.xls")
}


getwd()

list.files()


filelist <- list.files(path = ".", pattern = "\\.csv$")


i <- filelist[1]

data1 <- read.csv(i, header = TRUE, fill = TRUE)

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



