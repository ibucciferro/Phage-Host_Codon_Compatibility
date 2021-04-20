# read in the file (phage host codon correlation) to a data table (and remove the first four lines)
datatable <- read.table(file = 'phage_host_codon_correlation.txt', header = FALSE, sep = '\t', skip = 4)

#rename the columns in the data table
colnames(datatable) <- c("Gene", "Corr")

#use ggplot to create the box plot
library(ggplot2)
dev <- ggplot(datatable, aes(y = Corr, x = "")) + geom_boxplot() + coord_flip() + ylab("Codon Correlation") + xlab('')

#save the file as a PNG
dev.print(file="corrPlot.png", device=png, width=800)
