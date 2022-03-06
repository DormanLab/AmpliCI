#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressMessages(library(dada2,quietly = TRUE))
file_name<-args[1]
output_name<-args[2]


R2.Seq.273.clean <- file.path("./", file_name)
derepRs <- derepFastq(R2.Seq.273.clean, verbose=TRUE)
table_abun<-table(derepRs[3]$map)

data<-as.data.frame(table(table_abun))
colnames(data)<-c("Abun","Freq")
data$Abun<-as.numeric(levels(data$Abun[data$Abun]))


y<-NULL
for(i in 1:length(data$Abun)){
  y[data$Abun[i]]<-data$Freq[i]
}
y[is.na(y)]<-0


write(y, file = output_name,sep = "\n")