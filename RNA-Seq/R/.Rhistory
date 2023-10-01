read.delim('C:/Sharing Files/RNA-Seq/readcounts_control_1.txt',header = TRUE,sep = '\t') -> control1
read.delim('C:/Sharing Files/RNA-Seq/readcounts_control_2.txt',header = TRUE,sep = '\t') -> control2
read.delim('C:/Sharing Files/RNA-Seq/readcounts_e2_1.txt',header = TRUE,sep = '\t') -> e21
read.delim('C:/Sharing Files/RNA-Seq/readcounts_e2_2.txt',header = TRUE,sep = '\t') -> e22
rep1 <- merge(e21,control1)
rep2 <- merge(e22,control2)
rep2_log <- data.frame(Genid = rep2$Geneid, control = log(rep2$sorted_MCF7_control_rep2.bam), e2 = log(rep2$sorted_MCF7_E2_rep2.bam))
rep1_log <- data.frame(Genid = rep1$Geneid, control = log(rep1$sorted_MCF7_control_rep1.bam), e2 = log(rep1$sorted_MCF7_E2_rep1.bam))
plot(rep1$sorted_MCF7_E2_rep1.bam,rep1$sorted_MCF7_control_rep1.bam,pch = 20, col = "blue", main = "Replikat 1",xlab = "E2",ylab = "Control", log = "xy")
plot(rep2$sorted_MCF7_E2_rep2.bam,rep2$sorted_MCF7_control_rep2.bam,main = "Replikat 2", xlab ="E2",ylab="Control",pch=20,col ="red",log="xy")
library(DESeq2)
library(ggplot2)
read.delim('C:/Sharing Files/RNA-Seq/R/coldata',header = TRUE, sep = " ") -> coldata
View(control1)
countdata <- merge(control1,control2,e21,e22)
View(coldata)
countdata <- c(control1,control2)
countdata <- merge(control1,control2)
countdata <- merge(countdata,e21)
countdata <- merge(countdata,e22)
View(countdata)
colnames(countdata) <- c("Gen_id","control_1","control_2","e2_1","e2_2")
rownames(countdata) <- countdata$Gen_id
countdata$Gen_id <- NULL
rownames(coldata) %in% colnames(coldata)
coldata$Name <- NULL
coldata[order(coldata$Group),] -> coldata
colnames(coldata) <- "group"
data.frame(coldata)
data.frame(Group = coldata) -> coldata
colname(countdata) -> rownames(coldata)
colnames(countdata) -> rownames(coldata)
library(DESeq2)
library(ggplot2)
dds <- DESeqDataSetFromMatrix(countdata,coldata,design=~Group)
dds <- DESeq(dds)
res <- results(dds)
head(results(dds,tidy = TRUE))
ressig <- res[which(res$padj < 0.1),]
write.table(as.data.frame(ressig[order(-ressig$log2FoldChange),]),file="hoch- und runterreguliert.txt")
summary(res)
