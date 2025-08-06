brary(DESeq2)

# 读取输入文件
counts <- read.table(snakemake@input$counts, header=TRUE, row.names=1, comment.char="#")
samples <- read.csv(snakemake@input$samples, row.names="sample")

# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
			        countData = counts[, rownames(samples)],
				  colData = samples,
				  design = ~ group
				  )

# 差异分析
dds <- DESeq(dds)
res <- results(dds, contrast=c("group", "treatment", "control"))

# 输出结果
write.csv(as.data.frame(res), file=snakemake@output[[1]])
