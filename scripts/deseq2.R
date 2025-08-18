library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(fgsea)
library(msigdbr)
library(tidyverse)
library(ashr)
options(timeout=300) 
# 读取snakemake传入的config路径
config <- jsonlite::fromJSON(snakemake@input$config)

# 使用绝对路径处理所有文件
gsea_path <- config$gsea
base_path <- config$path_r
setwd(base_path)
# 读取输入文件
count_data <- read.table(snakemake@input$counts, header=TRUE, row.names=1, comment.char="#")
samples <- read.csv(snakemake@input$samples,  header=TRUE)


# 清理count数据（去掉注释列）
count_data <- count_data[, 6:ncol(count_data)]
colnames(count_data) <- gsub("results.align.|\\.bam", "", colnames(count_data))

# 确保样本顺序一致
samples <- samples[match(colnames(count_data), samples$sample), ]
rownames(samples) <- samples$sample

# 2. 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = round(count_data),
  colData = samples,
  design = ~ group
)

# 3. 预处理和PCA分析
vsd <- vst(dds, blind=FALSE)

# PCA 图
pca_data <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  geom_text_repel(aes(label=name), show.legend=FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA Plot")
  
ggsave("results/deseq2/PCA_plot.pdf", pca_plot, width=8, height=6)

# 4. 样本聚类热图
sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)


cluster_plot <- pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  main = "Sample Clustering Heatmap")


pdf("results/deseq2/cluster_heatmap.pdf", width=8, height=6)
print(cluster_plot)
dev.off()

# 5. 差异表达分析
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("group", "treatment", "control"))
res <- lfcShrink(dds, contrast=c("group", "treatment", "control"), res=res, type="ashr")

# 获取原始counts（未标准化）
raw_counts <- counts(dds, normalized=FALSE)

# 确保基因顺序一致
raw_counts <- raw_counts[match(rownames(res), rownames(raw_counts)), ]

# 合并结果
res_df <- cbind(as.data.frame(res), raw_counts)

# 保存完整结果
write.csv(res_df, "results/deseq2/DE_results.csv")

# 筛选差异基因（同样保留原始counts）
res_sig <- subset(res_df, abs(log2FoldChange) > 1 & padj < 0.05)
write.csv(res_sig, "results/deseq2/DE_results_sig.csv")

res_sig_pvalue <- subset(res_df, abs(log2FoldChange) > 1 & pvalue < 0.05)
write.csv(res_sig_pvalue, "results/deseq2/DE_results_sig_pvalue.csv")

# 6. 火山图
volcano_data <- as.data.frame(res)
volcano_data$gene <- rownames(volcano_data)
volcano_data$diffexpressed <- ifelse(
  volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) > 1,
  ifelse(volcano_data$log2FoldChange > 1, "Up", "Down"),
  "No"
)

volcano_plot <- ggplot(volcano_data, aes(log2FoldChange, -log10(pvalue), color=diffexpressed)) +
  geom_point(aes(color=diffexpressed), alpha=0.5) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_bw() +
  labs(color="Expression") +
  ggtitle("Volcano Plot") +
  theme(legend.position="top")

ggsave("results/deseq2/volcano_plot.pdf", volcano_plot, width=8, height=6)


# 7. 差异基因热图
sig_genes <- rownames(res_sig)
norm_counts <- counts(dds, normalized=TRUE)
sig_norm_counts <- norm_counts[sig_genes, ]

heatmap_plot <- pheatmap(
  log2(sig_norm_counts + 1),
  scale="row",
  clustering_distance_rows="euclidean",
  clustering_distance_cols="euclidean",
  show_rownames=FALSE,
  annotation_col=data.frame(group=samples$group, row.names=colnames(sig_norm_counts)),
  main="Differentially Expressed Genes Heatmap"
)

pdf("results/deseq2/DE_heatmap.pdf", width=8, height=10)
print(heatmap_plot)
dev.off()

# 8. GSEA分析
res$ENSEMBL=str_split(rownames(res),'[.]',simplify = T)[,1] #把行名.和以后得字符去掉形成为名为ENSEMBL的一列
dat1=res[,c(6,2)] #取出res中ENSGID和logFC两列
gene=bitr(dat1[,1],fromType="ENSEMBL",toType="SYMBOL",OrgDb="org.Hs.eg.db")#开始ID转换，会有丢失，gene中为ENSGID何SYMBOL两列

gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)##去重

data_all <-merge(as.data.frame(dat1),gene,by="ENSEMBL",all=F)
data_all_sort <- data_all %>%
 arrange(desc(log2FoldChange))
head(data_all_sort)
data_all_sort=na.omit(data_all_sort)
geneList = data_all_sort$log2FoldChange #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$SYMBOL

kegg_gmt <- read.gmt(gsea_path) #读gmt文件
gsea <- GSEA(geneList,
            TERM2GENE = kegg_gmt,
			pvalueCutoff=0.05) #GSEA分析
			
pdf("results/deseq2/GSEA_ridgeplot.pdf", width=8, height=10)			
ridgeplot(gsea,10)	
dev.off()

write.csv(gsea, "results/deseq2/gsea_results_full.csv")


#9.kegg和GO
# 定义绘图函数
plot_enrich <- function(enrich_res, title) {
  if (nrow(enrich_res) == 0) return(NULL)  # 处理空结果
  
  # 提取前10个显著条目
  df <- head(enrich_res[order(enrich_res$p.adjust), ], 10)
  
  ggplot(df, aes(x=-log10(p.adjust), 
             y=reorder(Description, -log10(p.adjust)))) +
    geom_bar(stat="identity", fill="steelblue", width=0.7) +
    theme_bw() +
    labs(x="-log10(adjusted p-value)", y="", title=title) +
    theme(axis.text.y = element_text(size=10))
}
# 去除ENSEMBL ID的版本号（关键步骤！）
ensembl_ids <- gsub("\\..*", "", rownames(res_sig))

# 转换为ENTREZID
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)


# 命名向量：名称=无版本号的ENSEMBL ID，值=ENTREZID
names(entrez_ids) <- gsub("\\..*", "", names(entrez_ids))  # 确保名称无版本号
entrez_ids <- na.omit(entrez_ids)

# 定义阈值
fc_threshold <- 1      # |log2FC| > 1
p_threshold <- 0.05    # padj < 0.05

# 上调基因（log2FC > 1且padj < 0.05）
up_genes <- rownames(res_sig)[
  res_sig$log2FoldChange > fc_threshold & 
  res_sig$padj < p_threshold
]
up_genes <- gsub("\\..*", "", up_genes)  # 去除版本号
up_entrez <- entrez_ids[up_genes]        # 匹配ENTREZID

# 下调基因（log2FC < -1且padj < 0.05）
down_genes <- rownames(res_sig)[
  res_sig$log2FoldChange < -fc_threshold & 
  res_sig$padj < p_threshold
]
down_genes <- gsub("\\..*", "", down_genes)  # 去除版本号
down_entrez <- entrez_ids[down_genes]        # 匹配ENTREZID

# GO分析
go_up <- enrichGO(
  gene = up_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  readable = TRUE
)

# KEGG分析
kegg_up <- enrichKEGG(
  gene = up_entrez,
  organism = "hsa",
  pAdjustMethod = "BH"
)

# GO分析
go_down <- enrichGO(
  gene = down_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  readable = TRUE
)

# KEGG分析
kegg_down <- enrichKEGG(
  gene = down_entrez,
  organism = "hsa",
  pAdjustMethod = "BH"
)


# 保存结果
write.csv(as.data.frame(go_up), "results/deseq2/GO_up.csv")
write.csv(as.data.frame(go_down), "results/deseq2/GO_down.csv")

write.csv(as.data.frame(kegg_up), "results/deseq2/KEGG_up.csv")
write.csv(as.data.frame(kegg_down), "results/deseq2/KEGG_down.csv")

# 绘图
pdf("results/deseq2/GO_enrichment_plots_up.pdf", width=10, height=6)
if (length(up_entrez) > 0) print(plot_enrich(go_up, "Upregulated GO Terms"))
dev.off()

pdf("results/deseq2/GO_enrichment_plots_down.pdf", width=10, height=6)
if (length(down_entrez) > 0) print(plot_enrich(go_down, "downregulated GO Terms"))
dev.off()

pdf("results/deseq2/KEGGE_enrichment_plots_up.pdf", width=10, height=6)
if (length(up_entrez) > 0) print(plot_enrich(kegg_up, "Upregulated KEGG Terms"))
dev.off()

pdf("results/deseq2/KEGG_enrichment_plots_down.pdf", width=10, height=6)
if (length(down_entrez) > 0) print(plot_enrich(kegg_down, "downregulated KEGG Terms"))
dev.off()



#p-value
# 上调基因（log2FC > 1且pvalue < 0.05）
up_genes <- rownames(res_sig_pvalue)[
  res_sig_pvalue$log2FoldChange > fc_threshold & 
  res_sig_pvalue$pvalue < p_threshold
]
up_genes <- gsub("\\..*", "", up_genes)  # 去除版本号
up_entrez <- entrez_ids[up_genes]        # 匹配ENTREZID

# 下调基因（log2FC < -1且padj < 0.05）
down_genes <- rownames(res_sig_pvalue)[
  res_sig_pvalue$log2FoldChange < -fc_threshold & 
  res_sig_pvalue$pvalue < p_threshold
]
down_genes <- gsub("\\..*", "", down_genes)  # 去除版本号
down_entrez <- entrez_ids[down_genes]        # 匹配ENTREZID

# GO分析
go_up <- enrichGO(
  gene = up_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  readable = TRUE
)

# KEGG分析
kegg_up <- enrichKEGG(
  gene = up_entrez,
  organism = "hsa",
  pAdjustMethod = "BH"
)

# GO分析
go_down <- enrichGO(
  gene = down_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  readable = TRUE
)

# KEGG分析
kegg_down <- enrichKEGG(
  gene = down_entrez,
  organism = "hsa",
  pAdjustMethod = "BH"
)


# 保存结果
write.csv(as.data.frame(go_up), "results/deseq2/GO_up_pvalue.csv")
write.csv(as.data.frame(go_down), "results/deseq2/GO_down_pvalue.csv")

write.csv(as.data.frame(kegg_up), "results/deseq2/KEGG_up_pvalue.csv")
write.csv(as.data.frame(kegg_down), "results/deseq2/KEGG_down_pvalue.csv")

# 绘图
pdf("results/deseq2/GO_enrichment_plots_up_pvalue.pdf", width=10, height=6)
if (length(up_entrez) > 0) print(plot_enrich(go_up, "Upregulated GO Terms"))
dev.off()

pdf("results/deseq2/GO_enrichment_plots_down_pvalue.pdf", width=10, height=6)
if (length(down_entrez) > 0) print(plot_enrich(go_down, "downregulated GO Terms"))
dev.off()

pdf("results/deseq2/KEGGE_enrichment_plots_up_pvalue.pdf", width=10, height=6)
if (length(up_entrez) > 0) print(plot_enrich(kegg_up, "Upregulated KEGG Terms"))
dev.off()

pdf("results/deseq2/KEGG_enrichment_plots_down_pvalue.pdf", width=10, height=6)
if (length(down_entrez) > 0) print(plot_enrich(kegg_down, "downregulated KEGG Terms"))
dev.off()