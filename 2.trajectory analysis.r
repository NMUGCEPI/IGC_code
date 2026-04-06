##Cell differentiation trajectory was inferred using Monocle.
library(Seurat)
library(monocle)
library(dplyr)
exp<-as(as.matrix(epi_IGC@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(exp),gene_short_name=rownames(exp))
rownames(feature_ann)<-rownames(exp)
fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-epi_IGC@meta.data
rownames(sample_ann)<-colnames(exp)
pd<-new("AnnotatedDataFrame", data =sample_ann)
obj<-newCellDataSet(exp,phenoData =pd,featureData =fd,expressionFamily=negbinomial.size())
obj <- estimateSizeFactors(obj)
obj <- estimateDispersions(obj)
obj <- detectGenes(obj, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(obj),
    num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(obj[expressed_genes,],
              fullModelFormulaStr = "~type",core=1)
top800 <- diff_test_res  %>% top_n(n = -800, wt = qval)
ordering_genes <- top800$gene_id
monocle<- setOrderingFilter(obj, ordering_genes)
monocle<- reduceDimension(monocle, max_components=2, reduction_method= 'DDRTree',norm_method = 'log')
monocle<- orderCells(monocle)