library(DESeq2)
args = commandArgs(trailingOnly=TRUE)
print(args)

data_path <- args[1]
meta_path <- args[2]
condition <- args[3]
outpath <- args[4]

## from https://github.com/hbc/NGS_Data_Analysis_Course/blob/master/sessionIII/lessons/03_DEG_getting_started.md

## Load in data
data <- read.csv(data_path, header=T, row.names=1) 

meta <- read.csv(meta_path, header=T, row.names=1)

### Check classes of the data we just brought in
class(data)
class(meta)

### Check that sample names match in both files
all(names(data) %in% rownames(meta))
all(names(data) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = as.formula(paste("~", condition)))


## from https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html

## Run analysis
dds <- DESeq(dds)

## Save normalized counts
normalized_counts <- counts(dds, normalized=TRUE)

# Convert to a data frame for easier manipulation if needed
normalized_counts_df <- as.data.frame(normalized_counts)

# Write the normalized counts to a CSV file
write.csv(normalized_counts_df, paste0(outpath,"_deseq2_normalized_counts.csv"))

## from https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
res <- results(dds)
res <- res[order(res$padj), ]

write.csv(res, paste0(outpath, "_", condition, "_deseq2_results.csv"))
