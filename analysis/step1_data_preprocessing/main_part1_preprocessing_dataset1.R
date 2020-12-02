# read data
annotations = read.delim("data/GCF_009858895.2_ASM985889v3_genomic.gff",header=FALSE,comment.char="#",check.names=FALSE,stringsAsFactors=FALSE)
human_cu = read.table("data/HumanCodonUsage.txt",header=FALSE,check.names=FALSE,stringsAsFactors=FALSE)
sars_cov2_fasta = read.table("data/SARS-CoV2.fasta",header=FALSE,sep="\n",check.names=FALSE,stringsAsFactors=FALSE)
sars_cov2_fasta = toupper(strsplit(paste0(sars_cov2_fasta[-1,],collapse=""),split="")[[1]])
load(file="final_data/PRJEB37886/processed_variants.RData")

# consider processed variants
variants_data = NULL
for(i in 1:nrow(processed_variants)) {
    for(j in 1:ncol(processed_variants)) {
        if(!is.na(processed_variants[i,j])&&processed_variants[i,j]>0) {
            curr_split = strsplit(colnames(processed_variants)[j],"_")[[1]]
            variants_data = rbind(variants_data,c(curr_split,rownames(processed_variants)[i],processed_variants[i,j]))
        }
    }
    cat(i/nrow(processed_variants),"\n")
}
rownames(variants_data) = 1:nrow(variants_data)
colnames(variants_data) = c("POSITION","REF","ALT","PATIENT","FREQUENCY")

# save results
save(annotations,file="results/dataset1/annotations.RData")
save(human_cu,file="results/dataset1/human_cu.RData")
save(sars_cov2_fasta,file="results/dataset1/sars_cov2_fasta.RData")
save(variants_data,file="results/dataset1/variants_data.RData")
