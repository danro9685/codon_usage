# read data
annotations = read.delim("external_data/GCF_009858895.2_ASM985889v3_genomic.gff",header=FALSE,comment.char="#",check.names=FALSE,stringsAsFactors=FALSE)
human_cu = read.table("external_data/HumanCodonUsage.txt",header=FALSE,check.names=FALSE,stringsAsFactors=FALSE)
human_cu_values = read.table("external_data/HumanCodonUsage_Hive.txt",,sep="\t",header=FALSE,check.names=FALSE,stringsAsFactors=FALSE)
for(i in 1:nrow(human_cu)) {
    curr_value = NA
    if(length(grep(human_cu[i,1],human_cu_values[,1]))==1) {
        curr_value = human_cu_values[grep(human_cu[i,1],human_cu_values[,1]),2]/100
    }
    if(length(grep(human_cu[i,1],human_cu_values[,4]))==1) {
        curr_value = human_cu_values[grep(human_cu[i,1],human_cu_values[,4]),5]/100
    }
    if(length(grep(human_cu[i,1],human_cu_values[,7]))==1) {
        curr_value = human_cu_values[grep(human_cu[i,1],human_cu_values[,7]),8]/100
    }
    if(length(grep(human_cu[i,1],human_cu_values[,10]))==1) {
        curr_value = human_cu_values[grep(human_cu[i,1],human_cu_values[,10]),11]/100
    }
    human_cu[i,3] = curr_value
}
sars_cov2_fasta = read.table("external_data/SARS-CoV2.fasta",header=FALSE,sep="\n",check.names=FALSE,stringsAsFactors=FALSE)
sars_cov2_fasta = toupper(strsplit(paste0(sars_cov2_fasta[-1,],collapse=""),split="")[[1]])
load(file="final_data/PRJEB39849/processed_variants.RData")
data1 = processed_variants
load(file="final_data/PRJNA614995/processed_variants.RData")
data2 = processed_variants
load(file="final_data/PRJNA625551/processed_variants.RData")
data3 = processed_variants
load(file="final_data/PRJNA633948/processed_variants.RData")
data4 = processed_variants
load(file="final_data/PRJNA636748/processed_variants.RData")
data5 = processed_variants
load(file="final_data/PRJNA639066/processed_variants.RData")
data6 = processed_variants
load(file="final_data/PRJNA645906/processed_variants.RData")
data7 = processed_variants

# consider processed variants
variants_data = NULL
processed_variants = data1
for(i in 1:nrow(processed_variants)) {
    for(j in 1:ncol(processed_variants)) {
        if(!is.na(processed_variants[i,j])&&processed_variants[i,j]>0) {
            curr_split = strsplit(colnames(processed_variants)[j],"_")[[1]]
            variants_data = rbind(variants_data,c(curr_split,rownames(processed_variants)[i],processed_variants[i,j]))
        }
    }
    cat(i/nrow(processed_variants),"\n")
}
processed_variants = data2
for(i in 1:nrow(processed_variants)) {
    for(j in 1:ncol(processed_variants)) {
        if(!is.na(processed_variants[i,j])&&processed_variants[i,j]>0) {
            curr_split = strsplit(colnames(processed_variants)[j],"_")[[1]]
            variants_data = rbind(variants_data,c(curr_split,rownames(processed_variants)[i],processed_variants[i,j]))
        }
    }
    cat(i/nrow(processed_variants),"\n")
}
processed_variants = data3
for(i in 1:nrow(processed_variants)) {
    for(j in 1:ncol(processed_variants)) {
        if(!is.na(processed_variants[i,j])&&processed_variants[i,j]>0) {
            curr_split = strsplit(colnames(processed_variants)[j],"_")[[1]]
            variants_data = rbind(variants_data,c(curr_split,rownames(processed_variants)[i],processed_variants[i,j]))
        }
    }
    cat(i/nrow(processed_variants),"\n")
}
processed_variants = data4
for(i in 1:nrow(processed_variants)) {
    for(j in 1:ncol(processed_variants)) {
        if(!is.na(processed_variants[i,j])&&processed_variants[i,j]>0) {
            curr_split = strsplit(colnames(processed_variants)[j],"_")[[1]]
            variants_data = rbind(variants_data,c(curr_split,rownames(processed_variants)[i],processed_variants[i,j]))
        }
    }
    cat(i/nrow(processed_variants),"\n")
}
processed_variants = data5
for(i in 1:nrow(processed_variants)) {
    for(j in 1:ncol(processed_variants)) {
        if(!is.na(processed_variants[i,j])&&processed_variants[i,j]>0) {
            curr_split = strsplit(colnames(processed_variants)[j],"_")[[1]]
            variants_data = rbind(variants_data,c(curr_split,rownames(processed_variants)[i],processed_variants[i,j]))
        }
    }
    cat(i/nrow(processed_variants),"\n")
}
processed_variants = data6
for(i in 1:nrow(processed_variants)) {
    for(j in 1:ncol(processed_variants)) {
        if(!is.na(processed_variants[i,j])&&processed_variants[i,j]>0) {
            curr_split = strsplit(colnames(processed_variants)[j],"_")[[1]]
            variants_data = rbind(variants_data,c(curr_split,rownames(processed_variants)[i],processed_variants[i,j]))
        }
    }
    cat(i/nrow(processed_variants),"\n")
}
processed_variants = data7
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
save(annotations,file="results/dataset2/annotations.RData")
save(human_cu,file="results/dataset2/human_cu.RData")
save(sars_cov2_fasta,file="results/dataset2/sars_cov2_fasta.RData")
save(variants_data,file="results/dataset2/variants_data.RData")
