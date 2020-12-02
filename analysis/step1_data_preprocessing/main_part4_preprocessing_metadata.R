# load data
load(file="final_data/PRJEB37886/statistics.RData")
metadata = statistics
load(file="final_data/PRJEB39849/statistics.RData")
metadata = rbind(metadata,statistics)
load(file="final_data/PRJNA613958/statistics.RData")
metadata = rbind(metadata,statistics)
load(file="final_data/PRJNA614995/statistics.RData")
metadata = rbind(metadata,statistics)
load(file="final_data/PRJNA625551/statistics.RData")
metadata = rbind(metadata,statistics)
load(file="final_data/PRJNA633948/statistics.RData")
metadata = rbind(metadata,statistics)
load(file="final_data/PRJNA636748/statistics.RData")
metadata = rbind(metadata,statistics)
load(file="final_data/PRJNA639066/statistics.RData")
metadata = rbind(metadata,statistics)
load(file="final_data/PRJNA645906/statistics.RData")
metadata = rbind(metadata,statistics)
rownames(metadata) = 1:nrow(metadata)
all_metadata = metadata

# dataset 1
load("results/dataset1/variants_data.RData")
metadata = all_metadata[which(all_metadata$Run%in%variants_data[,"PATIENT"]),]
rownames(metadata) = 1:nrow(metadata)
save(metadata,file="results/dataset1/metadata.RData")

# dataset 2
load("results/dataset2/variants_data.RData")
metadata = all_metadata[which(all_metadata$Run%in%variants_data[,"PATIENT"]),]
rownames(metadata) = 1:nrow(metadata)
save(metadata,file="results/dataset2/metadata.RData")

# dataset 3
load("results/dataset3/variants_data.RData")
metadata = all_metadata[which(all_metadata$Run%in%variants_data[,"PATIENT"]),]
rownames(metadata) = 1:nrow(metadata)
save(metadata,file="results/dataset3/metadata.RData")
