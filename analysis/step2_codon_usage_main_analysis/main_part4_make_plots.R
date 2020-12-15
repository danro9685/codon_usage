# load required libraries
library("ggplot2")

# load data
load("final_dataset/dataset1_unique_variants.RData")
load("final_dataset/dataset1_variants.RData")
load("final_dataset/dataset2_unique_variants.RData")
load("final_dataset/dataset2_variants.RData")
load("final_dataset/dataset3_unique_variants.RData")
load("final_dataset/dataset3_variants.RData")

# DATASET 1

# PLOT 1-2: consider only unique variants
raw_data = dataset1_unique_variants
data = raw_data
data = data.frame(RCU=as.numeric(data[,"RCU"]),VARIANT=factor(as.character(data[,"VARIANT"]),levels=c("A>C","T>G","A>G","T>C","A>T","T>A","C>A","G>T","C>G","G>C","C>T","G>A")),VARIANT_REDUCED=as.character(data[,"VARIANT_REDUCED"]))
ggplot(data,aes(x=VARIANT,y=RCU,fill=VARIANT_REDUCED)) + geom_boxplot() + ggtitle("Boxplot - Unique variants (Dataset 1)")
slices = as.numeric(table(data$VARIANT))
labels = NULL
for(i in  1:length(table(data$VARIANT))) {
    labels = c(labels,paste0(names(table(data$VARIANT))[i]," (",as.numeric(table(data$VARIANT))[i],")"))
}
pie(slices,labels,main="Pie Chart - Unique variants (Dataset 1)")

# PLOT 3-4: consider all variants
raw_data = dataset1_variants
data = raw_data
data = data.frame(RCU=as.numeric(data[,"RCU"]),VARIANT=factor(as.character(data[,"VARIANT"]),levels=c("A>C","T>G","A>G","T>C","A>T","T>A","C>A","G>T","C>G","G>C","C>T","G>A")),VARIANT_REDUCED=as.character(data[,"VARIANT_REDUCED"]))
ggplot(data,aes(x=VARIANT,y=RCU,fill=VARIANT_REDUCED)) + geom_boxplot() + ggtitle("Boxplot - All variants (Dataset 1)")
slices = as.numeric(table(data$VARIANT))
labels = NULL
for(i in  1:length(table(data$VARIANT))) {
    labels = c(labels,paste0(names(table(data$VARIANT))[i]," (",as.numeric(table(data$VARIANT))[i],")"))
}
pie(slices,labels,main="Pie Chart - All variants (Dataset 1)")

# DATASET 2

# PLOT 5-6: consider only unique variants
raw_data = dataset2_unique_variants
data = raw_data
data = data.frame(RCU=as.numeric(data[,"RCU"]),VARIANT=factor(as.character(data[,"VARIANT"]),levels=c("A>C","T>G","A>G","T>C","A>T","T>A","C>A","G>T","C>G","G>C","C>T","G>A")),VARIANT_REDUCED=as.character(data[,"VARIANT_REDUCED"]))
ggplot(data,aes(x=VARIANT,y=RCU,fill=VARIANT_REDUCED)) + geom_boxplot() + ggtitle("Boxplot - Unique variants (Dataset 2)")
slices = as.numeric(table(data$VARIANT))
labels = NULL
for(i in  1:length(table(data$VARIANT))) {
    labels = c(labels,paste0(names(table(data$VARIANT))[i]," (",as.numeric(table(data$VARIANT))[i],")"))
}
pie(slices,labels,main="Pie Chart - Unique variants (Dataset 2)")

# PLOT 7-8: consider all variants
raw_data = dataset2_variants
data = raw_data
data = data.frame(RCU=as.numeric(data[,"RCU"]),VARIANT=factor(as.character(data[,"VARIANT"]),levels=c("A>C","T>G","A>G","T>C","A>T","T>A","C>A","G>T","C>G","G>C","C>T","G>A")),VARIANT_REDUCED=as.character(data[,"VARIANT_REDUCED"]))
ggplot(data,aes(x=VARIANT,y=RCU,fill=VARIANT_REDUCED)) + geom_boxplot() + ggtitle("Boxplot - All variants (Dataset 2)")
slices = as.numeric(table(data$VARIANT))
labels = NULL
for(i in  1:length(table(data$VARIANT))) {
    labels = c(labels,paste0(names(table(data$VARIANT))[i]," (",as.numeric(table(data$VARIANT))[i],")"))
}
pie(slices,labels,main="Pie Chart - All variants (Dataset 2)")

# DATASET 3

# PLOT 9-10: consider only unique variants
raw_data = dataset3_unique_variants
data = raw_data
data = data.frame(RCU=as.numeric(data[,"RCU"]),VARIANT=factor(as.character(data[,"VARIANT"]),levels=c("A>C","T>G","A>G","T>C","A>T","T>A","C>A","G>T","C>G","G>C","C>T","G>A")),VARIANT_REDUCED=as.character(data[,"VARIANT_REDUCED"]))
ggplot(data,aes(x=VARIANT,y=RCU,fill=VARIANT_REDUCED)) + geom_boxplot() + ggtitle("Boxplot - Unique variants (Dataset 3)")
slices = as.numeric(table(data$VARIANT))
labels = NULL
for(i in  1:length(table(data$VARIANT))) {
    labels = c(labels,paste0(names(table(data$VARIANT))[i]," (",as.numeric(table(data$VARIANT))[i],")"))
}
pie(slices,labels,main="Pie Chart - Unique variants (Dataset 3)")

# PLOT 11-12: consider all variants
raw_data = dataset3_variants
data = raw_data
data = data.frame(RCU=as.numeric(data[,"RCU"]),VARIANT=factor(as.character(data[,"VARIANT"]),levels=c("A>C","T>G","A>G","T>C","A>T","T>A","C>A","G>T","C>G","G>C","C>T","G>A")),VARIANT_REDUCED=as.character(data[,"VARIANT_REDUCED"]))
ggplot(data,aes(x=VARIANT,y=RCU,fill=VARIANT_REDUCED)) + geom_boxplot() + ggtitle("Boxplot - All variants (Dataset 3)")
slices = as.numeric(table(data$VARIANT))
labels = NULL
for(i in  1:length(table(data$VARIANT))) {
    labels = c(labels,paste0(names(table(data$VARIANT))[i]," (",as.numeric(table(data$VARIANT))[i],")"))
}
pie(slices,labels,main="Pie Chart - All variants (Dataset 3)")
