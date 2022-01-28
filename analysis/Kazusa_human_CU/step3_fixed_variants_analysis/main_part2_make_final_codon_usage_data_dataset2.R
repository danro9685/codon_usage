# load data
load("data/dataset2/annotations.RData")
load("results/dataset2/codon_usage_information.RData")
load("data/dataset2/human_cu.RData")
load("data/dataset2/sars_cov2_fasta.RData")

# We consider two periods of 10 months the first half and 9 months the second half, starting 2020/01/01

# DATA 1: consider only unique variants for first half of timeline
codon_usage_information_selected = codon_usage_information[which(as.Date(paste0(as.character(substr(codon_usage_information[,"DATE"],1,7)),"-01"),format="%Y-%m-%d")<=as.Date("2020-10-31",format="%Y-%m-%d")),]
codon_usage_information_selected = codon_usage_information_selected[which(as.numeric(codon_usage_information_selected[,"VF"])>0.90),]
data = unique(codon_usage_information_selected[,c("POS","REF","ALT","VWT","VMUT")])
selected_variants = paste0(data[,"POS"],"_",data[,"REF"],"_",data[,"ALT"])
RCU_val = as.numeric(log(as.numeric(data[,"VMUT"])/as.numeric(data[,"VWT"]),base=2))
unique_variants = paste0(data[,"POS"],"_",data[,"REF"],"_",data[,"ALT"])
all_variants = paste0(codon_usage_information[,"POS"],"_",codon_usage_information[,"REF"],"_",codon_usage_information[,"ALT"])
mean_vf = NULL
for(i in unique_variants) {
    mean_vf = c(mean_vf,mean(as.numeric(codon_usage_information[which(all_variants==i),"VF"])))
}
VF = mean_vf
VARIANT = paste0(data[,"REF"],">",data[,"ALT"])
reduced_variants = VARIANT
reduced_variants[which(reduced_variants=="C>A")] = "C>A/G>T"
reduced_variants[which(reduced_variants=="G>T")] = "C>A/G>T"
reduced_variants[which(reduced_variants=="C>G")] = "C>G/G>C"
reduced_variants[which(reduced_variants=="G>C")] = "C>G/G>C"
reduced_variants[which(reduced_variants=="C>T")] = "C>T/G>A"
reduced_variants[which(reduced_variants=="G>A")] = "C>T/G>A"
reduced_variants[which(reduced_variants=="A>T")] = "A>T/T>A"
reduced_variants[which(reduced_variants=="T>A")] = "A>T/T>A"
reduced_variants[which(reduced_variants=="A>G")] = "A>G/T>C"
reduced_variants[which(reduced_variants=="T>C")] = "A>G/T>C"
reduced_variants[which(reduced_variants=="A>C")] = "A>C/T>G"
reduced_variants[which(reduced_variants=="T>G")] = "A>C/T>G"
VARIANT_REDUCED = reduced_variants
dataset2_unique_variants_first_half = data.frame(VARIANT_ID=selected_variants,RCU=RCU_val,VF=VF,VARIANT=VARIANT,VARIANT_REDUCED=VARIANT_REDUCED)

# DATA 2: consider only unique variants for second half of timeline
codon_usage_information_selected = codon_usage_information[which(as.Date(paste0(as.character(substr(codon_usage_information[,"DATE"],1,7)),"-01"),format="%Y-%m-%d")>as.Date("2020-10-31",format="%Y-%m-%d")),]
codon_usage_information_selected = codon_usage_information_selected[which(as.numeric(codon_usage_information_selected[,"VF"])>0.90),]
data = unique(codon_usage_information_selected[,c("POS","REF","ALT","VWT","VMUT")])
selected_variants = paste0(data[,"POS"],"_",data[,"REF"],"_",data[,"ALT"])
RCU_val = as.numeric(log(as.numeric(data[,"VMUT"])/as.numeric(data[,"VWT"]),base=2))
unique_variants = paste0(data[,"POS"],"_",data[,"REF"],"_",data[,"ALT"])
all_variants = paste0(codon_usage_information[,"POS"],"_",codon_usage_information[,"REF"],"_",codon_usage_information[,"ALT"])
mean_vf = NULL
for(i in unique_variants) {
    mean_vf = c(mean_vf,mean(as.numeric(codon_usage_information[which(all_variants==i),"VF"])))
}
VF = mean_vf
VARIANT = paste0(data[,"REF"],">",data[,"ALT"])
reduced_variants = VARIANT
reduced_variants[which(reduced_variants=="C>A")] = "C>A/G>T"
reduced_variants[which(reduced_variants=="G>T")] = "C>A/G>T"
reduced_variants[which(reduced_variants=="C>G")] = "C>G/G>C"
reduced_variants[which(reduced_variants=="G>C")] = "C>G/G>C"
reduced_variants[which(reduced_variants=="C>T")] = "C>T/G>A"
reduced_variants[which(reduced_variants=="G>A")] = "C>T/G>A"
reduced_variants[which(reduced_variants=="A>T")] = "A>T/T>A"
reduced_variants[which(reduced_variants=="T>A")] = "A>T/T>A"
reduced_variants[which(reduced_variants=="A>G")] = "A>G/T>C"
reduced_variants[which(reduced_variants=="T>C")] = "A>G/T>C"
reduced_variants[which(reduced_variants=="A>C")] = "A>C/T>G"
reduced_variants[which(reduced_variants=="T>G")] = "A>C/T>G"
VARIANT_REDUCED = reduced_variants
dataset2_unique_variants_second_half = data.frame(VARIANT_ID=selected_variants,RCU=RCU_val,VF=VF,VARIANT=VARIANT,VARIANT_REDUCED=VARIANT_REDUCED)
dataset2_unique_variants_second_half = dataset2_unique_variants_second_half[which(!dataset2_unique_variants_second_half$VARIANT_ID%in%dataset2_unique_variants_first_half$VARIANT_ID),]
rownames(dataset2_unique_variants_second_half) = 1:nrow(dataset2_unique_variants_second_half)

# save results
save(dataset2_unique_variants_first_half,file="final_dataset/dataset2_unique_variants_first_half.RData")
save(dataset2_unique_variants_second_half,file="final_dataset/dataset2_unique_variants_second_half.RData")
