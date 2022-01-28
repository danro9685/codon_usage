# load data
load("data/dataset2/annotations.RData")
load("data/dataset2/human_cu.RData")
load("data/dataset2/metadata.RData")
load("data/dataset2/sars_cov2_fasta.RData")
load("data/dataset2/variants.RData")
variants_data = variants[,c("POSITION","REFERENCE_ALLELE","VARIANT_ALLELE","PATIENT_ID","VARIANT_FREQUENCY")]
colnames(variants_data) = c("POSITION","REF","ALT","PATIENT","FREQUENCY")

# process annotations
cds_start = as.numeric(annotations[which(annotations[,3]=="CDS"),4])
cds_end = as.numeric(annotations[which(annotations[,3]=="CDS"),5])

# consider only variants spanning CDS regions
cds_positions = NULL
for(i in 1:length(cds_start)) {
    cds_positions = c(cds_positions,cds_start[i]:cds_end[i])
}
variants_data = variants_data[which(as.numeric(variants_data[,"POSITION"])%in%cds_positions),]

# remove variants in uncertain positions (different between REF-ANC and Wuhan1)
uncertain_positions = NULL
uncertain_variant = 8782
uncertain_variant_offset = uncertain_variant-cds_start[which(cds_start<=uncertain_variant&cds_end>=uncertain_variant)[1]]+1
if(uncertain_variant_offset%%3==1) { # first position of triplet
    uncertain_positions = c(uncertain_positions,c(uncertain_variant,(uncertain_variant+1),(uncertain_variant+2)))
}
if(uncertain_variant_offset%%3==2) { # second position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-1),uncertain_variant,c(uncertain_variant+1)))
}
if(uncertain_variant_offset%%3==0) { # third position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-2),(uncertain_variant-1),uncertain_variant))
}
uncertain_variant = 28144
uncertain_variant_offset = uncertain_variant-cds_start[which(cds_start<=uncertain_variant&cds_end>=uncertain_variant)[1]]+1
if(uncertain_variant_offset%%3==1) { # first position of triplet
    uncertain_positions = c(uncertain_positions,c(uncertain_variant,(uncertain_variant+1),(uncertain_variant+2)))
}
if(uncertain_variant_offset%%3==2) { # second position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-1),uncertain_variant,c(uncertain_variant+1)))
}
if(uncertain_variant_offset%%3==0) { # third position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-2),(uncertain_variant-1),uncertain_variant))
}
uncertain_variant = 29095
uncertain_variant_offset = uncertain_variant-cds_start[which(cds_start<=uncertain_variant&cds_end>=uncertain_variant)[1]]+1
if(uncertain_variant_offset%%3==1) { # first position of triplet
    uncertain_positions = c(uncertain_positions,c(uncertain_variant,(uncertain_variant+1),(uncertain_variant+2)))
}
if(uncertain_variant_offset%%3==2) { # second position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-1),uncertain_variant,c(uncertain_variant+1)))
}
if(uncertain_variant_offset%%3==0) { # third position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-2),(uncertain_variant-1),uncertain_variant))
}
variants_data = variants_data[which(!as.numeric(variants_data[,"POSITION"])%in%uncertain_positions),]

# perform analysis
codon_usage_information = NULL
for(i in 1:nrow(variants_data)) {
    curr_variants_data = variants_data[i,,drop=TRUE]
    curr_position = which(cds_start<=as.numeric(curr_variants_data["POSITION"])&cds_end>=as.numeric(curr_variants_data["POSITION"]))
    curr_position = curr_position[1] # take first if I have multiple matches to CDS
    curr_offset = as.numeric(curr_variants_data["POSITION"])-cds_start[curr_position]+1
    if(curr_offset%%3==1) { # first position of triplet
        triplet_wt = sars_cov2_fasta[c(as.numeric(curr_variants_data["POSITION"]),(as.numeric(curr_variants_data["POSITION"])+1),(as.numeric(curr_variants_data["POSITION"])+2))]
        triplet_mut = c(as.character(curr_variants_data["ALT"]),triplet_wt[2],triplet_wt[3])
    }
    if(curr_offset%%3==2) { # second position of triplet
        triplet_wt = sars_cov2_fasta[c((as.numeric(curr_variants_data["POSITION"])-1),as.numeric(curr_variants_data["POSITION"]),(as.numeric(curr_variants_data["POSITION"])+1))]
        triplet_mut = c(triplet_wt[1],as.character(curr_variants_data["ALT"]),triplet_wt[3])
    }
    if(curr_offset%%3==0) { # third position of triplet
        triplet_wt = sars_cov2_fasta[c((as.numeric(curr_variants_data["POSITION"])-2),(as.numeric(curr_variants_data["POSITION"])-1),as.numeric(curr_variants_data["POSITION"]))]
        triplet_mut = c(triplet_wt[1],triplet_wt[2],as.character(curr_variants_data["ALT"]))
    }
    triplet_wt = paste0(triplet_wt,collapse="")
    triplet_wt = as.character(human_cu[which(human_cu[,1]==triplet_wt),])
    triplet_mut = paste0(triplet_mut,collapse="")
    triplet_mut = as.character(human_cu[which(human_cu[,1]==triplet_mut),])
    # if the variant is synonymous then it is valid
    if(triplet_wt[2]==triplet_mut[2]) {
        curr_res = as.character(c(curr_variants_data,triplet_wt,triplet_mut))
        curr_res = c(curr_res[c(1,2,3,5,6,7,8,9,10,11,4)],metadata$CollectionDate[which(metadata$Run==curr_res[4])])
        codon_usage_information = rbind(codon_usage_information,curr_res)
    }
    cat(i/nrow(variants_data),"\n")
}
rownames(codon_usage_information) = 1:nrow(codon_usage_information)
colnames(codon_usage_information) = c("POS","REF","ALT","VF","CWT","PWT","VWT","CMUT","PMUT","VMUT","PATIENT","DATE")

# save results
save(codon_usage_information,file="results/dataset2/codon_usage_information.RData")
