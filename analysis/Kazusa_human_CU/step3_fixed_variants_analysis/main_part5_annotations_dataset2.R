# load data
load("data/dataset2/annotations.RData")
load("data/dataset2/human_cu.RData")
load("data/dataset2/sars_cov2_fasta.RData")
load("data/dataset2/variants.RData")
variants_data = variants[,c("POSITION","REFERENCE_ALLELE","VARIANT_ALLELE","PATIENT_ID","VARIANT_FREQUENCY")]
colnames(variants_data) = c("POSITION","REF","ALT","PATIENT","FREQUENCY")

# process annotations
cds_start = as.numeric(annotations[which(annotations[,3]=="CDS"),4])
cds_end = as.numeric(annotations[which(annotations[,3]=="CDS"),5])
cds_name = c("ORF1ab","ORF1ab","ORF1ab","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")

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
dataset = NULL
cont = 0
for(i in sample(1:nrow(variants_data),size=100000)) {
    cont = cont + 1
    curr_variants_data = variants_data[i,,drop=TRUE]
    curr_var_pos = as.numeric(curr_variants_data["POSITION"])
    curr_var_ref = as.character(curr_variants_data["REF"])
    curr_var_alt = as.character(curr_variants_data["ALT"])
    curr_position = which(cds_start<=curr_var_pos&cds_end>=curr_var_pos)
    curr_position = curr_position[1] # take first if I have multiple matches to CDS
    curr_cds = cds_name[curr_position]
    curr_offset = curr_var_pos-cds_start[curr_position]+1
    curr_cds_codon = ceiling(curr_offset/3)
    if(curr_offset%%3==1) { # first position of triplet
        triplet_wt = sars_cov2_fasta[c(curr_var_pos,(curr_var_pos+1),(curr_var_pos+2))]
        triplet_mut = c(curr_var_alt,triplet_wt[2],triplet_wt[3])
        triplet_mut_alt1 = c(c("A","C","G","T")[which(!c("A","C","G","T")%in%c(curr_var_ref,curr_var_alt))[1]],triplet_wt[2],triplet_wt[3])
        triplet_mut_alt2 = c(c("A","C","G","T")[which(!c("A","C","G","T")%in%c(curr_var_ref,curr_var_alt))[2]],triplet_wt[2],triplet_wt[3])
    }
    if(curr_offset%%3==2) { # second position of triplet
        triplet_wt = sars_cov2_fasta[c((curr_var_pos-1),curr_var_pos,(curr_var_pos+1))]
        triplet_mut = c(triplet_wt[1],curr_var_alt,triplet_wt[3])
        triplet_mut_alt1 = c(triplet_wt[1],c("A","C","G","T")[which(!c("A","C","G","T")%in%c(curr_var_ref,curr_var_alt))[1]],triplet_wt[3])
        triplet_mut_alt2 = c(triplet_wt[1],c("A","C","G","T")[which(!c("A","C","G","T")%in%c(curr_var_ref,curr_var_alt))[2]],triplet_wt[3])
    }
    if(curr_offset%%3==0) { # third position of triplet
        triplet_wt = sars_cov2_fasta[c((curr_var_pos-2),(curr_var_pos-1),curr_var_pos)]
        triplet_mut = c(triplet_wt[1],triplet_wt[2],curr_var_alt)
        triplet_mut_alt1 = c(triplet_wt[1],triplet_wt[2],c("A","C","G","T")[which(!c("A","C","G","T")%in%c(curr_var_ref,curr_var_alt))[1]])
        triplet_mut_alt2 = c(triplet_wt[1],triplet_wt[2],c("A","C","G","T")[which(!c("A","C","G","T")%in%c(curr_var_ref,curr_var_alt))[2]])
    }
    triplet_wt = paste0(triplet_wt,collapse="")
    triplet_wt = as.character(human_cu[which(human_cu[,1]==triplet_wt),])
    triplet_mut = paste0(triplet_mut,collapse="")
    triplet_mut = as.character(human_cu[which(human_cu[,1]==triplet_mut),])
    triplet_mut_alt1 = paste0(triplet_mut_alt1,collapse="")
    triplet_mut_alt1 = as.character(human_cu[which(human_cu[,1]==triplet_mut_alt1),])
    triplet_mut_alt2 = paste0(triplet_mut_alt2,collapse="")
    triplet_mut_alt2 = as.character(human_cu[which(human_cu[,1]==triplet_mut_alt2),])
    # save the results
    curr_res = c(paste0(curr_var_pos,"_",curr_var_ref,"_",curr_var_alt),curr_cds,curr_cds_codon,as.character(curr_variants_data["PATIENT"]),triplet_wt,triplet_mut,triplet_mut_alt1,triplet_mut_alt2)
    dataset = rbind(dataset,curr_res)
    cat(cont/100000,"\n")
}
rownames(dataset) = 1:nrow(dataset)
colnames(dataset) = c("VARIANT_ID","CDS_NAME","CDS_POS","PATIENT_ID","CODON_WT","PROTEIN_WT","VALUE_WT","CODON_MUT","PROTEIN_MUT","VALUE_MUT","CODON_ALT1","PROTEIN_ALT1","VALUE_ALT1","CODON_ALT2","PROTEIN_ALT2","VALUE_ALT2")
dataset2_all_variants = dataset

# save results
save(dataset2_all_variants,file="final_dataset/dataset2_all_variants.RData")
