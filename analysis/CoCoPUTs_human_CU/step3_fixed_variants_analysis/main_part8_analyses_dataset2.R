# load data
load("final_dataset/dataset2_unique_variants_first_half.RData")
load("final_dataset/dataset2_unique_variants_second_half.RData")
load("final_dataset/dataset2_all_variants.RData")

# test proportion of variants vs time (perform z-test of two proportions)
APOBEC_group1 = as.numeric(table(dataset2_unique_variants_first_half$VARIANT_REDUCED)["C>T/G>A"])
APOBEC_group2 = as.numeric(table(dataset2_unique_variants_second_half$VARIANT_REDUCED)["C>T/G>A"])
total_group1 = sum(table(dataset2_unique_variants_first_half$VARIANT_REDUCED))
total_group2 = sum(table(dataset2_unique_variants_second_half$VARIANT_REDUCED))
print(prop.test(x=c(APOBEC_group1,APOBEC_group2),n=c(total_group1,total_group2))) # 0.4234234 vs 0.1517503
print(prop.test(x=c(APOBEC_group1,APOBEC_group2),n=c(total_group1,total_group2))$p.value) # p-value < 2.2e-16

# now evaluate codon usage of selected mutations vs all other possible mutations (wild type is not considered)
dataset_mut = NULL
dataset_alt = NULL
for(i in 1:nrow(dataset2_all_variants)) {
    curr_valid = NULL
    if(dataset2_all_variants[i,"PROTEIN_MUT"]==dataset2_all_variants[i,"PROTEIN_ALT1"]) {
        curr_valid = c(curr_valid,as.numeric(dataset2_all_variants[i,"VALUE_ALT1"]))
    }
    if(dataset2_all_variants[i,"PROTEIN_MUT"]==dataset2_all_variants[i,"PROTEIN_ALT2"]) {
        curr_valid = c(curr_valid,as.numeric(dataset2_all_variants[i,"VALUE_ALT2"]))
    }
    if(length(curr_valid)>0) {
        dataset_mut = c(dataset_mut,as.numeric(dataset2_all_variants[i,"VALUE_MUT"]))
        dataset_alt = c(dataset_alt,mean(curr_valid))
    }
    cat(i/nrow(dataset2_all_variants),"\n")
}
print(length(dataset_mut))
print(mean(dataset_mut))
print(mean(dataset_alt))
print(median(dataset_mut))
print(median(dataset_alt))
print(fivenum(dataset_mut)) # 0.0035 0.1180 0.1696 0.1915 0.3940
print(fivenum(dataset_alt)) # 0.0035 0.1108 0.1247 0.1673 0.3940
print(t.test(dataset_mut,dataset_alt)) # p-value < 2.2e-16
