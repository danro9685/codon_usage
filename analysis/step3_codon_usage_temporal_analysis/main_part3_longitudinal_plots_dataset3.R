# load required libraries
library("ggplot2")
library("trend")

# load data
load("results/dataset3/codon_usage_information.RData")
load("final_dataset/dataset3_unique_variants.RData")
load("final_dataset/dataset3_variants.RData")

# preprocessing dates
tmp_dates = as.character(codon_usage_information[,"DATE"])
codon_usage_information[which(tmp_dates<"2020-02-16"),"DATE"] = "T01"
codon_usage_information[which(tmp_dates>="2020-02-16"&tmp_dates<"2020-03-01"),"DATE"] = "T02"
codon_usage_information[which(tmp_dates>="2020-03-01"&tmp_dates<"2020-03-16"),"DATE"] = "T03"
codon_usage_information[which(tmp_dates>="2020-03-16"&tmp_dates<"2020-04-01"),"DATE"] = "T04"
codon_usage_information[which(tmp_dates>="2020-04-01"&tmp_dates<"2020-04-16"),"DATE"] = "T05"
codon_usage_information[which(tmp_dates>="2020-04-16"&tmp_dates<"2020-05-01"),"DATE"] = "T06"
codon_usage_information[which(tmp_dates>="2020-05-01"&tmp_dates<"2020-05-16"),"DATE"] = "T07"
codon_usage_information[which(tmp_dates>="2020-05-16"&tmp_dates<"2020-06-01"),"DATE"] = "T08"
codon_usage_information[which(tmp_dates>="2020-06-01"&tmp_dates<"2020-06-16"),"DATE"] = "T09"
codon_usage_information[which(tmp_dates>="2020-06-16"&tmp_dates<"2020-07-01"),"DATE"] = "T10"
codon_usage_information[which(tmp_dates>="2020-07-01"&tmp_dates<"2020-07-16"),"DATE"] = "T11"
codon_usage_information[which(tmp_dates>="2020-07-16"&tmp_dates<"2020-08-01"),"DATE"] = "T12"
codon_usage_information[which(tmp_dates>="2020-08-01"&tmp_dates<"2020-08-16"),"DATE"] = "T13"
codon_usage_information[which(tmp_dates>="2020-08-16"&tmp_dates<"2020-09-01"),"DATE"] = "T14"
codon_usage_information[which(tmp_dates>="2020-09-01"&tmp_dates<"2020-09-16"),"DATE"] = "T15"
codon_usage_information[which(tmp_dates>="2020-09-16"&tmp_dates<"2020-10-01"),"DATE"] = "T16"
codon_usage_information[which(tmp_dates>="2020-10-01"&tmp_dates<"2020-10-16"),"DATE"] = "T17"
codon_usage_information[which(tmp_dates>="2020-10-16"&tmp_dates<"2020-11-01"),"DATE"] = "T18"

# make plot
variants_id_codon = paste0(codon_usage_information[,"POS"],"_",codon_usage_information[,"REF"],"_",codon_usage_information[,"ALT"])
dates_codon = as.character(codon_usage_information[,"DATE"])
data = codon_usage_information[,c("POS","REF","ALT")]
variants_id_data = paste0(data[,"POS"],"_",data[,"REF"],"_",data[,"ALT"])
RCU = NULL
VF = NULL
VARIANT = NULL
DATE = NULL
for(i in sort(names(table(codon_usage_information[,"DATE"])))) {
    curr_variants = variants_id_codon[which(dates_codon==i)]
    curr_data = dataset3_variants[which(variants_id_data%in%curr_variants),]
    RCU = c(RCU,as.numeric(curr_data[,"RCU"]))
    VF = c(VF,as.numeric(curr_data[,"VF"]))
    VARIANT = c(VARIANT,as.character(curr_data[,"VARIANT"]))
    DATE = c(DATE,rep(i,nrow(curr_data)))
}
data = data.frame(RCU=RCU*VF,VARIANT=VARIANT,DATE=DATE)
data = data[which(data$VARIANT%in%c("C>T","G>A")),]
ggplot(data,aes(x=DATE,y=RCU,fill=DATE)) + geom_boxplot() + ggtitle("Timeline - All variants (APOBEC) considering VF (Dataset 3)")
data = data.frame(RCU=RCU*VF,VARIANT=VARIANT,DATE=DATE)
data = data[which(!data$VARIANT%in%c("C>T","G>A")),]
ggplot(data,aes(x=DATE,y=RCU,fill=DATE)) + geom_boxplot() + ggtitle("Timeline - All variants (NOT APOBEC) considering VF (Dataset 3)")

# perform Mann-Kendall Trend Test

# APOBEC
set.seed(55555)
mean_rcu = NULL
data = data.frame(RCU=RCU*VF,VARIANT=VARIANT,DATE=DATE)
data = data[which(data$VARIANT%in%c("C>T","G>A")),]
for(i in names(table(data$DATE))) {
    mean_rcu = c(mean_rcu,mean(data$RCU[which(data$DATE==i)]))
}
print(mk.test(mean_rcu,alternative="less")$p.value) # pvalue = 0.001724072

# Not APOBEC
set.seed(66666)
mean_rcu = NULL
data = data.frame(RCU=RCU*VF,VARIANT=VARIANT,DATE=DATE)
data = data[which(!data$VARIANT%in%c("C>T","G>A")),]
for(i in names(table(data$DATE))) {
    mean_rcu = c(mean_rcu,mean(data$RCU[which(data$DATE==i)]))
}
print(mk.test(mean_rcu,alternative="greater")$p.value) # pvalue = 0.04561951
