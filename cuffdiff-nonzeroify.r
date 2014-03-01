#!/usr/bin/env Rscript
# this is a script to convert the values in the gene_exp.diff file produced by
# cuffdiff.  zero values present a problem b/c when calculating the log of the
# ratio, zero's cause infinite values.  the replace is calcualte as follows:
# 1. temporary remove all values with zeros in either value1 or 2
# 2. fit a normal (just use mean() and s())
# 3. calculate the area to the left of the lowest values 
# 4. calculate the median of that area and use it as replacement for 0.
# 
# writtern for Xiao

cmd = commandArgs()
file = cmd[length(cmd)]
output = paste(file,".non-zeroified", sep="")

if (!file.exists(file)){
    cat("usage: ", cmd[0], " gene_exp.diff\n")
    quit(save = "no", status = 1, runLast = F)
}
gene = read.csv(file, header=T, sep="\t")
non_zero_values = gene[gene$value_1 > 0 & gene$value_2 > 0,]

log_value1 = log(non_zero_values$value_1)
log_value2 = log(non_zero_values$value_2)

mu1 = mean(log_value1)
mu2 = mean(log_value2)
sd1 = sd(log_value1)
sd2 = sd(log_value2)

min_percentile_log_value1 = pnorm(min(log_value1), mu1, sd1)
min_percentile_log_value2 = pnorm(min(log_value2), mu2, sd2)

replacement_log_value1 = qnorm(min_percentile_log_value1 / 2, mu1, sd1)
replacement_log_value2 = qnorm(min_percentile_log_value2 / 2, mu2, sd2)

replacement_value1 = exp(replacement_log_value1) 
replacement_value2 = exp(replacement_log_value2) 

# min_percentile_log_value1 
# min_percentile_log_value2 
# min(log_value1)
# min(log_value2)
# replacement_log_value1
# replacement_log_value2

# replacement_value1 
# replacement_value2 

indices_with_both_zero = gene$value_1 == 0 & gene$value_2 == 0
indices_with_zero_value1 = gene$value_1 == 0 & ! indices_with_both_zero
indices_with_zero_value2 = gene$value_2 == 0 & ! indices_with_both_zero
indices_with_any_zero = indices_with_zero_value1 | indices_with_zero_value2 

gene[indices_with_zero_value1,'value_1'] = replacement_value1
gene[indices_with_zero_value2,'value_2'] = replacement_value2
gene[indices_with_any_zero,'log2.fold_change.'] = log2( gene[indices_with_any_zero,'value_1'] / gene[indices_with_any_zero,'value_2'] )
gene[indices_with_any_zero,'p_value'] = 1
gene[indices_with_any_zero,'q_value'] = 1
gene[indices_with_any_zero,'significant'] = 'no'

write.table(gene, file=output, sep = "\t", quote=F)

