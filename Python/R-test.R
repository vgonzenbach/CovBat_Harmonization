# script to generate CovBat R output to compare with Python output

# BiocManager::install("bladderbatch")
# BiocManager::install("sva")
library(sva)
library(CovBat)
library(bladderbatch)
options(stringsAsFactors=FALSE)
data(bladderdata)

pheno = pData(bladderEset)
# add fake age variable for numeric
pheno$age = c(1:7, rep(1:10, 5))
write.table(data.frame(cel=rownames(pheno), pheno), row.names=F, quote=F, sep="\t", file="bladder-pheno.txt")

edata = exprs(bladderEset)
write.table(edata, row.names=T, quote=F, sep="\t", file="bladder-expr.txt")
# use dataframe instead of matrix
mod = model.matrix(~as.factor(cancer) + age, data=pheno)
t = Sys.time()
cov_out = covbat(edata, bat = as.factor(pheno$batch), mod = mod, 
                 std.var = TRUE, percent.var = 0.80)
covdata = cov_out$dat.covbat
print(Sys.time() - t)
print(covdata[1:5, 1:5])

t = Sys.time()
com_out = combat_modded(edata, bat = as.factor(pheno$batch), mod = mod)
comdata = com_out$dat.combat
print(Sys.time() - t)
print(comdata[1:5, 1:5])

write.table(covdata, "r-covbat.txt", sep="\t", quote=F)
write.table(comdata, "r-combat.txt", sep="\t", quote=F)
