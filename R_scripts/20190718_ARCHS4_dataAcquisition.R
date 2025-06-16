library(dplyr)
library(rhdf5)
library(readxl)
library(preprocessCore)
library(xlsx)


##### ARCHS4 processing -----------------------------------------
#
#
#
#
#
#

fname = readClipboard()
ict_list_fname = "C:\\Users\\antho\\Documents\\NYU\\NYU Langone\\PhD\\Feske Lab\\Experiments\\05.20.18_ORAI3 Analyses\\R\\2019.07.18_ICTlist.txt"


# Acquires genes & sample metadata
genes = h5read(fname, "meta/genes")
samp_geo = as.data.frame(h5read(fname,"meta/Sample_geo_accession"))
samp_series = as.data.frame(h5read(fname,"meta/Sample_series_id"))
samp_info1 = as.data.frame(h5read(fname, "meta/Sample_source_name_ch1"))
samp_info2 = as.data.frame(h5read(fname, "meta/Sample_characteristics_ch1"))
samp_info3 = as.data.frame(h5read(fname,"meta/Sample_title"))

samp = bind_cols(samp_geo, samp_series, samp_info1, samp_info2, samp_info3)
heading = c('GEO #', 'Series #', 'Source Name', 'Characteristics', 'Title')
colnames(samp) = heading

# Acquires indices for ICT genes
for (i in 1:length(genes)) {
  genes[i] = tolower(genes[i])
}

ict_list = read.table(ict_list_fname)
icts = c()
for (i in 1:length(ict_list[,1])) {
  ict = as.character(ict_list[i,1])
  loc = match(tolower(ict),genes)
  
  if (!(is.na(loc))) {
    icts = c(icts, loc)
  }
}
genes = toupper(genes)


# retrieves expression data defined by 'index' argument
expression = h5read(fname, "data/expression", index=list(icts, 1:dim(samp_info1)[1]))
genes = genes[icts]
H5close()

# normalize samples & correct for differences in gene count distribution
expression = log2(expression+1)
expression = normalize.quantiles(expression)

rownames(expression) = genes

write.table(expression, 'ARCHS_orai_expression.txt', sep='\t')
write.table(samp, 'ARCHS_orai_meta.txt', sep='\t')

