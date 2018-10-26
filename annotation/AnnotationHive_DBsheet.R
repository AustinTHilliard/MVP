# ------------------------------------------------------------------------------
# setup and load table of annotation databases
# ------------------------------------------------------------------------------
# --clear workspace, set options and working directory
rm(list=ls())
options(stringsAsFactors=FALSE)
setwd("~/Documents/")

# --set name of file containing information on annotation databases
dbs_file = "AnnotationHive List of Annotation DBs - results-20180802-095742.csv.tsv"

# --check how many rows and columns it has
# rows
system(paste0("wc -l \"", dbs_file, "\""))
# 7009 AnnotationHive List of Annotation DBs - results-20180802-095742.csv.tsv

# columns
system(paste0("awk 'NR==1{print NF}' \"", dbs_file, "\""))
# 8

# load file with annotation databases
#  it should have 7009 rows and 8 columns
# first define colClasses
dbs_colClasses = c("character","numeric",rep("character",4),"numeric","character")
dbs = read.table(dbs_file, header=TRUE, sep="\t", colClasses=dbs_colClasses)
dim(dbs)
# [1] 7009    8

# ------------------------------------------------------------------------------
# explore species represented in the databases
# ------------------------------------------------------------------------------
# both the Build column and the first part of AnnotationSetName seem to represent species

# split the AnnotationSetName column on "_" and grab the first element from each row
dbspl = strsplit(dbs$AnnotationSetName, "_")
dbspl1 = sapply(dbspl, function(f) f[1])

# do they match? no
all(dbspl1 == dbs$Build)
# [1] FALSE

# which values are different?
setdiff(unique(dbspl1), unique(dbs$Build))
# character(0)
setdiff(unique(dbs$Build), unique(dbspl1))
# [1] "GRCh37/hg19"

# examine the counts of different values in each
countsBuild = table(dbs$Build)
countsAnno1 = table(dbspl1)
cbind(Build=countsBuild, Anno1=countsAnno1[match(names(countsBuild), names(countsAnno1))])
#             Build Anno1
# ailMel1        27    27
# bosTau8        27    27
# canFam3        37    37
# cavPor3        38    38
# chlSab2        19    19
# eboVir3        20    20
# equCab2       211   211
# gorGor5        17    17
# GRCh37          2     2
# GRCh37/hg19     5    NA
# GRCh38          2     2
# hg19         5544  5549
# hg37            2     2
# hg38          413   413
# loxAfr3        21    21
# melGal5        40    40
# mm10          507   507
# panTro5        31    31
# papAnu2        27    27
# rhiRox1        19    19

subset(dbs, Build=="GRCh37/hg19")[, 1]
# [1] "hg19_wellderly_all" "hg19_gnomAD_genomes" "hg19_gnomAD_exomes" "hg19_CADD_InDels" "hg19_CADD_wholeGenomeSNVs"

# conclusions
# 1. don't worry about the Build column
# the only difference between the Build and AnnotationSetName species are the 5 databases above
#  Build refers to them as "GRCh37/hg19" and AnnotationSetName refers to them as "hg19"
#
# 2. ~80% of all databases are for the hg19 human genome build, another ~6% are hg37 or hg38
# mouse (mm10, ~7%) and horse (equCab2, ~3%) are the most well-represented other species

# make a vector distinguish human databases from other species
hsvec = ifelse(grepl("^hg|^GRCh", dbspl1), TRUE, FALSE)
table(hsvec)
# hsvec
# FALSE  TRUE 
#  1041  5968

# ------------------------------------------------------------------------------
# explore information in the AnnotationSetName column
# ------------------------------------------------------------------------------
# do all entries have same number of elements?
dbsplLen = sapply(dbspl, length)
table(dbsplLen)
# dbsplLen
# 2    3    4    5    6 
# 1 6614  268   62   64 

# those where dbsplLen==3 seem to be species_source_featureType
head(dbs[dbsplLen == 3, 1:4], 10)
#                AnnotationSetName AnnotationSetVersion AnnotationSetType           AnnotationSetFields
# 1          hg19_UCSC_laminB1Lads                    1           generic bin,chrom,chromStart,chromEnd
# 2  hg19_UCSC_hiSeqDepthTopPt5Pct                    1           generic bin,chrom,chromStart,chromEnd
# 3         hg19_UCSC_crisprRanges                    1           generic bin,chrom,chromStart,chromEnd
# 4    hg19_UCSC_hiSeqDepthTop5Pct                    1           generic bin,chrom,chromStart,chromEnd
# 5    hg19_UCSC_hiSeqDepthTop1Pct                    1           generic bin,chrom,chromStart,chromEnd
# 6   hg19_UCSC_hiSeqDepthTop10Pct                    1           generic bin,chrom,chromStart,chromEnd
# 7     hg19_UCSC_phastBiasTracts3                    1           generic bin,chrom,chromStart,chromEnd
# 8        hg19_UCSC_eioJcviNASNeg                    1           generic bin,chrom,chromStart,chromEnd
# 9        hg19_UCSC_eioJcviNASPos                    1           generic bin,chrom,chromStart,chromEnd
# 10 hg19_UCSC_hiSeqDepthTopPt1Pct                    1           generic bin,chrom,chromStart,chromEnd

# what are the databases where dbsplLen > 3?
table(dbspl1[dbsplLen > 3])
# ailMel1 bosTau8 canFam3 cavPor3 equCab2 gorGor5    hg19    hg38 melGal5    mm10 panTro5 papAnu2 
#       2       2       2       2     166       1      13       9       2     191       2       2 

# most are from non-human species, what do the human ones look like?


#  the additional elements mostly designate chromosome
#  but for each feature type there's also one database that specifies "all" 
#  e.g. for mrna in mouse:
dbs[dbsplLen>3 & !hsvec & grepl("mm.*mrna",dbs$AnnotationSetName), c(1,7)]
#                        AnnotationSetName AnnotationSetSize
# 5857       mm10_UCSC_chrUn_GL456368_mrna                 2
# 5861       mm10_UCSC_chrUn_GL456378_mrna                13
# 5865 mm10_UCSC_chr5_GL456354_random_mrna                28
# 5866                mm10_UCSC_chr10_mrna             32469
# 5867                mm10_UCSC_chr16_mrna             16245
# 5868 mm10_UCSC_chr4_GL456216_random_mrna                18
# 5869       mm10_UCSC_chrUn_GL456385_mrna                 2
# 5878       mm10_UCSC_chrUn_GL456366_mrna                 6
# 5882       mm10_UCSC_chrUn_GL456372_mrna                19
# 5884                 mm10_UCSC_chr9_mrna             32516
# 5886 mm10_UCSC_chr4_JH584293_random_mrna               615
# 5888 mm10_UCSC_chr4_JH584295_random_mrna                 1
# 5891 mm10_UCSC_chrY_JH584301_random_mrna                 9
# 5892       mm10_UCSC_chrUn_GL456394_mrna                 2
# 5893       mm10_UCSC_chrUn_GL456382_mrna                27
# 5896                mm10_UCSC_chr12_mrna             43882
# 5897 mm10_UCSC_chr5_JH584298_random_mrna                17
# 5903 mm10_UCSC_chr1_GL456210_random_mrna               119
# 5904 mm10_UCSC_chr1_GL456213_random_mrna                 2
# 5907                mm10_UCSC_chr14_mrna             24182
# 5909       mm10_UCSC_chrUn_GL456370_mrna                 8
# 5910 mm10_UCSC_chrY_JH584302_random_mrna                 6
# 5911       mm10_UCSC_chrUn_GL456359_mrna                 1
# 5914 mm10_UCSC_chr5_JH584299_random_mrna                71
# 5918                 mm10_UCSC_chr3_mrna             24066
# 5919                 mm10_UCSC_chr5_mrna             34920
# 5927       mm10_UCSC_chrUn_GL456393_mrna                13
# 5928                 mm10_UCSC_chr1_mrna             32923
# 5930                 mm10_UCSC_chr6_mrna             43669
# 5932                mm10_UCSC_chr17_mrna             28309
# 5935                 mm10_UCSC_chr7_mrna             39556
# 5936       mm10_UCSC_chrUn_GL456392_mrna                 1
# 5938                 mm10_UCSC_chr2_mrna             46510
# 5939 mm10_UCSC_chr7_GL456219_random_mrna                22
# 5941                 mm10_UCSC_chrX_mrna             24549
# 5942 mm10_UCSC_chr5_JH584297_random_mrna                25
# 5945                mm10_UCSC_chr13_mrna             21357
# 5952       mm10_UCSC_chrUn_JH584304_mrna                64
# 5953                  mm10_UCSC_all_mrna            610553
# 5954 mm10_UCSC_chr4_JH584294_random_mrna               311
# 5955 mm10_UCSC_chr1_GL456211_random_mrna               290
# 5956                mm10_UCSC_chr18_mrna             16805
# 5959                 mm10_UCSC_chrM_mrna               107
# 5963 mm10_UCSC_chr5_JH584296_random_mrna                20
# 5964 mm10_UCSC_chrX_GL456233_random_mrna                55
# 5966                mm10_UCSC_chr11_mrna             35665
# 5971 mm10_UCSC_chr1_GL456221_random_mrna               296
# 5972       mm10_UCSC_chrUn_GL456360_mrna                 1
# 5978       mm10_UCSC_chrUn_GL456367_mrna                 5
# 5981                 mm10_UCSC_chrY_mrna              2206
# 5986                 mm10_UCSC_chr4_mrna             35950
# 5999 mm10_UCSC_chr4_JH584292_random_mrna                 4
# 6003                mm10_UCSC_chr19_mrna             14786
# 6004       mm10_UCSC_chrUn_GL456379_mrna                 2
# 6011                 mm10_UCSC_chr8_mrna             32137
# 6013 mm10_UCSC_chrY_JH584300_random_mrna                 4
# 6014 mm10_UCSC_chrY_JH584303_random_mrna                 3
# 6018 mm10_UCSC_chr1_GL456212_random_mrna               203
# 6022 mm10_UCSC_chr4_GL456350_random_mrna               924
# 6561       mm10_UCSC_chrUn_GL456381_mrna                 0
# 6562       mm10_UCSC_chrUn_GL456387_mrna                 0
# 6569       mm10_UCSC_chrUn_GL456396_mrna                 0
# 6576       mm10_UCSC_chrUn_GL456390_mrna                 0
# 6588       mm10_UCSC_chrUn_GL456389_mrna                 0

# note that the size of mm10_UCSC_all_mrna is 610553
#  which is the same as the sum of the chromosome-specific databases
sum(dbs[dbsplLen>3 & !hsvec & grepl("mm.*mrna",dbs$AnnotationSetName), 7])
# [1] 1196571
sum(dbs[dbsplLen>3 & !hsvec & grepl("mm.*mrna",dbs$AnnotationSetName), 7]) - 610553
# [1] 586018
sum(dbs[dbsplLen>3 & !hsvec & grepl("mm.*mrna",dbs$AnnotationSetName) & !grepl("all",dbs$AnnotationSetName), 7])
# [1] 586018

# conclusion
# we could ignore chromosome-specific databases where dbsplLen>3 


# ------------------------------------------------------------------------------
# explore database sources
# ------------------------------------------------------------------------------
# the second part of AnnotationSetName describes source
dbspl2 = sapply(dbspl, function(f) f[2])

# all of the non-human databases are from UCSC, as are vast majority of human
tapply(dbspl2, hsvec, table)
# $`FALSE`
# 
# UCSC 
# 1041 
# 
# $`TRUE`
# 
# acembly      CADD   ClinVar    dbNSFP     dbSNP      ExAC    gnomAD      UCSC wellderly 
#       1         2         4         1         4         2         2      5951         1 

# look at non-UCSC human databases
hs_nonUCSC = grep("UCSC", names(tapply(dbspl2, hsvec, table)$"TRUE"), 
                  invert=TRUE, value=TRUE)

dbspl3 = sapply(dbspl, function(f) f[3])
dbspl4 = sapply(dbspl, function(f) f[4])
dbspl5 = sapply(dbspl, function(f) f[5])
dbspl6 = sapply(dbspl, function(f) f[6])

sort(table(dbspl1))

dbsplvecs = lapply(as.list(1:max(sapply(dbspl, length))), 
                   function(f) sapply(dbspl, 
                                      function(ff) ff[f]))

dbspl1_2 = sapply(dbspl, function(f) paste0(f[1],"_",f[2]))











spl_fields = strsplit(dbs$AnnotationSetFields,",")
all_fields = sort(unique(unlist(spl_fields)))

field_df = as.data.frame(matrix(0, nrow=nrow(dbs), ncol=length(all_fields),
                                dimnames=list(dbs$AnnotationSetName, all_fields)))

for (s in 1:nrow(field_df)) {
    field_df[s, match(spl_fields[[s]], names(field_df))] = 1
}; rm(s)

