# Oct. 26, 2018 - Austin Hilliard
# ------------------------------------------------------------------------------
# This is code and console output from an interactive R session examining the spreaedsheet here:
# https://docs.google.com/spreadsheets/d/1HRZGbS_LJt2LbrCXNJgFp2v9YgxCktTAgooc3ptXOwk/edit#gid=250069010
#
# The spreadsheet was shared with me by Amir Bahmani. 
# It contains information about annotation databases currently used by his AnnotationHive tool.
# There are 7009 databases represented in the table.
#
# My goal was to look for patterns that would facilitate some semi-automated filtering and reorganization, 
#  so that the information is better organized for future AnnotationHive users.
# It's probably easiest to use this file in RStudio, or any editor that will collapse lines starting with # ----
#
# The major unfinished element is the investigation of Encode databases towards the end.

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

head(dbs)
#               AnnotationSetName AnnotationSetVersion AnnotationSetType           AnnotationSetFields        CreationDate Build AnnotationSetSize                                                                              Info
# 1         hg19_UCSC_laminB1Lads                    1           generic bin,chrom,chromStart,chromEnd 2018/03/08 17:54:00  hg19              1302          Number of Instances: AutoScaling ExecutionTime (msec): 313710; Keywords=
# 2 hg19_UCSC_hiSeqDepthTopPt5Pct                    1           generic bin,chrom,chromStart,chromEnd 2018/03/08 17:41:44  hg19              1224                     Number of Instances: AutoScaling ExecutionTime (msec): 297817
# 3        hg19_UCSC_crisprRanges                    1           generic bin,chrom,chromStart,chromEnd 2018/03/08 16:55:50  hg19            189126                     Number of Instances: AutoScaling ExecutionTime (msec): 366013
# 4   hg19_UCSC_hiSeqDepthTop5Pct                    1           generic bin,chrom,chromStart,chromEnd 2018/03/08 17:40:45  hg19             16119                     Number of Instances: AutoScaling ExecutionTime (msec): 299294
# 5   hg19_UCSC_hiSeqDepthTop1Pct                    1           generic bin,chrom,chromStart,chromEnd 2018/03/08 17:40:25  hg19              2060                     Number of Instances: AutoScaling ExecutionTime (msec): 309624
# 6  hg19_UCSC_hiSeqDepthTop10Pct                    1           generic bin,chrom,chromStart,chromEnd 2018/03/08 17:40:09  hg19             30671                     Number of Instances: AutoScaling ExecutionTime (msec): 321241

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
# 2. ~80% of all databases are for the hg19 human genome build, another ~6% are hg37 or hg38.
# mouse (mm10, ~7%) and horse (equCab2, ~3%) are the most well-represented other species

# make a vector distinguish human databases from other species
hsvec = ifelse(grepl("^hg|^GRCh", dbspl1), TRUE, FALSE)
table(hsvec)
# hsvec
# FALSE  TRUE 
#  1041  5968

# ------------------------------------------------------------------------------
# explore databases with more/fewer elements in AnnotationSetName than others
# ------------------------------------------------------------------------------
# --do all databases have same number of elements in their names?
dbsplLen = sapply(dbspl, length)
table(dbsplLen)
# dbsplLen
# 2    3    4    5    6 
# 1 6614  268   62   64 

# those where dbsplLen==3 seem to be species_source_featureType
# for example:
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

# the one database with dbsplLen == 2 is apparently the entire hg19 assembly, but genes defined specifically by AceView?
#  https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/makeDb/trackDb/human/hg19/acembly.html
dbs[dbsplLen < 3, 1:4]
#      AnnotationSetName AnnotationSetVersion AnnotationSetType                                                               AnnotationSetFields
# 4562      hg19_acembly                    1           generic bin,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds

# --what species are the databases where dbsplLen > 3?
table(dbspl1[dbsplLen > 3])
# ailMel1 bosTau8 canFam3 cavPor3 equCab2 gorGor5    hg19    hg38 melGal5    mm10 panTro5 papAnu2 
#       2       2       2       2     166       1      13       9       2     191       2       2 

# --most are from non-human species, what do the human ones look like?
dbs[dbsplLen>3 & hsvec, c(1,7)][order(dbs[dbsplLen>3 & hsvec, 1]), ]
#                              AnnotationSetName AnnotationSetSize
# 5139                     hg19_UCSC_all_bacends           2289275
# 5140                         hg19_UCSC_all_est           8686047
# 5135                     hg19_UCSC_all_fosends           1504626
# 5131                        hg19_UCSC_all_mrna           8825218
# 147                hg19_UCSC_lincRNAsCTBrain_R             21626
# 136             hg19_UCSC_lincRNAsCTForeskin_R             21626
# 127                 hg19_UCSC_lincRNAsCThLF_r1             21626
# 139                 hg19_UCSC_lincRNAsCThLF_r2             21626
# 132             hg19_UCSC_lincRNAsCTPlacenta_R             21626
# 146               hg19_UCSC_lincRNAsCTTestes_R             21626
# 6923 hg19_UCSC_snpArrayIlluminaHumanCytoSNP_12            604412
# 6927 hg19_UCSC_snpArrayIlluminaHumanOmni1_Quad           2341005
# 5122                       hg19_UCSC_uniGene_3            124338
# 6969               hg38_dbNSFP_35a_allVariants          83306036
# 5534                         hg38_UCSC_all_est           9802373
# 5545                        hg38_UCSC_all_mrna           9013204
# 5247               hg38_UCSC_lincRNAsCTBrain_R             21532
# 5242            hg38_UCSC_lincRNAsCTForeskin_R             21532
# 5243                hg38_UCSC_lincRNAsCThLF_r1             21532
# 5251                hg38_UCSC_lincRNAsCThLF_r2             21532
# 5234            hg38_UCSC_lincRNAsCTPlacenta_R             21532
# 5253              hg38_UCSC_lincRNAsCTTestes_R             21532

# the extra elements designate something about the featureType 
#  ~half are lincRNAs in different tissues, appended with "R", "r1", or "r2"
#   these only represent 12/46 total human lincRNA databases in the table, e.g...
hslinc = hsvec & grepl("linc", dbs$AnnotationSetName, ignore.case=T)
head(dbs[hslinc, c(1, 7)])
#                      AnnotationSetName AnnotationSetSize
# 126          hg19_UCSC_lincRNAsCTOvary             21626
# 127         hg19_UCSC_lincRNAsCThLF_r1             21626
# 128 hg19_UCSC_lincRNAsCTSkeletalMuscle             21626
# 129         hg19_UCSC_lincRNAsCTBreast             21626
# 130          hg19_UCSC_lincRNAsCTHeart             21626
# 131          hg19_UCSC_lincRNAsCTBrain             21626

#   hg19 and hg38 have 23 entries each
table(dbspl1[hslinc])
# hg19 hg38 
#   23   23 

#   within each genome build there are databases for lincRNAs in different tissues, all the same size,
#    and one just called hg*_UCSC_lincRNAsTranscripts, which is a different size from the tissue-specific ones
table(dbs$AnnotationSetSize[hslinc])
# 21482 21532 21626 21630 
#     1    22    22     1 

# conclusion: the human lincRNA databases can probably be consolidated within each genome build,
#             since they have the same number of entries

# the other human databases where dbsplLen>3 are relatively large, especially dbNSFP
#  they also include Illumina SNP array definitions, all mRNAs, all ESTs, uniGene
#  except for the SNP arrays these seem like databases one would always want to query
dbs[dbsplLen>3 & hsvec & !hslinc, c(1,7)]
#                              AnnotationSetName AnnotationSetSize
# 5122                       hg19_UCSC_uniGene_3            124338
# 5131                        hg19_UCSC_all_mrna           8825218
# 5135                     hg19_UCSC_all_fosends           1504626
# 5139                     hg19_UCSC_all_bacends           2289275
# 5140                         hg19_UCSC_all_est           8686047
# 5534                         hg38_UCSC_all_est           9802373
# 5545                        hg38_UCSC_all_mrna           9013204
# 6923 hg19_UCSC_snpArrayIlluminaHumanCytoSNP_12            604412
# 6927 hg19_UCSC_snpArrayIlluminaHumanOmni1_Quad           2341005
# 6969               hg38_dbNSFP_35a_allVariants          83306036

# --what about non-human databases where dbsplLen > 3?
#   they contain same types of information as human non-lincRNA databases above,
#    except separated by chromosome.
#   but for each feature type there's also one database that specifies "all" 
#    e.g. for mrna in mouse:
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

# conclusion: we could ignore chromosome-specific databases where dbsplLen > 3 

# many of these non-human databases are also quite small, and in fact make up ~2/3 of all databases with <100 size
sum(!hsvec & dbsplLen>3 & dbs$AnnotationSetSize<100)
# [1] 106
sum(dbs$AnnotationSetSize<100)
# [1] 158

# conclusion: should look closer at small databases to see where they could be merged
# human databases with AnnotationSetSize < 100 are mostly Encode
subset(dbs, AnnotationSetSize < 100 & hsvec)[, c(1,7)]
#                                                AnnotationSetName AnnotationSetSize
# 31                                        hg19_UCSC_altLocations                 9
# 109                                 hg19_UCSC_altSeqHaplotypesP9                80
# 110                                hg19_UCSC_altSeqHaplotypesP10                80
# 1433     hg19_UCSC_wgEncodeSunyAlbanyGeneStK562SlbpRbpAssocRnaV2                25
# 1548    hg19_UCSC_wgEncodeSunySwitchgearHt1080Igf2bp1RbpAssocRna                24
# 1978 hg19_UCSC_wgEncodeSunyAlbanyGeneStGm12878T7tagRbpAssocRnaV2                 7
# 2003   hg19_UCSC_wgEncodeSunyAlbanyGeneStGm12878T7tagRbpAssocRna                 7
# 2044     hg19_UCSC_wgEncodeSunySwitchgearHt1080Elavl1RbpAssocRna                 5
# 2343       hg19_UCSC_wgEncodeSunyAlbanyGeneStK562SlbpRbpAssocRna                25
# 3234                     hg19_UCSC_wgEncodeSydhTfbsK562Brf1StdPk                74
# 3909                    hg19_UCSC_wgEncodeSydhTfbsK562Xrcc4StdPk                15
# 4005            hg19_UCSC_wgEncodeSydhTfbsHepg2Pgc1aForsklnStdPk                53
# 4086                     hg19_UCSC_wgEncodeSydhTfbsK562Brf2StdPk                19
# 4089                hg19_UCSC_wgEncodeSydhTfbsHelas3Brg1IggmusPk                65
# 4264                     hg19_UCSC_wgEncodeSydhTfbsK562Pol3StdPk                81
# 4467                hg19_UCSC_wgEncodeHaibGenotypeHcmRegionsRep1                81
# 4478              hg19_UCSC_wgEncodeHaibGenotypeImr90RegionsRep1                95
# 4494            hg19_UCSC_wgEncodeHaibGenotypeGm12892RegionsRep2                81
# 4511            hg19_UCSC_wgEncodeHaibGenotypeAg09309RegionsRep1                98
# 4514                hg19_UCSC_wgEncodeHaibGenotypeHcfRegionsRep2                95
# 4518               hg19_UCSC_wgEncodeHaibGenotypeHrpeRegionsRep2                96
# 4538               hg19_UCSC_wgEncodeHaibGenotypeHmecRegionsRep2                89
# 4544            hg19_UCSC_wgEncodeHaibGenotypeAg09319RegionsRep2                93
# 4576                                 hg19_UCSC_iscaCuratedBenign                28
# 4892                                        hg19_UCSC_ntOoaHaplo                13
# 4932                                    hg19_UCSC_numtSAssembled                78
# 5128                                            hg19_UCSC_ntMito                 1
# 5178                             hg38_UCSC_hg38Patch11Haplotypes                59
# 5182                                hg38_UCSC_hg38Patch11Patches                64
# 5412                                 hg38_UCSC_iscaCuratedBenign                36
# 5547                               hg38_UCSC_altSeqLiftOverPslP3                66
# 6954                                        hg19_UCSC_snp147Mult                56
# 6956                                        hg19_UCSC_snp146Mult                44
# 6960                                        hg19_UCSC_snp150Mult                 6
# 6985                                        hg38_UCSC_snp150Mult                28
# 7003                                         GRCh37_ClinVar_papu                67
# 7005                                         GRCh38_ClinVar_papu                78

# many of these don't seem very useful, and are maybe more important for technical considerations
# for example, 
#  ClinVar_papu databases contain human variations found in the pseudoautosomal region (PAR), alternate loci, 
#   patch sequences, and unlocalized/unplaced contigs (papu)
#  https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/
#
#  hg19_UCSC_snp146Mult contains only SNPs that have been mapped to multiple locations in the reference genome assembly
#  http://www.noncode.org/cgi-bin/hgTables?hgsid=1197790&hgta_doSchemaDb=hg38&hgta_doSchemaTable=snp146Mult
#
#  Some of the Encode databases are hyper-specific and related to transcription factor binding or specific cell lines
#  e.g. hg19_UCSC_wgEncodeHaibGenotype has CNVs and SNPs from a neuroblastoma cell line, 
#   where cells were differentiated with retinoic acid
#  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM999325

# ------------------------------------------------------------------------------
# explore database sources and AnnotationSetType
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

# --what are the non-UCSC human databases?
hs_nonUCSC = grep("UCSC", names(tapply(dbspl2, hsvec, table)$"TRUE"), 
                  invert=TRUE, value=TRUE)

dbs[dbspl2 %in% hs_nonUCSC, c(1,3,7)]
#                             AnnotationSetName AnnotationSetType AnnotationSetSize
# 4562                             hg19_acembly           generic            259440
# 5026 hg19_ExAC_genericAnnotationForwebCleaned           generic             18225
# 6918                       hg19_wellderly_all           variant         185134780
# 6919                      hg19_gnomAD_genomes           variant         241056551
# 6920                       hg19_gnomAD_exomes           variant          15014473
# 6921                    hg19_ExAC_sitesVepVcf           variant             11781
# 6965                       hg37_dbSNP_allPAPU           variant           8743105
# 6966                           hg37_dbSNP_all           variant         323138224
# 6967                           hg38_dbSNP_all           variant         325174796
# 6968                       hg38_dbSNP_allPAPU           variant          12571905
# 6969              hg38_dbNSFP_35a_allVariants           variant          83306036
# 7002                   GRCh37_ClinVar_variant           variant            392746
# 7003                      GRCh37_ClinVar_papu           variant                67
# 7004                   GRCh38_ClinVar_variant           variant            393045
# 7005                      GRCh38_ClinVar_papu           variant                78
# 7008                         hg19_CADD_InDels           variant          19995943
# 7009                hg19_CADD_wholeGenomeSNVs           variant        8594355672

# these are mostly databases everyone would probably want to query

# --what is the significance of the AnnotationSetType column?
#  most of these non-UCSC human databases have AnnotationSetType == "variant", but only 92 do overall
table(dbs$AnnotationSetType)
# generic variant 
#    6917      92 

#   most of the 92 are human
table(dbspl1[dbs$AnnotationSetType=="variant"])
# bosTau8  GRCh37  GRCh38    hg19    hg37    hg38    mm10 
#       2       2       2      49       2      27       8 

# --what are the non-human ones?
#   different builds of dbSNP for mouse and bull
dbs[dbs$AnnotationSetType=="variant" & !hsvec, c(1,3,7)]
#            AnnotationSetName AnnotationSetType AnnotationSetSize
# 6994    mm10_UCSC_snp138Mult           variant           6156508
# 6995  mm10_UCSC_snp137Common           variant           5459905
# 6996    mm10_UCSC_snp137Mult           variant           3352724
# 6997        mm10_UCSC_snp138           variant         144689553
# 6998        mm10_UCSC_snp142           variant         165284897
# 6999        mm10_UCSC_snp137           variant         142650642
# 7000    mm10_UCSC_snp142Mult           variant           6592824
# 7001  mm10_UCSC_snp142Common           variant          16557030
# 7006     bosTau8_UCSC_snp148           variant         200690865
# 7007 bosTau8_UCSC_snp148Mult           variant            242285

# that's also what ~50% of the human ones are, in addition to Illumina and Affy SNP array definitions, hapmap, ClinVar,
#  CADD, ExAC, wellderly, and gnomAD
# there are multiple dbSNP build for each genome build
# conclusion: dbSNP databases could be merged, or take only newest for each genome build.
#  can hapmap databases also be merged?
dbs_hsvar = dbs[dbs$AnnotationSetType=="variant" & hsvec, c(1,3,7)]
dbs_hsvar[order(dbs_hsvar$AnnotationSetName),]
#                              AnnotationSetName AnnotationSetType AnnotationSetSize
# 7003                       GRCh37_ClinVar_papu           variant                67
# 7002                    GRCh37_ClinVar_variant           variant            392746
# 7005                       GRCh38_ClinVar_papu           variant                78
# 7004                    GRCh38_ClinVar_variant           variant            393045
# 7008                          hg19_CADD_InDels           variant          19995943
# 7009                 hg19_CADD_wholeGenomeSNVs           variant        8594355672
# 6921                     hg19_ExAC_sitesVepVcf           variant             11781
# 6920                        hg19_gnomAD_exomes           variant          15014473
# 6919                       hg19_gnomAD_genomes           variant         241056551
# 6943              hg19_UCSC_hapmapAllelesChimp           variant           4056757
# 6942            hg19_UCSC_hapmapAllelesMacaque           variant           3749168
# 6944            hg19_UCSC_hapmapAllelesSummary           variant           7933354
# 6945           hg19_UCSC_hapmapPhaseIIISummary           variant           8329934
# 6935                   hg19_UCSC_hapmapSnpsASW           variant           3122130
# 6931                   hg19_UCSC_hapmapSnpsCEU           variant           8059596
# 6934                   hg19_UCSC_hapmapSnpsCHB           variant           8102720
# 6932                   hg19_UCSC_hapmapSnpsCHD           variant           2611760
# 6936                   hg19_UCSC_hapmapSnpsGIH           variant           2815094
# 6933                   hg19_UCSC_hapmapSnpsJPT           variant           8102888
# 6941                   hg19_UCSC_hapmapSnpsLWK           variant           3058764
# 6939                   hg19_UCSC_hapmapSnpsMEX           variant           2819824
# 6937                   hg19_UCSC_hapmapSnpsMKK           variant           3074510
# 6938                   hg19_UCSC_hapmapSnpsTSI           variant           2839124
# 6940                   hg19_UCSC_hapmapSnpsYRI           variant           7966848
# 6950                          hg19_UCSC_snp138           variant         131920678
# 6961                    hg19_UCSC_snp138Common           variant          28293602
# 6953                   hg19_UCSC_snp138Flagged           variant            201278
# 6955                      hg19_UCSC_snp138Mult           variant           7301455
# 6946                    hg19_UCSC_snp141Common           variant          27706766
# 6959                   hg19_UCSC_snp141Flagged           variant            181072
# 6947                      hg19_UCSC_snp142Mult           variant               532
# 6949                          hg19_UCSC_snp144           variant         303845944
# 6958                    hg19_UCSC_snp144Common           variant          30106292
# 6948                   hg19_UCSC_snp144Flagged           variant            281824
# 6957                      hg19_UCSC_snp144Mult           variant               379
# 6964                    hg19_UCSC_snp146Common           variant          30036196
# 6952                   hg19_UCSC_snp146Flagged           variant            306104
# 6956                      hg19_UCSC_snp146Mult           variant                44
# 6951                    hg19_UCSC_snp147Common           variant          30240241
# 6962                   hg19_UCSC_snp147Flagged           variant            319852
# 6954                      hg19_UCSC_snp147Mult           variant                56
# 6963                    hg19_UCSC_snp150Common           variant          30766949
# 6960                      hg19_UCSC_snp150Mult           variant                 6
# 6930              hg19_UCSC_snpArrayAffy250Nsp           variant            514318
# 6928                   hg19_UCSC_snpArrayAffy5           variant            881276
# 6929                   hg19_UCSC_snpArrayAffy6           variant           1818594
# 6922              hg19_UCSC_snpArrayIllumina1M           variant           2436962
# 6925             hg19_UCSC_snpArrayIllumina300           variant            636092
# 6926             hg19_UCSC_snpArrayIllumina550           variant           1121944
# 6924             hg19_UCSC_snpArrayIllumina650           variant           1320776
# 6923 hg19_UCSC_snpArrayIlluminaHumanCytoSNP_12           variant            604412
# 6927 hg19_UCSC_snpArrayIlluminaHumanOmni1_Quad           variant           2341005
# 6918                        hg19_wellderly_all           variant         185134780
# 6966                            hg37_dbSNP_all           variant         323138224
# 6965                        hg37_dbSNP_allPAPU           variant           8743105
# 6969               hg38_dbNSFP_35a_allVariants           variant          83306036
# 6967                            hg38_dbSNP_all           variant         325174796
# 6968                        hg38_dbSNP_allPAPU           variant          12571905
# 6970                          hg38_UCSC_snp141           variant         128845406
# 6981                    hg38_UCSC_snp141Common           variant          29288285
# 6988                   hg38_UCSC_snp141Flagged           variant            273772
# 6982                      hg38_UCSC_snp141Mult           variant              7875
# 6990                          hg38_UCSC_snp142           variant         233053458
# 6991                    hg38_UCSC_snp142Common           variant          27645835
# 6973                   hg38_UCSC_snp142Flagged           variant            287078
# 6976                      hg38_UCSC_snp142Mult           variant               144
# 6986                          hg38_UCSC_snp144           variant         306403960
# 6983                    hg38_UCSC_snp144Common           variant          30620263
# 6989                   hg38_UCSC_snp144Flagged           variant            309669
# 6987                      hg38_UCSC_snp144Mult           variant               643
# 6977                          hg38_UCSC_snp146           variant         313906433
# 6974                    hg38_UCSC_snp146Common           variant          30840214
# 6993                   hg38_UCSC_snp146Flagged           variant            335148
# 6975                      hg38_UCSC_snp146Mult           variant               220
# 6972                          hg38_UCSC_snp147           variant         320763892
# 6980                    hg38_UCSC_snp147Common           variant          31009536
# 6979                   hg38_UCSC_snp147Flagged           variant            349344
# 6971                      hg38_UCSC_snp147Mult           variant               220
# 6978                          hg38_UCSC_snp150           variant         684763755
# 6984                    hg38_UCSC_snp150Common           variant          31405398
# 6992                   hg38_UCSC_snp150Flagged           variant            444734
# 6985                      hg38_UCSC_snp150Mult           variant                28

# ------------------------------------------------------------------------------
# explore database content
# ------------------------------------------------------------------------------
# element 3 of the AnnotationSetName holds the type of feature/content in the database
dbspl3 = sapply(dbspl, function(f) f[3])

# grab the first 4 characters for each element in dbspl3 and count them to get an idea of overall trends
sort(table(substring(dbspl3, 1, 4)))
#  35a acem alle anim asse bacE cent cgap cons darn deni evoC evof exom exon fact fosE hg19 hgdp InDe muPI nean ntHu ntMi ntOo  pdb pept pgAk pgAn 
#    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
# pgCh pgGi pgIr pgKr pgLu pgQu pgSj pgVe pgWa pgYh pgYo pseu rgdQ rgdR rnaC scaf sest site swit targ tfbs uniG vist whol allP altL ccds deco dgvM 
#    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    2    2    2    2    2 
# dgvS eioJ gnfA gwas hgIk illu jaxQ kgPr lami lrgT mafS mm10 netV nsca ntSs papu pgGa pgMd pgNb pgSa pgTk poly sibG sibT snpe vari vega wgRn cris 
#    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    3 
# hinv kgTa know locu mgcF mgcG netL newS numt qual uMas blas cnvD gc5B geno hg38 netB netN oreg pgAb pgKb pubs sgpG ucsf chr8 hiSe netA stra netD 
#    3    3    3    3    3    3    3    3    3    3    3    4    4    4    4    4    4    4    4    4    4    4    4    4    5    5    5    5    6 
# orfe pgHG tRNA intr netE netX ucsc chr6 chr9 chrM estO netF chr7 ensG agil  gap gene mrna ncbi refG snpA chrX gold nest netS refS simp wind altS 
#    6    6    6    7    7    7    7    8    8    8    8    8    9    9   10   10   10   10   10   10   10   11   11   11   11   11   11   11   12 
# augu netR phyl chrY netT rmsk hapm micr netO affy netG netH chr3 netC netP cpgI gens isca chr5 mult burg chr4 phas  all netM xeno clon linc gtex 
#   12   13   13   14   14   14   15   15   15   16   16   16   17   17   17   19   20   22   23   23   24   26   26   27   30   35   36   46   47 
# chr2 chrU snp1 chr1 pgNA chai wgEn 
#   57   67   81  101  137  388 5053 

# the vast majority are from Encode, e.g...
#   will come back to these
head(dbs[grep("^wgEn",dbspl3), c(1, 7)])
#                                                  AnnotationSetName AnnotationSetSize
# 11   hg19_UCSC_wgEncodeCshlShortRnaSeqK562ChromatinShortTransfrags           2090890
# 12  hg19_UCSC_wgEncodeCshlShortRnaSeqGm12878CytosolShortTransfrags             34050
# 13   hg19_UCSC_wgEncodeCshlShortRnaSeqK562NucleolusShortTransfrags            346845
# 14     hg19_UCSC_wgEncodeCshlShortRnaSeqK562CytosolShortTransfrags             71483
# 15 hg19_UCSC_wgEncodeCshlShortRnaSeqK562NucleoplasmShortTransfrags           4473338
# 16     hg19_UCSC_wgEncodeCshlShortRnaSeqK562NucleusShortTransfrags             47576

# --chai
#   388 databases start with "chain", e.g...
head(dbs[grep("^chai",dbspl3), c(1, 7)])
#             AnnotationSetName AnnotationSetSize
# 39 hg19_UCSC_chainGasAcu1Link           7447621
# 40 hg19_UCSC_chainCavPor3Link         240602108
# 41 hg19_UCSC_chainSusScr2Link          72687349
# 42 hg19_UCSC_chainOryCun2Link          75698210
# 43 hg19_UCSC_chainGeoFor1Link           5655878
# 44 hg19_UCSC_chainEquCab2Link         122119319

# these reflect information about cross-species alignments, e.g...
# http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?db=mm10&hgta_group=compGeno&hgta_track=vertebrateChainNet&hgta_table=chainGasAcu1&hgta_doSchema=describe+table+schema
# conclusion: maybe these could be filtered out as they're mostly technical information and not too useful to most biologists

# --pgNA
#   137 databases contain information about variants in personal genomes, e.g...
head(dbs[grep("^pgNA",dbspl3), c(1, 7)])
#            AnnotationSetName AnnotationSetSize
# 566 hg19_UCSC_pgNA12879indel            498439
# 568 hg19_UCSC_pgNA12880indel            490064
# 569 hg19_UCSC_pgNA18558indel            519860
# 571      hg19_UCSC_pgNA06994           3348889
# 572 hg19_UCSC_pgNA21732indel            546092
# 573      hg19_UCSC_pgNA12878           2764382

# see http://www.noncode.org/cgi-bin/hgTables?db=hg19&hgta_group=varRep&hgta_track=pgSnp&hgta_table=pgNA12878&hgta_doSchema=describe+table+schema

# --chr
#   354 databases are chromosome specific (at least based on the third element of AnnotationSetName)
head(dbs[grep("^chr",dbspl3), c(1, 7)])
#                             AnnotationSetName AnnotationSetSize
# 5855 mm10_UCSC_chrY_JH584300_random_intronEst                 5
# 5856 mm10_UCSC_chr1_GL456210_random_intronEst                31
# 5857            mm10_UCSC_chrUn_GL456368_mrna                 2
# 5858       mm10_UCSC_chrUn_GL456372_intronEst                 4
# 5859                       mm10_UCSC_chrM_est             41292
# 5860                mm10_UCSC_chr18_intronEst             44226

# they're only in mouse and horse, and are basically the same from above where dbsplLen > 3
table(dbspl1[grep("^chr",dbspl3)])
# equCab2    mm10 
#     165     189 

table(dbsplLen[grep("^chr",dbspl3)])
#   4   5   6 
# 228  62  64 

# --GTEx has 47 databases
dbs[grepl("gtex", tolower(dbspl3)), c(1,6,7)]
#                               AnnotationSetName Build AnnotationSetSize
# 4559                    hg19_UCSC_gtexGeneModel  hg19             56318
# 5075                  hg19_UCSC_gtexEqtlCluster  hg19           1622952
# 5077     hg19_UCSC_gtexEqtlTissueArteryCoronary  hg19             45901
# 5078    hg19_UCSC_gtexEqtlTissueBrainCerebellum  hg19             83526
# 5079       hg19_UCSC_gtexEqtlTissueAdrenalGland  hg19             64360
# 5080       hg19_UCSC_gtexEqtlTissueArteryTibial  hg19            162090
# 5081      hg19_UCSC_gtexEqtlTissueAdiposeSubcut  hg19            168468
# 5082           hg19_UCSC_gtexEqtlTissueProstate  hg19             26930
# 5083  hg19_UCSC_gtexEqtlTissueEsophagusJunction  hg19             58657
# 5084    hg19_UCSC_gtexEqtlTissueColonTransverse  hg19             90696
# 5085  hg19_UCSC_gtexEqtlTissueEsophagusMuscular  hg19            139428
# 5086        hg19_UCSC_gtexEqtlTissueNerveTibial  hg19            183830
# 5087             hg19_UCSC_gtexEqtlTissueUterus  hg19             15679
# 5088             hg19_UCSC_gtexEqtlTissueVagina  hg19             16223
# 5089    hg19_UCSC_gtexEqtlTissueEsophagusMucosa  hg19            150152
# 5090     hg19_UCSC_gtexEqtlTissueSmallIntestine  hg19             23975
# 5091              hg19_UCSC_gtexEqtlTissueOvary  hg19             24795
# 5092 hg19_UCSC_gtexEqtlTissueXformedfibroblasts  hg19            173370
# 5093  hg19_UCSC_gtexEqtlTissueBrainNucAccumbens  hg19             39369
# 5094             hg19_UCSC_gtexEqtlTissueTestis  hg19            173084
# 5095   hg19_UCSC_gtexEqtlTissueBrainHippocampus  hg19             20792
# 5096    hg19_UCSC_gtexEqtlTissueAdiposeVisceral  hg19             89417
# 5097 hg19_UCSC_gtexEqtlTissueXformedlymphocytes  hg19             58333
# 5098     hg19_UCSC_gtexEqtlTissueSkinNotExposed  hg19            104683
# 5099  hg19_UCSC_gtexEqtlTissueBrainHypothalamus  hg19             22706
# 5100           hg19_UCSC_gtexEqtlTissuePancreas  hg19             82019
# 5101            hg19_UCSC_gtexEqtlTissueStomach  hg19             72346
# 5102             hg19_UCSC_gtexEqtlTissueSpleen  hg19             49766
# 5103        hg19_UCSC_gtexEqtlTissueArteryAorta  hg19            124924
# 5104            hg19_UCSC_gtexEqtlTissueThyroid  hg19            199084
# 5105              hg19_UCSC_gtexEqtlTissueLiver  hg19             29804
# 5106       hg19_UCSC_gtexEqtlTissueColonSigmoid  hg19             59197
# 5107   hg19_UCSC_gtexEqtlTissueBrainFrontCortex  hg19             39336
# 5108        hg19_UCSC_gtexEqtlTissueBrainCortex  hg19             48499
# 5109       hg19_UCSC_gtexEqtlTissueBrainPutamen  hg19             28919
# 5110               hg19_UCSC_gtexEqtlTissueLung  hg19            149881
# 5111    hg19_UCSC_gtexEqtlTissueBreastMamTissue  hg19             83906
# 5112   hg19_UCSC_gtexEqtlTissueBrainAnCinCortex  hg19             22240
# 5113  hg19_UCSC_gtexEqtlTissueHeartLeftVentricl  hg19             92379
# 5114  hg19_UCSC_gtexEqtlTissueHeartAtrialAppend  hg19             79427
# 5115          hg19_UCSC_gtexEqtlTissuePituitary  hg19             42026
# 5116       hg19_UCSC_gtexEqtlTissueBrainCaudate  hg19             47653
# 5117         hg19_UCSC_gtexEqtlTissueWholeBlood  hg19            140265
# 5118     hg19_UCSC_gtexEqtlTissueMuscleSkeletal  hg19            144423
# 5119   hg19_UCSC_gtexEqtlTissueBrainCerebelHemi  hg19             62462
# 5120        hg19_UCSC_gtexEqtlTissueSkinExposed  hg19            176041
# 5394                    hg38_UCSC_gtexGeneModel  hg38             55855

# conclusion: aside from hg38_UCSC_gtexGeneModel, hg19_UCSC_gtexGeneModel, and hg19_UCSC_gtexEqtlCluster,
#  all GTEx data are tissue-specific eQTL data, probably the same as what I describe here:
#  https://docs.google.com/document/d/1pO3hXX6E8K3zZRJgoyqutCd6BLo4iBL3VUPW2lgVLcA/edit#heading=h.sjka4zvqbixr

# ------------------------------------------------------------------------------
# explore Encode databases --- Need to finish this section
# ------------------------------------------------------------------------------

# --what species are repesented in the Encode databases?
table(dbspl1[grepl("^wgEncode", dbspl3)])
# hg19 hg38 mm10 
# 4908  125   20

# make a vector of the Encode database names to explore
vEncode = gsub("wgEncode", "", dbspl3[grep("^wgEn", dbspl3)])

# look at the different 3 letter combinations that start the name
sort(table(substring(vEncode, 1, 3)))
# Dac  Uch  Unc  Uma  UwR  Gis  Aff  Sun  Rik  Gen  Reg  Duk  UwA  UwT  Ope  Bro  Csh  Syd  UwH  UwD  Awg  Hai 
#   1    6    7    8   31   50   51   64   75   80   99  104  173  203  207  286  305  384  409  557  826 1127

vEncspl = split(vEncode, substring(vEncode, 1, 3))

lapply(vEncspl, head, 20)
# these look like RNA-seq on different cellular compartments in different cell lines 
# $Aff
# [1] "AffyRnaChipTransfragsKeratinocyteNucleusLongpolya"        "AffyRnaChipFiltTransfragsKeratinocyteCytosolLongnonpolya"
# [3] "AffyRnaChipTransfragsHepg2NucleusLongpolya"               "AffyRnaChipTransfragsHepg2NucleolusTotal"                
# [5] "AffyRnaChipTransfragsK562CellTotal"                       "AffyRnaChipTransfragsGm12878NucleusLongnonpolya"         
# [7] "AffyRnaChipTransfragsK562CytosolLongpolya"                "AffyRnaChipTransfragsKeratinocyteCytosolLongpolya"       
# [9] "AffyRnaChipTransfragsGm12878CytosolLongpolya"             "AffyRnaChipFiltTransfragsHepg2NucleusLongnonpolya"       
# [11] "AffyRnaChipTransfragsGm12878CellTotal"                    "AffyRnaChipFiltTransfragsHepg2CytosolLongpolya"          
# [13] "AffyRnaChipFiltTransfragsGm12878CytosolLongpolya"         "AffyRnaChipFiltTransfragsHepg2NucleusLongpolya"          
# [15] "AffyRnaChipFiltTransfragsHepg2NucleolusTotal"             "AffyRnaChipTransfragsK562PolysomeLongnonpolya"           
# [17] "AffyRnaChipFiltTransfragsK562PolysomeLongnonpolya"        "AffyRnaChipTransfragsProstateCellLongpolya"              
# [19] "AffyRnaChipFiltTransfragsGm12878CytosolLongnonpolya"      "AffyRnaChipFiltTransfragsGm12878NucleusLongpolya"        
# 
# $Awg
# these look like information on transcription factor binding and DNase hyper-sensitivity sites
#   not sure about Segmentation
# [1] "AwgDnaseMasterSites"            "AwgSegmentationSegwayGm12878"   "AwgSegmentationChromhmmHepg2"   "AwgSegmentationChromhmmHelas3" 
# [5] "AwgSegmentationSegwayHepg2"     "AwgSegmentationSegwayK562"      "AwgSegmentationCombinedK562"    "AwgSegmentationChromhmmHuvec"  
# [9] "AwgSegmentationChromhmmK562"    "AwgSegmentationSegwayH1hesc"    "AwgSegmentationCombinedH1hesc"  "AwgSegmentationCombinedHuvec"  
# [13] "AwgSegmentationSegwayHuvec"     "AwgSegmentationChromhmmH1hesc"  "AwgSegmentationCombinedGm12878" "AwgSegmentationChromhmmGm12878"
# [17] "AwgSegmentationCombinedHelas3"  "AwgSegmentationCombinedHepg2"   "AwgSegmentationSegwayHelas3"    "AwgDnaseUwdukeHsmmtubeUniPk"   
table(substring(gsub("Awg","",vEncspl$Awg), 1, 4))
# Dnas Segm Tfbs 
# 126   18  682 
# 
# $Bro
# these look like cell line-specific histone marks generated by Broad
# [1] "BroadHistoneA549H3k79me2Dex100nmPk"   "BroadHistoneA549H3k04me3Etoh02Pk"     "BroadHistoneCd20CtcfPk"              
# [4] "BroadHistoneH1hescH4k20me1StdPk"      "BroadHistoneNhaH3k36me3StdPk"         "BroadHistoneGm12878H3k27acStdPk"     
# [7] "BroadHistoneK562Chd7a301223a1Pk"      "BroadHistoneNhaH3k79me2Pk"            "BroadHistoneK562Plu1StdPk"           
# [10] "BroadHistoneGm12878H3k27me3StdPkV2"   "BroadHistoneK562Chd4mi2Pk"            "BroadHistoneMonocd14ro1746H3k27acPk" 
# [13] "BroadHistoneMonocd14ro1746H3k09acPk"  "BroadHistoneOsteoblH2azStdPk"         "BroadHistoneA549H2azEtoh02Pk"        
# [16] "BroadHistoneHelas3H3k9acStdPk"        "BroadHistoneA549H3k04me3Dex100nmPk"   "BroadHistoneMonocd14ro1746H4k20me1Pk"
# [19] "BroadHistoneNhlfH3k4me3StdPk"         "BroadHistoneCd20H3k04me2Pk"          
# not sure about these, something to do with hidden markov models
tail(vEncspl$Bro)
# [1] "BroadHmmHepg2HMM"   "BroadHmmNhekHMM"    "BroadHmmHsmmHMM"    "BroadHmmH1hescHMM"  "BroadHmmGm12878HMM" "BroadHmmHmecHMM"  

# 
# $Csh
# [1] "CshlShortRnaSeqK562ChromatinShortTransfrags"     "CshlShortRnaSeqGm12878CytosolShortTransfrags"    "CshlShortRnaSeqK562NucleolusShortTransfrags"    
# [4] "CshlShortRnaSeqK562CytosolShortTransfrags"       "CshlShortRnaSeqK562NucleoplasmShortTransfrags"   "CshlShortRnaSeqK562NucleusShortTransfrags"      
# [7] "CshlShortRnaSeqProstateCellShortTransfrags"      "CshlShortRnaSeqGm12878NucleusShortTransfrags"    "CshlShortRnaSeqK562PolysomeShortTransfrags"     
# [10] "CshlShortRnaSeqK562CellShortTransfrags"          "CshlLongRnaSeqA549CytosolPapContigs"             "CshlLongRnaSeqK562NucleusPamJunctions"          
# [13] "CshlLongRnaSeqA549NucleusPapContigs"             "CshlLongRnaSeqA549CellPamContigs"                "CshlShortRnaSeqHuvecNucleusShorttotalTapContigs"
# [16] "CshlShortRnaSeqA549NucleusContigs"               "CshlLongRnaSeqHepg2CytosolPapJunctions"          "CshlLongRnaSeqBjCellPamJunctions"               
# [19] "CshlLongRnaSeqHepg2NucleusPapJunctions"          "CshlLongRnaSeqH1hescNucleusPapContigs"          
# 
# $Dac
# [1] "DacMapabilityConsensusExcludable"
# 
# $Duk
# [1] "DukeMapabilityRegionsExcludable"         "DukeAffyExonAosmcSimpleSignalRep1"       "DukeAffyExonMedulloSimpleSignalRep3"    
# [4] "DukeAffyExonAstrocySimpleSignalRep1"     "DukeAffyExonIpsSimpleSignalRep3"         "DukeAffyExonLncapAndroSimpleSignalRep2" 
# [7] "DukeAffyExonHsmmtSimpleSignalRep1"       "DukeAffyExonUrotsaUt189SimpleSignalRep1" "DukeAffyExonGm18507SimpleSignalRep3"    
# [10] "DukeAffyExonGm12892SimpleSignalRep2"     "DukeAffyExonStellateSimpleSignalRep1"    "DukeAffyExonH1hescSimpleSignalRep2"     
# [13] "DukeAffyExonGm18507SimpleSignalRep2"     "DukeAffyExonK562SimpleSignalRep2"        "DukeAffyExonK562SimpleSignalRep1"       
# [16] "DukeAffyExonOsteoblSimpleSignalRep2"     "DukeAffyExonLncapSimpleSignalRep2"       "DukeAffyExonMelanoSimpleSignalRep1"     
# [19] "DukeAffyExonProgfibSimpleSignalRep1"     "DukeAffyExonMelanoSimpleSignalRep2"     
# 
# $Gen
# [1] "Gencode2wayConsPseudoV14"   "Gencode2wayConsPseudoV19"   "Gencode2wayConsPseudoV17"   "Gencode2wayConsPseudoV7"    "GencodeCompV17"            
# [6] "GencodePseudoGeneV27lift37" "GencodePseudoGeneV7"        "GencodePolyaV19"            "GencodeBasicV7"             "GencodePolyaV14"           
# [11] "GencodeCompV14"             "GencodeBasicV14"            "GencodeManualV4"            "GencodePolyaV7"             "GencodePolyaV4"            
# [16] "GencodePseudoGeneV24lift37" "GencodePseudoGeneV14"       "Gencode2wayConsPseudoV4"    "GencodeBasicV27lift37"      "GencodeCompV24lift37"      
# 
# $Gis
# [1] "GisRnaPetK562PolysomePapClustersRep1"      "GisRnaPetHepg2NucleusPapClustersRep1"      "GisRnaPetMcf7NucleusPapClusters"          
# [4] "GisChiaPetK562Pol2InteractionsRep2"        "GisRnaPetHuvecCytosolPapClustersRep1"      "GisChiaPetMcf7Pol2InteractionsRep1"       
# [7] "GisChiaPetMcf7EraaInteractionsRep2"        "GisChiaPetMcf7CtcfInteractionsRep2"        "GisRnaPetGm12878CytosolPapClustersRep1"   
# [10] "GisRnaPetHepg2CytosolPapClustersRep1"      "GisRnaPetHuvecCytosolPapClustersRep1V2"    "GisRnaPetGm12878CytosolPapClustersRep1V2" 
# [13] "GisRnaPetH1hescCellPapClustersRep1"        "GisRnaPetHelas3CellPapClustersRep1"        "GisRnaPetA549CytosolPapClusters"          
# [16] "GisRnaPetK562ChromatinTotalClustersRep1"   "GisRnaPetGm12878CytosolPapClustersFree"    "GisChiaPetMcf7EraaInteractionsRep3"       
# [19] "GisChiaPetHelas3Pol2InteractionsRep1"      "GisRnaPetK562NucleoplasmTotalClustersRep1"
# 
# $Hai
# [1] "HaibTfbsGm12891Pax5c20V0416101PkRep2"    "HaibTfbsGm12878Ets1Pcr1xPkRep1V2"        "HaibTfbsMcf7Pmlsc71910V0422111PkRep1"   
# [4] "HaibTfbsEcc1Taf1V0422111PkRep2"          "HaibTfbsH1hescRxraV0416102PkRep2"        "HaibTfbsA549E2f6V0422111PkRep2"         
# [7] "HaibTfbsHepg2SrfV0416101PkRep1"          "HaibTfbsGm12878BatfPcr1xPkRep2"          "HaibTfbsEcc1Yy1sc281V0422111PkRep1"     
# [10] "HaibTfbsHepg2Fosl2V0416101PkRep2"        "HaibTfbsA549Tead4sc101184V0422111PkRep2" "HaibTfbsHct116Cebpbsc150V0422111PkRep2" 
# [13] "HaibTfbsH1hescSrfPcr1xPkRep2"            "HaibTfbsH1hescBcl11aPcr1xPkRep1"         "HaibTfbsGm12892Yy1V0416101PkRep1"       
# [16] "HaibTfbsK562Zbtb33Pcr1xPkRep2"           "HaibTfbsGm12878P300Pcr1xPkRep2"          "HaibTfbsHl60NrsfV0422111PkRep1"         
# [19] "HaibTfbsHepg2Sp2V0422111PkRep2"          "HaibTfbsK562CtcfcPcr1xPkRep1V2"         
# 
# $Ope
# [1] "OpenChromChipGm10248CtcfPkRep1"    "OpenChromFaireHtr8Pk"              "OpenChromDnaseIpsnihi11Pk"         "OpenChromDnaseAdultcd4th1Pk"      
# [5] "OpenChromDnaseCerebrumfrontalocPk" "OpenChromDnaseK562PkV2"            "OpenChromChipProgfibPol2PkRep1"    "OpenChromDnaseFibroblgm03348Pk"   
# [9] "OpenChromDnaseUrothelUt189PkV2"    "OpenChromDnaseHepatocytesPk"       "OpenChromChipMedulloCtcfPkRep1"    "OpenChromChipHuvecCmycPk"         
# [13] "OpenChromDnaseGm12892Pk"           "OpenChromFaireUrotsaPk"            "OpenChromFaireMrtttc549Pk"         "OpenChromDnaseMonocd14Pk"         
# [17] "OpenChromDnaseHsmmPk"              "OpenChromFaireK562NabutPk"         "OpenChromDnaseFibropag20443Pk"     "OpenChromDnaseCd20ro01794Pk"      
# 
# $Reg
# [1] "RegTfbsClusteredV3"               "RegDnaseClusteredV3"              "RegTfbsClusteredV2"               "RegTfbsClustered"                
# [5] "RegDnaseClustered"                "RegDnaseUwNhdfneoPeak"            "RegDnaseUwHepg2Peak"              "RegDnaseUwNb4Peak"               
# [9] "RegDnaseUwHsmmtubePeak"           "RegDnaseUwHmvecdbladPeak"         "RegDnaseUwGm06990Peak"            "RegDnaseUwRpmi7951Peak"          
# [13] "RegDnaseUwAg09309Peak"            "RegDnaseUwLhcnm2Diff4dPeak"       "RegDnaseUwHcmPeak"                "RegDnaseUwH7hescDiffprota14dPeak"
# [17] "RegDnaseUwTh2Peak"                "RegDnaseUwHcfaaPeak"              "RegDnaseUwNhberaPeak"             "RegDnaseUwHbmecPeak"             
# 
# $Rik
# [1] "RikenCageK562NucleusPapTssHmm"         "RikenCageHepg2NucleolusTotalTssHmmV3"  "RikenCageMcf7CytosolPapTssHmm"        
# [4] "RikenCageHobCellPapTssHmm"             "RikenCageHfdpcCellPapTssHmm"           "RikenCageHmscatCellPapTssHmm"         
# [7] "RikenCageSknshraCellPapTssHmm"         "RikenCageK562NucleoplasmTotalTssHmmV3" "RikenCageHchCellPapTssHmm"            
# [10] "RikenCageMcf7NucleusPapTssHmm"         "RikenCageHuvecNucleusPapTssHmm"        "RikenCageHuvecCytosolPapTssHmm"       
# [13] "RikenCageCd34mobilizedCellPapTssHmm"   "RikenCageCd20CellPapTssHmm"            "RikenCageHelas3CytosolPamTssHmmV2"    
# [16] "RikenCageHuvecNucleusPamTssHmm"        "RikenCageH1hescNucleusPapTssHmm"       "RikenCageGm12878NucleusPapTssHmm"     
# [19] "RikenCageHepg2CellPapTssHmmV2"         "RikenCageGm12878CytosolPamTssHmmV2"   
# 
# $Sun
# [1] "SunyAlbanyGeneStGm12878SlbpRbpAssocRna"     "SunyAlbanyGeneStGm12878RipinputRbpAssocRna" "SunyAlbanyGeneStH1hescRipinputRbpAssocRna" 
# [4] "SunyRipSeqGm12878Pabpc1Pk"                  "SunyAlbanyGeneStK562Pabpc1RbpAssocRnaV2"    "SunyRipSeqGm12878Elavl1Pk"                 
# [7] "SunyAlbanyGeneStK562Celf1RbpAssocRna"       "SunyAlbanyGeneStK562T7tagRbpAssocRnaV2"     "SunyAlbanyGeneStK562RipinputRbpAssocRnaV2" 
# [10] "SunyAlbanyGeneStK562T7tagRbpAssocRna"       "SunyRipSeqK562Elavl1Pk"                     "SunyAlbanyTilingGm12878T7tagRbpAssocRna"   
# [13] "SunyAlbanyGeneStK562RipinputRbpAssocRna"    "SunyAlbanyGeneStGm12878Pabpc1RbpAssocRnaV2" "SunyAlbanyGeneStGm12878Elavl1RbpAssocRna"  
# [16] "SunyAlbanyGeneStGm12878Celf1RbpAssocRnaV2"  "SunyAlbanyGeneStHepg2Elavl1RbpAssocRna"     "SunyAlbanyGeneStGm12878Pabpc1RbpAssocRna"  
# [19] "SunyAlbanyGeneStGm12878Igf2bp1RbpAssocRna"  "SunyAlbanyGeneStHelas3Pabpc1RbpAssocRna"   
# 
# $Syd
# [1] "SydhTfbsK562Corestab24166IggrabPk" "SydhTfbsSknshNrf1IggrabPk"         "SydhTfbsK562Stat1Ifng30StdPk"      "SydhTfbsHelas3Gcn5StdPk"          
# [5] "SydhTfbsHelas3Nrf1IggmusPk"        "SydhTfbsHuvecCjunStdPk"            "SydhTfbsHelas3Zzz3StdPk"           "SydhTfbsGm12878Pol2IggmusPk"      
# [9] "SydhTfbsH1hescJundIggrabPk"        "SydhTfbsK562Tf3c110StdPk"          "SydhTfbsH1hescSuz12UcdPk"          "SydhTfbsK562Pol2Ifng6hStdPk"      
# [13] "SydhTfbsHelas3Irf3IggrabPk"        "SydhTfbsH1hescMxi1IggrabPk"        "SydhTfbsHepg2Tcf7l2UcdPk"          "SydhTfbsGm12878Elk112771IggmusPk" 
# [17] "SydhHistoneK562H3k4me1UcdPk"       "SydhTfbsNb4CmycStdPk"              "SydhTfbsHelas3E2f1StdPk"           "SydhTfbsMcf7Hae2f1UcdPk"          
# 
# $Uch
# [1] "UchicagoTfbsK562Enr4a1ControlPk" "UchicagoTfbsK562EfosControlPk"   "UchicagoTfbsK562EjundControlPk"  "UchicagoTfbsK562Ehdac8ControlPk"
# [5] "UchicagoTfbsK562Egata2ControlPk" "UchicagoTfbsK562EjunbControlPk" 
# 
# $Uma
# [1] "UmassDekker5CGm12878PkV2" "UmassDekker5CK562PkV2"    "UmassDekker5CHelas3PkV2"  "UmassDekker5CK562Pk"      "UmassDekker5CH1hescPk"   
# [6] "UmassDekker5CH1hescPkV2"  "UmassDekker5CHelas3Pk"    "UmassDekker5CGm12878Pk"  
# 
# $Unc
# [1] "UncBsuProtGm12878MitoSig"             "UncBsuProtK562MitoSig"                "UncBsuProtGm12878NucleusSig"         
# [4] "UncBsuProtGm12878MembranefractionSig" "UncBsuProtK562NucleusSig"             "UncBsuProtK562MembranefractionSig"   
# [7] "UncBsuProtK562CytosolSig"            
# 
# $UwA
# [1] "UwAffyExonArrayK562SimpleSignalRep1"         "UwAffyExonArrayNhlfSimpleSignalRep2"         "UwAffyExonArrayAg09309SimpleSignalRep1"     
# [4] "UwAffyExonArraySaecSimpleSignalRep1"         "UwAffyExonArrayHepg2SimpleSignalRep2"        "UwAffyExonArrayNhdfneoSimpleSignalRep2"     
# [7] "UwAffyExonArraySknshraSimpleSignalRep1"      "UwAffyExonArrayHmecSimpleSignalRep1"         "UwAffyExonArrayHl60SimpleSignalRep1"        
# [10] "UwAffyExonArrayAg09319SimpleSignalRep1"      "UwAffyExonArrayTh1SimpleSignalRep1"          "UwAffyExonArrayRptecSimpleSignalRep1"       
# [13] "UwAffyExonArrayAg10803SimpleSignalRep2"      "UwAffyExonArrayHek293SimpleSignalRep2"       "UwAffyExonArrayWi38SimpleSignalRep1"        
# [16] "UwAffyExonArrayAg04449SimpleSignalRep2"      "UwAffyExonArrayHbmecSimpleSignalRep2"        "UwAffyExonArrayHrceSimpleSignalRep2"        
# [19] "UwAffyExonArrayK562Znff41b2SimpleSignalRep1" "UwAffyExonArrayHcfSimpleSignalRep2"         
# 
# $UwD
# [1] "UwDnaseGm04504HotspotsRep2"        "UwDgfSaecHotspots"                 "UwDgfGm06990Hotspots"              "UwDgfTh1wb33676984Hotspots"       
# [5] "UwDnaseTh2HotspotsRep1"            "UwDnaseMcf7Est100nm1hHotspotsRep2" "UwDnaseTh1HotspotsRep1"            "UwDnaseHffmycHotspotsRep2"        
# [9] "UwDnaseCmkHotspotsRep1"            "UwDnaseK562Znff41b2HotspotsRep2"   "UwDgfHmvechHotspots"               "UwDnaseMscHotspotsRep1"           
# [13] "UwDgfCd4naivewb11970640Hotspots"   "UwDnaseTh1wb54553204HotspotsRep2"  "UwDnaseK562Znff41b2HotspotsRep1"   "UwDnaseHreHotspotsRep2"           
# [17] "UwDnaseTh2wb54553204HotspotsRep1"  "UwDnaseNhberaHotspotsRep2"         "UwDnaseK562Znfp5HotspotsRep2"      "UwDgfHgfHotspots"                 
# 
# $UwH
# [1] "UwHistoneHacH3k04me3StdHotspotsRep1"            "UwHistoneHreH3k27me3StdHotspotsRep1"            "UwHistoneMonocd14ro1746H3k04me3StdHotspotsRep1"
# [4] "UwHistoneNhdfneoH3k4me3StdHotspotsRep2"         "UwHistoneHct116H3k4me3StdHotspotsRep2"          "UwHistoneHelas3H3k4me3StdHotspotsRep2"         
# [7] "UwHistoneSknmcH3k04me3StdHotspotsRep1"          "UwHistoneRptecH3k04me3StdHotspotsRep1"          "UwHistoneAg04449H3k4me3StdHotspotsRep2"        
# [10] "UwHistoneGm12878H3k36me3StdHotspotsRep2"        "UwHistoneSknshraH3k4me3StdHotspotsRep1"         "UwHistoneAg04450H3k09me3StdHotspotsRep1"       
# [13] "UwHistoneK562H3k04me3StdZnff41b2HotspotsRep1"   "UwHistoneH7esH3k04me3StdDiffa2dHotspotsRep1"    "UwHistoneCd20ro01778H3k04me3StdHotspotsRep2"   
# [16] "UwHistoneSaecH3k27me3StdHotspotsRep2"           "UwHistoneH7esH3k36me3StdDiffa5dHotspotsRep2"    "UwHistoneAg10803H3k4me3StdHotspotsRep1"        
# [19] "UwHistoneSkmcH3k04me3StdHotspotsRep2"           "UwHistoneSaecH3k4me3StdHotspotsRep2"           
# 
# $UwR
# [1] "UwRepliSeqHepg2ValleysRep1"   "UwRepliSeqBjPkRep1"           "UwRepliSeqImr90ValleysRep1"   "UwRepliSeqNhekPkRep1"        
# [5] "UwRepliSeqBjValleysRep2"      "UwRepliSeqGm12801ValleysRep1" "UwRepliSeqHuvecValleysRep1"   "UwRepliSeqBjPkRep2"          
# [9] "UwRepliSeqMcf7PkRep1"         "UwRepliSeqBjValleysRep1"      "UwRepliSeqGm12878PkRep1"      "UwRepliSeqSknshPkRep1"       
# [13] "UwRepliSeqImr90PkRep1"        "UwRepliSeqBg02esPkRep1"       "UwRepliSeqHelas3ValleysRep1"  "UwRepliSeqBg02esValleysRep1" 
# [17] "UwRepliSeqGm06990PkRep1"      "UwRepliSeqGm12878ValleysRep1" "UwRepliSeqGm12812PkRep1"      "UwRepliSeqSknshValleysRep1"  
# 
# $UwT
# [1] "UwTfbsA549CtcfStdHotspotsRep2"    "UwTfbsGm12864CtcfStdHotspotsRep2" "UwTfbsAg09319CtcfStdHotspotsRep1" "UwTfbsHbmecCtcfStdHotspotsRep1"  
# [5] "UwTfbsHelas3CtcfStdHotspotsRep1"  "UwTfbsHl60CtcfStdHotspotsRep1"    "UwTfbsHaspCtcfStdHotspotsRep2"    "UwTfbsHcpeCtcfStdHotspotsRep1"   
# [9] "UwTfbsGm12865CtcfStdHotspotsRep3" "UwTfbsHeeCtcfStdHotspotsRep1"     "UwTfbsHrpeCtcfStdHotspotsRep2"    "UwTfbsGm12801CtcfStdHotspotsRep1"
# [13] "UwTfbsWi38CtcfStdHotspotsRep2"    "UwTfbsK562CtcfStdHotspotsRep1"    "UwTfbsGm12873CtcfStdHotspotsRep1" "UwTfbsHcmCtcfStdHotspotsRep1"    
# [17] "UwTfbsHffCtcfStdHotspotsRep1"     "UwTfbsGm12865CtcfStdHotspotsRep2" "UwTfbsGm12865CtcfStdHotspotsRep1" "UwTfbsHmfCtcfStdHotspotsRep1"



# ------------------------------------------------------------------------------
# ignore below...
# ------------------------------------------------------------------------------
# dbsplvecs = lapply(as.list(1:max(sapply(dbspl, length))), 
#                    function(f) sapply(dbspl, 
#                                       function(ff) ff[f]))
# 
# spl_fields = strsplit(dbs$AnnotationSetFields,",")
# all_fields = sort(unique(unlist(spl_fields)))
# 
# field_df = as.data.frame(matrix(0, nrow=nrow(dbs), ncol=length(all_fields),
#                                 dimnames=list(dbs$AnnotationSetName, all_fields)))
# 
# for (s in 1:nrow(field_df)) {
#     field_df[s, match(spl_fields[[s]], names(field_df))] = 1
# }; rm(s)
# 