setwd("~/Documents/lipids/")
options(stringsAsFactors=FALSE)
rm(list=ls())

# ------------------------------------------------------------------------------
# -- load and prepare PAGE data
# ------------------------------------------------------------------------------
# load data and add column to match format of MVP snps
snps_page = read.csv("PAGE_materials/PAGE_Lipids_Replication_Table_updated.csv")
dim(snps_page)
# [1] 54  4

snps_page$MarkerName = paste0(snps_page$CHR, ":", snps_page$POS_hg19)

# load proxy snps and add columns to match format of MVP snps
snps_page_proxies = read.csv("PAGE_materials/PAGE_Lipids_Replication_Proxies_updated.csv")
dim(snps_page_proxies)
# [1] 716   8

# order by trait and index SNP position 
xord = order(snps_page_proxies$Trait.s., 
             snps_page_proxies$Chr_indexSNP,
             snps_page_proxies$Pos_indexSNP)
snps_page_proxies = snps_page_proxies[xord, ]
rm(xord)

snps_page_proxies$Index.MarkerName = paste0(snps_page_proxies$Chr_indexSNP, ":", 
                                            snps_page_proxies$Pos_indexSNP)
snps_page_proxies$Proxy.MarkerName = paste0(snps_page_proxies$Chr_proxy, ":", 
                                            snps_page_proxies$Pos_proxy)

# get names of all PAGE snps, including proxies
c1 = all(snps_page$MarkerName %in% 
         snps_page_proxies$Index.MarkerName)
c2 = all(snps_page_proxies$Index.MarkerName %in% 
         snps_page_proxies$Proxy.MarkerName)
if (c1 && c2) {
    all_snps_page = unique(c(snps_page_proxies$Proxy.MarkerName))
}
rm(c1, c2)

length(all_snps_page)
# [1] 700

# ------------------------------------------------------------------------------                                       
# -- MVP - load and filter ethnicity-specific data
# ------------------------------------------------------------------------------
# get names of files downloaded from Box:
#   MVP GWAS/Lipids GWAS/1000G_GWAS_OUT_GENESIS_KLARIN/
mvp_files = grep("^allchr.*gz$", list.files("MVP_materials/"), value=T)
mvp_files
# [1] "allchr..EUR.LDL.results_10_21_17.gz" "allchr.AFR.HDL.results_10_31_17.gz"  "allchr.AFR.LDL.results_10_31_17.gz"  "allchr.AFR.TC.results_10_31_17.gz"   "allchr.AFR.TG.results_10_31_17.gz"  
# [6] "allchr.EUR.HDL.results_10_31_17.gz"  "allchr.EUR.TC.results_10_31_17.gz"   "allchr.EUR.TG.results_10_31_17.gz"   "allchr.HIS.HDL.results_10_31_17.gz"  "allchr.HIS.LDL.results_10_31_17.gz" 
# [11] "allchr.HIS.TC.results_10_31_17.gz"   "allchr.HIS.TG.results_10_31_17.gz"

# loop through files and keep rows for SNPs in the PAGE results (incl. proxies)
#  will probably take 10-20min per file
for (file in mvp_files) {
    # list to hold results for current file
    xlist = list()
    
    # print time to console, initialize line counter for printing updates
    cat(paste0("---- ", Sys.time(), "\n", file), "...\n")
    line_num = 0
    
    # read current file one line at a time
    for (xline in readLines(paste0("MVP_materials/", file))) {
        # print line number every 1e6 lines
        if (line_num %% 1e6 == 0) 
            cat(line_num, "... ")
        
        # split the line based on spaces and grab SNP id
        xline_split = strsplit(xline, " ")[[1]]
        snp = xline_split[1]
        
        # if SNP is in the PAGE SNPs, add the line to the results list
        if (snp %in% all_snps_page) 
            xlist[[snp]] = xline_split
        
        # iterate line counter
        line_num = line_num + 1
    }
    # sanity check
    if (length(xlist) == 0) {
        warning(paste0("No PAGE SNPs found in ", file))
        next
    }
    
    # collapse to data frame to create filtered version of input file
    xdf = as.data.frame(do.call(rbind, xlist))
    
    # set dimnames, get colnames from first line (header) of input file
    rownames(xdf) = NULL
    names(xdf) = strsplit(readLines(paste0("MVP_materials/", file), n=1), " ")[[1]]
    
    # save filtered mvp results to rds file and clean up
    saveRDS(xdf, file=paste0("MVP_materials/", file, "_PAGE_updated.rds"))
    cat("\n")
    rm(line_num, xline, xline_split, snp, xdf)
    gc()
}
rm(file)
# 
# ------------------------------------------------------------------------------                                       
# ---- 2018-09-11 13:21:23
# allchr..EUR.LDL.results_10_21_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 
# ---- 2018-09-11 13:34:27
# allchr.AFR.HDL.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-11 13:54:08
# allchr.AFR.LDL.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-11 14:13:23
# allchr.AFR.TC.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-11 14:32:47
# allchr.AFR.TG.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-11 14:52:38
# allchr.EUR.HDL.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 
# ---- 2018-09-11 15:04:13
# allchr.EUR.TC.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 
# ---- 2018-09-11 15:15:58
# allchr.EUR.TG.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 
# ---- 2018-09-11 15:29:25
# allchr.HIS.HDL.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-11 15:49:00
# allchr.HIS.LDL.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-11 16:08:39
# allchr.HIS.TC.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-11 16:28:54
# allchr.HIS.TG.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# There were 14 warnings (use warnings() to see them)
# rm(file)
# warnings()
# Warning messages:
#     1: In readLines(file) : line 1799228 appears to contain an embedded nul
# 2: In readLines(file) : line 2768007 appears to contain an embedded nul
# 3: In readLines(file) : line 3364225 appears to contain an embedded nul
# 4: In readLines(file) : line 3501819 appears to contain an embedded nul
# 5: In readLines(file) : line 4298882 appears to contain an embedded nul
# 6: In readLines(file) : line 8072032 appears to contain an embedded nul
# 7: In readLines(file) : line 9902711 appears to contain an embedded nul
# 8: In readLines(file) : line 10729089 appears to contain an embedded nul
# 9: In readLines(file) : line 10947781 appears to contain an embedded nul
# 10: In readLines(file) : line 12443547 appears to contain an embedded nul
# 11: In readLines(file) : line 16879537 appears to contain an embedded nul
# 12: In readLines(file) : line 17015662 appears to contain an embedded nul
# 13: In readLines(file) : line 17986523 appears to contain an embedded nul
# 14: In readLines(file) : line 18418242 appears to contain an embedded nul

# ------------------------------------------------------------------------------                                       
# The warnings have to do with lines that are incomplete in "allchr.EUR.TC.results_10_31_17.gz"
# Nothing is affected since none of those lines represent PAGE snps
# ------------------------------------------------------------------------------                                       

# ------------------------------------------------------------------------------                                       
# -- MVP - match and merge with primary PAGE snps and group results by trait
# ------------------------------------------------------------------------------
# get filenames for filtered mvp files
mvp_page_files = paste0(mvp_files, "_PAGE_updated.rds")
cbind(mvp_files, mvp_page_files)
#       mvp_files                             mvp_page_files                                
#  [1,] "allchr..EUR.LDL.results_10_21_17.gz" "allchr..EUR.LDL.results_10_21_17.gz_PAGE_updated.rds"
#  [2,] "allchr.AFR.HDL.results_10_31_17.gz"  "allchr.AFR.HDL.results_10_31_17.gz_PAGE_updated.rds" 
#  [3,] "allchr.AFR.LDL.results_10_31_17.gz"  "allchr.AFR.LDL.results_10_31_17.gz_PAGE_updated.rds" 
#  [4,] "allchr.AFR.TC.results_10_31_17.gz"   "allchr.AFR.TC.results_10_31_17.gz_PAGE_updated.rds"  
#  [5,] "allchr.AFR.TG.results_10_31_17.gz"   "allchr.AFR.TG.results_10_31_17.gz_PAGE_updated.rds"  
#  [6,] "allchr.EUR.HDL.results_10_31_17.gz"  "allchr.EUR.HDL.results_10_31_17.gz_PAGE_updated.rds" 
#  [7,] "allchr.EUR.TC.results_10_31_17.gz"   "allchr.EUR.TC.results_10_31_17.gz_PAGE_updated.rds"  
#  [8,] "allchr.EUR.TG.results_10_31_17.gz"   "allchr.EUR.TG.results_10_31_17.gz_PAGE_updated.rds"  
#  [9,] "allchr.HIS.HDL.results_10_31_17.gz"  "allchr.HIS.HDL.results_10_31_17.gz_PAGE_updated.rds" 
# [10,] "allchr.HIS.LDL.results_10_31_17.gz"  "allchr.HIS.LDL.results_10_31_17.gz_PAGE_updated.rds" 
# [11,] "allchr.HIS.TC.results_10_31_17.gz"   "allchr.HIS.TC.results_10_31_17.gz_PAGE_updated.rds"  
# [12,] "allchr.HIS.TG.results_10_31_17.gz"   "allchr.HIS.TG.results_10_31_17.gz_PAGE_updated.rds"  

# reload MVP data matched to PAGE snps (including proxies)
if (all(mvp_page_files %in% list.files("MVP_materials/"))) {
    xlist = lapply(as.list(paste0("MVP_materials/", mvp_page_files)), readRDS)
    names(xlist) = mvp_page_files
}
t(sapply(xlist, dim))
#                                                      [,1] [,2]
# allchr..EUR.LDL.results_10_21_17.gz_PAGE_updated.rds  656   16
# allchr.AFR.HDL.results_10_31_17.gz_PAGE_updated.rds   690   16
# allchr.AFR.LDL.results_10_31_17.gz_PAGE_updated.rds   690   16
# allchr.AFR.TC.results_10_31_17.gz_PAGE_updated.rds    690   16
# allchr.AFR.TG.results_10_31_17.gz_PAGE_updated.rds    689   16
# allchr.EUR.HDL.results_10_31_17.gz_PAGE_updated.rds   656   16
# allchr.EUR.TC.results_10_31_17.gz_PAGE_updated.rds    647   16
# allchr.EUR.TG.results_10_31_17.gz_PAGE_updated.rds    655   16
# allchr.HIS.HDL.results_10_31_17.gz_PAGE_updated.rds   692   16
# allchr.HIS.LDL.results_10_31_17.gz_PAGE_updated.rds   696   16
# allchr.HIS.TC.results_10_31_17.gz_PAGE_updated.rds    696   16
# allchr.HIS.TG.results_10_31_17.gz_PAGE_updated.rds    696   16

# fix names with two dots
names(xlist) = gsub("..", ".", names(xlist), fixed=T)

# match MVP results to just the primary PAGE snps
.matchToPAGE = function (xdf=NULL, snps_page=NULL) {
    xdf_matched = xdf[match(snps_page$MarkerName, xdf$rsid), ]
    xdf_matched = as.data.frame(cbind(snps_page, xdf_matched[, -c(1:3)]))
    rownames(xdf_matched) = NULL
    return(xdf_matched)
}
xlist_matched = lapply(xlist, 
                       function(f) .matchToPAGE(xdf=f, snps_page=snps_page))
t(sapply(xlist_matched, dim))
#                                             [,1] [,2]
# allchr.EUR.LDL.results_10_21_17.gz_PAGE.rds   54   18
# allchr.AFR.HDL.results_10_31_17.gz_PAGE.rds   54   18
# allchr.AFR.LDL.results_10_31_17.gz_PAGE.rds   54   18
# allchr.AFR.TC.results_10_31_17.gz_PAGE.rds    54   18
# allchr.AFR.TG.results_10_31_17.gz_PAGE.rds    54   18
# allchr.EUR.HDL.results_10_31_17.gz_PAGE.rds   54   18
# allchr.EUR.TC.results_10_31_17.gz_PAGE.rds    54   18
# allchr.EUR.TG.results_10_31_17.gz_PAGE.rds    54   18
# allchr.HIS.HDL.results_10_31_17.gz_PAGE.rds   54   18
# allchr.HIS.LDL.results_10_31_17.gz_PAGE.rds   54   18
# allchr.HIS.TC.results_10_31_17.gz_PAGE.rds    54   18
# allchr.HIS.TG.results_10_31_17.gz_PAGE.rds    54   18

# define the trait represented by each dataset
# -- why not just strsplit to get traits from names(xlist)? --
traits = c("HDL", "LDL", "TC", "TG")
xtraits = sapply(as.list(names(xlist)), 
                 function(ff) traits[sapply(as.list(traits), 
                                            function(f) grepl(f, ff))])

cbind(names(xlist), xtraits)
#                                                     xtraits
#  [1,] "allchr.EUR.LDL.results_10_21_17.gz_PAGE.rds" "LDL"  
#  [2,] "allchr.AFR.HDL.results_10_31_17.gz_PAGE.rds" "HDL"  
#  [3,] "allchr.AFR.LDL.results_10_31_17.gz_PAGE.rds" "LDL"  
#  [4,] "allchr.AFR.TC.results_10_31_17.gz_PAGE.rds"  "TC"   
#  [5,] "allchr.AFR.TG.results_10_31_17.gz_PAGE.rds"  "TG"   
#  [6,] "allchr.EUR.HDL.results_10_31_17.gz_PAGE.rds" "HDL"  
#  [7,] "allchr.EUR.TC.results_10_31_17.gz_PAGE.rds"  "TC"   
#  [8,] "allchr.EUR.TG.results_10_31_17.gz_PAGE.rds"  "TG"   
#  [9,] "allchr.HIS.HDL.results_10_31_17.gz_PAGE.rds" "HDL"  
# [10,] "allchr.HIS.LDL.results_10_31_17.gz_PAGE.rds" "LDL"  
# [11,] "allchr.HIS.TC.results_10_31_17.gz_PAGE.rds"  "TC"   
# [12,] "allchr.HIS.TG.results_10_31_17.gz_PAGE.rds"  "TG" 

# group the matched PAGE snps by trait
.groupByTrait = function (trait=NULL, xlist_matched=NULL, xtraits=NULL) {
    xlist_matched_trait = xlist_matched[xtraits == trait]
    xlist_matched_trait = lapply(xlist_matched_trait,
                                 function(f) subset(f, Trait == trait))
    xdf_matched_trait = as.data.frame(do.call(rbind, xlist_matched_trait))
    split_rownames = strsplit(rownames(xdf_matched_trait), ".", fixed=T)
    xdf_matched_trait$race = sapply(split_rownames, function(f) f[2])
    names(xdf_matched_trait)[6:19] = paste0("MVP.",
                                            names(xdf_matched_trait)[6:19])
    xdf_matched_trait = xdf_matched_trait[order(xdf_matched_trait$rsID),
                                          c(1:5, 19, 6:18)]
    rownames(xdf_matched_trait) = NULL
    return(xdf_matched_trait)
}
xlist_matched_traits = lapply(as.list(traits), 
                              function(f) .groupByTrait(trait=f,
                                                        xlist_matched=xlist_matched,
                                                        xtraits=xtraits))
names(xlist_matched_traits) = traits

# check that there are the right number of rows
# should be 3 race/ethnicities per SNP
if (all(sapply(xlist_matched_traits, nrow) == table(snps_page$Trait) * 3)) {
    print("ok")
}

# collapse resulting list to a single data frame
xdf_matched_traits = as.data.frame(do.call(rbind, xlist_matched_traits))
rownames(xdf_matched_traits) = NULL

dim(xdf_matched_traits)
# [1] 162  19

# order based on trait, chromosome, and race
xord = order(xdf_matched_traits$Trait, xdf_matched_traits$MarkerName, xdf_matched_traits$MVP.race)
xdf_matched_traits = xdf_matched_traits[xord, ]
rm(xord)

# look at distributions of race and traits
table(subset(xdf_matched_traits, is.na(MVP.alleleA))$MVP.race)
# AFR EUR HIS 
#   2  17   3 

table(subset(xdf_matched_traits, is.na(MVP.alleleA))$Trait)
# HDL LDL  TC  TG 
#  14   2   3   3 

sort(table(apply(xdf_matched_traits[,c("Trait","MVP.race")], 
                 1, paste0, collapse=":")[is.na(xdf_matched_traits$MVP.alleleA)]))
# HDL:AFR LDL:EUR HDL:HIS  TC:EUR  TG:EUR HDL:EUR 
#       2       2       3       3       3       9

# plot histograms of MVP p-values for each trait
par(mfrow=c(2,2))
xcuts = c(0, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1)
for (tr in c("HDL","TC","LDL","TG")) {
    toplot = xdf_matched_traits$MVP.frequentist_add_pvalue[xdf_matched_traits$Trait == tr]
    barplot(table(cut(as.numeric(toplot), xcuts)),
            main=tr, col='lightgrey')
    abline(h=5, col='darkgrey')
}

# ------------------------------------------------------------------------------
# -- get proxies for PAGE snps missing from MVP
# ------------------------------------------------------------------------------

.getProxyForSnp = function (missing_snp=NULL, miss_snps_proxies=NULL, 
                            xtraits=NULL, xlist=NULL) {
    # get names of proxy snps for the missing snp
    snp_rows = miss_snps_proxies$Index.MarkerName == missing_snp
    prox = subset(miss_snps_proxies, snp_rows)$Proxy.MarkerName
    
    # get rid of entry that's same as missing snp
    prox = prox[prox != missing_snp]
    
    # get the trait for missing snp
    trait = unique(subset(miss_snps_proxies, snp_rows)$Trait.s.)
    
    # check if missing snp was significant for >1 trait in PAGE
    if (grepl("_", trait)) {
        trait = unlist(strsplit(trait, "_"))
    }
    
    # get the MVP results for the trait and loop through races
    mvp_trait = xlist[xtraits %in% trait]
    
    # empty list to hold results from loop
    mvp_proxies = list()
    for (race_eth in 1:length(mvp_trait)) {
        
        # get results and name for this race
        xmvp = mvp_trait[[race_eth]]
        xmvp_name = names(mvp_trait)[race_eth]
        
        # subset to rows for proxy snps
        xmvp = xmvp[match(prox, xmvp$rsid), ]
        
        # add back names for any missing proxy snps
        xmvp$rsid = prox
        # if (!all(is.na(xmvp))) {
        #     NArows = apply(xmvp, 1, function(f) all(is.na(f)))
        #     xmvp = xmvp[!NArows, ]
        # }
        
        # add subsetted MVP results to list
        mvp_proxies[[xmvp_name]] = xmvp
    }
    # collapse list to one data frame and split rownames
    mvp_proxies = as.data.frame(do.call(rbind, mvp_proxies))
    spl_rownames = strsplit(rownames(mvp_proxies), "\\.")
    
    # add columns for trait, index snp, and race
    mvp_proxies$Trait = sapply(spl_rownames, function(f) f[3])
    mvp_proxies$Index.MarkerName = rep(missing_snp, nrow(mvp_proxies))
    mvp_proxies$MVP.race = sapply(spl_rownames, function(f) f[2])
    
    # reorder columns and rename those with results for proxies
    mvp_proxies = mvp_proxies[, c(17:19, 1:16)]
    names(mvp_proxies)[-c(1:3)] = paste0("Proxy.MVP.",
                                         names(mvp_proxies)[-c(1:3)])
    # clear the rownames
    rownames(mvp_proxies) = NULL
    
    # return list with proxy names, trait, and MVP results for proxies
    return(list(Proxy.MarkerName=prox, Trait=trait, MVP=mvp_proxies))
}

# get names of missing snps and their proxies
miss_snps = unique(subset(xdf_matched_traits, is.na(MVP.alleleA))$MarkerName)
miss_snps_proxies = subset(snps_page_proxies, Index.MarkerName %in% miss_snps)

# get MVP results for proxy snps
mvp_proxies = lapply(as.list(miss_snps), 
                     function(f) .getProxyForSnp(missing_snp=f, 
                                                 miss_snps_proxies=miss_snps_proxies,
                                                 xtraits=xtraits,
                                                 xlist=xlist))
# name with missing snps
names(mvp_proxies) = miss_snps

# loop through results 
for (m in 1:length(mvp_proxies)) {
    # check whether missing snp had any proxies besides itself
    if (nrow(mvp_proxies[[m]]$MVP) == 0) {
        
        # if yes create a 1 row data frame with NAs
        mvp_proxies[[m]]$MVP = mvp_proxies[[m]]$MVP[1, ]
        mvp_proxies[[m]]$MVP$Trait = mvp_proxies[[m]]$Trait
        mvp_proxies[[m]]$MVP$Index.MarkerName = names(mvp_proxies)[m]
        
        # get original results for missing snp
        x = subset(xdf_matched_traits, MarkerName == names(mvp_proxies)[m])
        
        # see which races it was missing from
        xNA = is.na(x$MVP.alleleA)
        
        # add race information
        if (sum(xNA) > 1) {
            mvp_proxies[[m]]$MVP = mvp_proxies[[m]]$MVP[1:sum(xNA), ]
            mvp_proxies[[m]]$MVP[, 1:3] = x[, c(1,5,6)]
        } else {
            mvp_proxies[[m]]$MVP$MVP.race = x$MVP.race[xNA]
        }
    }
}
rm(m, x, xNA)

# collapse to one data frame
mvp_proxies_df = as.data.frame(do.call(rbind, lapply(mvp_proxies,
                                                     function(f) f$MVP)))
rownames(mvp_proxies_df) = NULL

# make unique ids for each missing snp
flt = subset(xdf_matched_traits, MarkerName %in% miss_snps & is.na(MVP.alleleA))
flt = apply(flt[, c("Trait", "MarkerName", "MVP.race")], 1,
            paste0, collapse="_")

# create same ids for proxy results 
tocheck = apply(mvp_proxies_df[, c("Trait", "Index.MarkerName", "MVP.race")], 1,
                paste0, collapse="_")

# filter proxy results to only relevatn trait/race combinations
mvp_proxies_df = mvp_proxies_df[tocheck %in% flt, ]

# plot histograms of proxy MVP p-values for each trait
par(mfrow=c(2,2))
xcuts = c(0, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1)
for (tr in c("HDL","TC","LDL","TG")) {
    toplot = mvp_proxies_df$Proxy.MVP.frequentist_add_pvalue[mvp_proxies_df$Trait == tr]
    barplot(table(cut(as.numeric(toplot), xcuts)),
            main=tr, col='lightgrey')
    abline(h=5, col='darkgrey')
}
    
# ------------------------------------------------------------------------------                                       
# -- write output, then save workspace with session info
# ------------------------------------------------------------------------------
# MVP results for primary PAGE snps
write.csv(xdf_matched_traits, row.names=F,
          file="PAGE_Lipids_Replication_MVP.ethnic-specific_noProxies_updated_new.csv")

# MVP results for proxies of missing PAGE snps
write.csv(mvp_proxies_df, row.names=F,
          file="PAGE_Lipids_Replication_MVP.ethnic-specific_Proxies_for_missing_updated_new.csv")


# get session info
session_info = sessionInfo()

# workspace 
save(list=ls(),
     file="lipids_replication_ethnic-specific_WORKSPACE_updated_new.RData")

