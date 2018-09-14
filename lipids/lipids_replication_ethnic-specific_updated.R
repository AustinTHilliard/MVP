setwd("~/Documents/lipids/")
options(stringsAsFactors=FALSE)
rm(list=ls())

# ------------------------------------------------------------------------------
# -- load and prepare PAGE data
# ------------------------------------------------------------------------------
# load data and add column to match format of MVP snps
snps_page = read.csv("PAGE_Lipids_Replication_Table_updated.csv")
dim(snps_page)
# [1] 54  4
head(snps_page)
#   Trait        rsID CHR  POS_hg19
# 1   HDL rs539621506   1 109046381
# 2   HDL  rs61519011   2 126435003
# 3   HDL  rs17102282   5 144103408
# 4   HDL rs373140531   7  89604535
# 5   HDL  rs11782435   8  13536115
# 6   HDL rs145312881   9  89053469

snps_page$MarkerName = paste0(snps_page$CHR, ":", snps_page$POS_hg19)
head(snps_page)
#   Trait        rsID CHR  POS_hg19  MarkerName
# 1   HDL rs539621506   1 109046381 1:109046381
# 2   HDL  rs61519011   2 126435003 2:126435003
# 3   HDL  rs17102282   5 144103408 5:144103408
# 4   HDL rs373140531   7  89604535  7:89604535
# 5   HDL  rs11782435   8  13536115  8:13536115
# 6   HDL rs145312881   9  89053469  9:89053469

# load proxy snps and add columns to match format of MVP snps
snps_page_proxies = read.csv("PAGE_Lipids_Replication_Proxies_updated.csv")
dim(snps_page_proxies)
# [1] 716   8
head(snps_page_proxies)
#   Trait.s. Chr_indexSNP Pos_indexSNP    IndexSNP Chr_proxy Pos_proxy    ProxySNP R2_index_proxy_snps
# 1      HDL            1    109046381 rs539621506         1 109046381 rs539621506            1.000000
# 2      HDL            2    126435003  rs61519011         2 126421304   rs1427309            0.838408
# 3      HDL            2    126435003  rs61519011         2 126421423   rs5834119            0.838497
# 4      HDL            2    126435003  rs61519011         2 126422062  rs62158251            0.858461
# 5      HDL            2    126435003  rs61519011         2 126426024   rs2195073            0.861774
# 6      HDL            2    126435003  rs61519011         2 126429997   rs1427307            0.838749

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
head(snps_page_proxies)
#   Trait.s. Chr_indexSNP Pos_indexSNP    IndexSNP Chr_proxy Pos_proxy    ProxySNP R2_index_proxy_snps Index.MarkerName Proxy.MarkerName
# 1      HDL            1    109046381 rs539621506         1 109046381 rs539621506            1.000000      1:109046381      1:109046381
# 2      HDL            2    126435003  rs61519011         2 126421304   rs1427309            0.838408      2:126435003      2:126421304
# 3      HDL            2    126435003  rs61519011         2 126421423   rs5834119            0.838497      2:126435003      2:126421423
# 4      HDL            2    126435003  rs61519011         2 126422062  rs62158251            0.858461      2:126435003      2:126422062
# 5      HDL            2    126435003  rs61519011         2 126426024   rs2195073            0.861774      2:126435003      2:126426024
# 6      HDL            2    126435003  rs61519011         2 126429997   rs1427307            0.838749      2:126435003      2:126429997


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
head(all_snps_page)
# [1] "1:109046381" "2:126421304" "2:126421423" "2:126422062" "2:126426024" "2:126429997"

# ------------------------------------------------------------------------------                                       
# -- MVP - load and filter ethnicity-specific data
# ------------------------------------------------------------------------------
# get names of files downloaded from Box:
#   MVP GWAS/Lipids GWAS/1000G_GWAS_OUT_GENESIS_KLARIN/
mvp_files = grep("^allchr.*gz$", list.files(), value=T)
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
    for (xline in readLines(file)) {
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
    names(xdf) = strsplit(readLines(file, n=1), " ")[[1]]
    
    # save filtered mvp results to rds file and clean up
    saveRDS(xdf, file=paste0(file, "_PAGE_updated.rds"))
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
# -- MVP - match and merge with primary PAGE snps and group results by trait
# ------------------------------------------------------------------------------
# get filenames for filtered mvp files
mvp_page_files = paste0(mvp_files, "_PAGE.rds")
cbind(mvp_files, mvp_page_files)
#       mvp_files                             mvp_page_files                                
#  [1,] "allchr..EUR.LDL.results_10_21_17.gz" "allchr..EUR.LDL.results_10_21_17.gz_PAGE.rds"
#  [2,] "allchr.AFR.HDL.results_10_31_17.gz"  "allchr.AFR.HDL.results_10_31_17.gz_PAGE.rds" 
#  [3,] "allchr.AFR.LDL.results_10_31_17.gz"  "allchr.AFR.LDL.results_10_31_17.gz_PAGE.rds" 
#  [4,] "allchr.AFR.TC.results_10_31_17.gz"   "allchr.AFR.TC.results_10_31_17.gz_PAGE.rds"  
#  [5,] "allchr.AFR.TG.results_10_31_17.gz"   "allchr.AFR.TG.results_10_31_17.gz_PAGE.rds"  
#  [6,] "allchr.EUR.HDL.results_10_31_17.gz"  "allchr.EUR.HDL.results_10_31_17.gz_PAGE.rds" 
#  [7,] "allchr.EUR.TC.results_10_31_17.gz"   "allchr.EUR.TC.results_10_31_17.gz_PAGE.rds"  
#  [8,] "allchr.EUR.TG.results_10_31_17.gz"   "allchr.EUR.TG.results_10_31_17.gz_PAGE.rds"  
#  [9,] "allchr.HIS.HDL.results_10_31_17.gz"  "allchr.HIS.HDL.results_10_31_17.gz_PAGE.rds" 
# [10,] "allchr.HIS.LDL.results_10_31_17.gz"  "allchr.HIS.LDL.results_10_31_17.gz_PAGE.rds" 
# [11,] "allchr.HIS.TC.results_10_31_17.gz"   "allchr.HIS.TC.results_10_31_17.gz_PAGE.rds"  
# [12,] "allchr.HIS.TG.results_10_31_17.gz"   "allchr.HIS.TG.results_10_31_17.gz_PAGE.rds"  

# reload MVP data matched to PAGE snps (including proxies)
if (all(mvp_page_files %in% list.files())) {
    xlist = lapply(as.list(mvp_page_files), readRDS)
    names(xlist) = mvp_page_files
}
t(sapply(xlist, dim))
#                                              [,1] [,2]
# allchr..EUR.LDL.results_10_21_17.gz_PAGE.rds  554   16
# allchr.AFR.HDL.results_10_31_17.gz_PAGE.rds   569   16
# allchr.AFR.LDL.results_10_31_17.gz_PAGE.rds   569   16
# allchr.AFR.TC.results_10_31_17.gz_PAGE.rds    569   16
# allchr.AFR.TG.results_10_31_17.gz_PAGE.rds    568   16
# allchr.EUR.HDL.results_10_31_17.gz_PAGE.rds   554   16
# allchr.EUR.TC.results_10_31_17.gz_PAGE.rds    545   16
# allchr.EUR.TG.results_10_31_17.gz_PAGE.rds    554   16
# allchr.HIS.HDL.results_10_31_17.gz_PAGE.rds   572   16
# allchr.HIS.LDL.results_10_31_17.gz_PAGE.rds   572   16
# allchr.HIS.TC.results_10_31_17.gz_PAGE.rds    572   16
# allchr.HIS.TG.results_10_31_17.gz_PAGE.rds    572   16

head(xlist[[1]])
#         rsid chromosome position   alleleA alleleB average_maximum_posterior_call     info all_AA  all_AB   all_BB all_total missing_data_proportion cohort_1_hwe frequentist_add_pvalue frequentist_add_beta_1 frequentist_add_se_1
# 1 1:16072621          1 16072621 GGTTCAATC       G                       0.978903 0.954859 118609 82423.6    14163    215196             1.08195e-15     0.334828               0.349722            -0.00334047            0.0035722
# 2 1:16074494          1 16074494         C       T                        0.98257 0.962683 118818 82295.7    14082    215196            -5.08515e-14     0.295931               0.304151            -0.00365917           0.00356098
# 3 1:16075835          1 16075835         C       T                       0.983586 0.964652 118681 82378.7  14136.7    215196            -5.73432e-14     0.323369               0.339895            -0.00339287           0.00355506
# 4 1:16080501          1 16080501         G       A                        0.99853 0.996087 128183 75850.5  11162.2    215196             8.75701e-14     0.672743                  0.404            -0.00303878           0.00364144
# 5 1:16080618          1 16080618         G       A                       0.997748 0.993889 128109 75897.4  11189.3    215196              5.6667e-14     0.709041               0.418954            -0.00294518           0.00364395
# 6 2:17471785          2 17471785         G       A                       0.999779 0.923304 214668 527.752 0.015996    215196            -6.53902e-14            1               0.766482                0.01347            0.0453567

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

head(xlist_matched[[1]])
#   Trait        rsID CHR  POS_hg19  MarkerName alleleA alleleB average_maximum_posterior_call     info  all_AA  all_AB  all_BB all_total missing_data_proportion cohort_1_hwe frequentist_add_pvalue frequentist_add_beta_1 frequentist_add_se_1
# 1   HDL rs539621506   1 109046381 1:109046381    <NA>    <NA>                           <NA>     <NA>    <NA>    <NA>    <NA>      <NA>                    <NA>         <NA>                   <NA>                   <NA>                 <NA>
# 2   HDL  rs61519011   2 126435003 2:126435003       G      GT                       0.949273  0.90945 24256.3 95812.6 95127.1    215196             -4.0573e-15     0.578126               0.422256             0.00271522           0.00338338
# 3   HDL  rs17102282   5 144103408 5:144103408       G       A                       0.997997 0.991909  162430 49051.5 3714.24    215196            -9.00045e-14     0.872317               0.112884             0.00718486           0.00453198
# 4   HDL rs373140531   7  89604535  7:89604535    <NA>    <NA>                           <NA>     <NA>    <NA>    <NA>    <NA>      <NA>                    <NA>         <NA>                   <NA>                   <NA>                 <NA>
# 5   HDL  rs11782435   8  13536115  8:13536115       C       T                       0.954946 0.897483  131388 73596.2 10211.5    215196             1.17662e-14     0.467829               0.689281            -0.00155833           0.00389746
# 6   HDL rs145312881   9  89053469  9:89053469    <NA>    <NA>                           <NA>     <NA>    <NA>    <NA>    <NA>      <NA>                    <NA>         <NA>                   <NA>                   <NA>                 <NA>

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

head(xlist_matched_traits[[1]])
#   Trait        rsID CHR  POS_hg19   MarkerName MVP.race MVP.alleleA MVP.alleleB MVP.average_maximum_posterior_call MVP.info MVP.all_AA MVP.all_AB MVP.all_BB MVP.all_total MVP.missing_data_proportion MVP.cohort_1_hwe MVP.frequentist_add_pvalue MVP.frequentist_add_beta_1 MVP.frequentist_add_se_1
# 1   HDL  rs10430621  10 133667886 10:133667886      AFR        <NA>        <NA>                               <NA>     <NA>       <NA>       <NA>       <NA>          <NA>                        <NA>             <NA>                       <NA>                       <NA>                     <NA>
# 2   HDL  rs10430621  10 133667886 10:133667886      EUR        <NA>        <NA>                               <NA>     <NA>       <NA>       <NA>       <NA>          <NA>                        <NA>             <NA>                       <NA>                       <NA>                     <NA>
# 3   HDL  rs10430621  10 133667886 10:133667886      HIS        <NA>        <NA>                               <NA>     <NA>       <NA>       <NA>       <NA>          <NA>                        <NA>             <NA>                       <NA>                       <NA>                     <NA>
# 4   HDL rs115207757  17  79371257  17:79371257      AFR        <NA>        <NA>                               <NA>     <NA>       <NA>       <NA>       <NA>          <NA>                        <NA>             <NA>                       <NA>                       <NA>                     <NA>
# 5   HDL rs115207757  17  79371257  17:79371257      EUR        <NA>        <NA>                               <NA>     <NA>       <NA>       <NA>       <NA>          <NA>                        <NA>             <NA>                       <NA>                       <NA>                     <NA>
# 6   HDL rs115207757  17  79371257  17:79371257      HIS        <NA>        <NA>                               <NA>     <NA>       <NA>       <NA>       <NA>          <NA>                        <NA>             <NA>                       <NA>                       <NA>                     <NA>

# collapse resulting list to a single data frame
xdf_matched_traits = as.data.frame(do.call(rbind, xlist_matched_traits))
rownames(xdf_matched_traits) = NULL

dim(xdf_matched_traits)
# [1] 162  19

xord = order(xdf_matched_traits$Trait, xdf_matched_traits$CHR, xdf_matched_traits$MVP.race)
xdf_matched_traits = xdf_matched_traits[xord, ]
rm(xord)

table(subset(xdf_matched_traits, is.na(MVP.alleleA))$MVP.race)
# AFR EUR HIS 
#  19  27  19 

table(subset(xdf_matched_traits, is.na(MVP.alleleA))$Trait)
# HDL LDL  TC  TG 
#  21   5  14  25 

sort(table(apply(xdf_matched_traits[,c("Trait","MVP.race")], 
                 1, paste0, collapse=":")[is.na(xdf_matched_traits$MVP.alleleA)]))
# LDL:AFR LDL:HIS LDL:EUR  TC:AFR  TC:HIS HDL:AFR HDL:HIS  TC:EUR  TG:AFR  TG:HIS HDL:EUR  TG:EUR 
#       1       1       3       4       4       6       6       6       8       8       9       9 

xdf_m_spl = split(xdf_matched_traits[,c(1,2,5:7)], xdf_matched_traits$MarkerName)
# # ------------------------------------------------------------------------------                                       
# # -- get proxies for PAGE snps missing from MVP
# # ------------------------------------------------------------------------------
# # assumes miss_snps_proxies, xtraits, and xlist exist in workspace
# .getProxyForSnp = function (missing_snp) {
#     snp_rows = miss_snps_proxies$Index.MarkerName == missing_snp
#     prox = subset(miss_snps_proxies, snp_rows)$Proxy.MarkerName
#     prox = prox[prox != missing_snp]
#     trait = unique(subset(miss_snps_proxies, snp_rows)$Trait.s.)
#     if (grepl("_", trait)) {
#         trait = unlist(strsplit(trait, "_"))
#     }
#     mvp_trait = xlist[xtraits %in% trait]
#     mvp_proxies = list()
#     for (race_eth in 1:length(mvp_trait)) {
#         xmvp = mvp_trait[[race_eth]]
#         xmvp_name = names(mvp_trait)[race_eth]
#         xmvp = xmvp[match(prox, xmvp$rsid), ]
#         xmvp$rsid = prox
#         if (!all(is.na(xmvp))) {
#             NArows = apply(xmvp, 1, function(f) all(is.na(f)))
#             xmvp = xmvp[!NArows, ]
#         }
#         mvp_proxies[[xmvp_name]] = xmvp
#     }
#     mvp_proxies = as.data.frame(do.call(rbind, mvp_proxies))
#     spl_rownames = strsplit(rownames(mvp_proxies), "\\.")
#     
#     mvp_proxies$Trait = sapply(spl_rownames, function(f) f[3])
#     mvp_proxies$Index.MarkerName = rep(missing_snp, nrow(mvp_proxies))
#     mvp_proxies$MVP.race = sapply(spl_rownames, function(f) f[2])
#     
#     mvp_proxies = mvp_proxies[, c(17:19, 1:16)]
#     names(mvp_proxies)[-c(1:3)] = paste0("Proxy.MVP.", 
#                                          names(mvp_proxies)[-c(1:3)])
#     rownames(mvp_proxies) = NULL
#     return(list(Proxy.MarkerName=prox, Trait=trait, MVP=mvp_proxies))
# }
# # get names of missing snps and their proxies
# miss_snps = unique(subset(xdf_matched_traits, is.na(MVP.alleleA))$MarkerName)
# miss_snps_proxies = subset(snps_page_proxies, Index.MarkerName %in% miss_snps)
# 
# #
# mvp_proxies = lapply(as.list(miss_snps), .getProxyForSnp)
# names(mvp_proxies) = miss_snps
# 
# for (m in 1:length(mvp_proxies)) {
#     if (nrow(mvp_proxies[[m]]$MVP) == 0) {
#         mvp_proxies[[m]]$MVP = mvp_proxies[[m]]$MVP[1, ]
#         mvp_proxies[[m]]$MVP$Trait = mvp_proxies[[m]]$Trait
#         mvp_proxies[[m]]$MVP$Index.MarkerName = names(mvp_proxies)[m]
#         
#         x = subset(xdf_matched_traits, MarkerName == names(mvp_proxies)[m])
#         xNA = is.na(x$MVP.alleleA)
#         if (sum(xNA) > 1) {
#             mvp_proxies[[m]]$MVP = mvp_proxies[[m]]$MVP[1:sum(xNA), ]
#             mvp_proxies[[m]]$MVP[, 1:3] = x[, c(1,5,6)]
#         } else {
#             mvp_proxies[[m]]$MVP$MVP.race = x$MVP.race[xNA]
#         }
#     }
# }
# rm(m, x, xNA)
# 
# #
# mvp_proxies_df = as.data.frame(do.call(rbind, lapply(mvp_proxies, 
#                                                      function(f) f$MVP)))
# rownames(mvp_proxies_df) = NULL
# 
# #
# flt = subset(xdf_matched_traits, MarkerName %in% miss_snps & is.na(MVP.alleleA))
# flt = apply(flt[, c("Trait", "MarkerName", "MVP.race")], 1, 
#             paste0, collapse="_")
# 
# tocheck = apply(mvp_proxies_df[, c("Trait", "Index.MarkerName", "MVP.race")], 1,
#                 paste0, collapse="_")
# 
# mvp_proxies_df = mvp_proxies_df[tocheck %in% flt, ]
    
# ------------------------------------------------------------------------------                                       
# -- write output, then save workspace with session info
# ------------------------------------------------------------------------------
# MVP results for primary PAGE snps
write.csv(xdf_matched_traits, row.names=F,
          file="PAGE_Lipids_Replication_MVP.ethnic-specific_noProxies_updated.csv")

# get session info
session_info = sessionInfo()

# workspace 
save(list=ls(),
     file="lipids_replication_ethnic-specific_noProxiesWORKSPACE.RData")

# # MVP results for proxies of missing PAGE snps
# write.csv(mvp_proxies_df, row.names=F,
#           file="PAGE_Lipids_Replication_MVP.ethnic-specific_Proxies_for_missing_updated.csv")

