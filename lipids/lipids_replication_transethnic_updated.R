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

# split snps by trait
xpage = split(snps_page, snps_page$Trait)
names(xpage)
# [1] "HDL" "LDL" "TC"  "TG"
sapply(xpage, nrow)
# HDL LDL  TC  TG 
#  18   8  14  14

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

# ------------------------------------------------------------------------------                                       
# -- MVP - load transethnic data
# ------------------------------------------------------------------------------
# get names of files downloaded from Box:
#   MVP-CAP other local/Manuscripts/MVP lipids GWAS/Summary_Statistics/
mvp_te_files = grep("^MVP.*gz$", list.files(), value=T)
mvp_te_files
# [1] "MVP.te.HDL.gwas.tsv.gz" "MVP.te.LDL.gwas.tsv.gz" "MVP.te.TC.gwas.tsv.gz"  "MVP.te.TG.gwas.tsv.gz" 

# load the data 
#  will take a minute or two for each file since they're each ~570mb
#  if on a machine with not much memory should use readLines approach instead
#   e.g. as in lipids_replication_wthnic-specific.R
for (file in mvp_te_files) {
    cat(file, "...\n")
    assign(file, read.table(file, header=T, as.is=T))
}
rm(file)

gc()
ls()
# [1] "mvp_te_files" 
#     "MVP.te.HDL.gwas.tsv.gz" "MVP.te.LDL.gwas.tsv.gz" "MVP.te.TC.gwas.tsv.gz" "MVP.te.TG.gwas.tsv.gz"
#     "snps_page" "snps_page_proxies" "xpage" 

sapply(as.list(mvp_te_files), function(f) dim(get(f)))
#          [,1]     [,2]     [,3]     [,4]
# [1,] 28946834 29221975 29085659 29184545
# [2,]        9        9        9        9
head(MVP.te.HDL.gwas.tsv.gz)
#    MarkerName Allele1 Allele2  Freq1 FreqSE  Effect StdErr P.value Direction
# 1  5:85928892       t       c 0.1646 0.1075 -0.0021 0.0047 0.66520       -+-
# 2 12:48151416       a       g 0.9982 0.0004  0.1185 0.0784 0.13070       ?+-
# 3  3:26280776       c       g 0.0038 0.0008  0.0414 0.0264 0.11670       +-+
# 4 14:36082010       t       c 0.0553 0.0111 -0.0074 0.0132 0.57600       ?--
# 5 2:170966953       t       c 0.9835 0.0041  0.0127 0.0115 0.27020       +++
# 6  6:65890004       a       c 0.0023 0.0005 -0.1564 0.0874 0.07354       ?-+

# ------------------------------------------------------------------------------                                       
# -- MVP - work on transethnic data
# ------------------------------------------------------------------------------
# make a list to hold the results, then loop through the traits
mvp_te_page = list()
for (trait in names(xpage)) {
    cat(trait, "...\n")
    
    # get the appropriate MVP dataset and subset based on PAGE snps
    x = get(grep(paste0("^MVP.*", trait), ls(), value=T))
    x = subset(x, MarkerName %in% xpage[[trait]]$MarkerName)
    
    # use match to insert NA rows where MVP data is missing a PAGE snp
    x = x[match(xpage[[trait]]$MarkerName, x$MarkerName), ]
    
    # combine the MVP and PAGE information
    names(x) = paste0("MVP.te_", names(x))
    mvp_te_page[[trait]] = cbind(xpage[[trait]], x[, -1])
}
rm(trait, x)
gc()

sapply(mvp_te_page, nrow)
# HDL LDL  TC  TG 
#  18   8  14  14 

# collapse to one data frame
mvp_te_page = as.data.frame(do.call(rbind, mvp_te_page))
rownames(mvp_te_page) = NULL
mvp_te_page
#    Trait        rsID CHR  POS_hg19   MarkerName MVP.te_Allele1 MVP.te_Allele2 MVP.te_Freq1 MVP.te_FreqSE MVP.te_Effect MVP.te_StdErr MVP.te_P.value MVP.te_Direction
# 1    HDL rs539621506   1 109046381  1:109046381              a              t       0.0034        0.0006        0.0155        0.0581      7.901e-01              ?+-
# 2    HDL  rs61519011   2 126435003  2:126435003              d              i       0.3322        0.0169        0.0024        0.0029      4.040e-01              +--
# 3    HDL  rs17102282   5 144103408  5:144103408              a              g       0.2601        0.1632        0.0015        0.0034      6.624e-01              --+
# 4    HDL rs373140531   7  89604535   7:89604535              d              i       0.9977        0.0003       -0.0529        0.0662      4.239e-01              ?--
# 5    HDL  rs11782435   8  13536115   8:13536115              t              c       0.2015        0.0402        0.0052        0.0036      1.508e-01              +++
# 6    HDL rs145312881   9  89053469   9:89053469              a              g       0.0030        0.0006        0.0594        0.0605      3.261e-01              ?+-
# 7    HDL  rs10430621  10 133667886 10:133667886           <NA>           <NA>           NA            NA            NA            NA             NA             <NA>
# 8    HDL  rs73548234  11  44041402  11:44041402              a              t       0.9575        0.0085       -0.0011        0.0149      9.417e-01              ?+-
# 9    HDL  rs13328893  11  50272145  11:50272145              t              c       0.1553        0.0001        0.0442        0.0048      1.360e-20              +?+
# 10   HDL  rs73478487  11  56582285  11:56582285              a              g       0.0805        0.0035        0.0318        0.0053      2.077e-09              +-+
# 11   HDL    rs526276  12 123192454 12:123192454              t              c       0.2230        0.0597        0.0202        0.0033      6.049e-10              ++-
# 12   HDL  rs75405126  14  53797383  14:53797383           <NA>           <NA>           NA            NA            NA            NA             NA             <NA>
# 13   HDL  rs60255124  14 105898715 14:105898715           <NA>           <NA>           NA            NA            NA            NA             NA             <NA>
# 14   HDL rs151291132  15  44842210  15:44842210              a              g       0.9785        0.0014        0.1368        0.0158      5.185e-18              +?+
# 15   HDL rs564672295  16  68685128  16:68685128              a              g       0.0029        0.0005        0.5710        0.0607      5.256e-21              ?++
# 16   HDL  rs12940636  17  53400110  17:53400110              t              c       0.6875        0.0507       -0.0061        0.0030      4.198e-02              ---
# 17   HDL rs115207757  17  79371257  17:79371257              a              g       0.0048        0.0011        0.0128        0.0430      7.665e-01              ?++
# 18   HDL  rs59465323  19   4021239   19:4021239              t              c       0.1055        0.0109        0.0077        0.0043      7.473e-02              +++
# 19   LDL rs141269052   2  17471577   2:17471577              a              g       0.2677        0.1863       -0.0123        0.0037      8.932e-04              --?
# 20   LDL  rs57520889   2  17472023   2:17472023              c              g       0.6466        0.0530        0.0294        0.0061      1.593e-06              ++?
# 21   LDL   rs2637618   4 130597120  4:130597120              a              g       0.9347        0.0138        0.0004        0.0129      9.782e-01              ?+-
# 22   LDL rs534316532   6  30243235   6:30243235              t              c       0.6317        0.0644        0.0055        0.0027      4.455e-02              +++
# 23   LDL   rs2894475   7 107262558  7:107262558           <NA>           <NA>           NA            NA            NA            NA             NA             <NA>
# 24   LDL  rs73729083   7 137559799  7:137559799              t              c       0.8896        0.0256        0.0723        0.0100      4.161e-13              ?++
# 25   LDL  rs35882350  12    623129    12:623129              a              g       0.7587        0.0312       -0.0098        0.0032      2.599e-03              --+
# 26   LDL   rs3747910  20   5528518   20:5528518              a              g       0.7922        0.0218        0.0053        0.0034      1.175e-01              +++
# 27    TC  rs67484410   1  16075835   1:16075835              t              c       0.2540        0.0156       -0.0072        0.0031      1.837e-02              ---
# 28    TC  rs10204498   2  17469783   2:17469783              t              c       0.7602        0.1573        0.0102        0.0035      3.620e-03              +++
# 29    TC  rs57520889   2  17472023   2:17472023              c              g       0.6464        0.0529        0.0274        0.0061      7.250e-06              ++?
# 30    TC   rs2020214   3  73815697   3:73815697              a              t       0.8860        0.0756       -0.0032        0.0053      5.547e-01              -++
# 31    TC    rs903381   5  95399878   5:95399878              a              g       0.0366        0.0096        0.0054        0.0158      7.321e-01              -+-
# 32    TC rs188713108   5 175952733  5:175952733              a              g       0.0061        0.0011       -0.0016        0.0231      9.458e-01              ++-
# 33    TC  rs62621992   6  30459165   6:30459165              t              c       0.0190        0.0040        0.0100        0.0105      3.390e-01              +++
# 34    TC rs201061478   7  52477383   7:52477383              d              i       0.0616        0.0136        0.0291        0.0198      1.410e-01              ?-+
# 35    TC  rs73729083   7 137559799  7:137559799              t              c       0.8895        0.0256        0.0821        0.0100      1.670e-16              ?++
# 36    TC  rs73729087   7 137562668  7:137562668              t              c       0.8982        0.0231        0.0906        0.0102      8.189e-19              ?++
# 37    TC  rs10512369   9 110765359  9:110765359           <NA>           <NA>           NA            NA            NA            NA             NA             <NA>
# 38    TC  rs35882350  12    623129    12:623129              a              g       0.7587        0.0312       -0.0103        0.0032      1.508e-03              --+
# 39    TC  rs17532490  13  41683061  13:41683061              a              g       0.0801        0.0033       -0.0240        0.0050      1.686e-06              ---
# 40    TC rs199986018  20   5544985   20:5544985              a              c       0.8004        0.0209        0.0092        0.0034      7.689e-03              +++
# 41    TG rs200645940   1 186621613  1:186621613              t              c       0.9785        0.0019        0.0123        0.0094      1.904e-01              ++-
# 42    TG rs188096511   1 240387861  1:240387861              a              t       0.0011        0.0002        0.0509        0.0526      3.325e-01              +--
# 43    TG rs186278890   2  74516034   2:74516034              t              c       0.9959        0.0009       -0.0196        0.0482      6.848e-01              --+
# 44    TG rs537734545   2  74728135   2:74728135              c              g       0.0040        0.0008        0.0094        0.0476      8.444e-01              ++-
# 45    TG   rs7612825   3  57462856   3:57462856              t              c       0.5743        0.1311        0.0081        0.0028      3.288e-03              +++
# 46    TG  rs75630351   4 179223627  4:179223627              t              c       0.0160        0.0027        0.0173        0.0128      1.749e-01              +-+
# 47    TG  rs11285757   5 132406598  5:132406598              d              i       0.2417        0.0200       -0.0175        0.0032      5.914e-08              ---
# 48    TG  rs55782534   7   9381826    7:9381826              t              c       0.0361        0.0058       -0.0010        0.0083      9.000e-01              -+-
# 49    TG rs147392990   7 124424635  7:124424635              d              i       0.0500        0.0097        0.0325        0.0141      2.127e-02              ?++
# 50    TG rs569580765   7 152729451  7:152729451              d              i       0.0101        0.0021       -0.0443        0.0348      2.021e-01              ?--
# 51    TG   rs2410684   8  21141988   8:21141988           <NA>           <NA>           NA            NA            NA            NA             NA             <NA>
# 52    TG   rs2511520  11 113274771 11:113274771              t              c       0.1665        0.0500        0.0165        0.0038      1.265e-05              +++
# 53    TG   rs1509517  11 113473871 11:113473871              a              t       0.0975        0.0279        0.0515        0.0154      8.413e-04              ?++
# 54    TG   rs6589401  11 113811847 11:113811847              t              c       0.0782        0.0503        0.0164        0.0078      3.429e-02              -++

# ------------------------------------------------------------------------------                                       
# -- get proxies for missing snps
# ------------------------------------------------------------------------------
# get names of missing snps and their proxies
miss_snps = unique(subset(mvp_te_page, is.na(MVP.te_Allele1))$MarkerName)
miss_snps
# [1] "10:133667886" "14:53797383"  "14:105898715" "7:107262558"  "9:110765359"  "8:21141988" 

miss_snps_proxies = subset(snps_page_proxies, Index.MarkerName %in% miss_snps)
dim(miss_snps_proxies)
# [1] 23 10
miss_snps_proxies
#     Trait.s. Chr_indexSNP Pos_indexSNP   IndexSNP Chr_proxy Pos_proxy    ProxySNP R2_index_proxy_snps Index.MarkerName Proxy.MarkerName
# 627      HDL           10    133667886 rs10430621        10 133664144  rs75564023           0.9037938     10:133667886     10:133664144
# 628      HDL           10    133667886 rs10430621        10 133667886  rs10430621           1.0000000     10:133667886     10:133667886
# 629      HDL           10    133667886 rs10430621        10 133671165  rs11156495           0.9729202     10:133667886     10:133671165
# 630      HDL           10    133667886 rs10430621        10 133673432   rs9888122           0.9209565     10:133667886     10:133673432
# 313      HDL           14     53797383 rs75405126        14  53797383  rs75405126           1.0000000      14:53797383      14:53797383
# 314      HDL           14    105898715 rs60255124        14 105898715  rs60255124           1.0000000     14:105898715     14:105898715
# 315      HDL           14    105898715 rs60255124        14 105898720 rs587741838           1.0000000     14:105898715     14:105898720
# 316      HDL           14    105898715 rs60255124        14 105911275  rs73362846           0.9008720     14:105898715     14:105911275
# 344      LDL            7    107262558  rs2894475         7 107222725   rs2057837           0.8567130      7:107262558      7:107222725
# 345      LDL            7    107262558  rs2894475         7 107224667   rs7778270           0.8567400      7:107262558      7:107224667
# 346      LDL            7    107262558  rs2894475         7 107238447   rs3801948           0.8822910      7:107262558      7:107238447
# 347      LDL            7    107262558  rs2894475         7 107248929   rs7803102           0.8779420      7:107262558      7:107248929
# 348      LDL            7    107262558  rs2894475         7 107250484   rs7793613           0.8797350      7:107262558      7:107250484
# 349      LDL            7    107262558  rs2894475         7 107253304   rs2395907           0.8041000      7:107262558      7:107253304
# 350      LDL            7    107262558  rs2894475         7 107255548   rs3801944           0.9652150      7:107262558      7:107255548
# 351      LDL            7    107262558  rs2894475         7 107258121  rs10273733           0.8051800      7:107262558      7:107258121
# 352      LDL            7    107262558  rs2894475         7 107260856      rs2808           0.9995820      7:107262558      7:107260856
# 353      LDL            7    107262558  rs2894475         7 107261556  rs10274041           0.9998800      7:107262558      7:107261556
# 354      LDL            7    107262558  rs2894475         7 107262558   rs2894475           1.0000000      7:107262558      7:107262558
# 355      LDL            7    107262558  rs2894475         7 107265361 rs542203850           0.8380120      7:107262558      7:107265361
# 356      LDL            7    107262558  rs2894475         7 107265403   rs7806835           0.8048700      7:107262558      7:107265403
# 610       TC            9    110765359 rs10512369         9 110765359  rs10512369           1.0000000      9:110765359      9:110765359
# 514       TG            8     21141988  rs2410684         8  21141988   rs2410684           1.0000000       8:21141988       8:21141988

# some of the missing SNPS only have themselves listed as a proxy
sapply(split(miss_snps_proxies, miss_snps_proxies$Index.MarkerName), 
       function(f) 
           all(f$Index.MarkerName == f$Proxy.MarkerName))
# 10:133667886 14:105898715  14:53797383  7:107262558   8:21141988  9:110765359 
#        FALSE        FALSE         TRUE        FALSE         TRUE         TRUE

# ------------------------------------------------------------------------------                                       
# -- find proxies in MVP data
# ------------------------------------------------------------------------------
# assumes that MVP data frames exist in the workspace:
#  MVP.te.HDL.gwas.tsv.gz, MVP.te.LDL.gwas.tsv.gz, MVP.te.TC.gwas.tsv.gz, MVP.te.TG.gwas.tsv.gz
.getProxiesInMVPte = function (miss_snps=NULL, miss_snps_proxies=NULL) {
    # create list to store results and loop through missing SNPs
    mvp_proxies = list()
    for (s in miss_snps) {
        
        # get proxies and trait for this missing SNP
        xprox0 = subset(miss_snps_proxies, Index.MarkerName == s)
        xprox = xprox0$Proxy.MarkerName
        xtrait = unique(xprox0$Trait.s.)

        # get appropriate set of MVP results from the workspace
        xname = grep(paste0("^MVP.*", xtrait), ls(name=".GlobalEnv"), value=T)
        xmvp = get(xname)
        
        # print updates and sanity check
        cat(paste0("--", s, " missing from ", xtrait, " results\n",
                   "  checking proxy SNPs ", paste0(xprox, collapse=", "), "\n", 
                   "  in ", xname, "\n\n"))
        
        # subset MVP results to the proxy SNPs
        tmp_mvp = xmvp[match(xprox, xmvp$MarkerName), ]
        
        # make sure missing proxies have IDs in the output
        tmp_mvp$MarkerName = xprox
        
        # add columns for missing index SNP and trait
        tmp_mvp$Index.MarkerName = rep(s, nrow(tmp_mvp))
        tmp_mvp$Trait = rep(xtrait, nrow(tmp_mvp))
        
        # remove rows where missing SNP is listed as its own proxy
        #  unless only proxy is itself, i.e. nrow(tmp_mvp) == 1
        if (nrow(tmp_mvp) > 1) {
            sameRow = tmp_mvp$Index.MarkerName == tmp_mvp$MarkerName
            tmp_mvp = tmp_mvp[!sameRow, ]
        }
        
        mvp_proxies[[s]] = tmp_mvp
    }
    
    # collapse to single data frame, reorder columns, set dimnames
    mvp_proxies = as.data.frame(do.call(rbind, mvp_proxies))
    colNum = ncol(mvp_proxies)
    colOrd = c(colNum,   # should be Trait
               colNum-1, # should be Index.MarkerName
               1:(colNum-2))
    mvp_proxies = mvp_proxies[, colOrd]
    names(mvp_proxies)[-c(1:2)] = paste0("Proxy.", 
                                         names(mvp_proxies)[-c(1:2)])
    rownames(mvp_proxies) = NULL
    return(mvp_proxies)
}

mvp_proxies = .getProxiesInMVPte(miss_snps, miss_snps_proxies)
# --10:133667886 missing from HDL results
#   checking proxy SNPs 10:133664144, 10:133667886, 10:133671165, 10:133673432
#   in MVP.te.HDL.gwas.tsv.gz
# 
# --14:53797383 missing from HDL results
#   checking proxy SNPs 14:53797383
#   in MVP.te.HDL.gwas.tsv.gz
# 
# --14:105898715 missing from HDL results
#   checking proxy SNPs 14:105898715, 14:105898720, 14:105911275
#   in MVP.te.HDL.gwas.tsv.gz
# 
# --7:107262558 missing from LDL results
#   checking proxy SNPs 7:107222725, 7:107224667, 7:107238447, 7:107248929, 7:107250484, 7:107253304, 7:107255548, 7:107258121, 7:107260856, 7:107261556, 7:107262558, 7:107265361, 7:107265403
#   in MVP.te.LDL.gwas.tsv.gz
# 
# --9:110765359 missing from TC results
#   checking proxy SNPs 9:110765359
#   in MVP.te.TC.gwas.tsv.gz
# 
# --8:21141988 missing from TG results
#   checking proxy SNPs 8:21141988
#   in MVP.te.TG.gwas.tsv.gz

mvp_proxies
#    Trait Index.MarkerName Proxy.MarkerName Proxy.Allele1 Proxy.Allele2 Proxy.Freq1 Proxy.FreqSE Proxy.Effect Proxy.StdErr Proxy.P.value Proxy.Direction
# 1    HDL     10:133667886     10:133664144          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 2    HDL     10:133667886     10:133671165          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 3    HDL     10:133667886     10:133673432          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 4    HDL      14:53797383      14:53797383          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 5    HDL     14:105898715     14:105898720          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 6    HDL     14:105898715     14:105911275             a             g      0.0114       0.0018       0.0066       0.0313     8.321e-01             ?+-
# 7    LDL      7:107262558      7:107222725             c             g      0.7284       0.0296      -0.0151       0.0034     9.004e-06             --?
# 8    LDL      7:107262558      7:107224667             a             g      0.2713       0.0294       0.0151       0.0034     9.750e-06             ++?
# 9    LDL      7:107262558      7:107238447          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 10   LDL      7:107262558      7:107248929          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 11   LDL      7:107262558      7:107250484          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 12   LDL      7:107262558      7:107253304          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 13   LDL      7:107262558      7:107255548          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 14   LDL      7:107262558      7:107258121          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 15   LDL      7:107262558      7:107260856          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 16   LDL      7:107262558      7:107261556          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 17   LDL      7:107262558      7:107265361          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 18   LDL      7:107262558      7:107265403          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 19    TC      9:110765359      9:110765359          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>
# 20    TG       8:21141988       8:21141988          <NA>          <NA>          NA           NA           NA           NA            NA            <NA>

# ------------------------------------------------------------------------------                                       
# -- write output, then save workspace with session info
# ------------------------------------------------------------------------------
# MVP results for PAGE snps
write.csv(mvp_te_page, row.names=F,
          file="PAGE_Lipids_Replication_MVP.te_noProxies_updated.csv")

# MVP results for proxies of PAGE snps missing from MVP
write.csv(mvp_proxies, row.names=F,
          file="PAGE_Lipids_Replication_MVP.te_Proxies_for_missing_updated.csv")

# csv files are combined in PAGE_Lipids_Replication_MVP.transethnic_updated.xlsx

# get session and workspace info
session_info = sessionInfo()
full_workspace = ls()

# workspace without full MVP results
save(list=ls()[-match(mvp_te_files, ls())],
     file="lipids_replication_transethnic_WORKSPACE-no_full_MVP_updated.RData")

# workspace including full MVP results
save(list=ls(),
     file="lipids_replication_transethnic_WORKSPACE-with_full_MVP_updated.RData")
