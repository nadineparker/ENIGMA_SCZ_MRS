
library(data.table); library(tidyverse)

# read in hm3plus reference
ref <- readRDS("Protocol/ref/map_hm3_plus.rds")

# read in sumstats
SCZ_EUR <- fread("Protocol/sumstats/PGC_SCZ_2022_EUR_no23andMe_GRCh38.sumstats.gz")
SCZ_EUR_T1 <- fread("Protocol/sumstats/PGC_SCZ_2022_EUR_Type1_GRCh38.sumstats.gz")
SCZ_EUR_T2 <- fread("Protocol/sumstats/PGC_SCZ_2022_EUR_Type2_GRCh38.sumstats.gz")
SCZ_Multi <- fread("Protocol/sumstats/PGC_SCZ_2022_Multi_no23andMe_GRCh38.sumstats.gz")

## Restict to hm3plus
SCZ_EUR <- SCZ_EUR[SCZ_EUR$RSID %in% ref$rsid,]
SCZ_EUR_T1 <- SCZ_EUR_T1[SCZ_EUR_T1$RSID %in% ref$rsid,]
SCZ_EUR_T2 <- SCZ_EUR_T2[SCZ_EUR_T2$RSID %in% ref$rsid,]
SCZ_Multi <- SCZ_Multi[SCZ_Multi$RSID %in% ref$rsid,]

# Write out hm3plus restricted sumstats
# fwrite(
#   SCZ_EUR, quote=F, row.names = F, sep = "\t",
#   file = "Protocol/sumstats/PGC_SCZ_2022_EUR_no23andMe_hm3plus_GRCh38.sumstats.gz")
# 
# fwrite(
#   SCZ_EUR_T1, quote=F, row.names = F, sep = "\t",
#   file = "Protocol/sumstats/PGC_SCZ_2022_EUR_Type1_hm3plus_GRCh38.sumstats.gz")
# 
# fwrite(
#   SCZ_EUR_T2, quote=F, row.names = F, sep = "\t",
#   file = "Protocol/sumstats/PGC_SCZ_2022_EUR_Type2_hm3plus_GRCh38.sumstats.gz")
# 
# fwrite(
#   SCZ_Multi, quote=F, row.names = F, sep = "\t",
#   file = "Protocol/sumstats/PGC_SCZ_2022_Multi_no23andMe_hm3plus_GRCh38.sumstats.gz")
