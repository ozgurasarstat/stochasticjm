
library(joineRML)

head(renal$prot)

renal$prot$id <- as.character(renal$prot$id)
renal$haem$id <- as.character(renal$haem$id)
renal$gfr$id  <- as.character(renal$gfr$id)
renal$surv$id <- as.character(renal$surv$id)

renal$prot$gender <- as.character(renal$prot$gender)
renal$haem$gender <- as.character(renal$haem$gender)
renal$gfr$gender  <- as.character(renal$gfr$gender)
renal$surv$gender <- as.character(renal$surv$gender)

# ids which has data only before 1 year after transplantation
id_prot <- names(which(tapply(renal$prot$years, renal$prot$id, function(x) max(x) >= 1) == FALSE))
id_haem <- names(which(tapply(renal$haem$years, renal$haem$id, function(x) max(x) >= 1) == FALSE))
id_gfr  <- names(which(tapply(renal$gfr$years, renal$gfr$id, function(x) max(x) >= 1) == FALSE))

# data for prot
prot_surv_short <- list(prot_short = renal$prot[!(renal$prot$id %in% id_prot), ],
                        surv_short = renal$surv[!(renal$surv$id %in% id_prot), ])
prot_surv_short$prot_short <- prot_surv_short$prot_short[prot_surv_short$prot_short$years >= 1, ]
prot_surv_short$prot_short$years_after1year   <- prot_surv_short$prot_short$years - 1
prot_surv_short$prot_short$fuyears_after1year <- prot_surv_short$prot_short$fuyears - 1
prot_surv_short$surv_short$fuyears_after1year <- prot_surv_short$surv_short$fuyears - 1

# data for haem
haem_surv_short <- list(haem_short = renal$haem[!(renal$haem$id %in% id_haem), ],
                        surv_short = renal$surv[!(renal$surv$id %in% id_haem), ])
haem_surv_short$haem_short <- haem_surv_short$haem_short[haem_surv_short$haem_short$years >= 1, ]
haem_surv_short$haem_short$years_after1year   <- haem_surv_short$haem_short$years - 1
haem_surv_short$haem_short$fuyears_after1year <- haem_surv_short$haem_short$fuyears - 1
haem_surv_short$surv_short$fuyears_after1year <- haem_surv_short$surv_short$fuyears - 1

# data for gfr
gfr_surv_short <- list(gfr_short  = renal$gfr[!(renal$gfr$id %in% id_gfr), ],
                       surv_short = renal$surv[!(renal$surv$id %in% id_gfr), ])
gfr_surv_short$gfr_short <- gfr_surv_short$gfr_short[gfr_surv_short$gfr_short$years >= 1, ]
gfr_surv_short$gfr_short$years_after1year    <- gfr_surv_short$gfr_short$years - 1
gfr_surv_short$gfr_short$fuyears_after1year  <- gfr_surv_short$gfr_short$fuyears - 1
gfr_surv_short$surv_short$fuyears_after1year <- gfr_surv_short$surv_short$fuyears - 1

renal_data <- list(renal_original = renal,
                   prot_surv_short = prot_surv_short,
                   haem_surv_short = haem_surv_short,
                   gfr_surv_short = gfr_surv_short)

usethis::use_data(renal_data, renal_data, overwrite = TRUE)
