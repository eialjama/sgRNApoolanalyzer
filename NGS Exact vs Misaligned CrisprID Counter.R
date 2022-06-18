

pat <- "(.*)-(.*)" #.* reads up to last desired sting, "-" in this case, and second .* reads everything after encountered string
load("NAR.Rda")
#MCC <- data.frame(str_split_fixed(NAR,"-",3))
MCC <- data.frame(NAR$name)

MCC['CountN'] <- c(1)
MCC <- ddply(MCC,"NAR.name", numcolwise(sum))
colnames(MCC) <- c('Crispr_ID', 'CountN')




pat <- "(.*)-(.*)" #.* reads up to last desired sting, "-" in this case, and second .* reads everything after encountered string
load("EAR.Rda")
EAC <- data.frame(EAR$name)

EAC['CountE'] <- c(1)
EAC <- ddply(EAC,"EAR.name", numcolwise(sum))
colnames(EAC) <- c('Crispr_ID', 'CountE')

EA_Final <- unlist(EAC$Crispr_ID)
MC_Final <- unlist(MCC$Crispr_ID)
Exact_versus_Non_Exact <- data.frame(MC_Final %in% EA_Final)
Exact_versus_Non_Exact <- data.frame(Exact_versus_Non_Exact,
                                     MCC)

NFEC <- Exact_versus_Non_Exact %>% filter(Exact_versus_Non_Exact$MC_Final..in..EA_Final == FALSE)
save(NFEC, file="NFEC.Rda")

if(nrow(MCC) >= nrow(EAC)) {
  a <- MCC
  b <- EAC
}else{
  a <- EAC
  b <- MCC
}


EVMAC <- do.call(rbind, lapply(1:nrow(a), function(x)  {
  print(x)
  matchIX <- grepl(a$Crispr_ID[x],
                   data.frame(t(b), stringsAsFactors = FALSE),
                   ignore.case = TRUE)
  if(any(matchIX)){
    cbind(a[x, ], b[matchIX, ])
  }
}))


save(EVMAC, file="EVMAC.Rda")

load("EVMAC.Rda")
load("NFEC.Rda")
Reference_Sheet <- read_excel("Reference.xlsx")

Check <- rbind.fill(EVMAC,NFEC)
Final <- rbind.fill(EVMAC)

Final_unlist <- unlist(Final$Crispr_ID)
Reference_Sheet_unlist <- unlist(Reference_Sheet$`Crispr ID`)
Check_unlist <- unlist(Check$Crispr_ID)

Final_missing <- Reference_Sheet_unlist %in% Check_unlist
total_CRISPR_Check <- data.frame(Reference_Sheet,
                                 Final_missing)

Filtered_List <- total_CRISPR_Check %>% filter(total_CRISPR_Check$Final_missing == FALSE)


Final$ratio <- Final$CountE / Final$CountN
Final$Base_2e_log <- log(Final$ratio, base = exp(2))

colnames(Final) <- c('Crispr_IDN',  'CountN', 'Crispr_IDE', 'CountE', 'ratio', 'Base_2Log')

Final <- Final %>% arrange(ratio)

save(Final, file="NGS_CRISPR.RDA")

################################################################################

load("NGS_CRISPR.RDA")

Final_Plot <- Final$Base_2Log
barplot(Final_Plot,
        ylim = c(-3,3),
        ylab = 'Read depth of Log Base 2',
        xlab = '3320 targets ')


plot(x=Final$CountE, y=Final$CountN, ylab = 'Number of Misaligned Sequences', xlab = 'Number of Perfectly Aligned Sequences')
Threshold_CountN <- Final %>% filter(Final$CountN >= 2500)
for(i in 1:nrow(Threshold_CountN))  {
  text(Threshold_CountN$CountE[i], Threshold_CountN$CountN[i]-1, labels=Threshold_CountN$Crispr_IDN[i], cex =.6)
}


plot(x=Final$CountN, y=Final$CountE)
Threshold_CountE <- Final %>% filter(Final$CountE >= 20000)
for(i in 1:nrow(Threshold_CountE))  {
  text(Threshold_CountE$CountN[i], Threshold_CountE$CountE[i]-1, labels=Threshold_CountE$Crispr_IDN[i], cex =.6)
}

plot(x=Final$ratio, y=Final$Base_2Log)
Threshold_ratio <- Final %>% filter(Final$ratio >= 130)
Threshold_Base_2e_log <- Final

for(i in 1:nrow(Threshold_ratio))  {
  text(Threshold_ratio$ratio[i], Threshold_ratio$Base_2Log[i]-1, labels=Threshold_ratio$Crispr_IDN[i], cex =.6)
}


for(i in 1:nrow(Threshold_Base_2e_log)) {
  text(Threshold_Base_2e_log$Base_2Log[i], Threshold_Base_2e_log$ratio[i]-1, labels=Threshold_Base_2e_log$Crispr_IDN[i], cex =.6)
}


pairs(~CountE+ratio+Base_2Log+CountN, data = Final,
      main = "Count Matrix")



