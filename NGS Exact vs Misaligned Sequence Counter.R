pat <- "(.*)-(.*)" #.* reads up to last desired sting, "-" in this case, and second .* reads everything after encountered string
load("NAR.Rda")
NAR$full_name <- str_c(NAR$name,
             '-',
             NAR$matching_Sequence)

NAR['CountN'] <- c(1)

NAR_Count <- ddply(NAR,"full_name", numcolwise(sum))
NAR_Count_Sequence <- transform(NAR_Count, full_name = sub(pat, "\\2", full_name))
MC <- data.frame(str_split_fixed(NAR_Count$full_name,"-",3))
MC <- data.frame(str_c(MC$X1,
                        '-',
                        MC$X2))
MC <- cbind(MC$str_c.MC.X1.......MC.X2., NAR_Count_Sequence$full_name)
MC <- data.frame(MC)
MC <- cbind(MC, NAR_Count$CountN)
colnames(MC) <- c('Crispr_ID', 'Sequence', 'CountN')



pat <- "(.*)-(.*)"
load("EAR.Rda")
EAR$full_name <- str_c(EAR$name,
             '-',
             EAR$matching_Sequence)

EAR['CountE'] <- c(1)

EAR_Count <- ddply(EAR, "full_name", numcolwise(sum))
EAR_Count_Sequence <- transform(EAR_Count, full_name = sub(pat, "\\2", full_name))
EA <- data.frame(str_split_fixed(EAR_Count$full_name,"-",3))
EA <- data.frame(str_c(EA$X1,
                        '-',
                        EA$X2))
EA <- cbind(EA$str_c.EA.X1.......EA.X2., EAR_Count_Sequence$full_name)
EA <- data.frame(EA)
EA <- cbind(EA,EAR_Count$CountE)
colnames(EA) <- c( 'Crispr_ID', 'Sequence', 'CountE')



EA_Final <- unlist(EA$Crispr_ID)
MC_Final <- unlist(MC$Crispr_ID)
Exact_versus_Non_Exact <- data.frame(MC_Final %in% EA_Final)
Exact_versus_Non_Exact <- data.frame(Exact_versus_Non_Exact,
                                     MC)

NFE <- Exact_versus_Non_Exact %>% filter(Exact_versus_Non_Exact$MC_Final..in..EA_Final == TRUE)
ers <- Exact_versus_Non_Exact %>% filter(Exact_versus_Non_Exact$MC_Final..in..EA_Final == FALSE)

save(NFE, file="NFE.Rda")



################################################################################

reference <- read_excel("Reference.xlsx")



EA_test <- data.frame(str_c(EA$Crispr_ID,
                       '-',
                       EA$Sequence),
                      EA$CountE)


REF_test <- data.frame(str_c(reference$`Crispr ID`,
                             '-',
                             reference$Sequence))

misaligned_Counts <- data.frame(EA_test$str_c.EA.Crispr_ID.......EA.Sequence. %in% REF_test$str_c.reference..Crispr.ID........reference.Sequence.)

EA_test <- data.frame(EA_test,
                      misaligned_Counts)


EA <- EA_test %>% filter(EA_test$EA_test.str_c.EA.Crispr_ID.......EA.Sequence...in..REF_test.str_c.reference..Crispr.ID........reference.Sequence. == TRUE)
INCORRECTLY_ALIGNED <- EA_test %>% filter(EA_test$EA_test.str_c.EA.Crispr_ID.......EA.Sequence...in..REF_test.str_c.reference..Crispr.ID........reference.Sequence. == FALSE)


names(EA) <- c("full_name", "CountE", 'boolean_column')

EA_CrisprID <- transform(EA, full_name = sub(pat, "\\1", full_name))
EA_Sequence <- transform(EA, full_name = sub(pat, "\\2", full_name))

EA <- data.frame(EA_CrisprID$full_name,
                 EA_Sequence$full_name,
                 EA$CountE)

names(EA) <- c("Crispr_ID", "Sequence", 'CountE')


################################################################################

if(nrow(MC) >= nrow(EA)) {
  a <- MC
  b <- EA
}else{
    a <- EA
    b <- MC
}



EVMA <- do.call(rbind, lapply(1:nrow(a), function(x)  {
                                           print(x)
                                           matchIX <- grepl(a$Crispr_ID[x],
                                                            data.frame(t(b), stringsAsFactors = FALSE),
                                                            ignore.case = TRUE)
                                           if(any(matchIX)){
                                             cbind(a[x, ], b[matchIX, ])
                                           }
                                           }))
save(EVMA, file="EVMA.Rda")



load("EVMA.Rda")
load("NFE.Rda")
Reference_Sheet <- read_excel("Reference.xlsx")

Check <- rbind.fill(EVMA,NFE)
Final <- rbind.fill(EVMA)

Final_unlist <- unlist(Final$Crispr_ID)
Reference_Sheet_unlist <- unlist(Reference_Sheet$`Crispr ID`)
Check_unlist <- unlist(Check$Crispr_ID)

Final_missing <- Reference_Sheet_unlist %in% Check_unlist
total_CRISPR_Check <- data.frame(Reference_Sheet,
                                 Final_missing)

Filtered_List <- total_CRISPR_Check %>% filter(total_CRISPR_Check$Final_missing == FALSE)
Final$ratio <- Final$CountE / Final$CountN
Final$Base_2e_log <- log(Final$ratio, base = exp(2))
colnames(Final) <- c('Crispr_IDN',  'SequenceN', 'CountN', 'Crispr_IDE', 'SequenceE', 'CountE', 'ratio', 'Base_2e_log')

Final<- Final %>%
           arrange(ratio)

save(Final,file="NGS_Sequence.RDA")




################################################################################




load("NGS_Sequence.RDA")


Final_Plot <- Final$Base_2e_log
barplot(Final_Plot,
        ylim = c(-7,7))


plot(x=Final$CountE, y=Final$CountN)
Threshold_CountN <- Final %>% filter(Final$CountN >= 2500)
for(i in 1:nrow(Threshold_CountN))  {
  text(Threshold_CountN$CountE[i], Threshold_CountN$CountN[i]-1, labels=Threshold_CountN$full_name[i], cex =.6)
}


plot(x=Final$CountN, y=Final$CountE)
Threshold_CountE <- Final %>% filter(Final$CountE >= 20000)
for(i in 1:nrow(Threshold_CountE))  {
  text(Threshold_CountE$CountN[i], Threshold_CountE$CountE[i]-1, labels=Threshold_CountE$full_name[i], cex =.6)
}

plot(x=Final$ratio, y=Final$Base_2e_log)
Threshold_ratio <- Final %>% filter(Final$ratio >= 130)
Threshold_Base_2e_log <- Final

for(i in 1:nrow(Threshold_ratio))  {
  text(Threshold_ratio$ratio[i], Threshold_ratio$Base_2e_log[i]-1, labels=Threshold_ratio$full_name[i], cex =.6)
}
for(i in 1:nrow(Threshold_Base_2e_log)) {
  text(Threshold_Base_2e_log$Base_2e_log[i], Threshold_Base_2e_log$ratio[i]-1, labels=Threshold_Base_2e_log$full_name[i], cex =.6)
}


pairs(~CountE+ratio+Base_2e_log+CountN, data = Final,
      main = "Count Matrix")









       
