load("NGS_Crispr.RDA")
NGS_Crispr <- Final

load("NGS_Sequence.RDA")
NGS_Sequence <- Final

NGS_Crispr_Threshold <- NGS_Crispr



NGS_Sequence_Condesnsed <- data.frame(NGS_Sequence$Crispr_IDN,
                                      NGS_Sequence$SequenceN,
                                      NGS_Sequence$SequenceE,
                                      NGS_Sequence$CountN,
                                      NGS_Sequence$CountE,
                                      NGS_Sequence$ratio,
                                      NGS_Sequence$Base_2e_log)

colnames(NGS_Sequence_Condesnsed) <- c('Crispr_IDN', 
                                       'SequenceN', 
                                       'SequenceE', 
                                       'CountN', 
                                       'CountE', 
                                       'ratio', 
                                       'Base_2e_log')

if(nrow(NGS_Crispr_Threshold) >= nrow(NGS_Sequence_Condesnsed)) {
  a <- NGS_Crispr_Threshold
  b <- NGS_Sequence_Condesnsed
}else{
  a <- NGS_Sequence_Condesnsed
  b <- NGS_Crispr_Threshold
}


THRESHOLD_0 <- do.call(rbind, lapply(1:nrow(a), function(x)  {
  print(x)
  matchIX <- grepl(a$Crispr_IDN[x],
                   data.frame(t(b), stringsAsFactors = FALSE),
                   ignore.case = TRUE)
  if(any(matchIX)){
    cbind(a[x, ], b[matchIX, ])
  }}))


colnames(THRESHOLD_0) <- c('Crispr_ID', 
                           'Sequence_Misaligned', 
                           'Sequence_Exact', 
                           'Count_Misaligned_For_Sequence', 
                           'Count_Exact_For_Sequence', 
                           'Ratio_For_Sequence', 
                           'Base_2e_Log_For_Sequence', 
                           'x', 
                           'Count_Misaligned_For_Crispr_ID', 
                           'y', 
                           'Count_Exact_For_Crispr_ID', 
                           'Ratio_For_Crispr_ID', 
                           'Base_2e_Log_For_Crispr_ID')

Final = subset(THRESHOLD_0, select = -c(x,y) )
Final <- Final %>% arrange(Crispr_ID)
Final$AlignmentError <- Final$Count_Exact_For_Crispr_ID - Final$Count_Exact_For_Sequence

write_xlsx(Final,"Exact vs Misaligned Sequence and Gene comparison.xlsx")
write_xlsx(NGS_Crispr_Threshold,"Exact vs Misaligned Gene comparison.xlsx")
