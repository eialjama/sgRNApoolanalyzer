library("dplyr")
library("Rsamtools")
library("readxl")
library("writexl")
library("stringr")
library("plyr")

setwd("C:/Users/matthew.woods/Desktop/Ion Torrent Read Counter and Aligner/Ion Read Counter and Aligner")

bam <- scanBam("Rerun.bam")


dataset <- data.frame(Sequence = bam[[1]][["seq"]],
                      name = bam[[1]][["rname"]])


save(dataset, file="SequenceData.Rda")


load("SequenceData.Rda")

REGEXMatch <- 'ACCG(.{20})GTTT' #Same REGEX used by other groups
reference_file <- read_excel("Reference.xlsx")


#dataset <- filter(dataset, grepl(REGEXMatch, Sequence))



Exact_Alignment <- filter(dataset, grepl(REGEXMatch, Sequence))

save(Exact_Alignment, file="Exact_Alignment.Rda")

load("Exact_Alignment.Rda")

Exact_Alignment$matching_Sequence <- str_extract(Exact_Alignment$Sequence, REGEXMatch)
Exact_Alignment$matching_Sequence <- str_sub(Exact_Alignment$matching_Sequence, 5,-5) #Reduce integer by 1 to get actual string slice


Exact_Alignment_unlist <- unlist(Exact_Alignment$matching_Sequence)
reference_file_unlist <- unlist(reference_file$Sequence)

Final_Count <- data.frame(Exact_Alignment,
                          int <-Exact_Alignment_unlist %in% reference_file_unlist)

NA_List <-Final_Count[is.na(Final_Count$name),]

Final_Count <- na.omit(Final_Count)


EAR <- Final_Count %>% filter(Final_Count$int == TRUE)
NAR <- Final_Count %>% filter(Final_Count$int != TRUE)


save(EAR, file='EAR.Rda')
save(NAR, file='NAR.Rda')


################################################################################
reference_file <- read_excel("Reference.xlsx")


load('NAR.Rda')
load("EAR.Rda")

EAR['Count'] <- c(1)


EAR_Count <- ddply(EAR, "matching_Sequence", numcolwise(sum))

NAR <- str_c(NAR$name,
                  '-',
            NAR$matching_Sequence)

NAR <- data.frame(NAR)

NAR['Count'] <- c(1)

NAR_Count <- ddply(NAR, "NAR", numcolwise(sum))




Final_output <- do.call(rbind, lapply(1:nrow(EAR_Count),
                                            function(x){
                                              print(x)
                                              matchIX <- grepl(EAR_Count$matching_Sequence[x],
                                                               data.frame(t(reference_file), stringsAsFactors = FALSE),
                                                               ignore.case = TRUE)
                                              if(any(matchIX))
                                                cbind(EAR_Count[x, ], reference_file[matchIX, ])
                                            }))  




Final_output_unlist <- unlist(Final_output$`Crispr ID`)
reference_file_unlist <- unlist(reference_file$`Crispr ID`)



Final_Count <- data.frame(reference_file_unlist,
                          int <-reference_file_unlist %in% Final_output_unlist)

Final_Count_to_Ref <- Final_Count %>% filter(Final_Count$int....reference_file_unlist..in..Final_output_unlist == FALSE)


write_xlsx(Final_Count_to_Ref, "Missing Guides.xlsx")
write_xlsx(Final_output, "Exact alignment Counts.xlsx")
write_xlsx(NAR_Count, 'Non-Perfect alignment Counts.xlsx')







