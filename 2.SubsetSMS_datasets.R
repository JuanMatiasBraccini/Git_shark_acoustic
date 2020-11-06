
########################- SCRIPT FOR SUBSETTING SMS DATASETS  -##########################

#notes: This script subset the data set downloaded from the SMS website

library(chron) #for extracting time
library(data.table)
library(dplyr)
library(stringr)


# Data section ------------------------------------------------------------
setwd("C:/Users/myb/Desktop/SMN_downloads")
DetectionsComplete=fread("DetectionsComplete.csv",data.table=FALSE)
TagsCompleteReport=fread("TagsCompleteReport.csv",data.table=FALSE)

UTC.WA=8

# Procedure section ------------------------------------------------------------
DetectionsComplete=DetectionsComplete%>%
                    mutate(DateTime=as.POSIXlt(Datetime,  format="%d/%m/%Y %H:%M",tz="GMT"),
                           DateTime.local=DateTime+(UTC.WA*3600),
                           Date.local=trunc(DateTime.local,"days"),
                           Time.local=times(format(DateTime.local, "%H:%M:%S")),
                           Latitude=-abs(Latitude)) 

TagsCompleteReport=TagsCompleteReport%>%
                    mutate(ReleaseDate=as.POSIXlt(ReleaseDate,  format="%d/%m/%Y %H:%M",tz="GMT"),
                           ReleaseDate.local=ReleaseDate+(UTC.WA*3600),
                           ReleaseLatitude=-abs(ReleaseLatitude))

ALL.tags=TagsCompleteReport%>%
  mutate(TagCode=ifelse(nchar(Transmitter)>10,as.numeric(gsub("-", "", str_match(Transmitter, '(?:-[^-]+){1}$')[,1])),
                        ifelse(nchar(Transmitter)<10,as.numeric(gsub("-", "", Transmitter)),
                               NA)))

TagsCompleteReport=TagsCompleteReport%>%
                    filter(Transmitter%in%unique(DetectionsComplete$Transmitter))

# duplicates=with(TagsCompleteReport,table(Transmitter,Species))
# duplicates=rowSums(duplicates)
# duplicates=duplicates[which(duplicates>1)]
# no.duplicates=TagsCompleteReport%>%
#                 filter(!Transmitter%in%names(duplicates))%>%
#                 pull(Transmitter)
# duplicates=TagsCompleteReport%>%
#         filter(Transmitter%in%names(duplicates))%>%
#         arrange(Transmitter)%>%
#   dplyr::select(Transmitter,TagType,FinTag,Species,Sex,ReleaseDate,ReleaseLength)



DATA=left_join(DetectionsComplete,TagsCompleteReport,by=c('Transmitter','Species'))%>%
        mutate(TagCode=ifelse(nchar(Transmitter)>10,as.numeric(gsub("-", "", str_match(Transmitter, '(?:-[^-]+){1}$')[,1])),
                       ifelse(nchar(Transmitter)<10,as.numeric(gsub("-", "", Transmitter)),
                       NA)),
               SerialNumber=sub(".*-", "", Receiver),
               StationName2=StationArray,
               Depth=NA)%>%
  dplyr::select(c(TagCode,Species,Sex,Name,SerialNumber,ReleaseDate,ReleaseLatitude,
                  ReleaseLongitude,Latitude,Longitude,Depth,StationName2,
                  DateTime.local,Date.local,Time.local,
                  FinTag,Sex,ReleaseLength,ReleaseSite,ReleaseState))%>%
  mutate(Project.rel=ifelse(ReleaseState=='SA',"South.Australia",
                            ifelse(ReleaseState=='WA',"SMN",       
                            ReleaseState)))




#Export data
setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/SMS_donwloads_processed")
write.csv(DATA,file=paste(paste(unique(DATA$Species),collapse='_'),".",Sys.Date(),".csv",sep=''))


setwd("C:/Matias/Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/SMS_tags")
write.csv(ALL.tags,file=paste('ALL.tags',".",Sys.Date(),".csv",sep=''))