library(tidyverse)
library(strex)
library(lubridate)
library(data.table)

library(actel)  #acoustic telemetry packages
library(RSP)

# Input data --------------------------------------------------------------------
hndl=function(x) paste('M:/Fisheries Research/FinFish/Shark/Matias/Commercial_species_10-years/',x,sep='')
Alldetections=read.csv(hndl('Alldetections.csv'))
Stations=read.csv(hndl('Stations.csv'))
Taglist=read.csv(hndl('Taglist.csv'))


# Manipulate data --------------------------------------------------------------------

# Fill in NA TagCodes
Alldetections=Alldetections%>%
                mutate(No.TagCode=ifelse(is.na(TagCode),'Yes','No'),
                       TagCode.original=TagCode,
                       TagCode=ifelse(is.na(TagCode),as.numeric(str_after_nth(Transmitter, "-", 2)),TagCode))

# Remove duplicates (i.e. same date-time-location-species)
Alldetections=Alldetections%>%
  mutate(dupli=paste(Date,Time,StationLatitude,StationLongitude,Species,TagCode))

Duplis=table(Alldetections$dupli)
Duplis=Duplis[Duplis>1]
Alldetections.duplicated=Alldetections%>%
  filter(dupli%in%names(Duplis))

Alldetections=Alldetections%>%
  distinct(dupli,.keep_all=TRUE)%>%
  dplyr::select(-dupli)

# Add tag details
Alldetections=Alldetections%>%
  left_join(Taglist%>%
              dplyr::select(Transmitter,Species,TagType,FinTag,
                            Sex,ReleaseDate,ReleaseLength,ReleaseLatitude,
                            ReleaseLongitude,ReleaseSite,RelatedTags),
            by=c('Transmitter','Species'))

# Format date and time
Alldetections=Alldetections%>%
  mutate(Timestamp=ymd_hm(paste(as.POSIXlt(Date,format='%d/%m/%Y'),Time)),
         Date.old=Date,
         Time.old=Time,
         ReleaseDate.old=ReleaseDate,
         ReleaseDate=as.POSIXlt(ReleaseDate,format='%d/%m/%Y'),
         Date=as.IDate(Timestamp),
         Time=as.ITime(Timestamp,units=c("hours", "minutes")))

#Fix dodgy longitude (112.3168 when it should be 115.3168)
Alldetections=Alldetections%>%
              mutate(StationLongitude=ifelse(StationLongitude<113 &StationName=='OTN 40',115.3168,StationLongitude))

# Add station region
Alldetections=Alldetections%>%
                mutate(Region=case_when(StationLatitude>(-24)~'Ningaloo',
                                        StationLatitude<(-24) & StationLatitude>(-33) & StationLongitude<116 ~'Perth',
                                        StationLatitude<(-33) ~'South'))


# Plot timeline --------------------------------------------------------------------
Sp=unique(Alldetections$Species)
Spatial.range=data.frame(Lat=range(Alldetections$StationLatitude),
                         Long=range(Alldetections$StationLongitude))
fn.timeline=function(d)
{
  Levels=d%>%distinct(TagCode,ReleaseDate)%>%arrange(ReleaseDate)%>%mutate(Levels=row_number())%>%dplyr::select(-ReleaseDate)
  LBLS=Levels$TagCode
  names(LBLS)=Levels$Levels
  
   p=d%>%distinct(StationLatitude, StationLongitude,Region)%>%
     ggplot(aes(StationLongitude,StationLatitude))+
     geom_point(aes(color=Region))+
     ggtitle(unique(d$Species))+
     xlim(Spatial.range$Long)+
     ylim(Spatial.range$Lat)
   print(p)
  
  p=d%>%
    left_join(Levels,by='TagCode')%>%
    distinct(Levels,TagCode,Date,Region)%>%
    ggplot(aes(Date, as.factor(Levels)))+
    geom_point(aes(color=Region))+
    scale_y_discrete("TagCode",labels=LBLS)+
    ggtitle(unique(d$Species))+
    theme(axis.text=element_text(size=8))
  
  print(p)
}
pdf(hndl("Detections thru time.pdf"))
for(s in 1:length(Sp)) fn.timeline(d=Alldetections%>%filter(Species==Sp[s]))
dev.off()
