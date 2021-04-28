#SOURCE ALL ACOUSTIC TRANSMITTES RELEVANT TO STUDY
library(lubridate)
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

WA.Fisheries=read.csv(file=handl_OneDrive("Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/TaglistReport.1.11.2018.csv"),stringsAsFactors =F)
South.Africa=read.csv(file=handl_OneDrive("Data/Tagging/Acoustic_tagging/Other researcher's tags/ATAP_South african tags_Nov2013.csv"),stringsAsFactors =F)
South.Oz=read.csv(file=handl_OneDrive("Data/Tagging/Acoustic_tagging/Other researcher's tags/Charlie's tags.csv"),stringsAsFactors =F)

South.Africa=South.Africa[,match(names(WA.Fisheries),names(South.Africa))]
South.Oz=South.Oz[,match(names(WA.Fisheries),names(South.Oz))]

WA.Fisheries$Project.rel=ifelse(WA.Fisheries$Code2%in%South.Oz$Code2,"South.Australia","SMN")
South.Africa$Project.rel="South.Africa"



#convert to date
WA.Fisheries$ReleaseDate2=as.POSIXct.Date(as.Date(WA.Fisheries$ReleaseDate2, "%d-%b-%y"))
South.Africa$ReleaseDate2=as.POSIXct.Date(as.Date(South.Africa$ReleaseDate2, "%d-%b-%y"))
South.Oz$ReleaseDate2=as.POSIXct.Date(as.Date(South.Oz$ReleaseDate2, "%d/%m/%Y"))  

TAGS=rbind(WA.Fisheries,South.Africa)

TAGS$ReleaseLatitude2=as.numeric(TAGS$ReleaseLatitude2)
TAGS$ReleaseLongitude2=as.numeric(TAGS$ReleaseLongitude2)
TAGS$ReleaseLength=as.numeric(TAGS$ReleaseLength)

rm(WA.Fisheries,South.Africa)



