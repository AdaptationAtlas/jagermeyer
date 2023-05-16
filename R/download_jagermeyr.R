DataDir<-"/home/jovyan/common_data"

# Create intermediate directory
SaveDir<-paste0(DataDir,"/jagermeyr/raw")

if(!dir.exists(SaveDir)){
    dir.create(SaveDir,recursive=T)
    }

url <-"http://www.pik-potsdam.de/~jonasjae/GGCMI_Phase3_by_GCM_GGCM_mean/"
doc <- XML::htmlParse(url)
links <- XML::xpathSApply(doc, "//a/@href")
links<-as.vector(grep(".nc4",links,value=T))

URLS<-paste0(url,links)

destfiles<-paste0(SaveDir,"/",links)

options(timeout=480)

for(i in 1:length(URLS)){
    URL<-URLS[i]
    destfile<-destfiles[i]
    # Display progress
    cat('\r                                                                                     ')
    cat('\r',paste0("Downloading file: ",URL))
    flush.console()
    
    if(!file.exists(destfile)){
            download.file(URL, destfile)
        }  
    
}