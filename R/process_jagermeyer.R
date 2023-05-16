require(data.table)
require(terra)
require(ggplot2)

options(scipen=999)

DataDir<-"/home/jovyan/common_data/jagermeyer/raw"

plotdir<-"/home/jovyan/common_data/jagermeyer/plots"
if(!dir.exists(plotdir)){dir.create(plotdir)}

Files<-list.files(DataDir,full.names=T)

# Load mapspam
mapspam<-terra::rast(c("/home/jovyan/common_data/mapspam_2010/raw/spam2010V2r0_global_A_MAIZ_A.tif",
             "/home/jovyan/common_data/mapspam_2010/raw/spam2010V2r0_global_A_SOYB_A.tif",
             "/home/jovyan/common_data/mapspam_2010/raw/spam2010V2r0_global_A_WHEA_A.tif",
             "/home/jovyan/common_data/mapspam_2010/raw/spam2010V2r0_global_A_RICE_A.tif"))

names(mapspam)<-c("Maize","Soybean","Wheat","Rice")

mapspam[mapspam<100]<-NA
mapspam[mapspam>=100]<-1

# Load regions layer
admin1<-terra::vect("/home/jovyan/common_data/cg_regions/raw/CGIAR_countries_simplified.shp")
cg_regions<-terra::aggregate(admin1,by="CG_REG",dissolve=T)

# https://www.nature.com/articles/s43016-021-00400-y
Files<-list.files(DataDir)
Files<-data.frame(crop=unlist(data.table::tstrsplit(Files,"_",keep=4)),
           scenario=unlist(data.table::tstrsplit(Files,"_",keep=5)),
           model=unlist(data.table::tstrsplit(Files,"_",keep=6)),
           type=unlist(data.table::tstrsplit(Files,"_",keep=7)),
           path=list.files(DataDir,full.names=T))

area_km<-terra::cellSize(terra::rast(Files$path[1])[[1]],unit="km")
regions<-data.frame(region=cg_regions$CG_REG,ID=1:6)

cuts<-c(5,10,30,50)
years<-2000:2060
years2<-2045:2055
Files<-Files[Files$type=="default",]

reg_dat_all<-rbindlist(lapply(1:nrow(Files),FUN=function(i){
    print(i)
    mask<-mapspam[[Files$crop[i]]]
    dat<-terra::rast(Files$path[i])
    mask<-terra::resample(mask,dat,method="near")
    dat<-dat*mask
    dat<-c(dat[[grep(paste(years,collapse="|"),names(dat))]],area_km)
    names(dat)<-c(paste0("y",years),"area")
    reg_dat<-terra::extract(dat,cg_regions)
    
        reg_gtlt<-do.call("cbind",lapply(cuts,FUN=function(CUT){
            dat_gt<-dat[[paste0("y",years2)]]
            dat_lt<-dat[[paste0("y",years2)]]
            dat_gt[dat_gt<CUT]<-NA
            dat_gt[!is.na(dat_gt)]<-1

            dat_lt[dat_lt>(-CUT)]<-NA
            dat_lt[!is.na(dat_lt)]<-1

            dat_gt<-sum(dat_gt,na.rm=T)
            dat_lt<-sum(dat_lt,na.rm=T)
            
            # This step would be faster with a zonal extraction from cg_regions raster
            reg_dat_gt<-terra::extract(dat_gt,cg_regions)
            reg_dat_lt<-terra::extract(dat_lt,cg_regions)

            X<-data.frame(reg_dat_gt=reg_dat_gt$sum,reg_dat_lt=reg_dat_lt$sum)
            colnames(X)<-paste0(colnames(X),"_",CUT)
            X
        }))
    
    reg_dat<-cbind(reg_dat,reg_gtlt)

    reg_dat$region<-regions$region[match(reg_dat$ID,regions$ID)]
    reg_dat$ID<-NULL
    reg_dat<-data.table(reg_dat)
    reg_dat$crop<-Files$crop[i]
    reg_dat$scenario<-Files$scenario[i]
    reg_dat$model<-Files$model[i]
    reg_dat$type<-Files$type[i]
    reg_dat
}))


reg_dat<-melt(reg_dat_all,id.vars=c("region","area","crop","scenario","model","type"))

reg_dat<-reg_dat[!is.na(area)]
setnames(reg_dat,"variable","year")
reg_dat[,year:=as.numeric(gsub("y","",year))]
reg_dat<-reg_dat[!is.na(year)]
summary<-reg_dat[,list(mean=mean(value,na.rm=T),
                       mean_dec=mean(value[value<0],na.rm=T),
                       mean_inc=mean(value[value>0],na.rm=T),
                       median=median(value,na.rm=T),
                       median_dec=median(value[value<0],na.rm=T),
                       median_inc=median(value[value>0],na.rm=T),
                       min=min(value,na.rm=T),
                       q25=quantile(value,na.rm=T)[2],
                       q50=quantile(value,na.rm=T)[3],
                       q75=quantile(value,na.rm=T)[4],
                       max=max(value,na.rm=T),
                       area_tot=sum(area,na.rm=T),
                       area_dec=sum(area[value<0],na.rm=T),
                       area_dec_5=sum(area[value<(-5)],na.rm=T),
                       area_dec_10=sum(area[value<(-10)],na.rm=T),
                       area_dec_25=sum(area[value<(-25)],na.rm=T),
                       area_dec_50=sum(area[value<(-50)],na.rm=T),
                       area_inc=sum(area[value>0],na.rm=T),
                       area_inc_5=sum(area[value>5],na.rm=T),
                       area_inc_10=sum(area[value>10],na.rm=T),
                       area_inc_25=sum(area[value>25],na.rm=T),
                       area_inc_50=sum(area[value>50],na.rm=T)),
                 by=list(region,year,crop,scenario,model,type)]

summary2<-summary[,list(mean=mean(mean,na.rm=T),
                       mean_dec=mean(mean_dec,na.rm=T),
                       mean_inc=mean(mean_inc,na.rm=T),
                       median=mean(median,na.rm=T),
                       median_dec=mean(median_dec,na.rm=T),
                       median_inc=mean(median_inc,na.rm=T),
                       min_mean=mean(min,na.rm=T),
                       min_min=min(min,na.rm=T),
                       q25=mean(q25),
                       q50=mean(q50),
                       q75=mean(q75),
                       max_max=max(max,na.rm=T),
                       max_mean=mean(max,na.rm=T),
                       area_tot=mean(area_tot,na.rm=T),
                       area_dec=mean(area_dec,na.rm=T),
                       area_dec_5=mean(area_dec_5,na.rm=T),
                       area_dec_10=mean(area_dec_10,na.rm=T),
                       area_dec_25=mean(area_dec_25,na.rm=T),
                       area_dec_50=mean(area_dec_50,na.rm=T),
                       area_inc=mean(area_inc,na.rm=T),
                       area_inc_5=mean(area_inc_5,na.rm=T),
                       area_inc_10=mean(area_inc_10,na.rm=T),
                       area_inc_25=mean(area_inc_25,na.rm=T),
                       area_inc_50=mean(area_inc_50,na.rm=T)),
                 by=list(region,crop,scenario,model,type)]


CI_t <- function (x, ci = 0.95)
{
`%>%` <- magrittr::`%>%`
Margin_Error <- qt(ci + (1 - ci)/2, df = length(x) - 1) * sd(x)/sqrt(length(x))
df_out <- data.frame( sample_size=length(x), Mean=mean(x), sd=sd(x),
Margin_Error=Margin_Error,
'CI lower limit'=(mean(x) - Margin_Error),
'CI Upper limit'=(mean(x) + Margin_Error)) %>%
tidyr::pivot_longer(names_to = "Measurements", values_to ="values", 1:6 )
return(df_out)
}

CI_t(c(1,5,8,6,5,8,54,5,5,8,4,11,1,2,3,4))[5,2]

# Average and variance of models x scenarios by year and crop
summary3<-summary[,list(mean_sd=sd(mean,na.rm=T),
                       mean_low=as.numeric(CI_t(mean)[5,2]),
                        mean_high=as.numeric(CI_t(mean)[6,2]),
                        mean=mean(mean,na.rm=T),
                        median=mean(median,na.rm=T),
                        min_mean=mean(min,na.rm=T),
                        min_min=min(min,na.rm=T),
                        max_max=max(max,na.rm=T),
                        max_mean=mean(max,na.rm=T)),
                 by=list(region,crop,year,type,scenario)]

(g<-ggplot(summary3,aes(x=year,y=mean,ymin=mean_low,ymax=mean_high,col=crop,fill=crop))+
geom_line(size=0.25)+
geom_ribbon(alpha=0.2,col=NA)+
theme_bw()+
theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank())+
facet_grid(scenario~region,scales="free_y")+ 
scale_x_continuous(breaks =c(2020,2060))+
geom_hline(yintercept=0,lty="dashed")+
coord_cartesian(
  xlim = NULL,
  ylim = c(-40,40),
  expand = TRUE,
  default = FALSE,
  clip = "on"
)+
labs(y="Percent change from baseline"))
#facet_wrap(vars(region))

ggsave(filename = "lines.png",
       plot = g,
       path = plotdir,
       width= 189,
       height = 100,
       units = "mm",
       scale = 1.1,
       dpi = 600)


# GOOD/BAD YEARS

reg_dat2<-reg_dat_all[,list(crop,region,scenario,model,area,reg_dat_gt_5,reg_dat_lt_5,reg_dat_gt_10,reg_dat_lt_10,reg_dat_gt_30,reg_dat_lt_30,y2010)]
reg_dat2<-reg_dat2[is.na(y2010),area:=NA][!is.na(area)][,y2010:=NULL]

reg_dat2[,area_dec30:=area*reg_dat_lt_30
          ][,area_inc30:=area*reg_dat_gt_30
           ][,area_dec10:=(area*reg_dat_lt_10)-area_dec30
            ][,area_inc10:=(area*reg_dat_gt_10)-area_inc30
             ][,area_dec5:=(area*reg_dat_lt_5)-area_dec30-area_dec10
              ][,area_inc5:=(area*reg_dat_gt_5)-area_inc30-area_inc10]

reg_dat2<-reg_dat2[,list(area_tot=sum(area*length(years2),na.rm=T),
                         area_dec5=sum(area_dec5,na.rm=T),
                         area_inc5=sum(area_inc5,na.rm=T),
                         area_dec10=sum(area_dec10,na.rm=T),
                         area_inc10=sum(area_inc10,na.rm=T),
                         area_dec30=sum(area_dec30,na.rm=T),
                         area_inc30=sum(area_inc30,na.rm=T)),
                   by=list(crop,scenario,region,model)
                   ][,prop_dec:=round(100*(area_dec5+area_dec10+area_dec30)/area_tot,2)
                    ][,prop_inc:=round(100*(area_inc10+area_inc30)/area_tot,2)
                     ][,area_nochange:=area_tot-area_dec10-area_dec30-area_inc10-area_inc30-area_dec5-area_inc5]

reg_dat2<-reg_dat2[,list(area_dec5=mean(area_dec5,na.rm=T),
                         area_inc5=mean(area_inc5,na.rm=T),
                         area_dec10=mean(area_dec10,na.rm=T),
                         area_inc10=mean(area_inc10,na.rm=T),
                         area_dec30=mean(area_dec30,na.rm=T),
                         area_inc30=mean(area_inc30,na.rm=T),
                          prop_dec=mean(prop_dec,na.rm=T),
                         prop_inc=mean(prop_inc,na.rm=T)),by=list(crop,region,scenario,area_tot)][,area_nochange:=area_tot-area_dec5-area_dec10-area_dec30-area_inc5-area_inc10-area_inc30]


plotdat<-melt(reg_dat2[,list(crop,scenario,region,area_nochange,area_dec5,area_inc5,area_dec10,area_inc10,area_dec30,area_inc30)],
              id.vars=c("crop","scenario","region"),value.name="area")

vals<-c(area_dec30="<-30",area_dec10="-30 to -10",area_dec5="-10 to -5",area_nochange="-5 to 5",area_inc5="5 to 10",area_inc10="10 to 30",area_inc30=">30")

plotdat[,variable:=vals[match(plotdat$variable,names(vals))]]

plotdat[,variable:=factor(variable,levels=as.vector(vals))
       ][,area:=area/10^6]

require(viridis)

(g2<-ggplot(plotdat,aes(x=crop,y=area,fill=variable))+
geom_bar(stat="identity",col="grey20",size=0.2)+
theme_bw()+
theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey70", color = NA),
      panel.border=element_blank())+
facet_grid(region~scenario,scales="free_y")+
coord_flip()+
scale_fill_brewer(type="div")+
scale_y_continuous(expand=c(0,0))+
labs(y="Total area 2045-2055 (M km2)",fill="% yield change")
 )


ggsave(filename = "bars.png",
       plot = g2,
       path = plotdir,
       width= 189,
       height = 80,
       units = "mm",
       scale = 1.3,
       dpi = 600)