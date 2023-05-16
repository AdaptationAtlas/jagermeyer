library(data.table)
library(terra)
library(tidyr)
library(ggplot2)
library(viridis)

# Set options
options(scipen=999)

# Define directories
data_dir <- "/home/jovyan/common_data/jagermeyer/raw"
plot_dir <- "/home/jovyan/common_data/jagermeyer/plots"

# Focal regions
#focal_regions<-c("Ikungi","Simanjiro","Micheweni")
focal_regions<-c("Singida","Pemba North","Manyara")

# Admin Level
admin_level<-"_1"

# Create a directory if it does not exist
if(!dir.exists(plot_dir)){dir.create(plot_dir)}

# Load files from data_dir
Files <- list.files(data_dir,full.names=T)

# Mask data using mapspam
mapspam_mask<-F

if(mapspam_mask){
  # Load mapspam raster
  mapspam <- terra::rast(c("/home/jovyan/common_data/mapspam_2010/raw/spam2010V2r0_global_A_MAIZ_A.tif",
                           "/home/jovyan/common_data/mapspam_2010/raw/spam2010V2r0_global_A_SOYB_A.tif",
                           "/home/jovyan/common_data/mapspam_2010/raw/spam2010V2r0_global_A_WHEA_A.tif",
                           "/home/jovyan/common_data/mapspam_2010/raw/spam2010V2r0_global_A_RICE_A.tif"))

  # Give names to raster layers
  names(mapspam) <- c("Maize","Soybean","Wheat","Rice")

  # Replace all pixel values below 100 with 0, and equal and above 100 with 1
  mapspam[mapspam < 100] <- NA
  mapspam[mapspam >= 100] <- 1
}

# Load regions layer and unzip files of Tanzania shapefile
TZ <- "/home/jovyan/common_data/atlas_boundaries/raw/gadm41/TZA.zip"
Files_unzip <- unzip(TZ, list=T)$Name
Files_unzip <- grep(admin_level, Files_unzip, value=T)
unzip(TZ, files=Files_unzip, exdir="/home/jovyan/common_data/jagermeyer")

# Read Tanzania shapefile
tanzania_shp <- vect(paste0("/home/jovyan/common_data/jagermeyer/gadm41_TZA",admin_level,".shp"))

# Select all regions
focal_shp <- tanzania_shp
focal_shp$ID <- 1:nrow(focal_shp)

# Load crop files information
Files <- list.files(data_dir)
Files <- data.frame(crop=unlist(data.table::tstrsplit(Files,"_",keep=4)),
                    scenario=unlist(data.table::tstrsplit(Files,"_",keep=5)),
                    model=unlist(data.table::tstrsplit(Files,"_",keep=6)),
                    type=unlist(data.table::tstrsplit(Files,"_",keep=7)),
                    path=list.files(data_dir,full.names=T))

# Calculate the area of each cell
area_km <- cellSize(rast(Files$path[1])[[1]], unit="km")

# Specify regions and years of interest
regions <- focal_shp
if(admin_level=="_2"){
  regions$region <- regions$NAME_2
}else{
  regions$region <- regions$NAME_1
}
cuts <- c(5,10,30,50)
years <- 2000:2060
years2 <- 2045:2055

# Filter files dataframe by type, kept as 'default'
Files <- Files[Files$type=="default",]

# Process data for each file
reg_dat_all <- rbindlist(lapply(1:nrow(Files),FUN=function(i){
  # Display progress
  cat('\r                                                                                                                                          ')
  cat('\r',paste0(i,"/",nrow(Files)))
  flush.console()

  # Extract mask and data rasters and rescale them
  dat <- terra::rast(Files$path[i])

  if(mapspam_mask){
    mask <- mapspam[[Files$crop[i]]]
    mask <- terra::resample(mask, dat, method="near")
    dat <- dat*mask
  }

  dat <- c(dat[[grep(paste(years,collapse="|"), names(dat))]], area_km)
  names(dat) <- c(paste0("y",years),"area")
  reg_dat <- terra::extract(dat, regions)

  # Extract data for each cut
  reg_gtlt <- do.call("cbind", lapply(cuts, FUN=function(CUT){
    dat_gt <- dat[[paste0("y",years2)]]
    dat_lt <- dat[[paste0("y",years2)]]

    # Replace all pixel values above CUT with 1, and less than or equal to CUT with NA
    dat_gt[dat_gt<CUT] <- NA
    dat_gt[!is.na(dat_gt)] <- 1

    # Replace all pixel values below -CUT with 1, and greater than or equal to -CUT with NA
    dat_lt[dat_lt>(-CUT)] <- NA
    dat_lt[!is.na(dat_lt)] <- 1

    # Calculate the sum of positive pixel values after replacement
    dat_gt <- sum(dat_gt, na.rm=T)
    dat_lt <- sum(dat_lt, na.rm=T)

    # Extract data for each region
    # This step would be faster with a zonal extraction from cg_regions raster
    reg_dat_gt <- terra::extract(dat_gt, regions)
    reg_dat_lt <- terra::extract(dat_lt, regions)

    # Construct dataframe with the extracted data for each region
    X <- data.frame(reg_dat_gt=reg_dat_gt$sum, reg_dat_lt=reg_dat_lt$sum)
    # Give appropriate column names to the dataframe
    colnames(X) <- paste0(colnames(X),"_",CUT)
    X
  }))

  # Combine extracted data and cut information into a single dataframe
  reg_dat <- cbind(reg_dat, reg_gtlt)

  # Add names to columns
  reg_dat$region <- regions$region[reg_dat$ID]
  reg_dat$ID <- NULL
  reg_dat <- data.table(reg_dat)
  reg_dat$crop <- Files$crop[i]
  reg_dat$scenario <- Files$scenario[i]
  reg_dat$model <- Files$model[i]
  reg_dat$type <- Files$type[i]
  reg_dat
}))

# Transform the dataframe from wide to long format
reg_dat <- melt(reg_dat_all, id.vars=c("region","area","crop","scenario","model","type"))

# Remove all rows with NA in column 'area'
# Column 'area' is the area of the cell binding each row data
reg_dat <- reg_dat[!is.na(area)]

# Rename column 'variable' to 'year'
setnames(reg_dat, 'variable', 'year')

# Convert year column to numeric by removing the year prefix
reg_dat[, year := as.numeric(gsub("y","",year))]

# Remove all rows with NA in column 'year'
reg_dat <- reg_dat[!is.na(year)]

# Summarize extractions by region ####

# create summary data frame grouped by "region", "year", "crop", "scenario", "model", and "type"
summary <- reg_dat[, list(
  mean = mean(value, na.rm = TRUE),
  mean_dec = mean(value[value < 0], na.rm = TRUE),
  mean_inc = mean(value[value > 0], na.rm = TRUE),
  median = median(value, na.rm = TRUE),
  median_dec = median(value[value < 0], na.rm = TRUE),
  median_inc = median(value[value > 0], na.rm = TRUE),
  min = min(value, na.rm = TRUE),
  q25 = quantile(value, na.rm = TRUE)[2],
  q50 = quantile(value, na.rm = TRUE)[3],
  q75 = quantile(value, na.rm = TRUE)[4],
  max = max(value, na.rm = TRUE),
  area_tot = sum(area, na.rm = TRUE),
  area_dec = sum(area[value < 0], na.rm = TRUE),
  area_dec_5 = sum(area[value < (-5)], na.rm = TRUE),
  area_dec_10 = sum(area[value < (-10)], na.rm = TRUE),
  area_dec_25 = sum(area[value < (-25)], na.rm = TRUE),
  area_dec_50 = sum(area[value < (-50)], na.rm = TRUE),
  area_inc = sum(area[value > 0], na.rm = TRUE),
  area_inc_5 = sum(area[value > 5], na.rm = TRUE),
  area_inc_10 = sum(area[value > 10], na.rm = TRUE),
  area_inc_25 = sum(area[value > 25], na.rm = TRUE),
  area_inc_50 = sum(area[value > 50], na.rm = TRUE)),
  by = list(region, year, crop, scenario, model, type)]

fwrite(summary,"s1-region-year-crop-scen-model-type.csv")

# create summarized data frame grouped by "region", "crop", "scenario", "model", and "type"
summary2 <- summary[, list(
  mean = mean(mean, na.rm = TRUE),
  mean_dec = mean(mean_dec, na.rm = TRUE),
  mean_inc = mean(mean_inc, na.rm = TRUE),
  median = mean(median, na.rm = TRUE),
  median_dec = mean(median_dec, na.rm = TRUE),
  median_inc = mean(median_inc, na.rm = TRUE),
  min_mean = mean(min, na.rm = TRUE),
  min_min = min(min, na.rm = TRUE),
  q25 = mean(q25),
  q50 = mean(q50),
  q75 = mean(q75),
  max_max = max(max, na.rm = TRUE),
  max_mean = mean(max, na.rm = TRUE),
  area_tot = mean(area_tot, na.rm = TRUE),
  area_dec = mean(area_dec, na.rm = TRUE),
  area_dec_5 = mean(area_dec_5, na.rm = TRUE),
  area_dec_10 = mean(area_dec_10, na.rm = TRUE),
  area_dec_25 = mean(area_dec_25, na.rm = TRUE),
  area_dec_50 = mean(area_dec_50, na.rm = TRUE),
  area_inc = mean(area_inc, na.rm = TRUE),
  area_inc_5 = mean(area_inc_5, na.rm = TRUE),
  area_inc_10 = mean(area_inc_10, na.rm = TRUE),
  area_inc_25 = mean(area_inc_25, na.rm = TRUE),
  area_inc_50 = mean(area_inc_50, na.rm = TRUE)),
  by = list(region, crop, scenario, model, type)]

fwrite(summary2,"s2-region-crop-scen-model-type.csv")


# function to calculate the confidence interval
CI_t <- function (x, ci = 0.95)
{
  `%>%` <- magrittr::`%>%` # pipe operator
  Margin_Error <- qt(ci + (1 - ci)/2, df = length(x) - 1) * sd(x)/sqrt(length(x)) # calculate margin of error
  df_out <- data.frame( sample_size=length(x), Mean=mean(x), sd=sd(x), # create data frame
                        Margin_Error=Margin_Error, # add margin of error
                        'CI lower limit'=(mean(x) - Margin_Error), # add lower limit
                        'CI Upper limit'=(mean(x) + Margin_Error)) %>% # add upper limit
    tidyr::pivot_longer(names_to = "Measurements", values_to ="values", 1:6 ) # pivot data frame
  return(df_out) # return data frame
}


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

fwrite(summary3,"s3-region-crop-year-type-scenario.csv")

# Plot Results ####

plot_dat<-summary3[region %in% focal_regions]


(g<-ggplot(plot_dat,aes(x=year,y=mean,ymin=mean_low,ymax=mean_high,col=crop,fill=crop))+
geom_line(size=0.25)+
geom_ribbon(alpha=0.1,col=NA)+
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
       path = plot_dir,
       width= 189,
       height = 100,
       units = "mm",
       scale = 1.1,
       dpi = 600)


# Calculate good and bad years ####

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


# Plot Results ####
plot_dat<-plotdat[region %in% focal_regions]

(g2<-ggplot(plot_dat,aes(x=crop,y=area,fill=variable))+
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
       path = plot_dir,
       width= 189,
       height = 80,
       units = "mm",
       scale = 1.3,
       dpi = 600)
