{library(ggplot2)
library(dplyr)
library(tibble)
library(forcats)
library(ggridges)
library(readr)
library(plotly)
library(MASS)
library(ggpointdensity)
library(ggbeeswarm)
library(magick)
library(viridis)
library(shiny)}
source("https://raw.githubusercontent.com/loganslade/ImgCytometR/refs/heads/main/plotting_functions.R")

#####Import and filter data from cell profiler results######
mycolorv <- viridis(n=11)
mycolorm <- magma(n=11)
mycolorr <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

hci.edu <- read_csv(file.choose())

col <- colnames(hci.edu)
col

#####Config######
channels <- c("DNA", 
              "FITC",
              "mCherry",
              "Cy5")
DNA <- "dna"
Cy5 <- "cdt1_ha"
FITC <- "edu"
mCherry <- "cdt1"


cyto <- F
cyto_var <- "prb1"

#####Re-Name#####

hci.edu.2 <- hci.edu %>% dplyr::rename(Well = "Metadata_Well...6", col = "Metadata_col", row = "Metadata_row", scene = "Metadata_scene...11",
                                       area = "AreaShape_Area", convex="AreaShape_ConvexArea", eccen = "AreaShape_Eccentricity", 
                                       form = "AreaShape_FormFactor", length = "AreaShape_MajorAxisLength", 
                                       min_feret = "AreaShape_MinFeretDiameter", min_length = "AreaShape_MinorAxisLength",
                                       image = "ImageNumber", nneigh = "Neighbors_NumberOfNeighbors_80",
                                       obj_id="Number_Object_Number",
                                       x = "Location_Center_X", y = "Location_Center_Y",
                                       telophase="Math_telophase"
) %>%
  
  {if("DNA" %in% channels)  dplyr::rename(.,
                                          !!paste0("int_",DNA) := "Intensity_IntegratedIntensity_DNA",
                                          !!paste0("int_",DNA,"cor") := "Intensity_IntegratedIntensity_dnacor", 
                                          !!paste0("max_",DNA) := "Intensity_MaxIntensity_DNA",
                                         # !!paste0("max_",DNA,"cor") := "Intensity_MaxIntensity_dnacor",
                                          !!paste0("mean_",DNA) := "Intensity_MeanIntensity_DNA",
                                          !!paste0("mean_",DNA,"cor") := "Intensity_MeanIntensity_dnacor",
                                          !!paste0("file_",DNA) := "FileName_DNA", 
                                          !!paste0("path_",DNA) := "PathName_DNA" 
  ) else . 
  } %>%
  
  {if("FITC" %in% channels)  dplyr::rename(.,
                                           !!paste0("int_",FITC) := "Intensity_IntegratedIntensity_FITC", 
                                           !!paste0("int_",FITC,"cor") := "Intensity_IntegratedIntensity_fitccor", 
                                           !!paste0("mean_",FITC) := "Intensity_MeanIntensity_FITC",
                                           !!paste0("mean_",FITC,"cor") := "Intensity_MeanIntensity_fitccor",
                                           !!paste0("file_",FITC) := "FileName_FITC", 
                                           !!paste0("path_",FITC) := "PathName_FITC") else . 
  } %>%               
  
  {if("mCherry" %in% channels)  dplyr::rename(.,
                                               !!paste0("int_",mCherry) := "Intensity_IntegratedIntensity_mCherry", 
                                               !!paste0("int_",mCherry,"cor") := "Intensity_IntegratedIntensity_mcherrycor", 
                                               !!paste0("mean_",mCherry) := "Intensity_MeanIntensity_mCherry",
                                               !!paste0("mean_",mCherry,"cor") := "Intensity_MeanIntensity_mcherrycor",
                                               !!paste0("file_",mCherry) := "FileName_mCherry", 
                                               !!paste0("path_",mCherry) := "PathName_mCherry") else . 
  } %>%
  
  {if("Cy5" %in% channels)  dplyr::rename(.,
                                          !!paste0("int_",Cy5) := "Intensity_IntegratedIntensity_Cy5", 
                                          !!paste0("int_",Cy5,"cor") := "Intensity_IntegratedIntensity_cy5cor", 
                                          !!paste0("mean_",Cy5) := "Intensity_MeanIntensity_Cy5",
                                          !!paste0("mean_",Cy5,"cor") := "Intensity_MeanIntensity_cy5cor",
                                          !!paste0("file_",Cy5) := "FileName_Cy5", 
                                          !!paste0("path_",Cy5) := "PathName_Cy5") else . 
  } %>%  
  
  {if("Cy3" %in% channels) dplyr::rename(.,
                                         !!paste0("int_",Cy3) := "Intensity_IntegratedIntensity_Cy3", 
                                         !!paste0("int_",Cy3,"cor") := "Intensity_IntegratedIntensity_cy3cor", 
                                         !!paste0("mean_",Cy3) := "Intensity_MeanIntensity_Cy3",
                                         !!paste0("mean_",Cy3,"cor") := "Intensity_MeanIntensity_cy3cor",
                                         !!paste0("file_",Cy3) := "FileName_Cy3", 
                                         !!paste0("path_",Cy3) := "PathName_Cy3") else .
  } 


col2 <- colnames(hci.edu.2)
coltest <- data.frame(col,col2)
coltest

hci.edu.2 <- hci.edu.2 %>% mutate(i_id = paste(image, obj_id, sep="_"))
#####Cyto joining####
if(cyto == T){
 
  hci.cyto <- read_csv(file.choose()) 
  
  cyto.2 <- hci.cyto %>% dplyr::rename(Well = "Metadata_Well...6", col = "Metadata_col", row = "Metadata_row", scene = "Metadata_scene...11",
                                         image = "ImageNumber",
                                         obj_id="Parent_Nuclei",
  ) %>%
    
    {if("DNA" %in% channels)  dplyr::rename(.,
                                            !!paste0("int_",DNA) := "Intensity_IntegratedIntensity_DNA",
                                            !!paste0("int_",DNA,"cor") := "Intensity_IntegratedIntensity_dnacor", 
                                            !!paste0("mean_",DNA) := "Intensity_MeanIntensity_DNA",
                                            !!paste0("mean_",DNA,"cor") := "Intensity_MeanIntensity_dnacor",
                                            
    ) else . 
    } %>%
    
    {if("FITC" %in% channels)  dplyr::rename(.,
                                             !!paste0("int_",FITC) := "Intensity_IntegratedIntensity_FITC", 
                                             !!paste0("int_",FITC,"cor") := "Intensity_IntegratedIntensity_fitccor", 
                                             !!paste0("mean_",FITC) := "Intensity_MeanIntensity_FITC",
                                             !!paste0("mean_",FITC,"cor") := "Intensity_MeanIntensity_fitccor",
    ) else . 
    } %>%               
    
    {if("TexasRed" %in% channels)  dplyr::rename(.,
                                                 !!paste0("int_",TexasRed) := "Intensity_IntegratedIntensity_TexasRed", 
                                                 !!paste0("int_",TexasRed,"cor") := "Intensity_IntegratedIntensity_texasredcor", 
                                                 !!paste0("mean_",TexasRed) := "Intensity_MeanIntensity_TexasRed",
                                                 !!paste0("mean_",TexasRed,"cor") := "Intensity_MeanIntensity_texasredcor",
    ) else . 
    } %>%
    
    {if("Cy5" %in% channels)  dplyr::rename(.,
                                            !!paste0("int_",Cy5) := "Intensity_IntegratedIntensity_Cy5", 
                                            !!paste0("int_",Cy5,"cor") := "Intensity_IntegratedIntensity_cy5cor", 
                                            !!paste0("mean_",Cy5) := "Intensity_MeanIntensity_Cy5",
                                            !!paste0("mean_",Cy5,"cor") := "Intensity_MeanIntensity_cy5cor",
    ) else . 
    } %>%  
    
    {if("Cy3" %in% channels) dplyr::rename(.,
                                           !!paste0("int_",Cy3) := "Intensity_IntegratedIntensity_Cy3", 
                                           !!paste0("int_",Cy3,"cor") := "Intensity_IntegratedIntensity_cy3cor", 
                                           !!paste0("mean_",Cy3) := "Intensity_MeanIntensity_Cy3",
                                           !!paste0("mean_",Cy3,"cor") := "Intensity_MeanIntensity_cy3cor",
    ) else .
    } 
  
  cyto.2 <- cyto.2 %>% mutate(i_id = paste(image, obj_id, sep="_")) %>% dplyr::select(i_id, !!paste0("mean_",cyto_var,"cor"))
  
  hci.edu.cn <- inner_join(hci.edu.2, cyto.2, by = c("i_id" = "i_id"), suffix = c(x="",y="_cyto")) %>% mutate(cn = .data[[paste0("mean_",cyto_var,"cor","_cyto")]]/.data[[paste0("mean_",cyto_var,"cor")]])
  
  
}

#####Config 2#####

plate <- read_csv(file.choose())

hci.edu.p <- {if(cyto == T){inner_join(hci.edu.cn, plate, by = c("Well" = "Well"))}
  else{inner_join(hci.edu.2, plate, by = c("Well" = "Well"))}}

#from excel column names#
varNames <- c("treat_1", "treat_2", "treat_3")
varOrder <- c("order_1", "order_2", "order_3")


hci.edu.e <- hci.edu.p %>% filter(form > 0.77) #form less 0.77 was mostly segmentation errors

hci.edu.e <- hci.edu.e %>% FactorOrder(varNames, varOrder)

for(i in 1:length(varNames)) print(levels(hci.edu.e[[varNames[i]]])) #Check order is correct

hci.edu.e$scenetreat <- interaction(hci.edu.e$Well, hci.edu.e$scene)

if(length(varNames) > 1) {hci.edu.e <- hci.edu.e %>% combinetreat(varNames)
                          hci.edu.e$treat <- as.factor(hci.edu.e$treat)} else{hci.edu.e <- hci.edu.e %>% mutate(treat = .[[varNames[1]]])}
varNames1 <- "treat"

hci.edu.e <- hci.edu.e %>% FactorOrder(varNames1, varOrder)
print(levels(hci.edu.e$treat))

#hci.edu.e <- hci.edu.e %>% filter(ab == "MCM")
#######Normalization of intensity value across experiments######
hci.edu.a <- hci.edu.e %>% 
  corany("int_dnacor", "scenetreat", "quart") %>%  
  corany("int_dna", "scenetreat", "quart") %>% 
  {if("FITC" %in% channels) corany(.,paste0("mean_",FITC,"cor"), "scenetreat", "quant") else .} %>% 
  {if("FITC" %in% channels) corany(.,paste0("mean_",FITC), "scenetreat", "quant") else .} %>% 
  {if("mCherry" %in% channels) corany(.,paste0("mean_",mCherry,"cor"), "scenetreat", "quant") else .}  %>%
  {if("Cy5" %in% channels) corany(.,paste0("mean_",Cy5,"cor"), "scenetreat", "quant") else .} %>% 
  {if("Cy5" %in% channels) corany(.,paste0("mean_",Cy5), "scenetreat", "quant") else .} %>% 
  {if("Cy3" %in% channels) corany(.,paste0("mean_",Cy3,"cor"), "scenetreat", "quant") else .} %>%
  {if(cyto == T) corany(.,paste0("mean_",cyto_var,"cor_cyto"), "scenetreat", "quant") else .}

hci.edu.b <- hci.edu.a %>% align("scenetreat", "int_dnacor", 0.2)

######Telo stats for later#######
hci.edu.c <- hci.edu.b %>% corany("mean_dnacor", "scenetreat", "quant") %>%
  mutate(da = cor_mean_dnacor/area,
         scale_da = da/min(da)) %>% corany("scale_da", "scenetreat", "mean")

######Downsample######
hci.edu.c <- hci.edu.c %>% filter() 
hci.edu.count <- hci.edu.c %>% group_by(.[varNames]) %>% summarise(n())
hci.edu.count
min.count <- min(hci.edu.count$`n()`)
hci.edu.ds <- hci.edu.c %>% group_by(.[varNames]) %>% sample_n(min.count)
if(min.count > 3999) hci.edu.ds2 <- hci.edu.c %>% group_by(.[varNames]) %>% sample_n(4000)

#######QC plots########

al <- hci.edu.ds %>% filter() %>% ggplot(aes(x=image, y=cor_int_dnacor)) +
  geom_pointdensity(size = 0.4, adjust = 0.08) +
  scale_y_continuous()+
  scale_x_continuous()+
  scale_color_gradientn(colours = mycolorv)+
  labs(title = "", x="Image", y="Mean DNA", colour = "Density")+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(cols = vars(NULL))
al
ggplotly(al)

#badimage <- c(105:117, 210:215, 223, 224, 243, 244)
  

dna.dense <- hci.edu.c %>% filter() %>% ggplot(aes(int_dnacor, fill = treat)) 
dna.dense + geom_density(adjust=0.2, alpha=0.3) + 
  theme_classic()+
  labs(title="RPE1_hTERT", x="DNA intensity", y="Density", fill = "Treatment") +  
  scale_x_continuous(limits = c(0, 1000))+
  #scale_fill_manual( values = c("turquoise3", "red1"))+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_wrap(vars(treat))

#####Set Phase gating parameters#####
gg.pre <- hci.edu.ds  %>% ggplot(aes(x=cor_int_dnacor, y=cor_mean_educor))+
  geom_pointdensity(size = 0.1, adjust = 0.05) +
  scale_y_continuous(trans="log10", limits = c(NA,NA))+
  xlim(0, 280)+
  scale_color_gradientn(colours = mycolorv)+
  labs(title = "", x="DNA Integrated Intensity", y="MFI EdU")+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_wrap(vars(treat))
gg.pre
ggplotly(gg.pre)


x1 <- 25
x2 <- 98
x3 <- x2+(x2-x1)
x5 <- 81
x6 <- 120

y1 <- min(hci.edu.c$cor_mean_educor)
y2<- 0.026
y4 <- max(hci.edu.c$cor_mean_educor)

gg.boxes <- hci.edu.ds %>% ggplot(aes(x=cor_int_dnacor, y=cor_mean_educor))
gg.boxes + geom_pointdensity(size = 0.4, adjust = 0.04) +
  scale_y_continuous(trans = "log10")+
  xlim(0, 280)+
  scale_color_gradientn(colours = mycolorv)+
  labs(title = "", x="DNA Integrated Intensity", y="Mean EdU")+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  annotate("rect", xmin= x1, xmax= x2, ymin= y1, ymax= y2, 
           alpha=0, color = "red", linetype = "dotted")+
  annotate("rect", xmin= x1, xmax= x5, ymin= y2, ymax= y4, 
           alpha=0, color = "red", linetype = "dotted")+
  annotate("rect", xmin= x5, xmax= x6, ymin= y2, ymax= y4, 
           alpha=0, color = "red", linetype = "dotted")+
  annotate("rect", xmin= x6, xmax= x3, ymin= y2, ymax= y4, 
           alpha=0, color = "red", linetype = "dotted")+
  annotate("rect", xmin=x2, xmax=x3, ymin=y1, ymax=y2, 
           alpha=0, color = "red", linetype = "dotted")+
  facet_grid(cols = vars(NULL))

#####Gates to cell cycle phase#####
gated.edu.f <- gating(hci.edu.c, "cor_mean_educor", "cor_int_dnacor")

gated.edu.f$phase2 <- factor(gated.edu.f$phase2, 
                             levels = c("G1", "S", "G2"))

gated.edu.f$phase <- factor(gated.edu.f$phase, 
                            levels = c("G1", "ES", "LS", "G2"))

gated.edu.f$phase5 <- factor(gated.edu.f$phase5, 
                             levels = c("G1", "ES","MS", "LS", "G2"))

phase <- gated.edu.f %>% group_by(treat, phase2) %>% summarise(count = n()) %>% mutate(per= prop.table(count) * 100)
phase

g <- phase %>% filter() %>% ggplot(aes(x=treat,y=per))
g + geom_col(aes(fill=fct_rev(phase2)), color = "black", linewidth = 1, width = 0.75)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)))+
  labs(title = "", x="", y= "% in Phase", fill = "Phase")+
  theme_classic()+
  scale_fill_manual(values = c("darkgrey", "darkorange", "deepskyblue"))+
  theme(text = element_text(size=24),
        aspect.ratio=1,
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12, color = "black"),
        axis.text.x = element_text(size=14, angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank())+
  geom_text(
    aes(label = ifelse(round(per,1) >1.9, paste(round(per,1)), paste("")), y = per),
    position = position_stack(vjust = 0.2),
    vjust = 0,
    fontface='bold',
    size=4)+
  facet_grid(cols = vars(NULL))


#####Cell density calculation######
#Remove cells that have neighbor circle on edge#
cmax <- 1024
sub <- 80

gated.edu.s <- gated.edu.f %>% filter(x > 0+sub & y > 0+sub,
                                      x < cmax-sub & y < cmax-sub)

hci.edu.s <- hci.edu.c %>% filter(x > 0+sub & y > 0+sub,
                                  x < cmax-sub & y < cmax-sub)

gated.edu.s <- gated.edu.s %>% mutate(bin = cut_interval(nneigh, length=7.5))

phase <- gated.edu.s %>% filter() %>%
  group_by(bin, phase2) %>% summarise(count = n()) %>% mutate(per= prop.table(count) * 100)
phase

g <- phase %>% filter() %>% ggplot(aes(y=bin,x=per))
g + geom_col(aes(fill=fct_rev(phase2)), color = "black", linewidth = 1, width = 0.75)+
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)))+
  labs(title = "", x="% in Phase", y= "Number of Neighbors", fill = "Phase")+
  theme_classic()+
  scale_fill_manual(values = c("darkgrey", "darkorange", "deepskyblue"))+
  theme(text = element_text(size=24),
        aspect.ratio=1,
        legend.title = element_text(size = 10),
        axis.text = element_text(size=10),
        axis.title = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),)+
  facet_grid(cols = vars(treat_2))

gated.edu.f <- gated.edu.s


#####Telo ID#######
t1 <- 12
t2 <- 25
t3 <- max(hci.edu.c$cor_scale_da)

dd <- hci.edu.ds %>% filter() %>% 
    ggplot(aes(x=cor_int_dnacor, y=cor_scale_da))
dd +   geom_pointdensity(size = 0.6, adjust = 0.1) +
  scale_color_gradientn(colours = mycolorv)+
  scale_x_continuous(limits = c(0, 250))+
  scale_y_continuous()+
  labs(title = "RPE1-hTERT", x="DNA Content", y= "Telophase index")+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
        facet_grid(cols = vars(NULL))+
  annotate("rect", xmin= x1, xmax= x2, ymin= t1, ymax= t2, 
          alpha=0, color = "red", linetype = "dotted")+
  annotate("rect", xmin= x1, xmax= x2, ymin= t2, ymax= t3, 
          alpha=0, color = "red", linetype = "dotted")


late_telo <- 1100
dd <- hci.edu.ds %>% filter() %>% 
  ggplot(aes(x=cor_int_dnacor, y=telophase))
dd +   geom_pointdensity(size = 0.6, adjust = 0.1) +
  scale_color_gradientn(colours = mycolorv)+
  scale_x_continuous(limits = c(0, 250))+
  scale_y_continuous()+
  labs(title = "RPE1-hTERT", x="DNA Content", y= "Telophase index")+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(cols = vars(treat_1))+
  geom_hline(yintercept = late_telo)



telo_values <- hci.edu.c %>% mutate(telo_test = ifelse(cor_scale_da > t1 & between(cor_int_dnacor,x1,x2), "telophase", 
                                                     "not_telo")) 
telo <- telo_values %>% filter() %>% filter(telo_test == "telophase")

#####Get telophases with neighbors and find pairs#######
matches <- telo_id(telo)
      
test <- as.data.frame(matches)
test
pairs2 <- test %>% filter(match == T) %>% mutate(dups = (x*y)+(x_pair*y_pair))
pairs2 <-  pairs2[!duplicated(pairs2$dups),]
pairs2 <- pairs2 %>% dplyr::select(-c("dups"))
pairs2

pair_values <- inner_join(pairs2, telo, by = c("parent"="i_id"), suffix=c("","_pair1"))
pair_values <- inner_join(pair_values, telo, by = c("pair"="i_id"), suffix=c("","_pair2"))

pair.count <- pair_values %>% group_by(.[varNames]) %>% summarise(n())
pair.count



#####Verified Pairs##########
####un-pair####
pairs <- pair_values
obj_1 <- pairs$parent
obj_2 <- pairs$pair
objs <- c(obj_1,obj_2)
#v.unpair <- hci.edu.s %>% filter(i_id %in% objs)
v.unpair <- gated.edu.f %>% filter(i_id %in% objs)
late_telophase <- v.unpair %>% filter(telophase < late_telo)
####Plotting Apps####
shiny2D() #For 2D plots 
shinyViolin() #Violin plots
shinyPercent() #For Gated percent plots


####Manual Violins##### 

colors <- c( "grey77", "deepskyblue3")
%>% filter(telophase < late_telo, treat_2 %in% c("NT", "1uM Palbo"), treat_1 %in% c("EtOH" , "100ng/mL DOX")) 

v.unpair <- v.unpair %>%
  group_by(treat_2) %>%
  mutate(dens = density(cor_mean_cdt1cor))
                        
test <- density(v.unpair$cor_mean_cdt1cor)$y    

#Kernal Density Heatmap violins#
v.unpair <- v.unpair %>%
  group_by(treat) %>%
  mutate(
    dens = {
      d <- density(cor_mean_cdt1cor, na.rm = TRUE)
      approx(d$x, d$y, xout = cor_mean_cdt1cor)$y
    }
  ) %>%
  ungroup()


vio <- v.unpair %>% filter() %>%  ggplot(aes(x=treat_1,y=cor_mean_cdt1cor))+
  geom_violin(linewidth = 1, scale = "width", fill="white", alpha = 1, draw_quantiles = c(0.5))+
  #geom_bar(linewidth = 0.5, aes(fill = treat_2), alpha = 0.5, stat = "summary", fun = "mean")+
  #geom_boxplot(linewidth = 0.5, aes(fill = treat_2), alpha = 0.2, fatten = 1, width = 0.5, outlier.shape = NA)+
  geom_quasirandom(size = 1.1, alpha=0.5, aes(color = dens))+
  scale_color_gradientn(colours = mycolorv)+
  scale_y_continuous(trans = "log10")+
  labs(x="", y="Mean CB-MCM2", title = "Late Telophase") +  
  #scale_fill_manual(values = colors)+
  theme_classic()+
  theme(text = element_text(size=40),
        aspect.ratio = 3,
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size=25),
        axis.text.x = element_text(size=30, angle = 45, hjust=1, color = "black"),
        plot.title = element_text(size=39, hjust = 0.5),
        #axis.line = element_blank(),
        strip.background = element_blank(),
        #panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        legend.position = "none")+
  facet_grid(cols=vars(treat_2))
vio 


vio <- v.unpair %>% filter() %>%  ggplot(aes(treat_1,cor_mean_cdt1cor))+
  geom_violin(linewidth = 1, scale = "width", aes(fill = treat_2), alpha = 0.3, draw_quantiles = c(0.5))+
  ggplot2::stat_density_1d(
    aes(x = cor_mean_cdt1cor, color = after_stat(density)),
    geom = "point",
    size = 2.1,
    alpha = 0.8
  ) +
  scale_color_gradientn(colors = mycolorv)+
  scale_y_continuous(trans = "log10")+
  labs(x="", y="Mean CB-MCM2", title = "Late Telophase") +  
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(text = element_text(size=40),
        aspect.ratio = 3,
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size=25),
        axis.text.x = element_text(size=30, angle = 45, hjust=1, color = "black"),
        plot.title = element_text(size=39, hjust = 0.5),
        strip.background = element_blank(),
        legend.position = "none")+
  facet_grid(cols=vars(treat_2))
vio 

vio <- gated.edu.f %>% filter(phase5 == "G1") %>% ggplot(aes(treat_3,cor_mean_mcmcor))+
  #geom_violin(linewidth = 0.5, scale = "width", aes(fill = treat_2), alpha = 0.3, draw_quantiles = c(0.5))+
  #geom_bar(linewidth = 0.5, aes(fill = treat_2), alpha = 0.5, stat = "summary", fun = "mean")+
  geom_boxplot(linewidth = 0.5, aes(fill = treat_2), alpha = 0.2, fatten = 1, width = 0.5, outlier.shape = NA)+
  geom_quasirandom(size = 0.1, alpha=0.01)+
  scale_y_continuous(trans = "log10")+
  labs(x="", y="Mean Nuclear CDT1", title = "Metaphase") +  
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(text = element_text(size=23),
        aspect.ratio = 3,
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16, angle = 45, hjust=1, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(cols=vars(treat_1))
vio 

ggplotly(vio)

vio <- gated.edu.f %>% filter(phase5 == "G1") %>% ggplot(aes(treat_3,int_mcmcor))+
  geom_violin(linewidth = 0.5, scale = "width", aes(fill = treat_2), alpha = 0.3, draw_quantiles = c(0.5))+
  #geom_bar(linewidth = 0.5, aes(fill = treat_2), alpha = 0.3, stat = "summary", fun = "mean")+
  geom_quasirandom(size = 0.1, alpha=0.1)+
  scale_y_continuous(trans = "identity")+
  labs(x="", y="Mean Nuclear CDT1", title = "Metaphase") +  
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(text = element_text(size=23),
        aspect.ratio = 3,
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16, angle = 45, hjust=1, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(cols=vars(treat_1))
vio 
ggplotly(vio)

means <- gated.edu.f %>% group_by(treat,phase5) %>% summarise(mean = mean(cor_mean_mcmcor-mean_mcmcor_cyto))
mcmmax <- means$mean[2]


vio <- gated.edu.f %>% filter(treat_2 == "NT") %>% ggplot(aes(phase5,cor_mean_mcmcor/mcmmax*100))+
  geom_violin(linewidth = 0.5, scale = "width", aes(fill = treat_2), alpha = 0.3, draw_quantiles = c(0.5))+
  geom_quasirandom(size = 0.5, alpha=0.1)+
  stat_summary(fun=mean, geom="point", shape=21, size=2, color="black", fill="red") +
  scale_y_continuous(trans = "identity")+
  labs(x="", y="Mean Nuclear CDT1", title = "Metaphase") +  
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(text = element_text(size=23),
        aspect.ratio = 1.5,
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16, angle = 45, hjust=1, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(cols=vars(NULL))
vio 


#######Metaphase#########
dd <- gated.edu.s %>% filter() %>% 
  ggplot(aes(x=aint_dnacor, y=cor_scale_da))
dd +   geom_pointdensity(size = 1, adjust = 1) +
  scale_x_continuous(limits = c(10, 110))+
  scale_y_continuous(trans="log2")+
  scale_color_gradientn(colours = mycolorv)+
  labs(title = "RPE1-hTERT", x="DNA Content", y= "Telophase index")+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  annotate("rect", xmin= x2, xmax= x3, ymin= t1, ymax= t2, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin= x2, xmax= x3, ymin= t2, ymax= t3, 
           alpha=0, color = "black", linetype = "dotted")+
  facet_grid(cols = vars(treat))



meta_values <- gated.edu.f %>% mutate(telo_test = ifelse(cor_scale_da > t1 & between(cor_int_dnacor,x2,x3), "telophase", 
                                                         "not_telo")) 
meta <- meta_values %>% filter() %>% filter(telo_test == "telophase")


metavio <- meta %>% filter(!treat_1 %in% c("DMSO", "mTORi")) %>% ggplot(aes(treat,cor_mean_cdt1cor))+
  geom_violin(linewidth = 1, scale = "width", fill = "deepskyblue4", alpha = 0.2, draw_quantiles = c(0.5))+
  geom_quasirandom(size = 2.1, alpha=0.1)+
  scale_y_continuous(trans = "log10")+
  labs(x="", y="Mean Nuclear CDT1", title = "Metaphase") +  
  theme_classic()+
  theme(text = element_text(size=23),
        aspect.ratio = 2,
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16, angle = 45, hjust=1, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(cols=vars(NULL))
metavio  







#################################################
####Everything past this point is for testing####
#################################################
v <- aov(foci ~ treat, data = dplyr::filter(gated.edu.f,phase2=="G2"))
summary(v)
TukeyHSD(v)
#######Manual Pair verification#######
pair.veri.flag <- F
if(pair.veri.flag == T){
  dir <- gsub("\\\\", "/", readClipboard())
  dir
  verified_pairs <- pair_values %>% mutate(path = paste(dir,"/",file_dna,sep="")) %>%
    mutate(x_center=((x+x_pair)/2),
           y_center=(y+y_pair)/2) %>% 
    mutate(cropx = x_center-45,
           cropy = y_center-45)
  
  
  save <- gsub("\\\\", "/", readClipboard())
  setwd(save)
  
  #TeloPhase crops with object annotation 
  for(i in 1:nrow(verified_pairs)){
    
    image <- image_read(verified_pairs$path[i])
    crop <- image_annotate(image, paste0(verified_pairs$parent[i]), size = 5, color = "red", degrees = 0, location = paste0("+",verified_pairs$x[i],"+",verified_pairs$y[i])) 
    crop <- image_annotate(crop, paste0(verified_pairs$pair[i]), size = 5, color = "red", degrees = 0, location = paste0("+",verified_pairs$x_pair[i],"+",verified_pairs$y_pair[i]))
    crop <- image_crop(crop, paste0("90x90+",verified_pairs$cropx[i],"+",verified_pairs$cropy[i]))
    crop <- image_scale(crop, "180x180")
    print(crop)
    print(i)
    verified_pairs$veri_test[i] <- readline("Verified?:")
    if(verified_pairs$veri_test[i] == "no"){break}
  }
  
  #Pair verification 
  pairs <- verified_pairs %>% filter(veri_test %in% c("t"))
  write.csv(pairs, "pairs.csv")}

#####Verified Pair Visualization#######
if(pair.veri.flag == T){ setwd(gsub("\\\\", "/", readClipboard()))}

pairfoldtest <- function(a,b,c,d,e,f){
  a %>% mutate(
    "fold_{e}" := .data[[b]]/d,
    "fold_{e}_pair2" := .data[[c]]/d,
    "{e}_test" := ifelse(.data[[paste0("fold_",e)]] > f & .data[[paste0("fold_",e,"_pair2")]] > f, "plus", "minus"))
}

pairs <- pairs %>% pairfoldtest("cor_mean_cdt1cor", "cor_mean_cdt1cor_pair2",plus, "cdt1", 1)
bins <- dplyr::select(gated.edu.s, bin, i_id) 
pairs <- inner_join(pairs, bins, by=c("parent"="i_id"))

foldtest <- function(a,b,c,d,e){
  a %>% mutate(
    "fold_{d}" := .data[[b]]/c,
    "{d}_test" := ifelse(.data[[paste0("fold_",d)]] > e, "Pos+", "Neg-"))
}


v.unpair <- v.unpair %>% foldtest("cor_mean_mcmcor",plus, "mcm",1.25)


loaded <- v.unpair %>% filter(telophase < late_telo, exp == "a") %>% 
            group_by(bin, mcm_test) %>% summarise(count = n()) %>% mutate(per= prop.table(count) * 100)
loaded


g <- loaded %>% filter() %>% ggplot(aes(bin,y=per))+
  geom_col(aes(fill=fct_rev(mcm_test)), color = "black", linewidth = 1, width = 0.75)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)))+
  labs(title = "Late Telophase", x="", y= "Percent", fill = "CB-MCM2")+
  scale_fill_manual(values = c("darkorange1","deepskyblue"))+
  theme_classic()+
  theme(text = element_text(size=23),
        aspect.ratio=1.75,
        legend.title = element_text(size = 16),
        axis.text = element_text(size=18, color = "black"),
        axis.text.x = element_text(size=14, angle = 45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank())+
  geom_text(
    aes(label = ifelse(round(per,1) >1.9, paste(round(per,1)), paste("")), y = per),
    position = position_stack(vjust = 0.2),
    vjust = 0,
    fontface='bold',
    size=5)+
  facet_grid(cols = vars(NULL))
g 


dd <- v.unpair %>% filter(telophase < late_telo) %>%
  ggplot(aes(x=telophase, y=cor_mean_cdt1cor))+
  geom_jitter(aes(fill = telophase, size = area), color = "black", pch = 21, alpha=0.9)+
  scale_fill_gradient2(low = "white", mid = "white", high = "magenta", na.value = "magenta", 
                       midpoint = 200, limits = c(0,3000))+
  scale_size(
    range = c(1, 3),
  )+
  scale_x_continuous(limits = c(NA, NA), transform = "identity")+
  scale_y_continuous(limits = c(NA, NA), transform = "identity")+
  labs(title = "Telophase", x="Telophase index", y= "phos-RB Cyto/Nuc", fill="Telophase \nIndex", size="Nuclear \nArea")+
  theme_classic()+
  theme(text = element_text(size=23),
        aspect.ratio=1,
        legend.title = element_text(size = 14),
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size=16, color="black"),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  geom_smooth(method = "lm", se = T, color = "black", linetype = 2, linewidth =1)+
  facet_grid(cols = vars(treat), rows = vars(NULL))
dd

v.unpair <- v.unpair %>% mutate(bin_t = cut_interval(telophase, length=400))


telovio <- v.unpair %>% filter(telophase < late_telo) %>% ggplot(aes(treat,cor_mean_mcmcor))+
  geom_violin(linewidth = 1, scale = "width", fill = "deepskyblue4", alpha = 0.2, draw_quantiles = c(0.5))+
  geom_quasirandom(size = 1, alpha=0.1)+
  scale_y_continuous(trans = "identity")+
  labs(x="", y="Mean MCM2", title = "Late Telophase") +  
  theme_classic()+
  theme(text = element_text(size=23),
        aspect.ratio = 2,
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16, angle = 45, hjust=1, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(cols=vars(NULL))
telovio  


cdt1box <- v.unpair %>% ggplot(aes(treat, mean_cdt1cor))
cdt1box + geom_boxplot(fill = "grey65", colour = "black",  
                      outlier.shape = NA, alpha=0.75, 
                      fatten = 0.1, lwd = 0.5, width = 0.5, notch = TRUE)+
  labs(x="", y="MFI phos-RB1", title = "RPE-hTERT") +  
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))

cdt1dense <- v.unpair %>% filter(treat_2 %in% c("DMSO")) %>% ggplot(aes(cor_mean_cdt1cor_cyto, fill=treat_1)) 
cdt1dense +   geom_histogram(aes(y=after_stat(c(
  count[group==1 & PANEL == 1]/sum(count[group==1 & PANEL == 1]),
  count[group==2 & PANEL == 1]/sum(count[group==2 & PANEL == 1]))*100)), position = "dodge") + 
  theme_classic()+
  labs(title="RPE-hTERT:Telophases", x="MFI CDT1", y="Percent", fill = "Treatment") +  
  scale_x_continuous(trans="log10")+
  #scale_fill_manual( values = c("turquoise3", "red1"))+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))
  facet_grid(cols=vars(exp))
  

telohisto <- v.unpair %>% filter(telophase < 500) %>% ggplot(aes(cor_mean_cdt1cor, fill=treatment)) 
telohisto +   geom_histogram(aes(y=after_stat(c(
  count[group==1 & PANEL == 1]/sum(count[group==1  & PANEL == 1]),
  count[group==2  & PANEL == 1]/sum(count[group==2  & PANEL == 1]),
  count[group==3  & PANEL == 1]/sum(count[group==3  & PANEL == 1]),
  count[group==4 & PANEL == 1]/sum(count[group==4  & PANEL == 1]))*100)), position = "dodge", binwidth = 0.03) + 
  theme_classic()+
  labs(title="Late Telophase", x="Mean phos-RB", y="Percent", fill = "Treatment") +  
  scale_x_continuous()+
  scale_fill_manual( values = c("cyan", "magenta1", "cyan4", "magenta4"))+
  theme_classic()+
  theme(text = element_text(size=20),
        aspect.ratio = 1,
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))
  facet_grid(cols = vars(treatment))

#####Rank#####
v.unpair$rank <- rank(v.unpair$cor_mean_cdt1cor)
#####Cropping engine#######
dir <- gsub("\\\\", "/", readClipboard())
dir

crop_set <- function(df, file, cx, cy, mode){
crops <- df %>% mutate(.,path = paste(dir,"/",.data[[file]],sep=""))  %>%
          {if(mode == "single") 
          mutate(cropx = x-cx,
                 cropy = y-cy) else .} %>%
           { if(mode=="pair") .
              else .
            }

crops.f <- crops %>% filter(cropx > 0 & cropy >0,
                            cropx < 954 & cropy < 954)
} #add a crop center mode for telophase pairs, formula should be in prior versions#



crop.count <- pairs %>% group_by(bin) %>% summarise(n())
crop.count
crop.sample <- pairs %>% filter(bin != "(37.5,45]") %>%
                    group_by(bin) %>% sample_n(10)

crops.cdt1 <- crop.sample %>% filter() %>% 
  crop_set(., "file_cdt1", 35, 35, "single")

crops.cdna <- crop.sample %>% filter() %>% 
  crop_set(., "file_dna", 35, 35)

save <- gsub("\\\\", "/", readClipboard())
save
setwd(save)

for(i in 1:nrow(crop.sample)){
  
  image <- image_read(crop.sample$path[i])
  crop <- image_crop(image, paste0("70x70+",crop.sample$cropx[i],"+",crop.sample$cropy[i])) 
  crop <- image_scale(crop, "105x105")
  image_write(crop,paste0(crop.sample$bin[i],"_",crop.sample$i_id[i],".tif")) }




#Rectangles
images_seq <- unique(pairs$image)
for(i in images_seq){
  print(i)
  filt <- pairs %>% filter(image == i)
  image <- image_read(filt$path[1])
  image <- image_draw(image)
  
  for (o in 1:nrow(filt)){
    rect((filt$x_center[o]-100), (filt$y_center[o]-100), (filt$x_center[o]+100), (filt$y_center[o]+100), 
         border = "white", lty = "dashed", lwd = 2)
  }
  image_write(image,paste0(filt$sirna2[1],"_",filt$exp[1],"_",filt$scene[1],".tif"), compression="Lossless")
  remove(image,filt)
  graphics.off() 
}
#Montage#
save <- "C:/Users/lslad/Desktop/Research/UNC-Chapel Hill/2-26-cbcdt1-phospho-RB-Nutlin-3/LS_mVenus_GFP_6w_20x_092223_RPE_cdt1_phosRB_NL3_2/analysis/cdt1 crops"
setwd(save)

list <- list.files(save)
list_image <- image_read(list)
img <- image_scale(list_image, "160x160")
mont <- image_montage(list_image, "160x160+1+1", tile = '16x15', bg = "grey11")
image_write(mont, "montage1.png")


#######Metaphase#########
dd <- gated.edu.s %>% filter() %>% 
  ggplot(aes(x=aint_dnacor, y=cor_scale_da))
dd +   geom_pointdensity(size = 1, adjust = 1) +
  scale_x_continuous(limits = c(10, 110))+
  scale_y_continuous(trans="log2")+
  scale_color_gradientn(colours = mycolorv)+
  labs(title = "RPE1-hTERT", x="DNA Content", y= "Telophase index")+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  annotate("rect", xmin= x2, xmax= x3, ymin= t1, ymax= t2, 
           alpha=0, color = "black", linetype = "dotted")+
  annotate("rect", xmin= x2, xmax= x3, ymin= t2, ymax= t3, 
           alpha=0, color = "black", linetype = "dotted")+
  facet_grid(cols = vars(treat))



meta_values <- gated.edu.s %>% mutate(telo_test = ifelse(cor_scale_da > t1 & between(cor_int_dnacor,x2,x3), "telophase", 
                                                       "not_telo")) 
meta <- meta_values %>% filter() %>% filter(telo_test == "telophase")


metavio <- meta %>% filter() %>% ggplot(aes(bin,cor_mean_cdt1cor_cyto))+
  geom_violin(linewidth = 1, scale = "width", fill = "deepskyblue4", alpha = 0.2, draw_quantiles = c(0.5))+
  geom_quasirandom(size = 1, alpha=0.1)+
  scale_y_continuous(trans = "log10")+
  labs(x="", y="Mean Nuclear CDT1", title = "Metaphase") +  
  theme_classic()+
  theme(text = element_text(size=23),
        aspect.ratio = 1.5,
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16, angle = 45, hjust=1, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(cols=vars(NULL))
metavio  



metahisto <- meta %>% filter() %>% ggplot(aes(cor_mean_cdt1cor, fill=treat))+
  geom_histogram(aes(y=after_stat(c(
  count[group==1]/sum(count[group==1]),
  count[group==2]/sum(count[group==2]),
  count[group==3]/sum(count[group==3]),
  count[group==4]/sum(count[group==4]))*100)), position = "dodge") + 
  theme_classic()+
  labs(title="RPE-hTERT: Metaphases", x="MFI phos-RB (Ser807/811)", y="Percent", fill = "Treatment") +  
  scale_x_continuous(trans="identity")+
  #scale_fill_manual( values = c("turquoise3", "red1"))+
  theme_classic()+
  theme(text = element_text(size=18),
        aspect.ratio = 1,
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))
metahisto

t.test(mean_cdt1cor~treat, data=meta)

meta %>% group_by(treat) %>% summarise(n())

dir <- gsub("\\\\", "/", readClipboard())

verified_metas <- meta %>% mutate(path = paste(dir,"/",file_dna,sep="")) %>%
  mutate(cropx = x-40,
         cropy = y-40)


save <- gsub("\\\\", "/", readClipboard())
setwd(save)

#metaphase crops with object annotation 
for(i in 1:nrow(verified_metas)){
  
  image <- image_read(verified_metas$path[i])
  crop <- image_annotate(image, paste0(verified_metas$obj_id[i]), size = 5, color = "red", degrees = 0, location = paste0("+",verified_metas$x[i],"+",verified_metas$y[i])) 
  crop <- image_crop(crop, paste0("80x80+",verified_metas$cropx[i],"+",verified_metas$cropy[i]))
  crop <- image_scale(crop, "160x160")
  print(crop)
  print(i)
  verified_metas$veri_test[i] <- readline("Verified?:")
  if(verified_metas$veri_test[i] == "no"){break}
}

#meta verification 
metas <- verified_metas %>% filter()
write.csv(metas, "nl3_verified_metas.csv")
metas <- read_csv(file.choose())

meta_obj <- meta$i_id

metas <- meta %>% mutate(
  fold_cdt1 = mean_cdt1cor/cdt1_s,
  cdt1test = ifelse(fold_cdt1 > 2, "plus", "minus"))

loaded <- metas %>% group_by(treat, cdt1test) %>% summarise(count = n()) %>% mutate(per= prop.table(count) * 100)
loaded

g <- loaded %>% filter() %>% ggplot(aes(treat,y=per))
g + geom_col(aes(fill=fct_rev(cdt1test)), color = "black", linewidth = 1, width = 0.75)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)))+
  labs(title = "", x="", y= "Percent", fill = "Loaded")+
  theme_classic()+
  scale_fill_manual(values = c("darkorange1","blue1"))+
  theme(text = element_text(size=16),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black", linewidth =0.75),
        # strip.background = element_blank(),
        #panel.border = element_rect(linetype = "dashed", fill = NA)
  )


cdt1vio <- meta %>% filter() %>% ggplot(aes(treatment, mean_mcmcor))
cdt1vio + geom_violin(linewidth = 1, scale = "width", fill = "deepskyblue4", alpha = 0.2, draw_quantiles = c(0.5))+
  geom_quasirandom(size = 2, alpha=0.1)+
  scale_y_continuous(trans="log10")+
  labs(x="", y="Mean MCM2", title = "RPE-hTERT:Metaphase") +  
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
facet_grid(cols=vars(time))



metas$rank <- rank(metas$mean_cdt1cor)

#cdt1 crops
save <- gsub("\\\\", "/", readClipboard())
setwd(save)

cdt1dir <- gsub("\\\\", "/", readClipboard())
metas <- metas %>% mutate(path_cdt1 = paste(cdt1dir,"/",file_cdt1,sep=""))  %>%
  mutate(cropx = x-35,
         cropy = y-35)

for(i in 1:nrow(metas)){
  
  image <- image_read(metas$path[i])
  crop <- image_crop(image, paste0("70x70+",metas$cropx[i],"+",metas$cropy[i])) 
  crop <- image_scale(crop, "105x105")
  image_write(crop,paste0(metas$treat[i],"_",metas$rank[i],"_",metas$obj_id[i],".tif")) }



#####mitosis label in all######
meta_obj <- meta$i_id
hci.edu.c <- hci.edu.c %>% mutate(m_phase = ifelse(i_id %in% objs, "telophase",
                                                   ifelse(i_id %in% meta_obj, "metaphase", "interphase")))


hci.edu.c <- hci.edu.c %>% mutate(m_phase = ifelse(i_id %in% objs, "telophase",
                                                   "interphase"))


#####Remove cells that have neighbor circle on edge########
cmax <- 1024
sub <- 80

gated.edu.s <- gated.edu.f %>% filter(x > 0+sub & y > 0+sub,
                                      x < cmax-sub & y < cmax-sub)

hci.edu.s <- hci.edu.c %>% filter(x > 0+sub & y > 0+sub,
                                  x < cmax-sub & y < cmax-sub)

#####Cell density calculation######
gated.edu.s <- gated.edu.s %>% mutate(bin = cut_interval(nneigh, length=7.5))
gated.edu.f <- gated.edu.f %>% mutate(bin = cut_interval(nneigh, length=7.5))

phase <- gated.edu.s %>% filter() %>%
  group_by(bin,treat, phase2) %>% summarise(count = n()) %>% mutate(per= prop.table(count) * 100)
phase

g <- phase %>% filter() %>% ggplot(aes(y=bin,x=per))
g + geom_col(aes(fill=fct_rev(phase2)), color = "black", linewidth = 1, width = 0.75)+
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)))+
  labs(title = "", x="% in Phase", y= "Number of Neighbors", fill = "Phase")+
  theme_classic()+
  scale_fill_manual(values = c("darkgrey", "darkorange", "deepskyblue"))+
  theme(text = element_text(size=24),
        aspect.ratio=1,
        legend.title = element_text(size = 10),
        axis.text = element_text(size=10),
        axis.title = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),)+
  facet_grid(cols = vars(treat))

#####2d test plots#####
facet <- T
gate_line <- F
plus <- 0.03
max <- 0.13

fde <- gated.edu.f %>% filter(phase2 == "G1") %>% ggplot(aes(x=cor_mean_cdt1cor, y=cor_mean_educor)) +
  geom_pointdensity(size = 0.01, adjust = 0.01) +
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  #xlim(10,120)+
  scale_color_gradientn(colours = mycolorv)+
  labs(title = "RPE1-hTERT", x="DNA Content", y="Mean CDT1", colour = "Density")+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        aspect.ratio=1,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  {if(facet == T)  facet_wrap(vars(treat))}+
  {if(gate_line == T) geom_hline(yintercept = plus)}
fde


gatecol <- hci.edu.c %>% filter(treat == "50k") %>% ggplot(aes(x=cor_int_dnacor, y=mean_cdt1cor_cyto)) 
gatecol +  geom_point(aes(color = m_phase), size = 1, alpha=0.3) +
  scale_y_continuous(trans = "log10")+
  xlim(20, 180)+
  scale_color_manual(values = c("#E0E0E0", "magenta", "cyan1"))+
  labs(title = "RPE1-hTERT", x="DNA Content", y="phos-RB Cyto/Nuc", colour = "Phase")+
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        legend.key = element_rect(color = NA),
        aspect.ratio=NULL,
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(cols = vars(NULL))

#####Shading Correction Tester#####

gradcol<- gated.edu.f %>% filter(phase2 == "G1") %>% ggplot(aes(x=x, y=y))+
  geom_jitter(aes(color = int_dnacor), alpha=0.9, size = 1.6)+
  scale_y_continuous(trans = "identity")+
  scale_color_gradientn(colours = mycolorv)+
  #xlim(0, 1024)+
  labs()+
  theme_classic()+
  theme(text = element_text(size=18),
        aspect.ratio = 1,
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))
gradcol

ggplotly(gradcol)

####1D Plots####
dg <- gated.edu.f %>% filter(phase5 %in% c("G1", "ES")) %>% ggplot(aes(cor_mean_mcmcor, fill = treat))
dg  + geom_density(adjust=1.2, alpha=0.3) + 
  theme_classic()+
  labs(x="Mean MCM2", y="Density", fill = "Treatment") +  
  scale_x_continuous(trans="log10")+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(rows = vars(phase))

line <- means %>% filter(time != "NT") %>% ggplot(aes(time, cor_mean_mcmcor, group=treat_1))
line + geom_line(colour = "black", linewidth=0.7)+
  geom_point(aes(colour = time_pd), size = 5)+
  #scale_colour_manual(values = c("blue1","darkorange1","darkgrey"))+
  labs(x=" Time in MG132", y="MFI CDT1", title = "RPE-hTERT", colour="Time in Palbo") +  
  theme_classic()+
  theme(text = element_text(size=18),
        legend.title = element_text(size = 12),
        axis.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1))+
  facet_grid(cols = vars(phase2))


