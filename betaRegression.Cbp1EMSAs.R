library(dplyr)
require(ggplot2)
require(ggbeeswarm)
require(cowplot)
require(viridis)
require(PMCMRplus) #for Dunnett's test
require(betareg) #for beta regression to test interaction
require(lmtest) #Likelihood ratio test for beta regression

rm(list=ls())

# data import
# multiple replicates
data.WT <- read.csv("data/gel_quantifications/Cbp1_WT_EMSA.csv")
data.dHTH1 <- read.csv("data/gel_quantifications/dHTH1_EMSA.csv")
data.dHTH3 <- read.csv("data/gel_quantifications/dHTH3_EMSA.csv")

#change order for templates
data.WT$Template <- factor(data.WT$Template, 
                           levels=c("WT", "A6/7C","A13/14C", "A17/18C", "T19/20G", "A22/23C"))
data.dHTH1$Template <- factor(data.dHTH1$Template, 
                           levels=c("WT", "A6/7C","A13/14C", "A17/18C", "T19/20G", "A22/23C"))
data.dHTH3$Template <- factor(data.dHTH3$Template, 
                           levels=c("WT", "A6/7C","A13/14C", "A17/18C", "T19/20G", "A22/23C"))


#calculate mean and sd for each group, WT Cbp1
data.WT.m <- data.WT[data.WT$Conc == 20,] %>% group_by(Template) %>% dplyr::summarize(mean(Fraction.Bound))
colnames(data.WT.m)<-c("Template", "Mean")

data.WT.sd <- data.WT[data.WT$Conc == 20,] %>% group_by(Template) %>% dplyr::summarize(sd(Fraction.Bound))
colnames(data.WT.sd)<-c("Template", "SD")

data.WT.m$SD <- data.WT.sd$SD


#plot
gb <- ggplot() + geom_quasirandom(data=data.WT[data.WT$Conc == 20,], aes(x=Template, y=Fraction.Bound), 
                            colour="darkgrey", alpha = 0.4) 

gb <- gb + geom_errorbar(data=data.WT.m, aes(x=Template, ymin=Mean-SD, ymax= Mean+SD), width=0.1, col="darkgrey")


gb <- gb + geom_point(data=data.WT.m, aes(x=Template, y=Mean, fill=Template), 
                      shape=21, size=3, col="darkgrey", alpha=0.8)

gb <- gb + scale_fill_viridis_d(option="plasma")
gb <- gb + theme_bw() 
gb <- gb + ylab("Fraction Bound") + ggtitle("WT Cbp1") + scale_y_continuous(limits =c(0,1), expand = c(0,0), position = "right")
gb <- gb + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                 panel.border = element_rect(fill = NA, colour = "grey20"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 legend.position = "none"
                 )
gb


#for dHTH1
data.dHTH1.m <- data.dHTH1[data.dHTH1$Conc == 160,] %>% group_by(Template) %>% dplyr::summarize(mean(Fraction.Bound))
colnames(data.dHTH1.m)<-c("Template", "Mean")

data.dHTH1.sd <- data.dHTH1[data.dHTH1$Conc == 160,] %>% group_by(Template) %>% dplyr::summarize(sd(Fraction.Bound))
colnames(data.dHTH1.sd)<-c("Template", "SD")

data.dHTH1.m$SD <- data.dHTH1.sd$SD

#bargraph dHTH1
gb.1 <- ggplot() + geom_quasirandom(data=data.dHTH1[data.dHTH1$Conc == 160,], aes(x=Template, y=Fraction.Bound), 
                                  colour="darkgrey", alpha = 0.4) 
gb.1 <- gb.1 + geom_errorbar(data=data.dHTH1.m, aes(x=Template, ymin=ifelse(Mean-SD>0, Mean-SD, 0), ymax= Mean+SD), width=0.1, col="darkgrey")


gb.1 <- gb.1 + geom_point(data=data.dHTH1.m, aes(x=Template, y=Mean, fill=Template), 
                          shape=21, size=3, col="darkgrey", alpha=0.8)
gb.1 <- gb.1 + scale_fill_viridis_d(option="plasma")
gb.1 <- gb.1 + theme_bw() 
gb.1 <- gb.1 + ylab("Fraction Bound") + ggtitle("dHTH1") + scale_y_continuous(limits =c(0,1), expand = c(0,0), position = "right")
gb.1 <- gb.1 + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                 panel.border = element_rect(fill = NA, colour = "grey20"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 legend.position = "none"
)
gb.1



#for dHTH3
data.dHTH3.m <- data.dHTH3[data.dHTH3$Conc == 80,] %>% group_by(Template) %>% dplyr::summarize(mean(Fraction.Bound))
colnames(data.dHTH3.m)<-c("Template", "Mean")

data.dHTH3.sd <- data.dHTH3[data.dHTH3$Conc == 80,] %>% group_by(Template) %>% dplyr::summarize(sd(Fraction.Bound))
colnames(data.dHTH3.sd)<-c("Template", "SD")

data.dHTH3.m$SD <- data.dHTH3.sd$SD


#graph dHTH3
gb.3 <- ggplot() + geom_quasirandom(data=data.dHTH3[data.dHTH3$Conc == 80,], aes(x=Template, y=Fraction.Bound), 
                                    colour="darkgrey", alpha = 0.4) 
gb.3 <- gb.3 + geom_errorbar(data=data.dHTH3.m, aes(x=Template, ymin=Mean-SD, ymax= Mean+SD), width=0.1, col="darkgrey")


gb.3 <- gb.3 + geom_point(data=data.dHTH3.m, aes(x=Template, y=Mean, fill=Template), 
                          shape=21, size=3, col="darkgrey", alpha=0.8)
gb.3 <- gb.3 + scale_fill_viridis_d(option="plasma")
gb.3 <- gb.3 + theme_bw() 
gb.3 <- gb.3 + ylab("Fraction Bound") + ggtitle("dHTH3") + scale_y_continuous(limits =c(0,1), expand = c(0,0), position = "right")
gb.3 <- gb.3 + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                     panel.border = element_rect(fill = NA, colour = "grey20"),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     legend.position = "none"
)
gb.3


#align plots vertically

gb.grid <- plot_grid(gb, gb.1, gb.3, align="v", ncol=1)
gb.grid

save_plot("plots/EMSA_mutants.pdf", gb.grid, ncol=1, base_width=2.5, base_height=7.5)


#Testing effect of template mutations and interaction between template mutations and HTH deletions in Cbp1

#making data frame with selected protein concentrations yielding ~0.5 of DNA bound
data.WT.an <- data.WT[data.WT$Conc ==20,]
data.dHTH1.an <- data.dHTH1[data.dHTH1$Conc ==160,]
data.dHTH3.an <- data.dHTH3[data.dHTH3$Conc ==80,]
data.WT.an$prot <- "WT"
data.dHTH1.an$prot <- "dHTH1"
data.dHTH3.an$prot <- "dHTH3"
data.an <- rbind(data.WT.an, data.dHTH1.an, data.dHTH3.an)
data.an$prot <- as.factor(data.an$prot)

#changing order of proteins
data.an$prot <- factor(data.an$prot, levels=c("WT", "dHTH1","dHTH3")) #change order
levels(data.an$Template)


#Dunnett test for to test which templates are differently bound for WT
dunnett.output.WT <- PMCMRplus::dunnettTest(x = data.an$Frac[data.an$prot == "WT"],
                            g = data.an$Template[data.an$prot == "WT"],
                            alternative="less")
dunnett.output.WT
  

#beta regression with Likelihood ratio test for individual interaction terms

lrtest.beta <- data.frame(Protein.variant = vector(), Template.mutation = vector(), p.val.interaction = vector()) #data frame to collect p values for each combination

for(j in 2:3){
  for(i in 2:length(levels(data.an$Template))){
    d <- data.an[data.an$prot %in% levels(data.an$prot)[c(1,j)] & data.an$Template %in% levels(data.an$Template)[c(1,i)],]
    w <-betareg(data = d, 
            Fraction.Bound ~ prot+Template, link="log") #without interaction term
    x <-betareg(data = d, 
            Fraction.Bound ~ prot*Template, link="log") #with interaction term
    y <- lrtest(w,x) #extact the p-value for the interaction term
    z <- c(levels(data.an$prot)[j], levels(data.an$Template)[i], y$`Pr(>Chisq)`[2])
    
    lrtest.beta[(j-2)*5 + i-1,] <- z
  }
}
lrtest.beta$p.adj <- p.adjust(lrtest.beta$p.val.interaction, method= "bonferroni") #Bonferroni multiple testing correction
lrtest.beta
