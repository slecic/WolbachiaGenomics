

library(tidyverse)
library(ggplot2)

# FST per position
fst <- read_tsv("/home/slecic/PhD/Transect/fst/brno_vs_vienna_allsampels.50miss.mac3.q30.dp30-300.95miss.maf05.meandp20.weir.fst")
ggfst <- ggplot(fst, aes(POS, WEIR_AND_COCKERHAM_FST)) + geom_point()
quantile(fst$WEIR_AND_COCKERHAM_FST, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(fst$WEIR_AND_COCKERHAM_FST, 0.975, na.rm = T)
# make an outlier column in the data.frame
fst <- fst %>% mutate(outlier = ifelse(WEIR_AND_COCKERHAM_FST > my_threshold, "outlier", "background"))
fst %>% group_by(outlier) %>% tally()
ggplot(fst, aes(POS, WEIR_AND_COCKERHAM_FST, colour = outlier)) + geom_point()


myout95 <- subset(fst, outlier == "outlier")

# plot the results per site
plotfst<-function(x){
  require("ggplot2")
  dat<-read_tsv(x)
  print(head(dat))
  dat$POS<-1:length(dat$POS)
  print(head(dat))
  #print(dat)
  my_threshold <- quantile(dat$WEIR_AND_COCKERHAM_FST, 0.975, na.rm = T)
  print(my_threshold)
  dat <- dat %>% mutate(outlier = ifelse(dat$WEIR_AND_COCKERHAM_FST > my_threshold, "outlier", "background"))
  dat %>% group_by(outlier) %>% tally()
  print(head(dat))
  theplot<-ggplot(dat, aes(POS, WEIR_AND_COCKERHAM_FST , colour = outlier))+geom_point()+theme_grey(15)+labs(x="SNP index", y="FST")+
    scale_color_manual(values = c("grey50", "red"))
  pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
  ggsave(filename=pngName, width=20, height=5, units="in", theplot)
}
plotfst("/home/slecic/PhD/Transect/fst/brno_vs_vienna_allsampels.50miss.mac3.q30.dp30-300.95miss.maf05.meandp20.weir.fst")


# plot the Fst values over 10kb windows 
plotfst<-function(x){
  require("ggplot2")
  dat<-read_tsv(x)
  print(head(dat))
  dat$BIN_START<-1:length(dat$BIN_START)
  print(head(dat))
  #print(dat)
  my_threshold <- quantile(dat$WEIGHTED_FST, 0.975, na.rm = T)
  print(my_threshold)
  dat <- dat %>% mutate(outlier = ifelse(dat$WEIGHTED_FST > my_threshold, "outlier", "background"))
  dat %>% group_by(outlier) %>% tally()
  print(head(dat))
  theplot<-ggplot(dat, aes(BIN_START, WEIGHTED_FST , colour = outlier))+geom_point()+theme_grey(15)+labs(x="SNP index", y="FST")+
    scale_color_manual(values = c("grey40", "purple"))
  pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
  ggsave(filename=pngName, width=20, height=5, units="in", theplot)
}
plotfst("/home/slecic/PhD/Transect/fst/brno_vs_vienna_allsampels.50miss.mac3.q30.dp30-300.99miss.maf05.meandp20.win10kb.windowed.weir.fst")


### plot the distribution of mean allele frequencies
afs <-read_tsv("/home/slecic/PhD/Transect/fst/var.full-genome-0.1.NoWol.snp.50miss.mac3.q30.dp30-300.95miss.maf05.meandp20.vcf.recode.af.frq.tab") 
names(afs) <- c("contig", "pos", "Nalleles", "N_contig", "afMajor", "afMinor")
hist(afs$afMinor, xlim = c(0,1))
min(afs$afMinor) # 0.05
max(afs$afMinor) # 0.6

afs0109 <- subset(afs, afMinor >= 0.1)
subset(afs, afMinor < 0.1)

viennaaf <- read_tsv("/home/slecic/PhD/Transect/fst/plink2.95miss.vienna.freq.afreq")
brnoaf <- read_tsv("/home/slecic/PhD/Transect/fst/plink2.95miss.brno.freq.afreq")
dim(viennaaf)
dim(brnoaf)
dim(afs)

# make a table of overall mean af, and population specific af
allele <- data.frame(conting = rep(afs$contig,2),
                     pos = rep(afs$pos,2),
                     pop = c(rep("vienna", nrow(afs)), rep("brno", nrow(afs))),
                     meanMaf = rep(afs$afMinor, 2),
                     af = c(viennaaf$ALT_FREQS, brnoaf$ALT_FREQS),
                     dist = c(rep(0, nrow(afs)), rep(68, nrow(afs))))
allele %>%
  ggplot( aes(x=af, fill=pop)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 10) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_classic() +
  labs(fill="")

allele %>%
  ggplot( aes(x=af, fill=pop)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_classic() +
  labs(fill="") +
  facet_wrap(~pop)

afdist <- allele %>%
  ggplot( aes(x=af, fill=pop)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 10) +
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) +
  theme_classic() +
  labs(fill="") +
  facet_wrap(~pop)
ggsave("/home/slecic/PhD/Transect/fst/All317SNPsfull95miss_afdist.png", plot=afdist, width = 10, height = 10, dpi = 300, device = "png", units = "in")


library(lme4)

outs <- subset(allele, pos %in% myout95$POS)
outs0109 <- subset(outs, meanMaf >= 0.1)
allele0109 <- subset(allele, meanMaf >= 0.1)
#laf <- lm(af ~ dist, data=allele[c(25, 438),])
#summary(laf)
#gaf <- glm(af ~ dist, family = poisson(link="log"), data = allele[c(25, 438),])
#summary(gaf)

ggall95missAlltogether <- allele0109 %>%
  ggplot(aes(dist,af, color=pop)) +
  geom_point(aes(fill=pop),size=3) +
  theme_classic() +
  geom_line(aes(group = pos),
            color="grey") +
  #arrow = arrow(type = "closed",
  #              length=unit(0.075, "inches"))) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,70)) +
  labs(x = "Distance", y = "Allele frequency") +
  geom_smooth(formula = y ~ x, method = "lm", aes(fill=pos), lty=1) #+
  #annotate("text")
  #facet_wrap(~pos, nrow = 10)
ggsave("/home/slecic/PhD/Transect/fst/All317SNPsfull95misstogether.png", plot=ggall95missAlltogether, width = 10, height = 10, dpi = 300, device = "png", units = "in")


ggall95missAll <- allele0109 %>%
  ggplot(aes(dist,af, color=pop)) +
  geom_point(aes(fill=pop),size=3) +
  theme_classic() +
  geom_line(aes(group = pos),
            color="grey") +
            #arrow = arrow(type = "closed",
            #              length=unit(0.075, "inches"))) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,70)) +
  labs(x = "Distance", y = "Allele frequency") +
  geom_smooth(formula = y ~ x, method = "lm", aes(fill=pos), lty=1) +
  #annotate("text")
 facet_wrap(~pos, nrow = 10)
ggsave("/home/slecic/PhD/Transect/fst/All317SNPsfull95miss.png", plot=ggall95missAll, width = 25, height = 25, dpi = 300, device = "png", units = "in")

ggout95missAll <- outs0109 %>%
  ggplot(aes(dist,af, color=pop)) +
  geom_point(aes(fill=pop),size=3) +
  theme_classic() +
  geom_line(aes(group = pos),
            color="grey") +
  #arrow = arrow(type = "closed",
  #              length=unit(0.075, "inches"))) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,70)) +
  labs(x = "Distance", y = "Allele frequency") +
  geom_smooth(formula = y ~ x, method = "lm", aes(fill=pos), lty=1) +
  #annotate("text")
  facet_wrap(~pos, nrow = 2)
ggsave("/home/slecic/PhD/Transect/fst/All8outSNPsfull95miss.png", plot=ggout95missAll, width = 10, height = 10, dpi = 300, device = "png", units = "in")


modelglm <- ggplot(data=allele, aes(x=dist, y=af)) +
  theme_classic() +
  #geom_line(col = "grey20", alpha = 0.5) +
  stat_smooth(method = "glm", method.args = list(family=poisson(link = "log")), lty=1, se=T, formula = y ~ x) +
  scale_y_continuous(0, 1) +
  labs(x = "Distance", y = "Allele frequency", color = "sh")
ggsave("/home/slecic/PhD/Transect/fst/logreg_afdist.png", plot=modelglm, width = 5, height = 5, dpi = 300, device = "png", units = "in")



#################################################
#################################################
# quantify allele frequencies differences
#################################################
#################################################

plotPfst<-function(x){
  require("ggplot2")
  dat<-read.table( x, header=FALSE )
  dat$V2<-1:length(dat$V2)
  dat<-dat[dat$V3 < 0.5,]
  theplot<-ggplot(dat, aes(x=V2, y=-log10(V3)))+geom_point()+theme_grey(15)+labs(x="SNP index", y="-log10(pFst)")
  pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
  ggsave(filename=pngName, width=20, height=5, units="in", theplot)
}
plotPfst("/home/slecic/PhD/Transect/fst/brno_vs_vienna_fst_10kb.windowed.weir.fst")


p.adj <- p.adjust(pvalues, method = c("BH"),n = length(pvalues))
pFst$p.adj <- p.adj
subset(pFst, p.adj > 0.05)


pFst <- read.table("/home/slecic/PhD/Transect/fst/brnovsvienna.pfst01af.50miss.mac3.q30.dp30.95miss.maf5.meanddp20.all.counts", header=F)
head(pFst)
table(is.na(pFst))
names(pFst) <- c("contig", "pos", "pval")
library(qvalue)
pvalues <- pFst$pval
p.adj <- p.adjust(pvalues, method = c("BH"), n = length(pvalues))

qobj <- qvalue(p = pvalues)
qobj <- qvalue(p = pvalues, fdr.level=0.05, lambda=0.)
table(qobj$significant)
#plot(qobj)
#hist(qobj)

pFst$qval <- qobj$qvalues
pFst$cutoff <- qobj$significant
pFst$log10qval <- -log10(pFst$qval)
head(pFst)
cutoff <- -log10(pFst[which(pFst$cutoff == "TRUE"),]$qval)

cutoff <- -log10(pFst[which(pFst$pval < 0.05),]$pval)

# plot with abline cutoff
plotPfst<-function(x){
  require("ggplot2")
  dat<-x
  dat<-dat[dat$pval < 0.1,]
  dat$pos<-1:length(dat$pos)
  theplot<-ggplot(dat, aes(x=pos, y=-log10(pval)))+geom_point()+theme_grey(15)+labs(x="SNP index", y="-log10(p)")
  theplot <- theplot + geom_hline(yintercept = min(cutoff), colour = "purple")
  #pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
  ggsave(filename="/home/slecic/PhD/Transect/fst/pfst01af_BrnoVsVienna.50miss.mac3.q30.dp30.95miss.maf5.meanddp20.png", width=20, height=5, units="in", theplot)
}
plotPfst(pFst)


mypfstouts <- subset(pFst, pval < 0.05)
pfstouts <- subset(allele, pos %in% mypfstouts$pos)


###### SUBSETTING SIGNIFICANT POSITIONS ####
### subset significant positions with values based on the COUNTS and save in a table
sig <- subset(pFst, cutoff == "TRUE")
dim(sig)
length(unique(sig$contig))
write.table(sig, "fdrPASScontigs.txt", row.names = F, col.names = T)

sigFemMale <- subset(pFst, cutoff == "TRUE")
dim(sigFemMale)
length(unique(sigFemMale$contig))
write.table(sigFemMale, "fdrPASScontigsFEMale.txt", row.names = F, col.names = T)

# check the number of overlapping contigs and Snps between Brnovs.Vienna and Femalesvs Males beaesd on the counts
dim(sig)
dim(sigFemMale)
dim(inner_join(sig, sigFemMale, by = c('contig','pos')))

### subset significant positions with values based on the GENOTYPES and save in a table
sig.geno <- subset(pFst, cutoff == "TRUE")
dim(sig.geno)
length(unique(sig.geno$contig))
write.table(sig.geno, "fdrPASScontigs_geno.txt", row.names = F, col.names = T)

sigFemMale.geno <- subset(pFst, cutoff == "TRUE")
dim(sigFemMale.geno)
length(unique(sigFemMale.geno$contig))
write.table(sigFemMale.geno, "fdrPASScontigsFEMale_geno.txt", row.names = F, col.names = T)

# check the number of overlapping contigs and Snps between Brnovs.Vienna and Femalesvs Males beaesd on the GENOTYPES
dim(sig.geno)
dim(sigFemMale.geno)
dim(inner_join(sig.geno, sigFemMale.geno, by = c('contig','pos')))


# with colored outliers above the cutoff
plotPfst<-function(x){
  require("ggplot2")
  dat<- as_tibble(x)
  print(head(dat))
  dat$pos<-1:length(dat$pos)
  #dat<-dat[dat$V3 < 0.9,]
  cutoff <- -log10(dat[which(dat$cutoff == "TRUE"),]$qval)
  dat <- dat %>% mutate(outlier = ifelse(dat$log10qval >= 1.6, "outlier", "background"))
  dat %>% group_by(outlier) %>% tally()
  print(head(dat))
  theplot<-ggplot(dat, aes(pos, -log10(qval), colour = outlier))+geom_point()+theme_grey(15)+labs(x="SNP index", y="-log10(p)")+
    scale_color_manual(values = c("grey40", "purple"))
  #theplot <- theplot + geom_hline(yintercept = min(cutoff), colour = "red")
  #pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
  ggsave(filename="fdrBrnoVsViennaColoredOutliers.png", width=20, height=5, units="in", theplot)
}
plotPfst(pFst)

##### to check overlapping contigs and Snps intersect the first two columns of Brno vs. Vienna and Females vs. Males
pFstBV <- read.table("/home/slecic/PhD/Transect/fst/pfst02af_BrnoVsVienna.counts", header=F)
dim(pFstBV)
pFst <- read.table("/home/slecic/PhD/Transect/fst/pfst02af_FemsvsMales.counts", header=F)
dim(pFst)
library(dplyr)
dim(inner_join(pFstBV, pFst, by = c('V1','V2')))


# plot the Fst values over 10kb windows 
plotfst<-function(x){
  require("ggplot2")
  dat<-read_tsv(x)
  print(head(dat))
  dat$BIN_START<-1:length(dat$BIN_START)
  print(head(dat))
  #print(dat)
  my_threshold <- quantile(dat$WEIGHTED_FST, 0.975, na.rm = T)
  print(my_threshold)
  dat <- dat %>% mutate(outlier = ifelse(dat$WEIGHTED_FST > my_threshold, "outlier", "background"))
  dat %>% group_by(outlier) %>% tally()
  print(head(dat))
  theplot<-ggplot(dat, aes(BIN_START, WEIGHTED_FST , colour = outlier))+geom_point()+theme_grey(15)+labs(x="SNP index", y="FST")+
    scale_color_manual(values = c("grey50", "red"))
  pngName<-paste(c(x, format(Sys.time(), "%a%b%d_%H_%M_%S.png")), collapse="_")
  ggsave(filename=pngName, width=20, height=5, units="in", theplot)
}
plotfst("/home/slecic/PhD/Transect/fst/brno_vs_vienna_fst_10kb.windowed.weir.fst")



########### plot contings in the manhatan plot   ###################

# data to plot in manhattahn plot
data <- data.frame(chr = pFst$contig,
                   win = pFst$pos,
                   pval = -log10(pFst$qval))
# vector of chromosome names
chromosomes <- as.character(unique(pFst$contig))
chrlen <- read.table("/home/slecic/PhD/Transect/fst/pilon_round4.nowol.bed")
chrlensub <- chrlen[as.character(chrlen$V1) %in% chromosomes, ]
#chrlen=chrlen[order(chrlen$V1),]
#chromosomes=unique(data$chr)
#chromosomes=as.character(chromosomes)
chromosomes_length=matrix(chrlensub$V3, nrow = 1, dimnames = list(NULL, chromosomes))
pval_lim <- max(na.omit(data$pval))
off <- 0

png("manhattanContig.png", width = 3000, height = 1200, units = "px", pointsize = 12)
for(i in seq_along(chromosomes)){
  chr <- chromosomes[i]
  if(i%%2){col <- "deepskyblue4"}else{col <- "darkorange1"}
  temp <- data[which(data$chr == chr), ]
  print(head(temp))
  if(chr == chromosomes[1]){
    plot(temp$win, temp$pval, xaxt = "n", xlim = c(0, sum(chromosomes_length[1, ])), ylim = c(0, pval_lim), col = col, pch=19,
         ylab = "-log10(p)", xlab = "chromosomes", cex.lab=1.5)
  }else if(chr == chromosomes[1]){
    temp2 = subset(temp, win %in% snpdata$win)
    plot(temp2$win, temp2$pval, xaxt = "n", xlim = c(0, sum(chromosomes_length[1, ])), ylim = c(0, pval_lim), col = "red", pch=19,
         ylab = "-log10(p)", xlab = "chromosomes", cex.lab=1.5)
  }else{
    points(temp$win+off, temp$pval, col = col, pch=20)
  }
  off <- off + chromosomes_length[1, chr]
  abline(h=min(cutoff), col = "red")
  chrl <- c(0, chrlensub$V3)
  v=c(0 + cumsum(chrl))
  axis(1, at=v[-length(v)] + diff(v) / 2, labels=chromosomes, las=2)
}
dev.off()


########################
## PLINK pca -> plink --vcf var.diplo-0.1.NoWol.snps.mis-0.1.vcf --double-id --allow-extra-chr --pca
#######################

library(tidyverse)
pca <- read_table2("/home/slecic/PhD/Transect/fst/plink.eigenvec", col_names = FALSE)
eigenval <- scan("/home/slecic/PhD/Transect/fst/plink.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$ind))
spp[grep("fe", pca$ind)] <- "female"
spp[grep("male", pca$ind)] <- "male"
spp <- c("female", "unknown", "female", "male", "female", "male", "male", "female", "male")
# location
loc <- rep(NA, length(pca$ind))
loc[grep("Vienna", pca$ind)] <- "vienna"
loc[grep("Brno", pca$ind)] <- "brno"
# combine - if you want to plot each in different colours
# host
host <- rep(NA, length(pca$ind))
host[grep("pru", pca$ind)] <- "prunus"
host[grep("lon", pca$ind)] <- "lonicera"
# combine - if you want to plot each in different colours
spp_loc <- paste0(spp, "_", loc, "_", host)
# remake data.frame
# first convert to percentage variance explained
pve <- data.frame(PC = 1:9, pve = eigenval/sum(eigenval)*100)
# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = spp, shape = loc)) + geom_point(size = 3)
#b + ggplot()
b <- b + scale_colour_manual(values = c("red", "blue", "green"))
#b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))



########
# INDIVUAL BASED INBREEDING with vcftools -> vcftools --vcf var.diplo-0.1.NoWol.snps.mis-0.1.vcf --het --out output_het
########
library(ggpubr)
inbr <- read.table("/home/slecic/PhD/Transect/fst/output_het.het", header=T)
head(inbr)
loc <- rep(NA, length(inbr$INDV))
loc[grep("Vienna", inbr$INDV)] <- "vienna"
loc[grep("Brno", inbr$INDV)] <- "brno"
# combine - if you want to plot each in different colours
p <- ggplot(inbr, aes(x=loc, y=F, colour = loc)) + 
  theme_grey(15) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_colour_manual(values = c("darkblue", "orange")) +
  stat_compare_means(method = "t.test")
p


##################
# per SNP inbredding with vcflib -> 
# popStats --type GL --target 1,3,4,6 --file var.diplo-0.1.NoWol.snps.q20.nomiss.vcf > Brno.popstats
##################
brnoStat <- read.table("/home/slecic/PhD/Transect/fst/BrnoVienna_full.popstats", header=T)
head(brnoStat)
colnames(brnoStat) <- c("conting", "pos", "targetAF", "expHET", "obsHET", "Nhets", "homoREF", "homoALT", "targetFIS")
brnoStat$pop <- rep("Brno", nrow(brnoStat))
viennaStat <- read.table("/home/slecic/PhD/Transect/fst/Vienna_full.popstats", header=T)
head(viennaStat)
colnames(viennaStat) <- c("conting", "pos", "targetAF", "expHET", "obsHET", "Nhets", "homoREF", "homoALT", "targetFIS")
viennaStat$pop <- rep("Vienna", nrow(viennaStat))
max(brnoStat$targetFIS)
subset(brnoStat, targetFIS > 0.1)
max(viennaStat$targetFIS)
subset(viennaStat, targetFIS > 0.1)
hist(brnoStat$targetFIS, breaks = 5)
hist(viennaStat$targetFIS, breaks = 5)
comb <- rbind(brnoStat, viennaStat)
comb2 <- subset(comb, targetFIS > 0.01 & targetFIS < 1)

library(ggpubr)
# plot the whole range
p <- ggplot(comb, aes(x=pop, y=targetFIS, colour = pop)) + 
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=15,face="bold"),
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=15),
        legend.key.size = unit(0.5, 'in'), #change legend key size
        legend.key.height = unit(0.5, 'in'), #change legend key height
        legend.key.width = unit(0.5, 'in')) +
  geom_violin() +
  #geom_boxplot()
  #geom_jitter(size = 0.001) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4")) +
  stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 0.8, size = 4) +
  labs(x="Population", y="Fis")
#p
ggsave(p, file="/home/slecic/PhD/Transect/fst/ClineFisWholeRange.png", width=10, height=10, units = "in", dpi = 300, limitsize = F)

res <- wilcox.test(brnoStat$targetFIS, viennaStat$targetFIS)
welch <- t.test(brnoStat$targetFIS, viennaStat$targetFIS)

# plot only Fis values higher than 0 and lower then 1
p <- ggplot(comb2, aes(x=pop, y=targetFIS, colour = pop)) + 
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6,face="bold"),
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=5),
        legend.key.size = unit(.5, 'in'), #change legend key size
        legend.key.height = unit(0.5, 'in'), #change legend key height
        legend.key.width = unit(0.5, 'in')) +
  geom_boxplot(outlier.size = NA) +
  geom_jitter(size = 0.001) +
  scale_colour_manual(values = c("darkblue", "orange")) +
  stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 0.53, size = 2)
#p
ggsave(p, file="/home/slecic/PhD/Transect/fst/Fis.png", width=50, height=50, units = "in", dpi = 300, limitsize = F)

# plot the whole range
p <- ggplot(comb2, aes(x=pop, y=targetFIS, colour = pop)) + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=15,face="bold"),
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=15),
        legend.key.size = unit(0.5, 'in'), #change legend key size
        legend.key.height = unit(0.5, 'in'), #change legend key height
        legend.key.width = unit(0.5, 'in')) +
  geom_boxplot() +
  #geom_jitter(size = 0.001) +
  scale_colour_manual(values = c("darkblue", "orange")) +
  stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 0.53, size = 5)
#p
ggsave(p, file="/home/slecic/PhD/Transect/fst/FisBox.png", width=10, height=10, units = "in", dpi = 300, limitsize = F)



brnoStat <- read.table("/home/slecic/PhD/Transect/fst/Wcer2neg.popstats", header=T)
head(brnoStat)
colnames(brnoStat) <- c("conting", "pos", "targetAF", "expHET", "obsHET", "Nhets", "homoREF", "homoALT", "targetFIS")
brnoStat$pop <- rep("wCer2neg", nrow(brnoStat))
viennaStat <- read.table("/home/slecic/PhD/Transect/fst/Wcer2pos.popstats", header=T)
head(viennaStat)
colnames(viennaStat) <- c("conting", "pos", "targetAF", "expHET", "obsHET", "Nhets", "homoREF", "homoALT", "targetFIS")
viennaStat$pop <- rep("wCer2pos", nrow(viennaStat))
max(brnoStat$targetFIS)
subset(brnoStat, targetFIS > 0.1)
max(viennaStat$targetFIS)
subset(viennaStat, targetFIS > 0.1)
hist(brnoStat$targetFIS, breaks = 5)
hist(viennaStat$targetFIS, breaks = 5)
comb <- rbind(brnoStat, viennaStat)
comb2 <- subset(comb, targetFIS > 0.01 & targetFIS < 1)

library(ggpubr)
# plot the whole range
p <- ggplot(comb, aes(x=pop, y=targetFIS, colour = pop)) + 
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=15,face="bold"),
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=15),
        legend.key.size = unit(0.5, 'in'), #change legend key size
        legend.key.height = unit(0.5, 'in'), #change legend key height
        legend.key.width = unit(0.5, 'in')) +
  geom_violin() +
  #geom_boxplot()
  #geom_jitter(size = 0.001) +
  scale_colour_manual(values = c("darkblue", "orange")) +
  stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 0.8, size = 4) +
  labs(x="Population", y="Fis")
#p
ggsave(p, file="/home/slecic/PhD/Transect/fst/AllFisWholeRange.png", width=10, height=10, units = "in", dpi = 300, limitsize = F)


############################
# Nucleotide diversity
###########################
brnoPi <- read.table("/home/slecic/PhD/Transect/fst/var.diplo-0.1.NoWol.snps.q20.nomiss.brno.1kbpi.diversity.windowed.pi", header=T)
head(brnoPi)
brnoPi$pop <- rep("Brno", nrow(brnoPi))
hist(brnoPi$PI)
viennaPi <- read.table("/home/slecic/PhD/Transect/fst/var.diplo-0.1.NoWol.snps.q20.nomiss.vienna.1kbpi.diversity.windowed.pi", header=T)
head(viennaPi)
viennaPi$pop <- rep("Vienna", nrow(viennaPi))
hist(viennaPi$PI)
median(brnoPi$PI)
median(viennaPi$PI)

combPi <- rbind(brnoPi, viennaPi)
options(scipen = 999)
# plot
p <- ggplot(combPi, aes(x=pop, y=PI, colour = pop)) + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=15,face="bold"),
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=15),
        legend.key.size = unit(0.5, 'in'), #change legend key size
        legend.key.height = unit(0.5, 'in'), #change legend key height
        legend.key.width = unit(0.5, 'in')) +
  geom_boxplot(outlier.size = NA) +
  #geom_jitter(size = 0.1) +
  scale_colour_manual(values = c("darkblue", "orange")) +
  stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = 0.03, size = 4)
#p
ggsave(p, file="/home/slecic/PhD/Transect/fst/Pi10kb.png", width=10, height=10, units = "in", dpi = 300, limitsize = F)

pi.contigBrno <- subset(brnoPi, CHROM %in% sig$contig)
pi.contigBrno$genome <- rep(1: nrow(pi.contigBrno))
plot(pi.contigBrno$genome,pi.contigBrno$PI,xlab="position",ylab="diversity", type = "l", 
     ylim = c(-0.01, 0.04))
pi.contigVienna <- subset(viennaPi, CHROM %in% sig$contig)
pi.contigVienna$genome <- rep(1: nrow(pi.contigVienna))
plot(pi.contigVienna$genome, pi.contigVienna$PI,xlab="position",ylab="diversity", type = "l",
     ylim = c(-0.01, 0.04))



#--------------------------- InbreedR -------------------------------#
library(inbreedR)
library(vcfR)
library(reshape2)

# write a function to caculate the population statistic g2 - Identity disequilibrium calculated across individuals and divided by the square of the mean heterozygosity
# Identity disequilibrium which is a covarinace between loci in heterozygous state. Therefore g2 is a measure of variance in identity by descent (IBD) amongst individuals
myg2 <- function(vcffile){
  vcf_file <- vcffile
  # read vcf
  vcf <- read.vcfR(vcf_file, verbose = FALSE )
  # extract genotypes
  gt <- extract.gt(vcf)
  # transpose and data.frame
  gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
  # NA handling
  gt[gt == "."] <- NA
  # split columns
  snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
  # convert
  brno_snp_genotypes <- inbreedR::convert_raw(snp_geno)
  # check data
  check_data(brno_snp_genotypes)
  # identity disequilibrium
  g2_brno_snps <- g2_snps(brno_snp_genotypes, nperm = 100, nboot = 10, 
                          CI = 0.95, parallel = FALSE, ncores = NULL)
  g2plot <- plot(g2_brno_snps, main = "SNPs",
       col = "blue", cex.axis=0.85)
  # heterozygosity-heterozygosity correlation coefficients (HHCs) 
  hhcplot <- HHC_brno_snps <- HHC(brno_snp_genotypes, reps = 100)
  plot(HHC_brno_snps, main = "SNPs",
       col = "blue")
}
myg2("/home/slecic/PhD/Transect/fst/var.full-genome-0.1.NoWol.snp.50miss.mac3.q30.dp30-300.95miss.maf05.meandp20.vcf.recode.wcer2pos.vcf")



vcf_file <- "/home/slecic/PhD/Transect/fst/var.full-genome-0.1.NoWol.snp.50miss.mac3.q30.dp30-300.95miss.maf05.meandp20.vcf.recode.Vienna.vcf"
  # read vcf
vcf <- read.vcfR(vcf_file, verbose = FALSE )
  # extract genotypes
gt <- extract.gt(vcf)
  # transpose and data.frame
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
  # NA handling
gt[gt == "."] <- NA
  # split columns
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
  # convert
brno_snp_genotypes <- inbreedR::convert_raw(snp_geno)
  # check data
check_data(brno_snp_genotypes)
  # identity disequilibrium
g2_brno_snps <- g2_snps(brno_snp_genotypes, nperm = 100, nboot = 10, 
                          CI = 0.95, parallel = FALSE, ncores = NULL)

png("/home/slecic/PhD/Transect/fst/All_g2_hhc.png", width = 1000, height = 1000, units="px", pointsize = 20)
par(mfrow=c(2,2))
plot(g2_brno_snps, main = "",
                 col = "#00BFC4", cex.axis=0.85, xlim=c(-0.02, 0.025))
  # heterozygosity-heterozygosity correlation coefficients (HHCs) 
hhcplot <- HHC_brno_snps <- HHC(brno_snp_genotypes, reps = 100)
plot(HHC_brno_snps, main = "",
      col ="#00BFC4", xlim=c(-1.2, 1))
#myg2("/home/slecic/PhD/Transect/fst/var.full-genome-0.1.NoWol.snp.50miss.mac3.q30.dp30-300.95miss.maf05.meandp20.vcf.recode.wcer2pos.vcf")



vcf_file <- "/home/slecic/PhD/Transect/fst/var.full-genome-0.1.NoWol.snp.50miss.mac3.q30.dp30-300.95miss.maf05.meandp20.vcf.recode.Brno.vcf"
# read vcf
vcf <- read.vcfR(vcf_file, verbose = FALSE )
# extract genotypes
gt <- extract.gt(vcf)
# transpose and data.frame
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
# NA handling
gt[gt == "."] <- NA
# split columns
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
# convert
brno_snp_genotypes <- inbreedR::convert_raw(snp_geno)
# check data
check_data(brno_snp_genotypes)
# identity disequilibrium
g2_brno_snps <- g2_snps(brno_snp_genotypes, nperm = 100, nboot = 10, 
                        CI = 0.95, parallel = FALSE, ncores = NULL)
plot(g2_brno_snps, main = "",
     col = "#F8766D", cex.axis=0.85, xlim=c(-0.02, 0.025))
# heterozygosity-heterozygosity correlation coefficients (HHCs) 
HHC_brno_snps <- HHC(brno_snp_genotypes, reps = 100)
plot(HHC_brno_snps, main = "",
     col = "#F8766D", xlim=c(-1.2 ,1))
dev.off()



############## ALL POPULATIONS ########################
vcf_file <- "/home/slecic/PhD/Transect/fst/var.full-genome-0.1.NoWol.snp.50miss.mac3.q30.dp30-300.95miss.maf05.meandp20.vcf.recode.wcer2pos.vcf"
# read vcf
vcf <- read.vcfR(vcf_file, verbose = FALSE )
# extract genotypes
gt <- extract.gt(vcf)
# transpose and data.frame
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
# NA handling
gt[gt == "."] <- NA
# split columns
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
# convert
brno_snp_genotypes <- inbreedR::convert_raw(snp_geno)
# check data
check_data(brno_snp_genotypes)
# identity disequilibrium
g2_brno_snps <- g2_snps(brno_snp_genotypes, nperm = 100, nboot = 10, 
                        CI = 0.95, parallel = FALSE, ncores = NULL)

png("/home/slecic/PhD/Transect/fst/All_g2_hhc.png", width = 1000, height = 1000, units="px", pointsize = 20)
par(mfrow=c(2,2))
plot(g2_brno_snps, main = "",
     col = "orange", cex.axis=0.85, xlim=c(-0.005, 0.025))
# heterozygosity-heterozygosity correlation coefficients (HHCs) 
hhcplot <- HHC_brno_snps <- HHC(brno_snp_genotypes, reps = 100)
plot(HHC_brno_snps, main = "",
     col ="orange", xlim=c(-1.2, 1))
#myg2("/home/slecic/PhD/Transect/fst/var.full-genome-0.1.NoWol.snp.50miss.mac3.q30.dp30-300.95miss.maf05.meandp20.vcf.recode.wcer2pos.vcf")



vcf_file <- "/home/slecic/PhD/Transect/fst/var.full-genome-0.1.NoWol.snp.50miss.mac3.q30.dp30-300.95miss.maf05.meandp20.vcf.recode.wcer2neg.vcf"
# read vcf
vcf <- read.vcfR(vcf_file, verbose = FALSE )
# extract genotypes
gt <- extract.gt(vcf)
# transpose and data.frame
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
# NA handling
gt[gt == "."] <- NA
# split columns
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
# convert
brno_snp_genotypes <- inbreedR::convert_raw(snp_geno)
# check data
check_data(brno_snp_genotypes)
# identity disequilibrium
g2_brno_snps <- g2_snps(brno_snp_genotypes, nperm = 100, nboot = 10, 
                        CI = 0.95, parallel = FALSE, ncores = NULL)
plot(g2_brno_snps, main = "",
     col = "blue", cex.axis=0.85, xlim=c(-0.005, 0.025))
# heterozygosity-heterozygosity correlation coefficients (HHCs) 
HHC_brno_snps <- HHC(brno_snp_genotypes, reps = 100)
plot(HHC_brno_snps, main = "",
     col = "blue", xlim=c(-1.2 ,1))
dev.off()


