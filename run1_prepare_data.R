# 2015 11 04 I.Zliobaite
# this script takes occuence matrix, dental traits and climate variables as inputs, and prepares a summary dataset with sites in the rows and site characteristics in the columns

input_file_occurence <- 'input_data/data_occurence.csv'
input_file_traits <- 'input_data/data_traits.csv'
input_file_climate <- 'input_data/data_climate.csv'

output_file_data <- 'working_data/data_all.csv'
output_file_table1 <- 'results/table1.txt'
output_file_tableA1 <- 'results/tableA1.txt'
output_file_tableA2 <- 'results/tableA2.txt'
output_file_tableA3 <- 'results/tableA3.txt'
output_file_tableSF2 <- 'results/tableSF2.txt'

plot_name_cor1 <- 'results/figure1a.pdf'
plot_name_cor2 <- 'results/figure1b.pdf'

plot_name_wetlands <- 'results/figure_wetlands.pdf'
plot_name_wetlands2 <- 'results/figure_wetlands2.pdf'
plot_name_wetlands3 <- 'results/figure_wetlands3.pdf'

fet_extract <- c('HYP','HOD','AL','OL','SF','OT','CM') #features to extract
fet_mass <- 'MASS_KG' #feature for mass


#read data
data_occurence <- read.csv(input_file_occurence, header = TRUE)
data_traits <- read.csv(input_file_traits, header = TRUE,sep = '\t')
data_climate <- read.csv(input_file_climate, header = TRUE,sep = '\t')


#extract trait features
from_species_to_sites <- function(occurence_matrix,data_traits,fet_extract,fet_mass)
{
  no_species_fact <- apply(occurence_matrix,1,sum)
  #print(no_species_fact)
  
  extracted_mass <- c()
  extracted_features_means <- c()
  extracted_features_proportions <-c()
  
  # in the loop (not matrix multiplication) in order to double-check match of taxon
  for (sk in 1:dim(occurence_matrix)[1]){
    ind <- which(occurence_matrix[sk,]==1)
    ind_traits <- match(colnames(occurence_matrix[,ind]),data_traits[,'TAXON'])
    #print(length(ind_traits))
    #print(occurence_matrix[sk,])
    #print(data_traits[ind_traits,'TAXON'])
    # mass
    extracted_mass <- rbind(extracted_mass,mean(log(data_traits[ind_traits,fet_mass])))
    #dental means and proportions
    means_temp <- c()
    props_temp <- c()
    for (sk2 in 1:length(fet_extract)){ #for all dental traits
      trait_now <- fet_extract[sk2]
      means_temp <- c(means_temp,mean(data_traits[ind_traits,trait_now]))
      # extracting proportions
      unique_values <- unique(data_traits[,trait_now])
      fet_names_prop <- c()
      props_temp_for_one_feature <- c()
      for (sk3 in 1:length(unique_values)){
        un_value_now <- unique_values[sk3]
        fet_names_prop <- c(fet_names_prop,paste(trait_now,'_',as.character(un_value_now),sep=''))
        props_temp_for_one_feature <- c(props_temp_for_one_feature,sum(data_traits[ind_traits,trait_now]==un_value_now)/no_species_fact[sk])
      }
      names(props_temp_for_one_feature) <- fet_names_prop #assigning feature names
      props_temp <- c(props_temp,props_temp_for_one_feature)
    }
    extracted_features_means <- rbind(extracted_features_means,means_temp)
    extracted_features_proportions <- rbind(extracted_features_proportions,props_temp)
  }
  extracted_mass <- round(extracted_mass,digits = 3)
  extracted_features_means <- round(extracted_features_means,digits = 3)
  extracted_features_proportions <- round(extracted_features_proportions,digits = 3)
  #print(extracted_mass)
  species_count <- no_species_fact
  data_all <- cbind(species_count,extracted_mass,extracted_features_means)
  colnames(data_all)[2] <- 'MASS_log_mean'
  colnames(data_all) <- c(colnames(data_all)[1:2],fet_extract)
  data_all <- cbind(data_all,extracted_features_proportions)
  return(data_all)
}


rn <- data_occurence[,1]
occurence_matrix <- data_occurence[,c(2:dim(data_occurence)[2])]
occurence_matrix <- as.matrix(occurence_matrix)
rownames(occurence_matrix) <- rn
colnames(occurence_matrix) <- gsub('\\.',' ',colnames(occurence_matrix))
data_features <- from_species_to_sites(occurence_matrix,data_traits,fet_extract,fet_mass)


if (sum(rownames(data_features) != data_climate[,'SITE'])>0){
  print('sites do not match, check site row order and occurence matrix row order')
}else{
  data_all <- cbind(data_climate,data_features)
}


#compute NPP
NPPt <- 3000 / (1 + exp(1.315 - 0.119 * data_climate$TEMP))
NPPp <- 3000 * (1 - exp(-0.000664*data_climate$PREC))
NPP <- round(apply(cbind(NPPt,NPPp),1,min))

NPPt_low <- 3000 / (1 + exp(1.315 - 0.119 * data_climate$TEMP_low))
NPPt_MIN <- 3000 / (1 + exp(1.315 - 0.119 * data_climate$TEMP_MIN))
NPPp_low <- 3000 * (1 - exp(-0.000664*data_climate$PREC_low))
NPPp_MIN <- 3000 * (1 - exp(-0.000664*data_climate$PREC_MIN))

NPP_low_MIN <- round(apply(cbind(NPPt_low,NPPp_MIN),1,min))
NPP_low_low <- round(apply(cbind(NPPt_low,NPPp_low),1,min))
NPP_MIN_MIN <- round(apply(cbind(NPPt_MIN,NPPp_MIN),1,min))
NPP_MIN_low <- round(apply(cbind(NPPt_MIN,NPPp_low),1,min))

data_all <- cbind(data_all,NPP,NPP_low_MIN,NPP_low_low,NPP_MIN_MIN,NPP_MIN_low)


#write file
write.table(data_all, file = output_file_data, quote = FALSE, row.names = FALSE,col.names = TRUE,sep=',')


#write table A2
write.table(t(data_occurence),file = output_file_tableA2,quote = FALSE,row.names = TRUE,col.names = FALSE,sep='\t')

#write table A3

fet_targets <- c('ELEV','TEMP','TEMP_MIN','TEMP_low','TEMP_low_MIN','TEMP_MAX','TEMP_high','TEMP_high_MAX','PREC','PREC_MIN','PREC_low','PREC_MAX','PREC_high','NPP','NPP_MIN_MIN','NPP_low_low','NPP_low_MIN','NPP_MIN_low','NDVI','NDVI_MIN9y','NDVI_low1y','NDVI_low1y_MIN','NDVI_low','NDVI_low_MIN','species_count')

data_A3 <- data_all[,c('SITE',fet_targets)]

write.table(t(data_A3),file = output_file_tableA3,quote = FALSE,row.names = TRUE,col.names = FALSE,sep='\t')



#make Table 1
pred_tab <- data_all[,c('SITE','SITE_name','AREA','ELEV','TEMP','PREC','NPP','NDVI','species_count')]
pred_tab[,'AREA'] <- round(pred_tab[,'AREA'])
pred_tab[,'ELEV'] <- round(pred_tab[,'ELEV'])
pred_tab[,'TEMP'] <- round(pred_tab[,'TEMP'],digits = 1)
pred_tab[,'PREC'] <- round(pred_tab[,'PREC'])
pred_tab[,'NDVI'] <- round(pred_tab[,'NDVI'],digits = 2)
colnames(pred_tab) <- c('Site','Site name','Area km2','Elev., m','Av. temp., C','Ann. precip., mm','NPP, gC/m2year','NDVI','No. species')
write.table(pred_tab,file = output_file_table1,quote = FALSE,row.names = FALSE,sep='\t')

#make Table A1
pred_tab <- data_all[,c('SITE','SITE_name','TEMP_MIN','PREC_MIN','NPP_MIN_MIN','NDVI_MIN9y')]
pred_tab[,'TEMP_MIN'] <- round(pred_tab[,'TEMP_MIN'],digits = 1)
pred_tab[,'PREC_MIN'] <- round(pred_tab[,'PREC_MIN'])
pred_tab[,'NDVI_MIN9y'] <- round(pred_tab[,'NDVI_MIN9y'],digits = 2)
colnames(pred_tab) <- c('Site','Site name','Min. temp., C','Min. precip., mm','Min. NPP, gC/m2year','Min. NDVI')
write.table(pred_tab,file = output_file_tableA1,quote = FALSE,row.names = FALSE,sep='\t')


# corr plot 1, 63 samples
library(corrplot)
fet_cor <- c('HYP','HOD','AL','OL','SF','OT','CM')
cor_traits <- round(cor(data_traits[,fet_cor]),digits = 2)
cor_sites <- round(cor(data_all[,fet_cor]),digits=2)
colnames(cor_sites) <- paste('mean(',colnames(cor_sites),')',sep='')
rownames(cor_sites) <- paste('mean(',rownames(cor_sites),')',sep='')

pdf(plot_name_cor1,width = 3.4,height = 3.4)
#png(plot_name_cor,width = 100, res = 100)
corrplot(as.matrix(cor_traits),method="square",addCoef.col="black", addCoefasPercent = TRUE,cl.pos = "n")
dev.off()

pdf(plot_name_cor2,width = 4,height = 4)
#png(plot_name_cor,width = 100, res = 100)
corrplot(as.matrix(cor_sites),method="square",addCoef.col="black", addCoefasPercent = TRUE,cl.pos = "n")
dev.off()

#make Table SF=2
ind_spec_SF2 <- which(data_traits[,'SF']==1)
species_SF2 <- as.vector(data_traits[ind_spec_SF2,'TAXON'])
species_SF2 <- gsub(" ",".",species_SF2)
sp_parks <- as.vector(data_occurence[,1])
sp_sum <- apply(data_occurence[,c(2:64)],1,sum)
species_tab <- data_occurence[,species_SF2]
sp_sum_SF_2 <- apply(species_tab,1,sum)
sp_div <- round(sp_sum_SF_2/sp_sum,digits = 3)
species_tab <- t(species_tab)
species_tab <- rbind(species_tab,sp_sum,sp_div)
colnames(species_tab) <- sp_parks
species_tab <- cbind(rownames(species_tab),species_tab)
species_tab[5,1] <- 'No. of species'
species_tab[6,1] <- 'prop(SF=2)'
species_tab[c(1:4),c(2:14)] <- gsub("0"," ",species_tab[c(1:4),c(2:14)])
species_tab[c(1:4),c(2:14)] <- gsub("1","v",species_tab[c(1:4),c(2:14)])
#species_tab[c(1:4),1] <- gsub("."," ",species_tab[c(1:4),1])
write.table(species_tab,file = output_file_tableSF2,quote = FALSE,row.names = FALSE,sep='\t')

#occurence sum

occ_sum <- apply(data_occurence[,2:dim(data_occurence)[2]],2,sum)
tax_SF1_hyp <- c('Hippotragus equinus','Hippotragus niger','Kobus ellipsiprymnus','Oryx beisa','Redunca fulvorufula','Redunca redunca')
tax_SF1_mes <- c('Cephalophus adersi','Cephalophus harveyi','Cephalophus nigrifrons','Cephalophus silvicultor','Cephalophus weynsi','Hippopotamus amphibius','Hylochoerus meinertzhageni','Potamochoerus larvatus')
#average number of sites is about the same
mn_occ_hyp <- c()
mn_hyp_hyp <- c()
cnt_hyp <- c()
for (sk in 1:length(tax_SF1_hyp)){
  tax_now <- gsub(' ','.',tax_SF1_hyp[sk])
  ind <- which(names(occ_sum)==tax_now)
  mn_occ_hyp <- c(mn_occ_hyp,occ_sum[ind])
  ind <- which(data_occurence[,tax_now]==1)
  mn_hyp_hyp <- c(mn_hyp_hyp,mean(data_all[ind,'HYP']))
  cnt_hyp <- c(cnt_hyp,mean(data_all[ind,'species_count']))
}
mn_occ_mes <- c()
mn_hyp_mes <- c()
cnt_mes <- c()
for (sk in 1:length(tax_SF1_mes)){
  tax_now <- gsub(' ','.',tax_SF1_mes[sk])
  ind <- which(names(occ_sum)==tax_now)
  mn_occ_mes <- c(mn_occ_mes,occ_sum[ind])
  ind <- which(data_occurence[,tax_now]==1)
  mn_hyp_mes <- c(mn_hyp_mes,mean(data_all[ind,'HYP']))
  cnt_mes <- c(cnt_mes,mean(data_all[ind,'species_count']))
}

#hippo occurs or not
ind_occ_hip <- which(data_occurence[,'Hippopotamus.amphibius']==1)
ind_not_hip <- which(data_occurence[,'Hippopotamus.amphibius']==0)
ind_occ_red <- which(data_occurence[,'Redunca.redunca']==1)
ind_not_red <- which(data_occurence[,'Redunca.redunca']==0)
fet_interest <- c('HYP','HOD','AL','OL','SF','OT','CM','PREC','TEMP','NPP','NDVI','ELEV','species_count','PREC_low','NPP_low_low','NDVI_low_MIN','TEMP_high')
mn <- apply(data_all[,fet_interest],2,mean)
mn_occ_hip <- apply(data_all[ind_occ_hip,fet_interest],2,mean)
mn_not_hip <- apply(data_all[ind_not_hip,fet_interest],2,mean)
mn_occ_red <- apply(data_all[ind_occ_red,fet_interest],2,mean)
mn_not_red <- apply(data_all[ind_not_red,fet_interest],2,mean)
mn_all <- rbind(mn,mn_occ_hip,mn_not_hip,mn_occ_red,mn_not_red)
print(mn_all)
data_wet <- cbind(abs(data_occurence[,'Hippopotamus.amphibius']-1),abs(data_occurence[,'Redunca.redunca']-1),data_all[,fet_interest])
colnames(data_wet)[1:2] <- c('Hippo','Redunca')
col1 <- "#e79f00"
col2 <- "#9ad0f3"
col3 <- "#009E73"
pdf(plot_name_wetlands,width = 4.6, height = 11)
par(mfrow = c(6, 3),mai = c(0.4, 0.55, 0.3, 0.01))
#HYP
lmts <- range(data_wet$HYP)
boxplot(HYP~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'HYP',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(HYP~Hippo,data=data_wet,mean)
points(1:2, mn$HYP, col = "black",pch='-')
boxplot(HYP~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(HYP~Hippo,data=data_wet,mean)
points(1:2, mn$HYP, col = "black",pch='-')
boxplot(data_wet$HYP,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$HYP), col = "black",pch='-')
#HOD
lmts <- range(data_wet$HOD)
boxplot(HOD~Hippo,data=data_wet, main="Hippo", ylim=lmts, ylab = 'HOD',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(HOD~Hippo,data=data_wet,mean)
points(1:2, mn$HOD, col = "black",pch='-')
boxplot(HOD~Redunca,data=data_wet, main="Redunca", ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(HOD~Hippo,data=data_wet,mean)
points(1:2, mn$HOD, col = "black",pch='-')
boxplot(data_wet$HOD,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$HOD), col = "black",pch='-')
#AL
lmts <- range(data_wet$AL)
boxplot(AL~Hippo,data=data_wet, main="Hippo", ylim=lmts, ylab = 'AL',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(AL~Hippo,data=data_wet,mean)
points(1:2, mn$AL, col = "black",pch='-')
boxplot(AL~Redunca,data=data_wet, main="Redunca", ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(AL~Hippo,data=data_wet,mean)
points(1:2, mn$AL, col = "black",pch='-')
boxplot(data_wet$AL,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$AL), col = "black",pch='-')
#OL
lmts <- range(data_wet$OL)
boxplot(OL~Hippo,data=data_wet, main="Hippo", ylim=lmts, ylab = 'OL',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(OL~Hippo,data=data_wet,mean)
points(1:2, mn$OL, col = "black",pch='-')
boxplot(OL~Redunca,data=data_wet, main="Redunca", ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(OL~Hippo,data=data_wet,mean)
points(1:2, mn$OL, col = "black",pch='-')
boxplot(data_wet$OL,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$OL), col = "black",pch='-')
#SF
lmts <- range(data_wet$SF)
boxplot(SF~Hippo,data=data_wet, main="Hippo", ylim=lmts, ylab = 'SF',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(SF~Hippo,data=data_wet,mean)
points(1:2, mn$SF, col = "black",pch='-')
boxplot(SF~Redunca,data=data_wet, main="Redunca", ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(SF~Hippo,data=data_wet,mean)
points(1:2, mn$SF, col = "black",pch='-')
boxplot(data_wet$SF,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$SF), col = "black",pch='-')
#CM
lmts <- range(data_wet$CM)
boxplot(CM~Hippo,data=data_wet, main="Hippo", ylim=lmts, ylab = 'CM',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(CM~Hippo,data=data_wet,mean)
points(1:2, mn$CM, col = "black",pch='-')
boxplot(CM~Redunca,data=data_wet, main="Redunca", ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(CM~Hippo,data=data_wet,mean)
points(1:2, mn$CM, col = "black",pch='-')
boxplot(data_wet$CM,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$CM), col = "black",pch='-')
dev.off()

pdf(plot_name_wetlands2,width = 4.6, height = 11)
par(mfrow = c(6, 3),mai = c(0.4, 0.55, 0.3, 0.01))
#PREC
lmts <- range(data_wet$PREC)
boxplot(PREC~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'PREC',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(PREC~Hippo,data=data_wet,mean)
points(1:2, mn$PREC, col = "black",pch='-')
boxplot(PREC~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(PREC~Hippo,data=data_wet,mean)
points(1:2, mn$PREC, col = "black",pch='-')
boxplot(data_wet$PREC,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$PREC), col = "black",pch='-')
#NPP
lmts <- range(data_wet$NPP)
boxplot(NPP~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'NPP',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(NPP~Hippo,data=data_wet,mean)
points(1:2, mn$NPP, col = "black",pch='-')
boxplot(NPP~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(NPP~Hippo,data=data_wet,mean)
points(1:2, mn$NPP, col = "black",pch='-')
boxplot(data_wet$NPP,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$NPP), col = "black",pch='-')
#NDVI
lmts <- range(data_wet$NDVI)
boxplot(NDVI~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'NDVI',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(NDVI~Hippo,data=data_wet,mean)
points(1:2, mn$NDVI, col = "black",pch='-')
boxplot(NDVI~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(NDVI~Hippo,data=data_wet,mean)
points(1:2, mn$NDVI, col = "black",pch='-')
boxplot(data_wet$NDVI,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$NDVI), col = "black",pch='-')
#TEMP
lmts <- range(data_wet$TEMP)
boxplot(TEMP~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'TEMP',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(TEMP~Hippo,data=data_wet,mean)
points(1:2, mn$TEMP, col = "black",pch='-')
boxplot(TEMP~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(TEMP~Hippo,data=data_wet,mean)
points(1:2, mn$TEMP, col = "black",pch='-')
boxplot(data_wet$TEMP,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$TEMP), col = "black",pch='-')
#ELEV
lmts <- range(data_wet$ELEV)
boxplot(ELEV~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'ELEV',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(ELEV~Hippo,data=data_wet,mean)
points(1:2, mn$ELEV, col = "black",pch='-')
boxplot(ELEV~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(ELEV~Hippo,data=data_wet,mean)
points(1:2, mn$ELEV, col = "black",pch='-')
boxplot(data_wet$ELEV,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$ELEV), col = "black",pch='-')
#species_count
lmts <- range(data_wet$species_count)
boxplot(species_count~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'species_count',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(species_count~Hippo,data=data_wet,mean)
points(1:2, mn$species_count, col = "black",pch='-')
boxplot(species_count~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(species_count~Hippo,data=data_wet,mean)
points(1:2, mn$species_count, col = "black",pch='-')
boxplot(data_wet$species_count,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$species_count), col = "black",pch='-')
dev.off()

pdf(plot_name_wetlands3,width = 4.6, height = 7.5)
par(mfrow = c(4, 3),mai = c(0.4, 0.55, 0.3, 0.01))
#PREC_low
lmts <- range(data_wet$PREC_low)
boxplot(PREC_low~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'PREC_low (spatial)',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(PREC_low~Hippo,data=data_wet,mean)
points(1:2, mn$PREC_low, col = "black",pch='-')
boxplot(PREC_low~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(PREC_low~Hippo,data=data_wet,mean)
points(1:2, mn$PREC_low, col = "black",pch='-')
boxplot(data_wet$PREC_low,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$PREC_low), col = "black",pch='-')
#NPP_low_low
lmts <- range(data_wet$NPP_low_low)
boxplot(NPP_low_low~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'NPP_low_low',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(NPP_low_low~Hippo,data=data_wet,mean)
points(1:2, mn$NPP_low_low, col = "black",pch='-')
boxplot(NPP_low_low~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(NPP_low_low~Hippo,data=data_wet,mean)
points(1:2, mn$NPP_low_low, col = "black",pch='-')
boxplot(data_wet$NPP_low_low,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$NPP_low_low), col = "black",pch='-')
#NDVI_low_MIN
lmts <- range(data_wet$NDVI_low_MIN)
boxplot(NDVI_low_MIN~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'NDVI_low_MIN',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(NDVI_low_MIN~Hippo,data=data_wet,mean)
points(1:2, mn$NDVI_low_MIN, col = "black",pch='-')
boxplot(NDVI_low_MIN~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(NDVI_low_MIN~Hippo,data=data_wet,mean)
points(1:2, mn$NDVI_low_MIN, col = "black",pch='-')
boxplot(data_wet$NDVI_low_MIN,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$NDVI_low_MIN), col = "black",pch='-')
#TEMP_high
lmts <- range(data_wet$TEMP_high)
boxplot(TEMP_high~Hippo,data=data_wet, main="Hippo",ylim=lmts, ylab = 'TEMP_high',names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(TEMP_high~Hippo,data=data_wet,mean)
points(1:2, mn$TEMP_high, col = "black",pch='-')
boxplot(TEMP_high~Redunca,data=data_wet, main="Redunca",ylim=lmts,names = c('present','absent'),col = c(col1,col2))
mn <- aggregate(TEMP_high~Hippo,data=data_wet,mean)
points(1:2, mn$TEMP_high, col = "black",pch='-')
boxplot(data_wet$TEMP_high,data=data_wet, main="All sites",ylim=lmts,col = col3)
points(1, mean(data_wet$TEMP_high), col = "black",pch='-')
dev.off()