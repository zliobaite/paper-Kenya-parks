# 22.1.2015 I.Zliobaite
# plot target variables against traits

library(corrplot)

#read
input_file_name <- 'working_data/data_all.csv'
input_file2 <- 'input_data/data_occurence.csv'

plot_name_fig1 <- 'results/figure2.pdf'
plot_name_fig1_png <- 'results/figure2.png'
plot_name_figA2 <- 'results/figureA2.pdf'

fet_targets <- c('ELEV','TEMP','TEMP_MIN','TEMP_low','TEMP_low_MIN','TEMP_MAX','TEMP_high','TEMP_high_MAX','PREC','PREC_MIN','PREC_low','PREC_MAX','PREC_high','NPP','NPP_MIN_MIN','NPP_low_low','NPP_low_MIN','NPP_MIN_low','NDVI','NDVI_MIN9y','NDVI_low1y','NDVI_low1y_MIN','NDVI_low','NDVI_low_MIN','species_count')

fet_inputs <- c('HYP','HOD','AL','OL','SF','OT','CM','HYP_1', 'HYP_2','HYP_3','HOD_1','HOD_2','HOD_3','MASS_log_mean','species_count','ELEV')

#read data
data_all <- read.csv(input_file_name, header = TRUE,sep = ',')
data_occ <- as.matrix(read.csv(input_file2, header = TRUE,sep = ','))
if (sum(data_occ[,1]!=data_all[,'SITE'])>0) {print('site mismatch')}
data_occ <- data_occ[,2:dim(data_occ)[2]]
sd <- as.vector(apply(data_occ,2,sd)) #standard deviations
ind_non_const <- which(sd>0)
data_occ <- data_occ[,ind_non_const]

#correlations traits
cor_all <- c()
for (sk2 in 1:length(fet_inputs))
{
  cor_row <- c()
  for (sk in 1:length(fet_targets))
  {
    cor_row <- cbind(cor_row, cor(data_all[,fet_targets[sk]],data_all[,fet_inputs[sk2]]))
  }
  cor_all <- rbind(cor_all,cor_row)
}
rownames(cor_all) <- fet_inputs
colnames(cor_all) <- fet_targets
cor_all <- t(cor_all)

process_feature_names <- function(Feature)
{
  Feature <- gsub('MASS_log_mean','log(MASS)',Feature)
  Feature <- gsub('KG_', 'KG', Feature)
  Feature <- gsub('_KG', '', Feature)
  Feature <- gsub('_', '=', Feature)
  for (sk in 1:length(Feature))
  {
    if (grepl('=',Feature[sk]))
    {
      #Feature[sk] <- paste('prop(',Feature[sk],')',sep='')
      Feature[sk] <- paste('p(',Feature[sk],')',sep='')
    }
    else
    {
      #Feature[sk] <- paste('mean(',Feature[sk],')',sep='')
      Feature[sk] <- paste('m(',Feature[sk],')',sep='')
    }
    if (grepl('med',Feature[sk]))
    {
      Feature[sk] <- gsub('med', '', Feature[sk])
      Feature[sk] <- gsub('mean', 'median', Feature[sk])
    }
  }
  return(Feature)
}

colnames(cor_all)[1:(length(fet_inputs)-2)] <- process_feature_names(colnames(cor_all)[1:(length(fet_inputs)-2)])

cor_all <- t(cor_all)
pdf(plot_name_fig1,width = 8, height = 9)
corrplot(as.matrix(cor_all),method="square",addCoef.col="black", addCoefasPercent = TRUE,cl.pos = "n",col=colorRampPalette(c("darkorange","white","steelblue"))(200))
dev.off()

png(plot_name_fig1_png,width = 8, height = 9, units = 'in', res = 1200)
corrplot(as.matrix(cor_all),method="square",addCoef.col="black", addCoefasPercent = TRUE,cl.pos = "n",col=colorRampPalette(c("darkorange","white","steelblue"))(200))
dev.off()


#correlations occurence

fet_species <- colnames(data_occ)
cor_all <- c()
for (sk2 in 1:length(fet_species)){
  cor_row <- c()
  for (sk in 1:length(fet_targets)){
    cor_row <- cbind(cor_row, cor(data_all[,fet_targets[sk]],as.numeric(data_occ[,fet_species[sk2]])))
  }
  cor_all <- rbind(cor_all,cor_row)
}
rownames(cor_all) <- gsub('\\.',' ',fet_species)
colnames(cor_all) <- fet_targets

pdf(plot_name_figA2,width = 7,height = 25)
#png(plot_name_cor,width = 100, res = 100)
corrplot(as.matrix(cor_all),method="square",addCoef.col="black", addCoefasPercent = TRUE,cl.pos = "n",col=colorRampPalette(c("darkorange","white","steelblue"))(200))
dev.off()
