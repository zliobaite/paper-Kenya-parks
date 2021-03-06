# 2015 11 04 I.Zliobaite
# this script takes occuence matrix, dental traits and climate variables as inputs, and prepares a summary dataset with sites in the rows and site characteristics in the columns

input_file_occurence <- 'input_data/data_occurence.csv'
input_file_traits <- 'input_data/data_traits.csv'
input_file_climate <- 'input_data/data_climate.csv'

output_file_data <- 'working_data/data_hyp4.csv'

plot_name_cor1 <- 'results/figure1a.pdf'
plot_name_cor2 <- 'results/figure1b.pdf'

fet_extract <- c('HYP','HOD','AL','OL','SF','OT','CM') #features to extract
fet_mass <- 'MASS_KG' #feature for mass


#read data
data_occurence <- read.csv(input_file_occurence, header = TRUE)
data_traits <- read.csv(input_file_traits, header = TRUE,sep = '\t')
data_climate <- read.csv(input_file_climate, header = TRUE,sep = '\t')

#impute equids HYP=4
ind <- which(data_traits[,'FAMILY'] == 'Equidae')
data_traits[ind,'HYP'] <- 4
#Alcelaphus buselaphus
#Beatragus hunteri
#Connochaetes taurinus
ind <- which(data_traits[,'TAXON'] == 'Alcelaphus buselaphus')
data_traits[ind,'HYP'] <- 4
ind <- which(data_traits[,'TAXON'] == 'Beatragus hunteri')
data_traits[ind,'HYP'] <- 4
ind <- which(data_traits[,'TAXON'] == 'Connochaetes taurinus')
data_traits[ind,'HYP'] <- 4
ind <- which(data_traits[,'TAXON'] == 'Loxodonta africana')
data_traits[ind,'HYP'] <- 4

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
