# 31.1.2015 I.Zliobaite (Ilaret)
# model elev, precip and temp targets

library(lars)

#Net primary productivity (NPP) is computed from temperature and precipitation using the classic formula from Leith (1972), as cited in Liu et al (2012). 

#read
input_file_name <- 'working_data/data_hyp4.csv'

output_file_R2_lars_fit <- 'working_data/out_R2_lars_fit.csv'
output_file_R2_lars_cv <- 'working_data/out_R2_lars_cv.csv'
output_file_predictions_t1_cv <- 'working_data/out_predictions_t1_cv.csv'
output_file_predictions_t1_fit <- 'working_data/out_predictions_t1_fit.csv'
output_file_models_lars <- 'working_data/out_models_lars.csv'
output_file_features_t1 <- 'working_data/out_feature_selection_t1.csv'
#output_file_fet_coefficients <- 'working_data/out_feature_coefficients.csv'
output_file_fet_coefficients_t1 <- 'working_data/out_feature_coefficients_t1.csv'
output_model1 <- 'working_data/out_model1_fit.txt'


input_file_features_t1 <- output_file_features_t1
input_file_fet_coefficients_t1 <- output_file_fet_coefficients_t1

output_file_table2 <- 'results/table6.txt'
output_file_table3 <- 'results/table7.txt'
output_file_table4 <- 'results/table8.txt'

plot_name_fig2 <- 'results/figure2.pdf'
plot_name_NPP <- 'results/figureNPP.pdf'

R2_files <- c(output_file_R2_lars_fit,output_file_R2_lars_cv)

fet_targets <- c('PREC','PREC_MIN','PREC_MAX','PRECsp_MIN','PRECsp_MAX','NPP','NPP_MIN','NDVI','NDVI_MIN')

#fet_targets <- c('ELEV','TEMP','TEMP_MAX','TEMP_MIN', 'TEMPmax_MAX','TEMPmin_MIN','TEMPmax','TEMPmin')
# TEMP2 is the most extreme (min of min)
# TEMP is mean of min
# TEMP3 is min of mean

data_sites_all <- read.csv(input_file_name, header = TRUE)
p <- dim(data_sites_all)[2]

#0s removed due to linear dependence
#fet_inputs <- c('HYP','HOR','AL','OL','SF','OT','CM','HYP_1',"HYP_2","HYP_3",'HOR_1',"HOR_2","HOR_3",'MASS_log_mean')
fet_inputs <- c('HYP','HOR','AL','OL','SF','OT','CM')
fet_inputs_all <- c('HYP','HOR','AL','OL','SF','OT','CM')
fet_inputs_one <- c('MASS_log_mean','no_species_fact')

param_steps <- 6 #buvo13
param_round_digits <- 3

#for model reporting
param_model_select_lars <- 6+1 #the first model is 0

#for recording models for analysis of variable importance
param_target1 <- 'TEMP_MAX'
param_select_t1 <- 1 #which model
dg1 <- 0


##############################################################
#functions
compute_R2_matrix <- function(predictions_fun)
{  
  R2_all <- c()
  for (sk in 2:dim(predictions_fun)[2])
  {
    R2_all <- c(R2_all,compute_R2(predictions_fun[,1],predictions_fun[,sk]))  
  }
  return(R2_all)
}

compute_R2 <- function(true_fun,pred_fun)
{
  R2 <- 1 - sum((true_fun - pred_fun)^2)/sum((true_fun - mean(true_fun))^2)
  #R2 <-sqrt(sum((true_fun - pred_fun)^2)/length(true_fun))/15
  return(R2)
}

convert_predictions <- function(true_fun,pred_fun)
{
  comp_fun <- dim(pred_fun)[3]
  predictions_fun <- true_fun
  for (sk in 1:comp_fun)
  {
    predictions_fun <- cbind(predictions_fun,pred_fun[,,sk])
  }
  return(predictions_fun)
}

make_write_down_model <- function(model_now,mnames)
{
  k <- length(mnames)
  mnames <- process_feature_names(mnames)
  fml <- as.character(model_now[k+1])
  for(sk in 1:k)
  {
    if (model_now[sk]!=0)
    {
      if (model_now[sk]>0)
      {
        sgn = '+'
      }else
      {
        sgn = '-'
      }
      fml = paste(fml,sgn,as.character(abs(model_now[sk])),mnames[sk],sep = ' ') 
    }
  }
  return(fml)
}

convert_models <- function(coef_fun)
{
  comp_fun <- dim(coef_fun)[3]
  coefficients_fun <- c()
  for (sk in 1:comp_fun)
  {
    coefficients_fun <- rbind(coefficients_fun,coef_fun[,,sk])
  }
  return(coefficients_fun)
}

process_feature_names <- function(Feature)
{
  Feature <- gsub('MASS_log_mean','log(MASS)',Feature)
  Feature <- gsub('_', '=', Feature)
  for (sk in 1:length(Feature))
  {
    if (grepl('=',Feature[sk]))
    {
      Feature[sk] <- paste('prop(',Feature[sk],')',sep='')
    }
    else
    {
      Feature[sk] <- paste('mean(',Feature[sk],')',sep='')
    }
    if (grepl('med',Feature[sk]))
    {
      Feature[sk] <- gsub('med', '', Feature[sk])
      Feature[sk] <- gsub('mean', 'median', Feature[sk])
    }
  }
  return(Feature)
}

##############################################################

# manual cross validation
N <- dim(data_sites_all)[1]
all_R2_lars_cv <- c(0:param_steps)
all_R2_lars_fit <- c(0:param_steps)
all_R2_ols_all_cv <- c()
all_R2_ols_mass_cv <- c()
all_R2_ols_nspec_cv <- c()

print('target now:')
for (cik in 1:length(fet_targets))
{
  target_now <- fet_targets[cik]
  print(target_now)
  
  data_sites <- data_sites_all[,c(target_now,union(union(fet_inputs,fet_inputs_one),fet_inputs_all))] #for pls filter data
  fml <- as.formula(paste(target_now,'  ~.')) #regression variables, formula for pls
  
  predictions_lars_cv <- c()
  predictions_lars_fit <- c()
  predictions_ols_all_cv <- c()
  predictions_ols_mass_cv <- c()
  predictions_ols_nspec_cv <- c()
  
  model_collection_lars <- c()
  
  fet_selection_lars <- data_sites*0
  fet_coefs_lars <- data_sites*0
  
  #cross validation
  for (ind_test in 1:N)
  {
    #training data
    ind_train <- setdiff(c(1:N),ind_test)
    
    #fit models
    fit_lars <- lars(as.matrix(data_sites[ind_train,fet_inputs]),as.matrix(data_sites[ind_train,target_now]), normalize = TRUE,max.steps = param_steps)
    
    fit_ols_all <- lm(as.formula(paste(target_now,'~.')), data=data_sites[ind_train,c(fet_inputs_all,target_now)])
    fit_ols_mass <- lm(as.formula(paste(target_now,'~ MASS_log_mean')), data=data_sites[ind_train,c(fet_inputs_one,target_now)])
    fit_ols_nspec <- lm(as.formula(paste(target_now,'~ no_species_fact')), data=data_sites[ind_train,c(fet_inputs_one,target_now)])
    
    #extract models (coefficients) plain regression
    model_lars <- coef(fit_lars, mode = 'step')

    intercept_lars <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs]*0, type="fit")
  
    model_collection_lars <- rbind(model_collection_lars,c(model_lars[param_model_select_lars,],intercept_lars$fit[param_model_select_lars]))
    
    #selection order
    for (sk6 in 1:param_steps)
    {
      if (length(fit_lars$actions[[sk6]])>1)
      {
        fta <- fit_lars$actions[[sk6]][1] 
        print(fit_lars$actions[[sk6]])
      }else
      {
        fta <- fit_lars$actions[[sk6]]
      }
      ind_fet_now <- which(colnames(fet_selection_lars)==names(fta))
      ind_coef_now <- which(as.data.frame(colnames(model_lars))==names(fta))
      if (fit_lars$actions[[sk6]][1]>0)
      {
        fet_selection_lars[ind_test,ind_fet_now] <- param_steps + 1 - sk6
        fet_coefs_lars[ind_test,ind_fet_now] <- sign(model_lars[sk6+1,ind_coef_now])
      }
    }
    
    #predictions together with ground truth in the first column
    pp_lars <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs], type="fit")
    predictions_lars_cv <- rbind(predictions_lars_cv,cbind(data_sites[ind_test,target_now],t(pp_lars$fit)))
    predictions_ols_all_cv <- rbind(predictions_ols_all_cv,cbind(data_sites[ind_test,target_now],predict(fit_ols_all,data_sites[ind_test,c(fet_inputs_all,target_now)])))
    predictions_ols_mass_cv <- rbind(predictions_ols_mass_cv,cbind(data_sites[ind_test,target_now],predict(fit_ols_mass,data_sites[ind_test,c(fet_inputs_one,target_now)])))
    predictions_ols_nspec_cv <- rbind(predictions_ols_nspec_cv,cbind(data_sites[ind_test,target_now],predict(fit_ols_nspec,data_sites[ind_test,c(fet_inputs_one,target_now)])))
    
    #manual prediction (model lars)
    #sum(model[13,]*data_sites[ind_test,fet_inputs]) + intercept$fit[13])
    #manual predictions pls
    #sum(fit_pls$coefficients[,,5]*data_sites[4,2:32]) + intercepts_pls[5]
  }
  
  R2_lars_cv <- compute_R2_matrix(predictions_lars_cv)
  R2_ols_all_cv <- compute_R2_matrix(predictions_ols_all_cv)
  R2_ols_mass_cv <- compute_R2_matrix(predictions_ols_mass_cv)
  R2_ols_nspec_cv <- compute_R2_matrix(predictions_ols_nspec_cv)
  
  all_R2_lars_cv <- cbind(all_R2_lars_cv,round(R2_lars_cv,digits = param_round_digits))
  all_R2_ols_all_cv <- cbind(all_R2_ols_all_cv,round(R2_ols_all_cv,digits = param_round_digits))
  all_R2_ols_mass_cv <- cbind(all_R2_ols_mass_cv,round(R2_ols_mass_cv,digits = param_round_digits))
  all_R2_ols_nspec_cv <- cbind(all_R2_ols_nspec_cv,round(R2_ols_nspec_cv,digits = param_round_digits))
  
  #full fit and add to models as the last row
  fit_lars <- lars(as.matrix(data_sites[,fet_inputs]),as.matrix(data_sites[,target_now]), normalize = TRUE,max.steps = param_steps)
 
  model_lars <- coef(fit_lars, mode = 'step')
 
  intercept_lars <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs]*0, type="fit")
  
  model_collection_lars <- rbind(model_collection_lars,c(model_lars[param_model_select_lars,],intercept_lars$fit[param_model_select_lars]))
  
  #predictions all
  pp_lars <- predict.lars(fit_lars, data_sites[,fet_inputs], type="fit")
  predictions_lars_fit <- cbind(data_sites[,target_now],pp_lars$fit)
  
  R2_lars_fit <- compute_R2_matrix(predictions_lars_fit)
  
  all_R2_lars_fit <- cbind(all_R2_lars_fit,round(R2_lars_fit,digits = param_round_digits))
  
  if (target_now == param_target1){
    mod <- coef(fit_lars, mode = 'step')
    int <- predict.lars(fit_lars, data_sites[ind_test,fet_inputs]*0, type="fit")
    model_now <- c(mod[param_select_t1+1,],int$fit[param_select_t1+1]) #+1 because 1st model 0
    mnames <- names(model_now)[1:(length(model_now)-1)]
    model_now <- as.vector(model_now)
    md <- make_write_down_model(round(model_now,digits = dg1),mnames)
    write.table(md,file = output_model1,quote = FALSE,row.names = FALSE,sep=',')
    # predictions_lars_cv <- rbind(predictions_lars_cv,cbind(data_sites[ind_test,target_now],t(pp_lars$fit)))
    predictions_t1_cv <- cbind(as.data.frame(round(predictions_lars_cv,digits = dg1)),data_sites_all[,c('SITE','LAT','LON')],round(compute_R2_matrix(predictions_lars_cv)[param_select_t1+1],digits = 3))
    predictions_t1_fit <- cbind(as.data.frame(round(predictions_lars_fit,digits = dg1)),data_sites_all[,c('SITE','LAT','LON')],round(compute_R2_matrix(predictions_lars_fit)[param_select_t1+1],digits = 3))
    colnames(predictions_t1_cv) <- c('true',rep('pr',(dim( predictions_t1_cv)[2]-5)),'SITE','LAT','LON','R2')
    colnames(predictions_t1_fit) <- c('true',rep('pr',(dim( predictions_t1_fit)[2]-5)),'SITE','LAT','LON','R2')
    fet_selection_t1 <- fet_selection_lars
    fet_coefs_t1 <- fet_coefs_lars
  }
  
  #jacknife
  #jack = sqrt(apply((model_collection[1:13,] - matrix(1,13,1)%*%model_collection[14,])^2,2,sum)*(N-1)/N)
}

colnames(all_R2_lars_cv) <- c('step',fet_targets)
colnames(all_R2_lars_fit) <- c('step',fet_targets)
colnames(all_R2_ols_all_cv) <- fet_targets
colnames(all_R2_ols_mass_cv) <- fet_targets
colnames(all_R2_ols_nspec_cv) <- fet_targets

write.table(all_R2_lars_cv,file = output_file_R2_lars_cv,quote = FALSE,row.names = FALSE,sep='\t')
write.table(all_R2_lars_fit,file = output_file_R2_lars_fit,quote = FALSE,row.names = FALSE,sep='\t')

write.table(predictions_t1_fit,file = output_file_predictions_t1_fit,quote = FALSE,row.names = FALSE,sep='\t')
write.table(predictions_t1_cv,file = output_file_predictions_t1_cv,quote = FALSE,row.names = FALSE,sep='\t')



#write LARS models
#write.table(model_collection_lars,file = output_file_models_lars,quote = FALSE,row.names = FALSE,sep='\t')
#write feature selections
write.table(fet_selection_t1,file = output_file_features_t1,quote = FALSE,row.names = FALSE,sep='\t')
write.table(fet_coefs_t1,file = output_file_fet_coefficients_t1,quote = FALSE,row.names = FALSE,sep='\t')


############################################################################################
#plotting 

xlabels <- c('selection step','selection step')
ylabels <- c('R2',NA)
naive <- c(0,-0.174)
titles <- c('LARS fit (modeling)','LARS CV (testing)')
#parula
mycolors <-c('#0072bd','#d95319','#edb120','#7e2f8e','#77ac30','#4dbeee','#a2142f','#808080','#000000')
mylwd <- 2.2

# figure 2
pdf(plot_name_fig2,width = 6, height = 3.15)
#png(plot_name)
par(mfrow = c(1,2), mar= c(4, 4, 1, 1) + 0.2) #matrix of subplots, mar is for margins
for (plt in 1:2)
{
  #read file
  R2 <- read.csv(R2_files[plt], header = TRUE, sep = '\t')
  n <- dim(R2)[1]
  p <- dim(R2)[2]
  #setting up plotting area
  plot(NA,xlab = xlabels[plt],ylab = ylabels[plt],main = titles[plt],xlim=c(R2[1,'step'],R2[n,'step']),ylim=c(-1,1))
  #,ylim=c(min(R2[,c(2:p)]),max(R2[,c(2:p)])))
  #xaxs="i" reduces white space within axes
  for (sk in 1:length(fet_targets))
  {
    points(R2[,'step'],R2[,fet_targets[sk]],type="o",col = mycolors[sk],lwd=mylwd,pch=sk,cex=0.6)
  }  
  lines(c(R2[1,'step'],R2[n,'step']),c(naive[plt],naive[plt]),col = 1,lwd=1, lty=2)
  if (plt == 1 )
  {
    legend("bottom",fet_targets,bty='n',pch = c(1:sk), col = mycolors,lwd=mylwd,cex = 0.6)
  }
}
#legend("right",fet_targets,bty='n',pch = c(1:sk), col = mycolors,lwd=mylwd,cex = 0.7)
dev.off()

#################################################
#feature selection tables

make_feature_table <- function(input_file_features,input_file_fet_coefficients,out_file,out_file_latex)
{
  feature_selection <-read.csv(input_file_features, header = TRUE, sep = '\t')
  feature_selection <- feature_selection[,c(2:dim(feature_selection)[2])]
  mean_rating <- apply(feature_selection,2,mean)
  ind <- order(-mean_rating)
  Feature <- names(mean_rating[ind])
  
  Feature <- process_feature_names(Feature)
  
  Importance <-  round(mean_rating[ind],digits = 1)
  Frequency <- round(apply(feature_selection[,ind]>0,2,sum)/dim(feature_selection)[1],digits=2)
  feature_table <- cbind(Feature,Importance,Frequency)
  
  feature_coefficients <-read.csv(input_file_fet_coefficients, header = TRUE, sep = '\t')
  #assume the same order as features
  feature_coefficients <- feature_coefficients[,c(2:dim(feature_coefficients)[2])]
  feature_coefficients <- feature_coefficients[,ind]
  Sign <- c()
  Consistency <- c()
  for (sk in 1:dim(feature_coefficients)[2])
  {
    ind_pl <- which(feature_coefficients[,sk]==1)
    ind_mn <- which(feature_coefficients[,sk]==-1)
    lpl <- length(ind_pl)
    lmn <- length(ind_mn)
    if (lpl>lmn)
    {
      Sign <- rbind(Sign,'+')
      Consistency <- rbind(Consistency,lpl/(lpl+lmn))
    }
    else
    {
      if (lpl<lmn)
      {
        Sign <- rbind(Sign,'-')
        Consistency <- rbind(Consistency,lmn/(lpl+lmn))
      }
      else
      {
        Sign <- rbind(Sign,'.')
        Consistency <- rbind(Consistency,NA)
      }
    }  
  }
  Consistency <- round(Consistency,digits = 2)
  colnames(Sign) <- 'Sign'
  colnames(Consistency) <- 'Consistency'
  
  feature_table <- cbind(feature_table,Sign,Consistency)
  
  write.table(feature_table,file = out_file,quote = FALSE,row.names = FALSE,sep='\t')
}

make_feature_table(input_file_features_t1,input_file_fet_coefficients_t1,output_file_table2)