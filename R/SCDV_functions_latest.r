#' QC check
#'
#' This function allows you to check the library size and detaction rate.
#' @param data_count Input count matrix.
#' @keywords qc
#' @export
#' @examples
#' \dontrun{
#'    qc_check(input_data)
#' }
qc_check <- function(data_count){
  library_size <- apply(data_count,2,sum)
  detected_gene_num <- apply(data_count,2,function(x)length(which(x>0)))
  return(list(library_size=library_size,detected_gene_num=detected_gene_num))
}

#' Get the nearest neighboring cells
#'
#' This function allows you to get the nunmber of nearest neighbor used for estimating dropout.
#' @param data_count Input count matrix.
#' @param gene_len Kilo-base-pair length of gene
#' @param max_num Maximum number of nearest neighboring cells
#' @param scale_factor Global scale factor for normalizing the nearest neighboring cells
#' @keywords nearest neighbor
#' @export
get_expected_cell <- function(data_count,gene_len,max_num=20,scale_factor=1e6){

  if(max_num > ncol(data_count)){
    max_num <- ncol(data_count)
  }

  gene_len_scale <- ceiling(gene_len/1000)
  data_FPKM <- t(t(data_count/gene_len_scale)*scale_factor/apply(data_count,2,sum))
  data_dist <- as.matrix(dist(t(data_FPKM)))
  data_expect <- matrix(data=NA,ncol=ncol(data_count),nrow=nrow(data_count))

  if(max_num==1){
    for(j in 1:ncol(data_FPKM)){
      nearest_cell_idx <- sort(data_dist[j,],decreasing=FALSE,index.return=TRUE)$ix[2]
      data_expect[,j] <- data_FPKM[,nearest_cell_idx]
    }
  }
  else{
    sim_stat <- list()
    sim_stat[[1]] <- rep(NA,ncol(data_FPKM))
    for(j in 1:ncol(data_FPKM)){
      nearest_cell_idx <- sort(data_dist[j,],decreasing=FALSE,index.return=TRUE)$ix[2]
      data_neighbor <- data_FPKM[,nearest_cell_idx]
      sim_stat[[1]][j] <- cor(data_FPKM[,j],data_neighbor)
    }

    for(i in 2:max_num){
      sim_stat[[i]] <- rep(NA,ncol(data_FPKM))
      for(j in 1:ncol(data_FPKM)){
        nearest_cell_idx <- sort(data_dist[j,],decreasing=FALSE,index.return=TRUE)$ix[2:(i+1)]
        data_neighbor <- apply(data_FPKM[,nearest_cell_idx],1,mean)
        sim_stat[[i]][j] <- cor(data_FPKM[,j],data_neighbor)
      }
    }

    sim_stat_combine <- do.call(cbind,sim_stat)
    sim_stat_mean <- apply(sim_stat_combine,2,mean)
    sim_stat_mean_diff <- diff(sim_stat_mean)

    x <- 1:length(sim_stat_mean_diff)
    neighbor_num <- which.min(sapply(x, function(k) {
      x2 <- pmax(0,x-k)
      sum(lm(sim_stat_mean_diff~x+x2)$residuals^2)
    }))
    neighbor_num <- neighbor_num + 1

    for(j in 1:ncol(data_FPKM)){
      nearest_cell_idx <- sort(data_dist[j,],decreasing=FALSE,index.return=TRUE)$ix[2:neighbor_num]
      data_expect[,j] <- apply(data_FPKM[,nearest_cell_idx],1,mean)
    }

  }

  return(list(data_expect=data_expect,neighbor_num=neighbor_num))

}

##estimate drop out
#' @export
estimate_drop_out <- function(sc_data,sc_data_expect,gene_len,per_tile_beta=4,per_tile_tau=4,alpha_init=c(1,-1),beta_init=c(0.1,0.1,0.1),tau_init=c(0.1,-0.1,-0.1),em_error_par=0.01,em_min_count=1,em_max_count=100,trace_flag=0){

  data_observe <- sc_data
	data_mui <- log(sc_data_expect+1)
	N_beta <- per_tile_beta
	N_tau <- per_tile_tau
	#data_mui_percentile_beta <- quantile(data_mui[which(data_mui>0)],seq(1/N_beta,(N_beta-1)/N_beta,1/N_beta))
	#data_mui_percentile_tau <- quantile(data_mui[which(data_mui>0)],seq(1/N_tau,(N_tau-1)/N_tau,1/N_tau))
	data_mui_percentile_beta <- (max(data_mui)/N_beta)*c(1:N_beta)
	data_mui_percentile_tau <- (max(data_mui)/N_tau)*c(1:N_tau)

	zero_idx <- which(data_observe==0)
	zero_indicator <- rep(0,length(data_observe))
	zero_indicator[zero_idx] <- 1

	##generate spline matrix for beta
	spline_mat_beta <- matrix(data=0,nrow=N_beta+3,ncol=length(data_observe))
	spline_mat_beta[1,] <- 1
	spline_mat_beta[2,] <- data_mui
	spline_mat_beta[3,] <- data_mui^2
	spline_mat_beta[4,] <- data_mui^3

	for(i in 5:(N_beta+3)){
		select_idx <- which(data_mui > data_mui_percentile_beta[i-4])
		spline_mat_beta[i,select_idx] <- (data_mui[select_idx]-data_mui_percentile_beta[i-4])^3
	}

	##generate spline matrix for tau
	spline_mat_tau <- matrix(data=0,nrow=N_tau+3,ncol=length(data_observe))
	spline_mat_tau[1,] <- 1
	spline_mat_tau[2,] <- data_mui
	spline_mat_tau[3,] <- data_mui^2
	spline_mat_tau[4,] <- data_mui^3

	for(i in 5:(N_tau+3)){
		select_idx <- which(data_mui > data_mui_percentile_tau[i-4])
		spline_mat_tau[i,select_idx] <- (data_mui[select_idx]-data_mui_percentile_tau[i-4])^3
	}

	##generate constrain matrix for beta
	data_mui_uni <- unique(data_mui)

	cons_mat <- matrix(data=0,nrow=N_beta+3,ncol=length(data_mui_uni))
	cons_mat[1,] <- 0
	cons_mat[2,] <- 1
	cons_mat[3,] <- 2*data_mui_uni
	cons_mat[4,] <- 3*data_mui_uni^2

	for(i in 5:(N_beta+3)){
		select_idx <- which(data_mui_uni > data_mui_percentile_beta[i-4])
		cons_mat[i,select_idx] <- 3*(data_mui_uni[select_idx]-data_mui_percentile_beta[i-4])^2
	}

	##initialize parameters
	beta_k <- c(beta_init,rep(0,N_beta))
	tau_k <- c(tau_init,rep(0,N_tau))
	alpha_k <- alpha_init
	lamda_k <- exp(colSums(beta_k*spline_mat_beta))
	phi_k <- exp(colSums(tau_k*spline_mat_tau))
	step_count <- 1
	flag <- 0
	logLik_trace <- matrix(data=NA,nrow=em_max_count,ncol=2)
  em_trace <- 1
  em_trace_count <- 1
	alpha_trace <- matrix(data=NA,nrow=em_max_count,ncol=2)
	beta_trace <- matrix(data=NA,nrow=em_max_count,ncol=length(beta_k))
	tau_trace <- matrix(data=NA,nrow=em_max_count,ncol=length(tau_k))

	while(flag==0 && step_count <= em_max_count){

		cat(step_count, " of maximum ", em_max_count,"EM steps\n")
		flush.console()
		##E step
		drop_out_z_kp1 <- 1/(1+exp(-alpha_k[1]-alpha_k[2]*data_mui)*(1/(1+gene_len*lamda_k*phi_k))^(1/phi_k))
		drop_out_z_kp1[which(data_observe > 0)] <- 0

    logLik_trace[step_count,2] <- sum(drop_out_z_kp1*((alpha_k[1]+alpha_k[2]*data_mui) - log(1+exp(alpha_k[1]+alpha_k[2]*data_mui)))*zero_indicator + (1-drop_out_z_kp1)*(-log(1+exp(alpha_k[1]+alpha_k[2]*data_mui)))) + sum((lgamma(data_observe+1/phi_k)-lgamma(1/phi_k)-(data_observe+1/phi_k)*log(1+gene_len*lamda_k*phi_k)+data_observe*log(lamda_k*phi_k))*(1-drop_out_z_kp1)) + sum((1-drop_out_z_kp1)*(data_observe*log(gene_len)-lgamma(data_observe+1)))

		##M step for alpha
		logi_fun <- function(alpha_est,drop_out_z_kp1_in,data_mui_in,zero_indicator_in){
			output <- sum(drop_out_z_kp1_in*((alpha_est[1]+alpha_est[2]*data_mui_in) - log(1+exp(alpha_est[1]+alpha_est[2]*data_mui_in)))*zero_indicator_in + (1-drop_out_z_kp1_in)*(-log(1+exp(alpha_est[1]+alpha_est[2]*data_mui_in))))
			return(output)
		}

		opt_result_logi <- optim(alpha_k, fn = logi_fun, gr = NULL, drop_out_z_kp1_in = drop_out_z_kp1, data_mui_in = data_mui, zero_indicator_in = zero_indicator, method = "Nelder-Mead",control = list(fnscale = -1))

		alpha_k <- opt_result_logi$par

		##M step for beta
		spline_fun <- function(param_est,beta_k_in,spline_mat_beta_in,spline_mat_tau_in,data_mui_in,drop_out_z_kp1_in,data_observe_in,gene_len_in){
			beta_est <- param_est[c(1:length(beta_k_in))]
			tau_est <- param_est[-c(1:length(beta_k_in))]

			sc_lamda <- exp(colSums(beta_est*spline_mat_beta_in))
			sc_phi <- exp(colSums(tau_est*spline_mat_tau_in))

			output <- sum((lgamma(data_observe_in+1/sc_phi)-lgamma(1/sc_phi)-(data_observe_in+1/sc_phi)*log(1+gene_len_in*sc_lamda*sc_phi)+data_observe_in*log(sc_lamda*sc_phi))*(1-drop_out_z_kp1_in))
			return(output)
		}

    if(em_trace == 0){
      ##monotonic contrain on lamda
      if(em_trace_count == 1){
        beta_k <- c(beta_init,rep(0,N_beta))
      }
		  opt_result <- constrOptim(c(beta_k,tau_k), f = spline_fun, grad = NULL, ui = t(rbind(cons_mat,matrix(data=0,nrow=N_tau+3,ncol=ncol(cons_mat)))), ci = rep(0,ncol(cons_mat)), control = list(fnscale = -1), beta_k_in = beta_k, spline_mat_beta_in = spline_mat_beta, spline_mat_tau_in = spline_mat_tau, data_mui_in = data_mui, drop_out_z_kp1_in = drop_out_z_kp1, data_observe_in = data_observe, gene_len_in = gene_len)
      em_trace_count <- em_trace_count + 1
    }else{
      ##without constrain on lamda
      opt_result <- optim(c(beta_k,tau_k), fn = spline_fun, gr = NULL, beta_k,spline_mat_beta,spline_mat_tau,data_mui,drop_out_z_kp1,data_observe,gene_len,method = "Nelder-Mead",control = list(fnscale = -1))
    }

    param_k <- opt_result$par
		beta_k <- param_k[c(1:length(beta_k))]
		tau_k <- param_k[-c(1:length(beta_k))]

		lamda_k <- exp(colSums(beta_k*spline_mat_beta))
		phi_k <- exp(colSums(tau_k*spline_mat_tau))

		logLik_trace[step_count,1] <- opt_result_logi$value + opt_result$value + sum((1-drop_out_z_kp1)*(data_observe*log(gene_len)-lgamma(data_observe+1)))

		alpha_trace[step_count,] <- alpha_k
	  beta_trace[step_count,] <- beta_k
		tau_trace[step_count,] <- tau_k

		##stop if converged
    if(step_count > em_min_count && em_trace_count > 2 && ((abs((logLik_trace[step_count,1] - logLik_trace[step_count-1,1])/logLik_trace[step_count-1,1]) < em_error_par) || (logLik_trace[step_count,1] < logLik_trace[step_count-1,1]))){
      flag <- 1
      message('EM converged in ',step_count,' steps')
      flush.console()
    }

    ##switch to constrOptim
		if(step_count > em_min_count && em_trace == 1 && ((abs((logLik_trace[step_count,1] - logLik_trace[step_count-1,1])/logLik_trace[step_count-1,1]) < em_error_par) || (logLik_trace[step_count,1] < logLik_trace[step_count-1,1]))){
			em_trace <- 0
    }

		step_count <- step_count + 1
	}

	alpha_out <- alpha_trace[step_count-1,]
	beta_out <- beta_trace[step_count-1,]
	tau_out <- tau_trace[step_count-1,]

	lamda_out <- exp(colSums(beta_out*spline_mat_beta))
	phi_out <- exp(colSums(tau_out*spline_mat_tau))

	drop_out_z_out <- 1/(1+exp(-alpha_out[1]-alpha_out[2]*data_mui)*(1/(1+gene_len*lamda_out*phi_out))^(1/phi_out))
	drop_out_z_out[which(data_observe > 0)] <- 0

	##get true expression
	data_true_out <- (data_observe*phi_out+1)/(gene_len*phi_out+1/lamda_out)
	data_true_out[which(data_observe == 0)] <- 0

	if(trace_flag==1){
		return(list(data_true=data_true_out,alpha_trace=alpha_trace,beta_trace=beta_trace,tau_trace=tau_trace,Loglik_trace=logLik_trace,
		spline_knot_beta_log=data_mui_percentile_beta,spline_knot_tau_log=data_mui_percentile_tau,post_weight=drop_out_z_out))
	}
	else{
		return(list(data_true=data_true_out,alpha_prior=alpha_out,beta_prior=beta_out,tau_prior=tau_out,spline_knot_beta_log=data_mui_percentile_beta,
		spline_knot_tau_log=data_mui_percentile_tau,post_weight=drop_out_z_out))
	}

}

##estimate drop out wrap up
#' @import parallel
#' @export
estimate_dropout_main <- function(sc_data_all,sc_data_expect_all,gene_len,ncore=1,per_tile_beta=4,per_tile_tau=4,alpha_init=c(1,-1),beta_init=c(0.1,0.1,0.1),tau_init=c(0.1,-0.1,-0.1),em_error_par=0.01,em_min_count=1,em_max_count=100,trace_flag=0){

  if(ncore > 1){
    result <- mclapply(c(1:ncol(sc_data_all)),function(i){estimate_drop_out(sc_data_all[,i],sc_data_expect_all[,i],gene_len,per_tile_beta,per_tile_tau,alpha_init,beta_init,tau_init,em_error_par,em_min_count,em_max_count,trace_flag)},mc.cores=ncore)
  }
  else{
    result <- lapply(c(1:ncol(sc_data_all)),function(i){estimate_drop_out(sc_data_all[,i],sc_data_expect_all[,i],gene_len,per_tile_beta,per_tile_tau,alpha_init,beta_init,tau_init,em_error_par,em_min_count,em_max_count,trace_flag)})
  }
  return(result)
}
##get weighted mean and variance
#' @export
get_weighted_stat <- function(data_in,weight_in){

	weight_norm <- weight_in/sum(weight_in)
	mean_weighted <- weighted.mean(data_in,weight_norm)
	var_weighted <- sum(weight_in*(data_in-mean_weighted)^2)/(sum(weight_in) - 1)

	return(list(mean_weighted=mean_weighted,var_weighted=var_weighted))

}

##adjust library size
#' @export
adjust_library_size <- function(dropout_est,hp_gene,gene_names){

  match_idx <- match(gene_names,hp_gene)
  data_hp <- sapply(dropout_est,function(x) x$data_true[!is.na(match_idx)])
  data_hp_weight <- sapply(dropout_est,function(x) x$post_weight[!is.na(match_idx)])

  row_median <- apply(data_hp,1,median)
  select_idx <- which(row_median>0)

  data_hp_sub <- data_hp[select_idx,]
  data_hp_weight_sub <- 1 - data_hp_weight[select_idx,]

  data_stat_weighted <- sapply(c(1:nrow(data_hp_sub)),function(x) get_weighted_stat(data_hp_sub[x,],data_hp_weight_sub[x,]))
  data_mean_weighted <- unlist(data_stat_weighted[1,])

  reference_mean <- mean(data_mean_weighted)
  library_scale <- sapply(1:ncol(data_hp_sub),function(x) reference_mean/get_weighted_stat(data_hp_sub[,x],data_hp_weight_sub[,x])$mean_weighted)

  return(library_scale)
}

##permutation test for differential expression
#' @export
permutation_test_mean <- function(data_1,weight_1,data_2,weight_2,num_permute=1000){

  stat_1 <- get_weighted_stat(data_1,weight_1)
  stat_2 <- get_weighted_stat(data_2,weight_2)
  n1 <- sum(weight_1)
  n2 <- sum(weight_2)
  group_id <- c(rep(1,length(data_1)),rep(2,length(data_2)))
  test_org <- (stat_1$mean_weighted - stat_2$mean_weighted)/sqrt((1/n1+1/n2)*((n1-1)*stat_1$var_weighted+(n2-1)*stat_2$var_weighted)/(n1+n2-2))
  data_combine <- c(data_1,data_2)
  weight_combine <- c(weight_1,weight_2)

  set.seed(12345)
  sample_mat <- t(sapply(1:num_permute,function(i) sample(group_id)))

  permut <- sapply(1:num_permute,function(id) {

    sample_group <- sample_mat[id,]
    sample_data_1 <- data_combine[sample_group==1]
    sample_data_2 <- data_combine[sample_group==2]
    sample_weight_1 <- weight_combine[sample_group==1]
    sample_weight_2 <- weight_combine[sample_group==2]

    stat_1 <- get_weighted_stat(sample_data_1,sample_weight_1)
    stat_2 <- get_weighted_stat(sample_data_2,sample_weight_2)
    n1 <- sum(sample_weight_1)
    n2 <- sum(sample_weight_2)
    (stat_1$mean_weighted - stat_2$mean_weighted)/sqrt((1/n1+1/n2)*((n1-1)*stat_1$var_weighted+(n2-1)*stat_2$var_weighted)/(n1+n2-2))

    })

  pval_greater <- mean(test_org < permut)
  pval_less <- mean(test_org > permut)
  pval_ts <- mean(abs(test_org) < abs(permut))

  return(list(pval_greater=pval_greater,pval_less=pval_less,pval_ts=pval_ts))
}

##permutation test for differential variance
#' @export
permutation_test_var <- function(data_1,weight_1,data_2,weight_2,num_permute=1000){

  stat_1 <- get_weighted_stat(data_1,weight_1)
  stat_2 <- get_weighted_stat(data_2,weight_2)
  n1 <- sum(weight_1)
  n2 <- sum(weight_2)
  group_id <- c(rep(1,length(data_1)),rep(2,length(data_2)))
  test_org <- stat_1$var_weighted/stat_2$var_weighted
  data_residual_combine <- c(data_1 - stat_1$mean_weighted,data_2 - stat_2$mean_weighted)
  weight_combine <- c(weight_1,weight_2)

  set.seed(12345)
  sample_mat <- t(sapply(1:num_permute,function(i) sample(group_id)))

  permut <- sapply(1:num_permute,function(id) {

    sample_group <- sample_mat[id,]
    sample_data_1 <- data_residual_combine[sample_group==1]
    sample_data_2 <- data_residual_combine[sample_group==2]
    sample_weight_1 <- weight_combine[sample_group==1]
    sample_weight_2 <- weight_combine[sample_group==2]
    n1 <- sum(sample_weight_1)
    n2 <- sum(sample_weight_2)
    (sum(sample_weight_1*(sample_data_1^2))/(n1 - 1)) / (sum(sample_weight_2*(sample_data_2^2))/(n2 - 1))

  })

  pval_greater <- mean(test_org < permut)
  pval_less <- mean(test_org > permut)
  pval_ts <- mean(abs(test_org) < abs(permut))

  return(list(pval_greater=pval_greater,pval_less=pval_less,pval_ts=pval_ts))
}

##get var, mean, and fitted var
#' @export
get_var_fit <- function(input_data,data_weight,span_param = 0.5){

	data_stat_weighted <- sapply(c(1:nrow(input_data)),function(x) get_weighted_stat(input_data[x,],data_weight[x,]))

	data_mean_weighted <- unlist(data_stat_weighted[1,])
	data_var_weighted <- unlist(data_stat_weighted[2,])

	data_reg <- data.frame(mean=log2(data_mean_weighted + 1),var=log2(data_var_weighted + 1))
	fitted_data <- loess(var ~ mean, data_reg,span=span_param)$fitted
	fitted_data[which(fitted_data < 0)] <- 0
	return(list(var_expect=fitted_data,mean=log2(data_mean_weighted + 1),var=log2(data_var_weighted + 1)))
}



##get estimates
#' @export
scdv_estimate <- function(treatment_data,treatment_data_weight,control_data,control_data_weight,span_param = 0.5){

	treatment_data <- as.matrix(treatment_data)
	control_data <- as.matrix(control_data)

	treatment_data_weight <- as.matrix(treatment_data_weight)
	control_data_weight <- as.matrix(control_data_weight)

	result_treatment <- get_var_fit(treatment_data,treatment_data_weight,span_param)
	var_expect_treatment <- result_treatment$var_expect
	mean_treatment <- result_treatment$mean
	var_treatment <- result_treatment$var

	result_control <- get_var_fit(control_data,control_data_weight,span_param)
	var_expect_control <- result_control$var_expect
	mean_control <- result_control$mean
	var_control <- result_control$var

	scale_factor_treatment <- var_treatment-var_expect_treatment
	scale_factor_control <- var_control-var_expect_control

	scale_factor_diff <- scale_factor_treatment - scale_factor_control

	return(list(scale_factor_treatment=scale_factor_treatment,scale_factor_control=scale_factor_control,
				var_expect_treatment=var_expect_treatment,var_expect_control=var_expect_control,
				mean_treatment=mean_treatment,mean_control=mean_control))
}


##check infinite data.frame
#' @export
is.infinite.data.frame <- function(obj){
    sapply(obj,FUN = function(x) all(is.infinite(x)))
}


##permute function
#' @export
scdv_permute <- function(treatment_data,treatment_data_weight,control_data,control_data_weight,var_expect_treatment,var_expect_control,per_time = 1000){

	df_treatment <- ncol(treatment_data)
	df_control <- ncol(control_data)

	treatmeat_data_mean <- unlist(sapply(c(1:nrow(treatment_data)),function(x) get_weighted_stat(treatment_data[x,],treatment_data_weight[x,]))[1,])
	treatment_data_residual <- (treatment_data - treatmeat_data_mean)/sqrt(2^var_expect_treatment-1)
	treatment_data_residual[is.na(treatment_data_residual)] <- 0

	for(k in 1:ncol(treatment_data_residual)){
		treatment_data_residual[is.infinite(treatment_data_residual[,k]),k] <- 0
	}

	control_data_mean <- unlist(sapply(c(1:nrow(control_data)),function(x) get_weighted_stat(control_data[x,],control_data_weight[x,]))[1,])
	control_data_residual <- (control_data - control_data_mean)/sqrt(2^var_expect_control-1)
	control_data_residual[is.na(control_data_residual)] <- 0

	for(k in 1:ncol(control_data_residual)){
		control_data_residual[is.infinite(control_data_residual[,k]),k] <- 0
	}

	combine_data <- cbind(treatment_data_residual,control_data_residual)
	combine_data <- combine_data^2

	combine_weight <- cbind(treatment_data_weight,control_data_weight)

	treatment_data_sf_per <- matrix(data=NA,nrow=nrow(treatment_data),ncol=per_time)
	control_data_sf_per <- matrix(data=NA,nrow=nrow(control_data),ncol=per_time)

	pb = txtProgressBar(min = 0, max = per_time, initial = 0, style = 3)

	set.seed(12345)
	for(i in 1:per_time){

		setTxtProgressBar(pb,i)

		per_idx <- sample(c(1:(df_treatment+df_control)),df_treatment)

		treatment_data_var_per <- sapply(1:nrow(combine_data),function(x) sum(combine_data[x,per_idx]*combine_weight[x,per_idx]/sum(combine_weight[x,per_idx])))
		control_data_var_per <- sapply(1:nrow(combine_data),function(x) sum(combine_data[x,-per_idx]*combine_weight[x,-per_idx]/sum(combine_weight[x,-per_idx])))

		treatment_data_sf_per[,i] <- log2(treatment_data_var_per+1)
		control_data_sf_per[,i] <- log2(control_data_var_per+1)
	}

	close(pb)
	return(list(treatment_data_sf_per=treatment_data_sf_per,control_data_sf_per=control_data_sf_per))
}


##permute function multi-core
#' @import parallel
#' @export
scdv_permute_mc <- function(treatment_data,treatment_data_weight,control_data,control_data_weight,var_expect_treatment,var_expect_control,per_time = 1000,ncore = 4){

	df_treatment <- ncol(treatment_data)
	df_control <- ncol(control_data)

	treatmeat_data_mean <- unlist(sapply(c(1:nrow(treatment_data)),function(x) get_weighted_stat(treatment_data[x,],treatment_data_weight[x,]))[1,])
	treatment_data_residual <- (treatment_data - treatmeat_data_mean)/sqrt(2^var_expect_treatment-1)
	treatment_data_residual[is.na(treatment_data_residual)] <- 0

	for(k in 1:ncol(treatment_data_residual)){
		treatment_data_residual[is.infinite(treatment_data_residual[,k]),k] <- 0
	}

	control_data_mean <- unlist(sapply(c(1:nrow(control_data)),function(x) get_weighted_stat(control_data[x,],control_data_weight[x,]))[1,])
	control_data_residual <- (control_data - control_data_mean)/sqrt(2^var_expect_control-1)
	control_data_residual[is.na(control_data_residual)] <- 0

	for(k in 1:ncol(control_data_residual)){
		control_data_residual[is.infinite(control_data_residual[,k]),k] <- 0
	}

	combine_data <- cbind(treatment_data_residual,control_data_residual)
	combine_data <- combine_data^2

	combine_weight <- cbind(treatment_data_weight,control_data_weight)

	treatment_data_sf_per <- matrix(data=NA,nrow=nrow(treatment_data),ncol=per_time)
	control_data_sf_per <- matrix(data=NA,nrow=nrow(control_data),ncol=per_time)

	set.seed(12345)

	permute_fun <- function(i){

		per_idx <- sample(c(1:(df_treatment+df_control)),df_treatment)

		treatment_data_var_per <- sapply(1:nrow(combine_data),function(x) sum(combine_data[x,per_idx]*combine_weight[x,per_idx]/sum(combine_weight[x,per_idx])))
		control_data_var_per <- sapply(1:nrow(combine_data),function(x) sum(combine_data[x,-per_idx]*combine_weight[x,-per_idx]/sum(combine_weight[x,-per_idx])))

		treatment_data_sf_per_temp <- log2(treatment_data_var_per+1)
		control_data_sf_per_temp <- log2(control_data_var_per+1)
		return(list(treatment_data_sf_per_temp=treatment_data_sf_per_temp,control_data_sf_per_temp=control_data_sf_per_temp))
	}

	output_list <- mclapply(c(1:per_time),permute_fun, mc.cores=ncore)

	for(i in 1:per_time){
		treatment_data_sf_per[,i] <- output_list[[i]]$treatment_data_sf_per_temp
		control_data_sf_per[,i] <- output_list[[i]]$control_data_sf_per_temp
	}

	return(list(treatment_data_sf_per=treatment_data_sf_per,control_data_sf_per=control_data_sf_per))
}


##main function
#' @import parallel
#' @export
scdv_main <- function(treatment_data,treatment_data_weight,control_data,control_data_weight,per_time=1000,span_param=0.5,ncore=1){

	message('Estimating variance scale factor')
	flush.console()
	result <- scdv_estimate(treatment_data,treatment_data_weight,control_data,control_data_weight,span_param)

	message('Permutation process to obtain empirical p-value')
	flush.console()
	if(ncore > 1){
		permute_result <- scdv_permute_mc(treatment_data,treatment_data_weight,control_data,control_data_weight,result$var_expect_treatment,result$var_expect_control,per_time,ncore)
	}
	else{
		permute_result <- scdv_permute(treatment_data,treatment_data_weight,control_data,control_data_weight,result$var_expect_treatment,result$var_expect_control,per_time)
	}

	sf_diff_pval <- rep(NA,nrow(treatment_data))
	sf_diff_pval_alt <- rep(NA,nrow(treatment_data))
	sf_diff_pval_ts <- rep(NA,nrow(treatment_data))

	sf_diff <- result$scale_factor_treatment - result$scale_factor_control

	for(i in 1:nrow(treatment_data)){
		sf_diff_null <- permute_result$treatment_data_sf_per[i,] - permute_result$control_data_sf_per[i,]

		sf_diff_pval[i] <- length(which(sf_diff_null >= sf_diff[i]))/per_time
		sf_diff_pval_alt[i] <- length(which(-sf_diff_null >= -sf_diff[i]))/per_time
		sf_diff_pval_ts[i] <- length(which(abs(sf_diff_null) >= abs(sf_diff[i])))/per_time
	}

	sf_diff_fdr <- p.adjust(sf_diff_pval,method="fdr")
	sf_diff_fdr_alt <- p.adjust(sf_diff_pval_alt,method="fdr")
	sf_diff_fdr_ts <- p.adjust(sf_diff_pval_ts,method="fdr")

	return(list(sf_treatment=result$scale_factor_treatment,sf_control=result$scale_factor_control,sf_diff=sf_diff,sf_diff_pval=sf_diff_pval,sf_diff_fdr=sf_diff_fdr,sf_diff_pval_alt=sf_diff_pval_alt,sf_diff_fdr_alt=sf_diff_fdr_alt,sf_diff_pval_ts=sf_diff_pval_ts,sf_diff_fdr_ts=sf_diff_fdr_ts))

}
