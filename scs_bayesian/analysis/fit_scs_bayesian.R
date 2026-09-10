# Reproducible Bayesian SCS analyses in R/CmdStan.
# The analysis is restricted to permanent implant recipients and contains no data.
# Runs a specified primary, sensitivity, predictive-check, or CV stage.
# Source data are de-identified and the primary estimand is conditional on permanent implantation.

suppressPackageStartupMessages({library(cmdstanr); library(posterior)})
set.seed(20260907)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- Sys.getenv("SCS_R_OUT", unset = file.path(root, "outputs"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN", unset = cmdstanr::cmdstan_path()))
stage <- Sys.getenv("SCS_STAGE", unset = "primary")
data_path <- Sys.getenv("SCS_DATA_CSV", unset = file.path(root, "data", "analysis_cohort.csv"))

source_data <- function() {
  d <- read.csv(data_path, check.names = FALSE)
  d <- d[d$Implantation == 1, , drop = FALSE]
  expected_n <- as.integer(Sys.getenv("SCS_EXPECTED_N", unset = "246"))
  expected_events <- as.integer(Sys.getenv("SCS_EXPECTED_EVENTS", unset = "121"))
  stopifnot(nrow(d) == expected_n, sum(d$`6M_Success`) == expected_events, !anyNA(d))
  d
}

make_x <- function(d, scaling = NULL, quadratic = FALSE, centre = FALSE) {
  nrs_col <- if ("NRS" %in% names(d)) "NRS" else if ("PreVAS" %in% names(d)) "PreVAS" else stop("Input requires NRS (or legacy PreVAS) column")
  raw <- data.frame(Age_sd=d$Age, BMI_sd=d$BMI, Duration_sd=d$Duration,
                    NRS_sd=d[[nrs_col]], log_MED_plus1_sd=log1p(d$MED))
  if (is.null(scaling)) scaling <- lapply(raw, function(x) c(mean=mean(x), sd=sqrt(mean((x-mean(x))^2))))
  x <- as.data.frame(Map(function(z,s) (z-s[["mean"]])/s[["sd"]], raw, scaling))
  if (quadratic) {
    z2 <- x$Duration_sd^2
    if (is.null(scaling$Duration_sq_sd)) scaling$Duration_sq_sd <- c(mean=mean(z2),sd=sqrt(mean((z2-mean(z2))^2)))
    x$Duration_sq_sd <- (z2-scaling$Duration_sq_sd[["mean"]])/scaling$Duration_sq_sd[["sd"]]
  }
  x$Sex <- d$Sex; x$Alcohol <- d$Alcohol; x$Smoking <- d$Smoking; x$Psy <- d$Psy
  for (lv in c(1,2,3)) x[[paste0("Pay_",lv,"_vs_0")]] <- as.integer(d$Pay == lv)
  for (lv in c(1,3,4)) x[[paste0("UpperLower_",lv,"_vs_2")]] <- as.integer(d$UpperLower == lv)
  if (centre) x$Hospital_2_vs_1 <- as.integer(d$Hospital == 2)
  list(x=as.matrix(x), scaling=scaling)
}

write_model <- function(kind) {
  f <- file.path(out, paste0("stan_",kind,".stan"))
  if (kind == "hierarchical") {
    lines <- c(
      "data { int<lower=1> N; int<lower=1> K; int<lower=1> G; matrix[N,K] X; array[N] int<lower=1,upper=G> diagnosis; array[N] int<lower=0,upper=1> y; real<lower=0> beta_sd; real<lower=0> sigma_sd; }",
      "parameters { real alpha; vector[K] beta; real<lower=0> sigma_dx; vector[G] z_dx; }",
      "transformed parameters { vector[G] dx_effect = z_dx * sigma_dx; }",
      "model { alpha ~ normal(0,2.5); beta ~ normal(0,beta_sd); sigma_dx ~ normal(0,sigma_sd); z_dx ~ normal(0,1); y ~ bernoulli_logit(alpha + X*beta + dx_effect[diagnosis]); }",
      "generated quantities { vector[N] p; for (n in 1:N) p[n] = inv_logit(alpha + X[n]*beta + dx_effect[diagnosis[n]]); }"
    )
  } else if (kind == "fixed") {
    lines <- c(
      "data { int<lower=1> N; int<lower=1> K; int<lower=1> G; matrix[N,K] X; matrix[N,G-1] D; array[N] int<lower=0,upper=1> y; real<lower=0> beta_sd; }",
      "parameters { real alpha; vector[K] beta; vector[G-1] gamma_dx; }",
      "model { alpha ~ normal(0,2.5); beta ~ normal(0,beta_sd); gamma_dx ~ normal(0,beta_sd); y ~ bernoulli_logit(alpha + X*beta + D*gamma_dx); }",
      "generated quantities { vector[N] p; for (n in 1:N) p[n] = inv_logit(alpha + X[n]*beta + D[n]*gamma_dx); }"
    )
  } else {
    lines <- c(
      "data { int<lower=1> N; int<lower=1> K; matrix[N,K] X; array[N] int<lower=0,upper=1> y; real<lower=0> beta_sd; }",
      "parameters { real alpha; vector[K] beta; }",
      "model { alpha ~ normal(0,2.5); beta ~ normal(0,beta_sd); y ~ bernoulli_logit(alpha + X*beta); }",
      "generated quantities { vector[N] p; for (n in 1:N) p[n] = inv_logit(alpha + X[n]*beta); }"
    )
  }
  writeLines(lines,f); cmdstan_model(f,quiet=TRUE)
}

fit_one <- function(d, kind="hierarchical", beta_sd=1, sigma_sd=.5, quadratic=FALSE, centre=FALSE,
                    seed=20260907, warmup=1500, sampling=1500, adapt_delta=.95, output_basename=NULL) {
  xx <- make_x(d, quadratic=quadratic, centre=centre)
  levels <- 1:5
  dx <- match(d$Dx,levels)
  stan_data <- list(N=nrow(d),K=ncol(xx$x),X=xx$x,y=as.integer(d$`6M_Success`),beta_sd=beta_sd)
  if(kind=="hierarchical") stan_data <- c(stan_data,list(G=5L,diagnosis=dx,sigma_sd=sigma_sd))
  if(kind=="fixed") {
    D <- sapply(2:5,function(g) as.integer(dx==g)); colnames(D) <- paste0("Dx_",2:5,"_vs_1")
    stan_data <- c(stan_data,list(G=5L,D=D))
  }
  model <- write_model(kind)
  fit <- model$sample(data=stan_data,seed=seed,chains=4,parallel_chains=2,
                      iter_warmup=warmup,iter_sampling=sampling,adapt_delta=adapt_delta,
                      refresh=250)
  list(fit=fit,x=xx$x,scaling=xx$scaling,levels=levels,kind=kind,dx=dx)
}

diagnostics <- function(fit) {
  sm <- fit$summary(variables=NULL, posterior::default_convergence_measures())
  sd <- fit$sampler_diagnostics(format="draws_array")
  energy <- sd[, , "energy__"]
  ebfmi <- apply(energy,2,function(e) mean(diff(e)^2)/var(e))
  list(max_rhat=max(sm$rhat,na.rm=TRUE),min_bulk_ess=min(sm$ess_bulk,na.rm=TRUE),min_tail_ess=min(sm$ess_tail,na.rm=TRUE),
       divergences=sum(sd[, , "divergent__"]),min_ebfmi=min(ebfmi),summary=sm)
}

save_fit <- function(obj,label,save_draws=TRUE) {
  fit <- obj$fit; dg <- diagnostics(fit)
  fit$save_object(file.path(out,paste0(label,"_fit.rds")))
  write.csv(dg$summary,file.path(out,paste0(label,"_diagnostics.csv")),row.names=FALSE)
  if(save_draws) {
    draw_variables <- switch(obj$kind,
      hierarchical=c("alpha","beta","sigma_dx","z_dx","dx_effect"),
      fixed=c("alpha","beta","gamma_dx"),
      none=c("alpha","beta"))
    write.csv(as.data.frame(fit$draws(variables=draw_variables,format="draws_df")),file.path(out,paste0(label,"_draws.csv")),row.names=FALSE)
  }
  write.csv(data.frame(feature=colnames(obj$x)),file.path(out,paste0(label,"_features.csv")),row.names=FALSE)
  as.data.frame(c(dg[1:5],list(kind=obj$kind,beta_sd=NA, sigma_sd=NA)))
}

coefficient_table <- function(fit,features,label) {
  d <- as.data.frame(fit$draws(variables="beta",format="draws_df"))
  ans <- do.call(rbind,lapply(seq_along(features),function(i){z<-d[[paste0("beta[",i,"]")]];data.frame(term=features[i],or_mean=mean(exp(z)),or_2.5=quantile(exp(z),.025),or_97.5=quantile(exp(z),.975),p_or_gt_1=mean(z>0))}))
  write.csv(ans,file.path(out,paste0(label,"_coefficients.csv")),row.names=FALSE); ans
}

posterior_probs <- function(fit) {
  d <- fit$draws(variables="p",format="draws_matrix")
  colMeans(d)
}

apparent_metrics <- function(fit,y,label) {
  d <- fit$draws(variables="p",format="draws_matrix")
  d <- d[seq(1,nrow(d),length.out=min(1000,nrow(d))),,drop=FALSE]
  auc <- apply(d,1,function(p) {r<-rank(p); (sum(r[y==1])-sum(seq_len(sum(y==1))))/(sum(y==1)*sum(y==0))})
  brier <- apply(d,1,function(p) mean((y-p)^2))
  pp <- colMeans(d); z <- qlogis(pmin(pmax(pp,1e-6),1-1e-6)); cc <- glm(y~z,family=binomial())
  ans<-data.frame(apparent_auc_mean=mean(auc),apparent_auc_2.5=quantile(auc,.025),apparent_auc_97.5=quantile(auc,.975),apparent_brier_mean=mean(brier),apparent_brier_2.5=quantile(brier,.025),apparent_brier_97.5=quantile(brier,.975),apparent_calibration_intercept=coef(cc)[1],apparent_calibration_slope=coef(cc)[2])
  write.csv(ans,file.path(out,paste0(label,"_apparent_performance.csv")),row.names=FALSE); ans
}

run_primary <- function() {
  d<-source_data(); obj<-fit_one(d,adapt_delta=.99,output_basename="primary")
  meta<-save_fit(obj,"primary"); coefficient_table(obj$fit,colnames(obj$x),"primary"); apparent_metrics(obj$fit,d$`6M_Success`,"primary")
  write.csv(d,file.path(out,"analysis_cohort.csv"),row.names=FALSE)
  write.csv(data.frame(diagnosis_code=1:5),file.path(out,"diagnosis_levels.csv"),row.names=FALSE)
  write.csv(data.frame(n=nrow(d),successes=sum(d$`6M_Success`),failures=sum(!d$`6M_Success`),missing_cells=sum(is.na(d))),file.path(out,"cohort_manifest.csv"),row.names=FALSE)
  write.csv(meta,file.path(out,"primary_convergence_summary.csv"),row.names=FALSE)
}

run_sensitivity <- function(name) {
  d<-source_data(); spec<-switch(name,
    sensitivity_prior_narrow=list(kind="hierarchical",beta_sd=.5,sigma_sd=.5,quadratic=FALSE,centre=FALSE,adapt_delta=.99),
    sensitivity_prior_wide=list(kind="hierarchical",beta_sd=2.5,sigma_sd=.5,quadratic=FALSE,centre=FALSE,adapt_delta=.99),
    sensitivity_diagnosis_fixed=list(kind="fixed",beta_sd=1,sigma_sd=.5,quadratic=FALSE,centre=FALSE,adapt_delta=.95),
    sensitivity_duration_quadratic=list(kind="hierarchical",beta_sd=1,sigma_sd=.5,quadratic=TRUE,centre=FALSE,adapt_delta=.99),
    sensitivity_institution_adjusted=list(kind="hierarchical",beta_sd=1,sigma_sd=.5,quadratic=FALSE,centre=TRUE,adapt_delta=.95),
    sensitivity_no_diagnosis=list(kind="none",beta_sd=1,sigma_sd=.5,quadratic=FALSE,centre=FALSE,adapt_delta=.95),
    stop("Unknown sensitivity stage"))
  obj<-do.call(fit_one,c(list(d=d,seed=20260907+nchar(name)),spec)); dg<-save_fit(obj,name); co<-coefficient_table(obj$fit,colnames(obj$x),name); ap<-apparent_metrics(obj$fit,d$`6M_Success`,name)
  focus<-co[co$term %in% c("Duration_sd","log_MED_plus1_sd","Sex"),]
  write.csv(cbind(data.frame(scenario=name,kind=spec$kind,beta_sd=spec$beta_sd,sigma_sd=spec$sigma_sd,quadratic=spec$quadratic,centre=spec$centre),dg,ap),file.path(out,paste0(name,"_summary.csv")),row.names=FALSE)
  write.csv(focus,file.path(out,paste0(name,"_focus.csv")),row.names=FALSE)
}

run_cv_fold <- function(k) {
  d<-source_data(); set.seed(20260907); fold<-sample(rep(1:5,length.out=nrow(d)))
  # Stratified assignment is recreated deterministically using outcome-specific shuffled allocation.
  fold<-integer(nrow(d)); for(y in 0:1){ii<-which(d$`6M_Success`==y); ii<-sample(ii); fold[ii]<-rep(1:5,length.out=length(ii))}
  tr<-which(fold!=k); te<-which(fold==k)
  obj<-fit_one(d[tr,,drop=FALSE],seed=20261000+k,warmup=1000,sampling=1000,adapt_delta=.99)
  dg<-save_fit(obj,paste0("cv_fold",k),save_draws=FALSE)
  xt<-make_x(d[te,,drop=FALSE],scaling=obj$scaling)$x
  xt<-xt[,colnames(obj$x),drop=FALSE]
  dr<-as.data.frame(obj$fit$draws(variables=c("alpha","beta","dx_effect"),format="draws_df")); beta<-as.matrix(dr[,paste0("beta[",seq_len(ncol(obj$x)),"]")]); alpha<-dr$alpha
  de<-as.matrix(dr[,paste0("dx_effect[",1:5,"]")]); dx<-match(d$Dx[te],1:5)
  pp<-colMeans(plogis(alpha + beta%*%t(xt) + de[,dx,drop=FALSE]))
  write.csv(data.frame(row=te,outcome=d$`6M_Success`[te],cv_predicted_probability=pp),file.path(out,paste0("cv_fold",k,"_predictions.csv")),row.names=FALSE)
  write.csv(dg,file.path(out,paste0("cv_fold",k,"_summary.csv")),row.names=FALSE)
}

run_cv_no_diagnosis_fold <- function(k) {
  d<-source_data(); set.seed(20260907); fold<-integer(nrow(d)); for(y in 0:1){ii<-which(d$`6M_Success`==y); ii<-sample(ii); fold[ii]<-rep(1:5,length.out=length(ii))}
  tr<-which(fold!=k); te<-which(fold==k); obj<-fit_one(d[tr,,drop=FALSE],kind="none",seed=20262000+k,warmup=1000,sampling=1000,adapt_delta=.95)
  dg<-diagnostics(obj$fit); write.csv(dg$summary,file.path(out,paste0("cv_no_diagnosis_fold",k,"_diagnostics.csv")),row.names=FALSE)
  xt<-make_x(d[te,,drop=FALSE],scaling=obj$scaling)$x; xt<-xt[,colnames(obj$x),drop=FALSE]
  dr<-as.data.frame(obj$fit$draws(variables=c("alpha","beta"),format="draws_df")); beta<-as.matrix(dr[,paste0("beta[",seq_len(ncol(obj$x)),"]")]); pp<-colMeans(plogis(dr$alpha + beta%*%t(xt)))
  write.csv(data.frame(row=te,outcome=d$`6M_Success`[te],cv_predicted_probability=pp),file.path(out,paste0("cv_no_diagnosis_fold",k,"_predictions.csv")),row.names=FALSE)
  write.csv(as.data.frame(dg[1:5]),file.path(out,paste0("cv_no_diagnosis_fold",k,"_summary.csv")),row.names=FALSE)
}

run_cv_performance <- function() {
  p<-do.call(rbind,lapply(1:5,function(k)read.csv(file.path(out,paste0("cv_fold",k,"_predictions.csv")))));p<-p[order(p$row),];y<-p$outcome; pr<-p$cv_predicted_probability
  auc<-function(y,p){r<-rank(p);(sum(r[y==1])-sum(seq_len(sum(y==1))))/(sum(y==1)*sum(y==0))}; calc<-function(ii){z<-qlogis(pmin(pmax(pr[ii],1e-6),1-1e-6)); cc<-glm(y[ii]~z,family=binomial());c(auc(y[ii],pr[ii]),mean((y[ii]-pr[ii])^2),coef(cc)[1],coef(cc)[2])}
  set.seed(20260907); boot<-t(replicate(2000,{ii<-sample(seq_along(y),replace=TRUE);if(length(unique(y[ii]))<2) return(rep(NA,4));calc(ii)}));boot<-boot[complete.cases(boot),]
  point<-calc(seq_along(y)); ans<-data.frame(cv_auc=point[1],cv_auc_2.5=quantile(boot[,1],.025),cv_auc_97.5=quantile(boot[,1],.975),cv_brier=point[2],cv_brier_2.5=quantile(boot[,2],.025),cv_brier_97.5=quantile(boot[,2],.975),cv_calibration_intercept=point[3],cv_calibration_intercept_2.5=quantile(boot[,3],.025),cv_calibration_intercept_97.5=quantile(boot[,3],.975),cv_calibration_slope=point[4],cv_calibration_slope_2.5=quantile(boot[,4],.025),cv_calibration_slope_97.5=quantile(boot[,4],.975))
  write.csv(ans,file.path(out,"cross_validated_performance.csv"),row.names=FALSE);write.csv(p,file.path(out,"cross_validated_predictions.csv"),row.names=FALSE)
  bins<-cut(pr,breaks=quantile(pr,seq(0,1,.1)),include.lowest=TRUE); tab<-aggregate(cbind(mean_pred=pr,observed=y),by=list(bin=bins),FUN=mean); tab$n<-as.integer(table(bins)); tab<-tab[,c("bin","n","mean_pred","observed")]; write.csv(tab,file.path(out,"cross_validated_calibration_deciles.csv"),row.names=FALSE)
}

run_cv_no_diagnosis_performance <- function() {
  p<-do.call(rbind,lapply(1:5,function(k)read.csv(file.path(out,paste0("cv_no_diagnosis_fold",k,"_predictions.csv")))));p<-p[order(p$row),];y<-p$outcome; pr<-p$cv_predicted_probability
  auc<-function(y,p){r<-rank(p);(sum(r[y==1])-sum(seq_len(sum(y==1))))/(sum(y==1)*sum(y==0))}; calc<-function(ii){z<-qlogis(pmin(pmax(pr[ii],1e-6),1-1e-6));cc<-glm(y[ii]~z,family=binomial());c(auc(y[ii],pr[ii]),mean((y[ii]-pr[ii])^2),coef(cc)[1],coef(cc)[2])}
  set.seed(20260907); boot<-t(replicate(2000,{ii<-sample(seq_along(y),replace=TRUE);if(length(unique(y[ii]))<2)return(rep(NA,4));calc(ii)}));boot<-boot[complete.cases(boot),];point<-calc(seq_along(y));ans<-data.frame(cv_auc=point[1],cv_auc_2.5=quantile(boot[,1],.025),cv_auc_97.5=quantile(boot[,1],.975),cv_brier=point[2],cv_brier_2.5=quantile(boot[,2],.025),cv_brier_97.5=quantile(boot[,2],.975),cv_calibration_intercept=point[3],cv_calibration_slope=point[4]);write.csv(ans,file.path(out,"cv_no_diagnosis_performance.csv"),row.names=FALSE)
}

run_ppc <- function() {
  d<-source_data(); fit<-readRDS(file.path(out,"primary_fit.rds")); p<-fit$draws(variables="p",format="draws_matrix");set.seed(20260907);rep<-rowSums(matrix(rbinom(length(p),1,as.vector(p)),nrow=nrow(p)));obs<-sum(d$`6M_Success`);write.csv(data.frame(statistic="replicated_success_count",observed=obs,posterior_predictive_mean=mean(rep),ppc_2.5=quantile(rep,.025),ppc_97.5=quantile(rep,.975),two_sided_tail_probability=min(1,2*min(mean(rep<=obs),mean(rep>=obs))),row.names=NULL),file.path(out,"posterior_predictive_check.csv"),row.names=FALSE)
}

run_prior <- function() {
  d<-source_data(); xx<-make_x(d)$x; set.seed(20260907); m<-5000L
  alpha<-rnorm(m,0,2.5); beta<-matrix(rnorm(m*ncol(xx),0,1),m,ncol(xx)); sigma<-abs(rnorm(m,0,.5)); z<-matrix(rnorm(m*5),m,5); dx<-match(d$Dx,1:5)
  p<-plogis(alpha + beta%*%t(xx) + (z*sigma)[,dx,drop=FALSE])
  write.csv(data.frame(quantile=c(.025,.5,.975),patient_probability=quantile(p,c(.025,.5,.975)),cohort_mean_probability=quantile(rowMeans(p),c(.025,.5,.975))),file.path(out,"prior_predictive_check.csv"),row.names=FALSE)
}

postprocess_sensitivities <- function() {
  d <- source_data()
  specs <- list(
    sensitivity_prior_narrow=list(kind="hierarchical",beta_sd=.5,sigma_sd=.5,quadratic=FALSE,centre=FALSE),
    sensitivity_prior_wide=list(kind="hierarchical",beta_sd=2.5,sigma_sd=.5,quadratic=FALSE,centre=FALSE),
    sensitivity_diagnosis_fixed=list(kind="fixed",beta_sd=1,sigma_sd=.5,quadratic=FALSE,centre=FALSE),
    sensitivity_duration_quadratic=list(kind="hierarchical",beta_sd=1,sigma_sd=.5,quadratic=TRUE,centre=FALSE),
    sensitivity_institution_adjusted=list(kind="hierarchical",beta_sd=1,sigma_sd=.5,quadratic=FALSE,centre=TRUE),
    sensitivity_no_diagnosis=list(kind="none",beta_sd=1,sigma_sd=.5,quadratic=FALSE,centre=FALSE)
  )
  for (name in names(specs)) {
    spec <- specs[[name]]; fit <- readRDS(file.path(out,paste0(name,"_fit.rds"))); x <- make_x(d,quadratic=spec$quadratic,centre=spec$centre)$x
    dg <- diagnostics(fit); write.csv(dg$summary,file.path(out,paste0(name,"_diagnostics.csv")),row.names=FALSE)
    co <- coefficient_table(fit,colnames(x),name); ap <- apparent_metrics(fit,d$`6M_Success`,name)
    focus <- co[co$term %in% c("Duration_sd","log_MED_plus1_sd","Sex"),]
    write.csv(focus,file.path(out,paste0(name,"_focus.csv")),row.names=FALSE)
    write.csv(cbind(data.frame(scenario=name,kind=spec$kind,beta_sd=spec$beta_sd,sigma_sd=spec$sigma_sd,quadratic=spec$quadratic,centre=spec$centre),as.data.frame(dg[1:5]),ap),file.path(out,paste0(name,"_summary.csv")),row.names=FALSE)
  }
}

if(stage=="primary") run_primary() else if(stage %in% c("sensitivity_prior_narrow","sensitivity_prior_wide","sensitivity_diagnosis_fixed","sensitivity_duration_quadratic","sensitivity_institution_adjusted","sensitivity_no_diagnosis")) run_sensitivity(stage) else if(stage=="sensitivity_all") for(s in c("sensitivity_prior_narrow","sensitivity_prior_wide","sensitivity_diagnosis_fixed","sensitivity_duration_quadratic","sensitivity_institution_adjusted","sensitivity_no_diagnosis")) run_sensitivity(s) else if(grepl("^cv_fold[1-5]$",stage)) run_cv_fold(as.integer(sub("cv_fold","",stage))) else if(stage=="cv_all") for(k in 1:5) run_cv_fold(k) else if(grepl("^cv_no_diagnosis_fold[1-5]$",stage)) run_cv_no_diagnosis_fold(as.integer(sub("cv_no_diagnosis_fold","",stage))) else if(stage=="cv_no_diagnosis_all") for(k in 1:5) run_cv_no_diagnosis_fold(k) else if(stage=="cv_performance") run_cv_performance() else if(stage=="cv_no_diagnosis_performance") run_cv_no_diagnosis_performance() else if(stage=="posterior_predictive_check") run_ppc() else if(stage=="prior_predictive_check") run_prior() else if(stage=="postprocess_sensitivities") postprocess_sensitivities() else stop("Unknown SCS_STAGE")


