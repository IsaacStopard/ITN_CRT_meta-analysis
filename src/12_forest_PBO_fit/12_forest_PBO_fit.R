library(rstan)
library(tidyverse)
library(loo)
library(patchwork)

ggplot2::theme_set(theme_bw() +
                     theme(text = element_text(size = 17),
                           legend.text = element_text(size = 11),
                           legend.title = element_text(size = 11),
                           panel.background = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text = element_text(colour = "black", size = 17),
                           strip.text = element_text(size = 17))
)

orderly2::orderly_dependency("1_data_cleaning", "latest()", c("stan_data.rds"))

orderly2::orderly_shared_resource("m1_re.stan")

orderly2::orderly_dependency("1_data_cleaning", "latest()",
                             c("stan_data.rds"))

list2env(readRDS("stan_data.rds"), envir = .GlobalEnv)

extract_data_in <- function(COMBO_stan, train_inds){

  u_ij <- unique(COMBO_stan[,c("i", "cluster", "BL_prev_num", "BL_prev_denom", "net_use_num", "net_use_denom")]) |> dplyr::mutate(ij = dplyr::row_number())
  COMBO_stan <- COMBO_stan |> dplyr::left_join(u_ij)

  # design matrices for dummy coding
  pmat_i <- fastDummies::dummy_cols(COMBO_stan$i)[,-1]
  pmat_l <- fastDummies::dummy_cols(COMBO_stan$l)[, -c(1, 2)] # data_1 is the pyrethroid only nets
  pmat_li <- fastDummies::dummy_cols(COMBO_stan$li)[,-c(1, 2)]
  pmat_ij <- fastDummies::dummy_cols(COMBO_stan$ij)[,-1]
  pmat_it <- fastDummies::dummy_cols(COMBO_stan$it)[,-1]

  pmat_i_pooled_gq <- matrix(0, nrow = nrow(pmat_i), ncol = ncol(pmat_i))
  pmat_li_pooled_gq <- matrix(0, nrow = nrow(pmat_li), ncol = ncol(pmat_li))
  pmat_it_pooled_gq <- matrix(0, nrow = nrow(pmat_it), ncol = ncol(pmat_it))

  N_ij <- as.integer(nrow(u_ij))
  N_ij_unq <- length(unique(u_ij$i))

  # calculating the indexing to match up pmat_ij with missing clusters and studies
  u_ij_train <- unique(COMBO_stan[train_inds, c("i", "cluster", "ij", "BL_prev_num", "BL_prev_denom", "net_use_num", "net_use_denom")]) |>
    dplyr::mutate(ij_train = dplyr::row_number())

  N_ij_train <- nrow(u_ij_train)
  N_ij_unq_train <- length(unique(u_ij_train$i))

  ij_unq_train_df <- data.frame(i = sort(unique(u_ij$i))) |> dplyr::mutate(ij_unq_train = i %in% unique(u_ij_train$i))

  x <- 0
  for(i in 1:nrow(ij_unq_train_df)){
    if(ij_unq_train_df[i, "ij_unq_train"]){
      x <- x+1
      ij_unq_train_df[i, "r_id_ij_unq_train"] <- x
    } else{
      ij_unq_train_df[i, "r_id_ij_unq_train"] <- -1 # index position that doesn't exist
    }
  }

  ij_unq_train <- ij_unq_train_df$ij_unq_train |> as.integer()
  r_id_ij_unq <- ij_unq_train_df$r_id_ij_unq_train # index of positions corresponding to fitted sigma

  #
  ij_train_df <- data.frame(ij = sort(unique(u_ij$ij))) |> dplyr::mutate(ij_train = ij %in% unique(u_ij_train$ij))
  x <- 0
  for(i in 1:nrow(ij_train_df)){
    if(ij_train_df[i, "ij_train"]){
      x <- x+1
      ij_train_df[i, "r_id_ij_train"] <- x
    } else{
      ij_train_df[i, "r_id_ij_train"] <- -1 # index position that doesn't exist
    }
  }

  ij_train <- ij_train_df$ij_train |> as.integer()
  r_id_ij <- ij_train_df$r_id_ij_train # index of positions corresponding to fitted cluster random effect

  m_prob_s_gq <- u_ij |> group_by(i) |> summarise(m = sum(BL_prev_num)/sum(BL_prev_denom)) |> select(m) |> unlist() |> as.vector()
  d_prob_s_gq <- as.vector(as.matrix(pmat_ij) %*% (u_ij$BL_prev_num / u_ij$BL_prev_denom - m_prob_s_gq[u_ij$i]))
  m_prob_s_gq <- as.vector((as.matrix(pmat_i) %*% m_prob_s_gq))

  N_gq <- nrow(COMBO_stan)

  m_prob_s_gq_pooled <- rep(u_ij |> summarise(m = sum(BL_prev_num)/sum(BL_prev_denom)) |> select(m) |> as.vector() |> unname() |> unlist(), N_gq)
  d_prob_s_gq_pooled <- rep(0, N_gq)

  m_net_s_gq <- u_ij |> group_by(i) |> summarise(m = sum(net_use_num)/sum(net_use_denom)) |> select(m) |> unlist() |> as.vector()
  d_net_s_gq <- as.vector(as.matrix(pmat_ij) %*% (u_ij$net_use_num / u_ij$net_use_denom - m_net_s_gq[u_ij$i]))
  m_net_s_gq <- as.vector((as.matrix(pmat_i) %*% m_net_s_gq))

  m_net_s_gq_pooled <- rep(u_ij |> summarise(m = sum(net_use_num)/sum(net_use_denom)) |> select(m) |> as.vector() |> unname() |> unlist(), N_gq)
  d_net_s_gq_pooled <- rep(0, N_gq)

  r_id <- u_ij$i
  N_i <- ncol(pmat_i)

  u_li <- unique(COMBO_stan[,c("l", "li", "i")]) |> arrange(li) |> filter(l!=1)
  #r_id_li <- u_li$l - 1
  N_i_pbo <- nrow(subset(u_li, l == 2))
  N_i_pp <- nrow(subset(u_li, l == 3))

  data_in <- list("N" = nrow(COMBO_stan),
                  "N_i" = N_i,
                  "N_l" = ncol(pmat_l),
                  "N_li" = ncol(pmat_li),
                  "N_ij" = N_ij,
                  "N_it" = ncol(pmat_it),
                  "pmat_i" = pmat_i,
                  "pmat_l" = pmat_l,
                  "pmat_li" = pmat_li,
                  "pmat_ij" = pmat_ij,
                  "pmat_it" = pmat_it,
                  #"r_id_li" = r_id_li,
                  "N_i_pbo" = N_i_pbo,
                  "N_i_pp" = N_i_pp,
                  "r_id" = r_id,
                  "N_ij_unq" = N_ij_unq, # cluster specific random effect for each study
                  "ij_train" = ij_train,
                  "r_id_ij" = r_id_ij,
                  "ij_unq_train" = ij_unq_train,
                  "r_id_ij_unq" = r_id_ij_unq,
                  "N_ij_train" = N_ij_train,
                  "N_ij_unq_train" = N_ij_unq_train,
                  "time" = COMBO_stan$time / 12,
                  "years" = COMBO_stan$time / 12,
                  "pos" = as.integer(COMBO_stan$prev_num),
                  "test" = as.integer(COMBO_stan$prev_denom),
                  "pos_s" = as.integer(u_ij$BL_prev_num),
                  "test_s" = as.integer(u_ij$BL_prev_denom),
                  "base_net_pos" = as.integer(u_ij$net_use_num),
                  "base_net_test" = as.integer(u_ij$net_use_denom),
                  "prior_sd" = 2.0,
                  "prior_sd_t" = 2.0,
                  "N_train" = length(train_inds),
                  "train_inds" = train_inds,
                  "N_gq" = N_gq,
                  "N_ij_gq" = N_ij,
                  "N_ij_unq_gq" = N_ij_unq,
                  "pmat_ij_gq" = pmat_ij,
                  "pmat_i_gq" = pmat_i,
                  "pmat_l_gq" = pmat_l,
                  "pmat_li_gq" = pmat_li,
                  "pmat_it_gq" = pmat_it,
                  "pmat_i_pooled_gq" = pmat_i_pooled_gq,
                  "pmat_li_pooled_gq" = pmat_li_pooled_gq,
                  "pmat_it_pooled_gq" = pmat_it_pooled_gq,
                  "r_id_gq" = u_ij$i,
                  "ij_unq_train_gq" = ij_unq_train,
                  "ij_train_gq" = ij_train,
                  "m_prob_s_gq" = m_prob_s_gq,
                  "d_prob_s_gq" = d_prob_s_gq,
                  "years_gq" = COMBO_stan$time / 12,
                  "gq" = 0,
                  "m_net_s_gq" = m_net_s_gq,
                  "d_net_s_gq" = d_net_s_gq,
                  "m_net_s_gq_pooled" = m_net_s_gq_pooled,
                  "d_net_s_gq_pooled" = d_net_s_gq_pooled,
                  "m_prob_s_gq_pooled" = m_prob_s_gq_pooled,
                  "d_prob_s_gq_pooled" = d_prob_s_gq_pooled
  )
  return(data_in)
}

COMBO_stan_subset <- subset(COMBO_stan, study != "Protopopoff et al") |>
  mutate(
  # study
  i = case_when(
    str_detect(author, "Staedke") ~ 1,
    str_detect(author, "Mosha") ~ 2,
    str_detect(author, "Accrombessi") ~ 3
  ),
  # net study
  li = case_when(
    l == 1 ~ 1, # if pyrethroid only nets then it is coded as 1, these should not be included in pmat, so just the intercept is used for these clusters
    i == 1 & l == 2 ~ 2,
    i == 2 & l == 2 ~ 3,
    i == 2 & l == 3 ~ 4,
    i == 3 & l == 3 ~ 5,
    .default = NA
  ))

data_in_subset <- extract_data_in(COMBO_stan_subset, train_inds = seq(1, nrow(COMBO_stan_subset)))

saveRDS(COMBO_stan_subset, file = "COMBO_stan_PBO_subset.rds")
saveRDS(data_in_subset, file = "data_in_PBO_subset.rds")

adapt <- 0.999
m_td <- 10.5
step_in <- 0.5

stan_model_in <- rstan::stan_model(file = "m1_re.stan")

fit_full <- rstan::sampling(stan_model_in,
                            data = data_in_subset,
                            iter = iter,
                            warmup = warmup,
                            control = list(max_treedepth = m_td, adapt_delta = adapt, stepsize = step_in),
                            cores = 4,
                            chains = 4,
                            seed = 123)

saveRDS(fit_full, file = "m1_re_PBO_subset_fit_full.rds")

