###############################################
##### TRIAL DATA CLEANING AND COMBINATION #####
###############################################

# setwd("/Users/ijs11/Documents/GitHub/trials/src/1_data_cleaning")

# packages
pacman::p_load(tidyverse,
               rio,
               fastDummies,
               lubridate,
               boot)

# main dataset
master <- import("HetMal_master_cluster_prev_data.csv") |>
  rename(cluster = Cluster_no)

#######################
##### PROTOPOPOFF #####
#######################

p_master <- master |> subset(Trial_no == 29)

##### net use
bl_net <- import("Baseline_Net_usage_per_Cluster_Muleba_JackieCook.dta") |>
  rename(net_use_num = sleepundernet,
         net_use_denom = totres)

##### trial net use
tr_net <- import("Post_Net_usage_per_Cluster_Muleba_JackieCook.dta") |>
  mutate(tr_net_use_num = sleepundernet,
         tr_net_use_denom = totres,
         time = case_when(
           XS_id == "A" ~ 4,
           XS_id == "B" ~ 9,
           XS_id == "C" ~ 16,
           XS_id == "D" ~ 21,
           XS_id == "E" ~ 28,
           XS_id == "F" ~ 33
         )) |>
  filter(cluster < 49) |> # one extra cluster not included
  select(cluster, time, tr_net_use_num, tr_net_use_denom)

##### post trial ITN distribution
prev <- p_master |>
  filter(Survey_no != 1) |>
  mutate(time = case_when(
    Survey_no == 2 ~ 4,
    Survey_no == 3 ~ 9,
    Survey_no == 4 ~ 16,
    Survey_no == 5 ~ 21,
    Survey_no == 6 ~ 28,
    Survey_no == 7 ~ 33
  ),
  date = lubridate::dmy(paste0("15/",Survey_month,"/",Survey_year)),
  seasonality = paste0("protopopoff_",lubridate::quarter(date)),
  ) |>
  select(cluster,
         time,
         prev_num = n_Pf_pos,
         prev_denom = N_tested,
         date,
         seasonality)

##### baseline prevalence
prev_BL <- p_master |>
  filter(Survey_no == 1) |>
  select(cluster,
         BL_prev_denom = N_tested,
         BL_prev_num = n_Pf_pos)

##### full dataset
PROTO <-  left_join(prev, prev_BL) |>
  left_join(bl_net) |>
  left_join(tr_net) |>
  mutate(Trial_no = 29) |>
  left_join(master |>
              select(Trial_no, cluster, Arm) |> distinct()) |> # adding the arm
  filter(!str_detect(Arm, "IRS")) # no arms with IRS

# checks
# tests to check clusters are all included

u_cl <- unique(p_master[, c("Trial_ID", "Trial_no", "Cluster_ID", "cluster", "Arm")]) |> filter(!str_detect(Arm, "IRS")) |> select(cluster) |> unlist() |> as.vector()

if(!(PROTO |> group_by(time) |> summarise(c = all(cluster %in% u_cl)) |> select(c) |> all())){
  warning("not all clusters identified in PROTO")
}

if(!(PROTO |> group_by(time) |> summarise(c = all(u_cl %in% cluster)) |> select(c) |> all())){
  warning("clusters are missing in PROTO")
}

rm(list = c("p_master", "prev", "prev_BL", "bl_net", "tr_net", "u_cl"))

#################
##### MOSHA #####
#################

m_master <- master |> subset(Trial_no == 12)

# checking the clusters and arms match
arm_check <- import("post_XS_ABCDE_cluster_level_int_netuse.dta") |> select(cluster = clusterno, intervention) |>
  left_join(m_master |> select(cluster, Arm) |> distinct())

if(any(arm_check$Arm != arm_check$intervention)){
  warning("Mosha cluster arm post codes do not match")
}

# trial net use
tr_net <- import("post_XS_ABCDE_cluster_level_int_netuse.dta") |>
  mutate(time = case_when(
    survey == "postA" ~ 12,
    survey == "postB" ~ 18,
    survey == "postC" ~ 24,
    survey == "postD" ~ 30,
    survey == "postE" ~ 36
  )
  ) |>
  select(cluster = clusterno,
         time,
         tr_net_use_num = netuse,
         tr_net_use_denom = sleepundernet)

# baseline net use
# arm level numbers of people tested for use of nets at baseline are from Table 1
# from https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext?rss=yes&sf163018506=1#supplementaryMaterial
# assumed to be equal among clusters

bl_net <- import("Baseline _cluster_itnuse.dta") |>
  left_join(m_master |> select(Trial_no, cluster, Arm) |> distinct()) |>
  left_join(m_master |> group_by(Arm) |> summarise(n = n()) |>
              mutate(net_use_denom_tot = case_when(Arm == "Std_llin" ~ 4962,
                                                   Arm == "PPF_Py" ~ 4520,
                                                   Arm == "CFP_Py" ~ 4803,
                                                   Arm == "PBO_Py" ~ 4369),
                     net_use_denom = round(net_use_denom_tot/n, digits = 0)
              ) |>
              select(net_use_denom, Arm)
  ) |>
  mutate(net_use_num = round(itnuse/100 * net_use_denom, digits = 0)) |>
  select(cluster, net_use_num, net_use_denom, Arm)

### prevalence
prev <- m_master |>
  filter(Survey_no != 1) |>
  mutate(Survey_year = ifelse(Survey_no == 2, 2020, Survey_year)) |>
  mutate(time = case_when(Survey_no == 2 ~ 12,
                          Survey_no == 3 ~ 18,
                          Survey_no == 4 ~ 24,
                          Survey_no == 5 ~ 30,
                          Survey_no == 6 ~ 36),
         date = lubridate::dmy(paste0("15/",Survey_month,"/",Survey_year)),
         seasonality = paste0("mosha_",lubridate::quarter(date))) |>
  select(cluster, time,
         prev_num = n_Pf_pos,
         prev_denom = N_tested,
         Pf_prev, Arm, date, seasonality)

# checking the prevalence estimates match
if(!all(round(prev$Pf_prev, digits = 7) == round(prev$prev_num/prev$prev_denom, digits = 7))){
  warning("not all prevalence estimates correspond - Mosha")
}

##### baseline prevalence
prev_BL <- m_master |>
  filter(Survey_no == 1) |>
  mutate(BL_prev_num = round(Pf_prev * N_tested, digits = 0)) |> # the numbers n_Pf_pos are not integers at baseline
  select(cluster,
         BL_prev_num,
         BL_prev_denom = N_tested,
         BL_Pf_prev = Pf_prev)

MOSHA <- left_join(prev |> select(-Pf_prev),
                   prev_BL |> select(-BL_Pf_prev)
) |>
  left_join(bl_net) |>
  left_join(tr_net) |>
  mutate(Trial_no = 12) |>
  left_join(master |> select(Trial_no, cluster, Arm) |> distinct())

# checks
# tests to check clusters are all included

u_cl <- unique(m_master[, c("cluster", "Arm")]) |> select(cluster) |> unlist() |> as.vector()

if(!(MOSHA |> group_by(time) |> summarise(c = all(cluster %in% u_cl)) |> select(c) |> all())){
  warning("not all clusters identified in MOSHA")
}

if(!(MOSHA |> group_by(time) |> summarise(c = all(u_cl %in% cluster)) |> select(c) |> all())){
  warning("clusters are missing in MOSHA")
}

rm(list = c("prev", "prev_BL", "bl_net", "tr_net", "u_cl", "m_master"))

###################
##### STAEDKE #####
###################

dt2 <- import("dt2.csv") |>
  mutate(distr_date = lubridate::dmy(paste0("15/", month,"/", year.1)))

s_master <- master |> subset(Trial_no == 14) |> mutate(cluster = parse_number(Cluster_ID))

# 3 clusters had no dominant net so were excluded
# https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30214-2/fulltext#seccestitle80
# removing these clusters

r_cl <- unique(s_master$cluster)[!(unique(s_master$cluster) %in% unique(dt2$cluster))]
u_cl <- unique(s_master$cluster)[(unique(s_master$cluster) %in% unique(dt2$cluster))]

s_master <- s_master |> subset(cluster %in% u_cl)

if(!all(unique(s_master$cluster) %in% unique(dt2$cluster)) |
   any(s_master$Arm == "")){
  warning("clusters not correct for Staedke")
}

### baseline net use
bl_net <- dt2 |>
  select(cluster,
         arm,
         net_use = ITN_Use_prior,
         net_use_denom = Number_tested_baseline) |>
  mutate(net_use_num = round(net_use*net_use_denom,0)) |>
  select(-net_use)

### trial net use
tr_net_use <- dt2 |>
  select(cluster,
         sleeping_under_net_6m,
         sleeping_under_net_12m,
         sleeping_under_net_18m,
         sleeping_under_net_25m) |>
  pivot_longer(values_to = "coverage", names_to = "month", cols = c(2:5)) |>
  mutate(time = parse_number(month)) |>
  select(-month)

tr_net_tested <- dt2 |>
  select(cluster,
         Number_tested_6m,
         Number_tested_12m,
         Number_tested_18m) |>
  mutate(Number_tested_25m = round((Number_tested_6m + Number_tested_12m + Number_tested_18m)/3)) |>
  pivot_longer(values_to = "tested", names_to = "month", cols = c(2:5)) |>
  mutate(time = parse_number(month))|>
  select(-month)

tr_net <- left_join(tr_net_use, tr_net_tested) |>
  mutate(tr_net_use_num = round(coverage*tested,0),
         tr_net_use_denom = tested) |>
  select(cluster, time, tr_net_use_num, tr_net_use_denom)

### prevalence
# numbers tested at 25m are missing - so assumed to be approximately the mean value for the other times for each clusters
prev_25m <- dt2 |>
  mutate(cluster = cluster,
         prev_denom = round((Number_tested_6m + Number_tested_12m + Number_tested_18m)/3),
         prev_num = round(Prevalence_25m*prev_denom,0),
         Pf_prev = Prevalence_25m,
         time = 25,
         #date = distr_date + months(time),
         .keep = "none")

prev <- s_master |>
  filter(Survey_no != 1) |>
  mutate(time = case_when(Survey_no == 2 ~ 6,
                          Survey_no == 3 ~ 12,
                          Survey_no == 4 ~ 18)
  ) |>
  select(cluster, time, prev_num = n_Pf_pos, prev_denom = N_tested, Pf_prev) |>
  bind_rows(prev_25m) |>
  left_join(dt2 |> select(distr_date, cluster)) |>
  mutate(date = distr_date + months(time),
         seasonality = paste0("staedke_",quarter(date))) |>
  select(-distr_date)

# baseline_prevalence
prev_BL <- s_master |>
  filter(Survey_no == 1) |>
  select(Pf_prev, cluster, BL_prev_num = n_Pf_pos, BL_prev_denom = N_tested)

if(!all(round(prev_BL$Pf_prev, digits = 5) == round(prev_BL$BL_prev_num/prev_BL$BL_prev_denom, digits = 5))){
  warning("Staedke baseline prevalence does not match")
}

STAEDKE <- left_join(prev, prev_BL |> select(-Pf_prev)) |>
  left_join(bl_net |> select(-arm)) |>
  left_join(tr_net) |>
  mutate(Trial_no = 14) |>
  left_join(s_master |> select(Trial_no, cluster, Arm) |> distinct())

# checks
if(!(STAEDKE |> group_by(time) |> summarise(c = all(cluster %in% u_cl)) |> select(c) |> all())){
  warning("not all clusters identified in STAEDKE")
}

if(!(STAEDKE |> group_by(time) |> summarise(c = all(u_cl %in% cluster)) |> select(c) |> all())){
  warning("clusters are missing in STAEDKE")
}

# checking the number of clusters are correct
if(any((STAEDKE |> group_by(Arm, time) |> summarise(n = n()) |> ungroup() |>filter(Arm == "LLIN") |> select(n)!=52))){
  warning("numbers of clusters in LLIN arm are not correct")
}

if(any((STAEDKE |> group_by(Arm, time) |> summarise(n = n()) |> ungroup() |>filter(Arm == "PBO_LLIN") |> select(n)!=49))){
  warning("numbers of clusters in PBO_LLIN arm are not correct")
}

# numbers of missing values in the prevalence
# 13 values are missing at 25 months
STAEDKE |> group_by(time) |> summarise(n = sum(is.na(prev_num))) |> filter(n > 0)

STAEDKE <- STAEDKE |> filter(!is.na(prev_num))

rm(list = c("bl_net", "tr_net_use", "tr_net_tested", "tr_net", "prev_25m", "prev", "prev_BL"))

#################
##### BENIN #####
#################

Benin <- import("Benin_cluster_prev_nets.xlsx") |> rename(Arm = arm)

# missing numbers of people asked about use of nets at baseline
# values obtained from Table 1
# https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)02319-4/fulltext#seccestitle150

bl_net <- Benin |>
  filter(survey == "1.XS Baseline") |>
  group_by(Arm) |>
  mutate(n = n(),
         useanynet_n = case_when(Arm == "Control" ~ round(1392/n, digits = 0),
                                 Arm == "CFP net" ~ round(1326/n, digits = 0),
                                 Arm == "PPF net" ~ round(1370/n, digits = 0)
         )
  ) |>
  mutate(cluster = cluster,
         net_use_num = round(useanynet*useanynet_n,0),
         net_use_denom = useanynet_n,
         .keep = "none")

tr_net <- Benin |>
  filter(survey != "1.XS Baseline") |>
  mutate(cluster = cluster,
         tr_net_use_num = useanynet_n,
         tr_net_use_denom = use,
         time = case_when(
           survey == "3.XS 6mo" ~ 6,
           survey == "5.XS 18mo" ~ 18,
           survey == "7.XS 30mo" ~ 30),
         .keep = "none")

prev <- Benin |>
  filter(survey != "1.XS Baseline") |>
  mutate(cluster = cluster,
         prev_num = pos,
         prev_denom = notested,
         time = case_when(
           survey == "3.XS 6mo" ~ 6,
           survey == "5.XS 18mo" ~ 18,
           survey == "7.XS 30mo" ~ 30),
         .keep = "none") |>
  mutate(
    date = lubridate::dmy("15/03/2020") + months(time),
    seasonality = paste0("accrombessi_",lubridate::quarter(date))
  )

prev_BL <- Benin |>
  filter(survey == "1.XS Baseline") |>
  select(cluster, BL_prev_num = pos, BL_prev_denom = notested)

ACCROMBESSI <- left_join(prev, prev_BL) |>
  left_join(bl_net) |>
  left_join(tr_net) |>
  mutate(Trial_no = 11) |>
  left_join(master |> select(Trial_no, cluster, Arm) |> distinct())

rm(list = c("prev_BL", "prev", "bl_net", "tr_net"))

# Combined
COMBO_stan <- bind_rows(
  PROTO |> mutate(author = "Protopopoff"),
  MOSHA |>  mutate(author = "Mosha"),
  STAEDKE |> select(-Pf_prev) |> mutate(author = "Staedke"),
  ACCROMBESSI |> mutate(author = "Accrombessi")) |>
  mutate(
    # net
    l = case_when(
      str_detect(Arm, "CFP") ~ 3,
      str_detect(Arm, "PPF") ~ 4,
      str_detect(Arm, "PBO") ~ 2,
      str_detect(Arm, "llin") | str_detect(Arm, "LLIN") | str_detect(Arm, "Control") ~ 1
    ),
    # study
    i = case_when(
      str_detect(author, "Protopopoff") ~ 1,
      str_detect(author, "Staedke") ~ 2,
      str_detect(author, "Mosha") ~ 3,
      str_detect(author, "Accrombessi") ~ 4
    ),
    # net study
    li = case_when(
      l == 1 ~ 1, # if pyrethroid only nets then it is coded as 1, these should not be included in pmat, so just the intercept is used for these clusters
      i == 1 & l == 2 ~ 2,
      i == 2 & l == 2 ~ 3,
      i == 3 & l == 2 ~ 4,
      i == 3 & l == 3 ~ 5,
      i == 3 & l == 4 ~ 6,
      i == 4 & l == 3 ~ 7,
      i == 4 & l == 4 ~ 8,
      .default = NA
    ),
    study = case_when(
      i == 1 ~ "Protopopoff et al",
      i == 2 ~ "Staedke et al",
      i == 3 ~ "Mosha et al",
      i == 4 ~ "Accrombessi et al"
    ),
    net = case_when(
      l == 1 ~ "Pyrethroid-only",
      l == 2 ~ "Pyrethroid-PBO",
      l == 3 ~ "Pyrethroid-pyrrole",
      l == 4 ~ "Pyrethroid-pyriproxyfen"
    ),
    prev = prev_num / prev_denom,
    prev_lower = qbeta(0.025, 1 + prev_num, 1 + prev_denom - prev_num),
    prev_upper = qbeta(0.975, 1 + prev_num, 1 + prev_denom - prev_num),
    it = seasonality) |>
  ######################################################
  ##### not including pyrethroid-pyriproxyfen nets #####
  #######################################################
  filter(l != 4)

# calculating the beta-binomial posteriors
# |>
#   group_by(i, cluster) |>
#   mutate(cluster_index = ifelse(l == 1,1,cur_group_id())) |>
#   group_by(cluster_index) |>
#   mutate(lij = cur_group_id()) |>
#   ungroup() |>
#   select(-c(cluster_index, Trial_no))


#######################################
##### checking for missing values #####
#######################################

check_na <- function(colname){
  c <- COMBO_stan[is.na(COMBO_stan[,colname]),]
  c |> group_by(author)
  a <- c[,"author"] |> unique()
  if(nrow(c) > 0){
    print(paste(colname, "has ",nrow(c), "NAs in", a))
  } else{
    print("no NAs")
  }
}

##### response variable
check_na("prev_num")
check_na("prev_denom")

##### explanatory variable
check_na("BL_prev_denom")
check_na("BL_prev_num")
check_na("net_use_denom")
check_na("tr_net_use_denom")
check_na("time")
check_na("seasonality")
check_na("l")
check_na("i")
check_na("li")
check_na("Arm")
check_na("date")

# variables with NAs
check_na("net_use_num")

# interpolating the values
na_pos <- which(is.na(COMBO_stan$net_use_num))

na_a <- unique(COMBO_stan[na_pos,"author"])[[1]]

set.seed(123)
fit <- lme4::glmer(cbind(net_use_num, net_use_denom - net_use_num) ~ (1|cluster), data = COMBO_stan[-na_pos,] |> subset(author == na_a),
                   family = binomial(link = "logit"))

# predicting just using the intercept
COMBO_stan[na_pos, "net_use_num"] <- rbinom(n = 1,
                                            prob = boot::inv.logit(summary(fit)$coefficients[1, "Estimate"]),
                                            size = COMBO_stan[na_pos, "net_use_denom"] |> unlist() |> unname() |> as.vector() |> unique())

# variables with NAs
check_na("tr_net_use_num")

na_pos_t <- which(is.na(COMBO_stan$tr_net_use_num))
na_a_t <- unique(COMBO_stan[na_pos_t,"author"])[[1]]

set.seed(123)
fit_t <- lme4::glmer(cbind(tr_net_use_num, tr_net_use_denom - tr_net_use_num) ~ time + (1|cluster), data = COMBO_stan[-na_pos_t,] |> subset(author == na_a_t),
                     family = binomial(link = "logit"))

# predicting just using the intercept
COMBO_stan[na_pos_t, "tr_net_use_num"] <- rbinom(n = length(na_pos_t),
                                                 prob = predict(fit_t, newdata = COMBO_stan[na_pos_t, ], type = "response"),
                                                 size = COMBO_stan[na_pos_t, "tr_net_use_denom"] |> unlist() |> unname() |> as.vector())

# checking all values are integers
int <- COMBO_stan |> select(prev_num, prev_denom, BL_prev_denom, BL_prev_num, net_use_denom, net_use_num, tr_net_use_denom) %%1

if(all(int!=0)){
  warning("some values in COMBO_stan are not integers")
}

######################################################
##### creating the unique and study combinations #####
######################################################

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

  # check
  # prob_s_train <- u_ij_train$BL_prev_num / u_ij_train$BL_prev_denom
  # count_m_prob <- rep(0, N_i)
  # total_m_prob <- rep(0, N_i)
  # prob_s <- rep(NA, N_ij)
  #
  # for(i in 1:N_ij){
  #   if(ij_train[i] == 1){
  #     prob_s[i] = prob_s_train[r_id_ij[i]];
  #     count_m_prob[r_id[i]] <- count_m_prob[r_id[i]] + 1;
  #     total_m_prob[r_id[i]] <- total_m_prob[r_id[i]] + prob_s[i];
  #   } else{
  #     prob_s[i] = 0;
  #   }
  # }
  #
  # m_prob_s_i <- rep(NA, N_i)
  # for(i in 1:N_i){
  #   if(ij_unq_train[i] == 1){
  #     m_prob_s_i[i] = total_m_prob[i] / count_m_prob[i];
  #   } else{
  #     m_prob_s_i[i] = 1;
  #   }
  # }
  #
  # d_prob_s_ij <- rep(NA, N_ij)
  # for(i in 1:N_ij){
  #   d_prob_s_ij[i] = prob_s[i] - m_prob_s_i[r_id[i]];
  # }
  #
  # d_prob_s = as.matrix(pmat_ij) %*% d_prob_s_ij;
  # m_prob_s = as.matrix(pmat_i) %*% m_prob_s_i;

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

data_in_full <- extract_data_in(COMBO_stan, train_inds = seq(1, nrow(COMBO_stan)))

iter <- 8000
warmup <- floor(iter/2)

lower <- 0.025
upper <- 0.975

##### interaction with baseline prevalence and and baseline net use
set.seed(123)

K_folds <- 5

COMBO_stan <- COMBO_stan |> dplyr::mutate(
  BL_prev = BL_prev_num/BL_prev_denom,
  BL_prev_lower = qbeta(0.025, 1 + BL_prev_num, 1 + BL_prev_denom - BL_prev_num),
  BL_prev_upper = qbeta(0.975, 1 + BL_prev_num, 1 + BL_prev_denom - BL_prev_num),
  net_use = net_use_num/net_use_denom,
  net_use_lower = qbeta(0.025, 1 + net_use_num, 1 + net_use_denom - net_use_num),
  net_use_upper = qbeta(0.975, 1 + net_use_num, 1 + net_use_denom - net_use_num),
  trial_net_use = tr_net_use_num/tr_net_use_denom,
  trial_net_use_lower = qbeta(0.025, 1 + tr_net_use_num, 1 + tr_net_use_denom - tr_net_use_num),
  trial_net_use_upper = qbeta(0.975, 1 + tr_net_use_num, 1 + tr_net_use_denom - tr_net_use_num),
  fold = loo::kfold_split_stratified(x = paste0(COMBO_stan[,"l"], "_", COMBO_stan[,"i"]), K = K_folds)
)

saveRDS(list("iter" = iter,
             "warmup" = warmup,
             "lower" = lower,
             "upper" = upper,
             "data_in_full" = data_in_full,
             "COMBO_stan" = COMBO_stan,
             "extract_data_in" = extract_data_in,
             "K_folds" = K_folds),
        file = "stan_data.rds"
)
