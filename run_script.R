rm(list = ls())

orderly2::orderly_run("1_data_cleaning")
orderly2::orderly_run("2_fits")

rm(list = ls())
orderly2::orderly_run("3_cv", list(model = "m0"))
rm(list = ls())
orderly2::orderly_run("3_cv", list(model = "m1"))
rm(list = ls())
orderly2::orderly_run("3_cv", list(model = "m2"))

rm(list = ls())
orderly2::orderly_run("4_time_m0_c_fits", list(model = "m0"))

orderly2::orderly_run("5_two_stage_meta_analysis")
orderly2::orderly_run("6_malariasim")
orderly2::orderly_run("7_malariasim_results")

orderly2::orderly_run("8_forest")

rm(list = ls())
orderly2::orderly_run("9_inv_logit")

rm(list = ls())
orderly2::orderly_run("10_cov_results")

orderly2::orderly_run("11_study_descriptions")

orderly2::orderly_run("12_forest_PBO_fit")

rm(list = ls())
orderly2::orderly_run("13_forest_PBO_plots")

orderly2::orderly_validate_archive(action = "orphan")

orderly2::orderly_prune_orphans()
