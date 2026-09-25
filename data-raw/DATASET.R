

# Age-specific hazard rate of death for gastric cancer patients
aff_death_haz = 0.536

# Age-specific hazard rate of gastric cancer for general population
dis_haz <- c(rep(0.05,15), rep(0.1,5),rep(0.3,5),rep(0.6,5),rep(1.2,5),rep(2.4,5),
             rep(4.0,5),rep(6.1,5),rep(8.8,5),rep(12.7,5),rep(17.7,5),rep(24.7,5),
             rep(31.0,5),rep(39.4,5),rep(44.8,5),rep(48.7,15))/100000

# Age-specific hazard rate of death in general population
unaff_death_haz <- c(558,rep(28,4),rep(15.3,10),rep(79.5,10), rep(163.4,10),
                     rep(255.4,10),rep(453.3,10),rep(992.1,10),rep(1978.7,10),
                     rep(4708.2,10),rep(14389.6,15))/100000

usethis::use_data(aff_death_haz,dis_haz,unaff_death_haz,overwrite = TRUE)
