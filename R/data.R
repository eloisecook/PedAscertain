#' Constant hazard of death for affected (diseased) individuals.
#'
#' A constant hazard rate of death for individuals affected by the disease that is
#' linked to the rare variant (e.g. gastric cancer).
#'
#' @format A numeric value.
#'
#' @source Schauer, M., Peiper, M., Theisen, J., and Knoefel, W. (2011). Prognostic factors in
#'    patients with diffuse type gastric cancer (linitis plastica) after operative treatment. European
#'    Journal of Medical Research, 16(1):29–33.
"aff_death_haz"


#' Age-specific hazard of disease (e.g. gastric cancer) in the general population.
#'
#' A numeric vector containing the age-specific hazard rates of
#' gastric cancer in the general population. The rates are expressed
#' per person-year.
#'
#' @format A numeric vector of length 100.
#'
#' @source SEER (2025).
#'    SEER*Explorer: An interactive website for SEER cancer statistics
#'    Surveillance Research Program, National Cancer Institute.
#'    https://seer.cancer.gov/statistics-network/explorer/. Data source(s): SEER Incidence
#'    Data, November 2024 Submission (1975-2022)
"dis_haz"


#' Age-specific hazard of death in the general population.
#'
#' A numeric vector containing the age-specific hazard rates of death
#' in the general population. The rates are expressed per person-year.
#'
#' @format A numeric vector of length 100.
#'
#' @source Xu, J., Murphy, S. L., Kochanek, K., and Arias, E. (2025). Deaths: Final data for 2022.
#'    National Vital Statistics Reports Vol. 74, No. 4, National Center for Health Statistics.
"unaff_death_haz"
