#' @title Multiple Pedigree Simulation
#' @description Simulates multiple pedigrees and groups them into one outputted dataset
#' of all individuals.
#' Based on pedigree simulation methodology described by Nieuwoudt et al. (2018).
#'
#' @param num Number of total pedigrees to simulate.
#' @param aff Number of disease cases among the carriers required for the pedigree to be ascertained.
#' @param frailty Indicator for whether include (1) or not include (0) a frailty term.
#' @param frail_var Value of the frailty variance, if frailty included.
#' @param dis_haz Hazard of the disease in the population of interest.
#' @param unaff_death_haz Hazard of death in the population of interest.
#' @param aff_death_haz Hazard of death in affected (diseased) population.
#' @param a_start Start year of the ascertainment interval (interval in which the proband develops the disease).
#' @param a_end Ending year of the ascertainment interval (interval in which the proband develops the disease).
#' @param f_start Start year of the founder birth interval (interval in which the founder is randomly born).
#' @param f_end Ending year of the founder birth interval (interval in which the founder is randomly born).
#' @param GRR Genetic relative risk i.e. the constant multiplier on the hazard for disease for carriers compared to non-carriers.
#' @param track_ped Indicator for whether the current pedigree number is printed (1) or not printed (0) after it is completed.
#' @param track_full Indicator for whether "full" is printed (1) or not (0) after one full cycle of pedigrees are completed.
#'
#' @details Pedigrees are simulated until one meets
#' the ascertainment requirements specified by the user, then this pedigree is added to
#' the dataset of specified number of pedigrees.
#' The founder is a carrier of a dominant rare genetic variant. Individuals in the pedigree can
#' have children, experience disease onset, and die. Waiting times for these events
#' will be generated when an individual is created, and will be in the form of their age when the event happens.
#' These events will continue for individuals until the end of the ascertainment interval,
#' or censoring from death. Death is based on unaff_death_haz and disease onset on dis_haz, and both
#' can only happen once per individual. A total number of children will be generated for an individual
#' and their births will be uniformly distributed between the ages of 18 and 45. If an individual has children,
#' they will be generated a non-carrier partner of the opposite sex who is over 18 and is within 5 years
#' of the age of the individual. The partner will not have children waiting times, but
#' will have disease onset and death waiting times. Children of carriers will have a 50%
#' probability of being a carrier. Children may not be born if the death waiting time
#' is earlier than the child's birth, or there is censoring from the end of the ascertainment
#' interval. After disease onset, a new waiting time will be generated using aff_death_haz.
#'
#' A proband is a carrier who has disease onset within the user specified ascertainment
#' interval. For the pedigree to be ascertained, this proband must exist, and there must
#' also be at least the number of user specified affected (diseased) carriers throughout
#' the study time. This is between the time the founder is born and the end of the ascertainment
#' interval. The founder is born uniformly within an interval specified by the user.
#'
#' A Gamma frailty term with a mean of 1 can be included to capture
#' between family heterogeneity in the general hazard of disease.
#'
#' @examples # Load example hazard datasets
#' data("dis_haz")
#' data("aff_death_haz")
#' data("unaff_death_haz")
#'
#' # Simulate multiple (e.g. 10) pedigrees
#' peds <- multiPed(num = 10, aff = 2, frailty = 1, frail_var = 4.519,
#' dis_haz = dis_haz, unaff_death_haz = unaff_death_haz, aff_death_haz = aff_death_haz,
#' a_start = 2005, a_end = 2025, f_start = 1885, f_end = 1935, GRR = 19,
#' track_ped = 1, track_full = 0)
#'
#' # View pedigrees dataset
#' peds
#'
#' @seealso [pedigree()]
#'
#' @references
#' Nieuwoudt, C., Jones, S. J., Brooks-Wilson, A., and Graham, J. (2018). Simulating pedigrees
#' ascertained for multiple disease-affected relatives. Source Code for Biology and Medicine,
#' 13(2). https://doi.org/10.1186/s13029-018-0069-6
#'
#' @export
### Simulate a specified number of pedigrees, output as one combined data frame
### Specify minimum number of cases per pedigree and frailty information

multiPed <- function(num, aff, frailty, frail_var, dis_haz, unaff_death_haz,
                     aff_death_haz, a_start, a_end, f_start, f_end, GRR, track_ped,
                     track_full){
  # Create first entry in data frames that can be added to
  full_df <- data.frame(matrix(0, nrow = 1, ncol = 15))
  names(full_df) <- c("FamID", "ID", "sex", "partner", "dadID", "mumID",
                      "affected","carrier", "birth", "onset", "death",
                      "proband", "event_age","start_age", "ind")
  # Simulate specified number of pedigrees
  for(j in 1:num){
    pedigree <- myPed(aff, frailty, frail_var, j, dis_haz, unaff_death_haz,
                      aff_death_haz, a_start, a_end, f_start, f_end, GRR)
    # Keep track of how many pedigrees have been simulated so far (optional)
    if(track_ped == 1){
      print(j)
    }
    # Collect all carriers for each method data frames
    full_df <- dplyr::bind_rows(full_df, full_c(pedigree, a_end))
  }
  # Remove the initial row using to first create the data frames
  final_full <- full_df[-1,]
  # Show that one full cycle is complete (optional)
  if(track_full == 1){
    print("full")
  }
  # Output the data frame of individuals
  return(final_full)
}
