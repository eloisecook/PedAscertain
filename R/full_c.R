#' @importFrom rlang .data

full_c <- function(allped,a_end)
{
  ped_df <- as.data.frame(allped)
  # Define event_age for individuals, which is age at onset or censoring from
  # death or the end of the study
  ped_df <- ped_df |>
    dplyr::mutate(event_age = pmin(.data$onset - .data$birth, .data$death - .data$birth,
                            a_end - .data$birth, na.rm = TRUE))
  # Define start_age for individuals, which is always 0 for the full group
  ped_df <- ped_df |> dplyr::mutate(start_age = 0)
  ped_df <- dplyr::filter(ped_df, .data$event_age > .data$start_age)
  return(ped_df)
}
