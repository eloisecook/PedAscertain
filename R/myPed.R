#' @importFrom rlang .data

myPed <- function(aff, frailty, frail_var, j, dis_haz, unaff_death_haz, aff_death_haz,
                  a_start, a_end, f_start, f_end, GRR){


  ########################### EVENT WAITING TIMES ###############################

  # Waiting times of founder and children in pedigree
  waiting_times <- function(personID,df,dis_haz_used){
    # Find the person in the pedigree using their ID
    person = df[df$ID == personID, ]
    # Birth year
    birth = person$birth
    # Carrier status
    c = person$carrier

    ##### Death waiting time
    wait_death <- min(100 + birth,PWEXP::rpwexp(1, rate = unaff_death_haz,
                                         breakpoint = 1:99) + birth, a_end)
    ##### Cancer Onset waiting time
    wait_onset <- min(PWEXP::rpwexp(1, rate = dis_haz_used*GRR*c + dis_haz_used*(1-c),
                             breakpoint = 1:99) + birth, a_end)
    ##### Children waiting time
    # Generate number of children
    num_offspring <- stats::rnbinom(1,2,4/7)

    # Uniformly distribute children between ages of 18 and 45
    wait_offspring <- sort(stats::runif(num_offspring, min = 18, max = 45))
    if (length(wait_offspring) > 0){
      for (i in 1:length(wait_offspring)){
        # If waiting time after stop year, make waiting time equal to stop year
        wait_offspring[i] = min(wait_offspring[i] + birth, a_end)
      }
    }
    # Return the waiting times
    return(c(wait_death, wait_onset, wait_offspring))
  }


  # Waiting times for partner, which exclude children waiting times
  # A person's partner will have children at the person's waiting times
  waiting_partner <- function(personID,df,dis_haz_used){
    # Find the partner in the pedigree using their ID
    person = df[df$ID == personID, ]
    # Birth year
    birth = person$birth
    # Carrier status
    c = person$carrier

    # Find birth of all children who has the partner as a parent
    child_na = c(df[df$dadID == personID | df$mumID == personID, "birth"])
    # Ignore children with NA as their parents
    child_no = child_na[!is.na(child_na)]

    ##### Death waiting time
    # Death must be after the latest birth of their children
    wait_death <- min(100 + birth,PWEXP::rpwexp(1, rate = unaff_death_haz,
                                         breakpoint = 1:99) + max(child_no), a_end)
    ##### Cancer onset waiting time
    wait_onset <- min(PWEXP::rpwexp(1, rate = dis_haz_used*GRR*c + dis_haz_used*(1-c),
                             breakpoint = 1:99) + birth, a_end)
    # Return the waiting times
    return(c(wait_death, wait_onset))
  }

  ################################# EVENTS ######################################

  # Next event
  events <- function(personID,df,dis_haz_used){
    # Find the person in the pedigree using their ID
    person = df[df$ID == personID, ]
    # Carrier status
    c = person$carrier
    # Generate waiting times for the person
    w = waiting_times(personID,df,dis_haz_used)
    birth = person$birth

    # Person has events until censoring occurs
    continue = TRUE
    while (continue == TRUE){
      # Smallest waiting time is the next event to occur
      n = which.min(w)
      # If smallest waiting time is stop year of study, no more events for person
      if (w[n] == a_end){
        continue = FALSE
      }
      else{
        # If next event is death, stop any further events
        if (n == 1){
          person$death = w[n]
          continue = FALSE
        }
        # If next event is cancer onset, generate new death waiting time
        if (n == 2){
          # Set the onset year as the current waiting time
          person$onset = w[n]
          # Set affected status to 1
          person$affected = 1
          w[1] <- min(stats::rexp(1, rate = aff_death_haz) + person$onset,
                      a_end, w[1])
          # Remove cancer onset waiting time by setting it to the stop year
          w[n] <- a_end
        }
        # If next event is a child, generate a child
        if (n >= 3){
          # Generate child in data frame
          df = child(personID,w[n],df)
          # Remove waiting time by setting it to the stop year
          w[n] = a_end
          # Update person in case they have had a partner ID added
          person = df[df$ID == personID, ]
        }
      }
    }
    # ind = 0 indicates that person has finished their events
    person$ind = 0
    # Put updated person back in the data frame
    df[df$ID == personID, ] <- person
    return(df)
  }

  # Different events for partner as they do not need children waiting times
  # Their children waiting times are the same as the person they are paired with

  events_partner <- function(personID,df,dis_haz_used){
    # Find the partner in the pedigree using their ID
    person = df[df$ID == personID, ]
    # Carrier status
    c = person$carrier
    # Generate waiting times for the partner
    w = waiting_partner(personID,df,dis_haz_used)

    # Partner has events until censoring occurs
    continue = TRUE
    while (continue == TRUE){
      # Smallest waiting time is the next event to occur
      n = which.min(w)
      # If smallest waiting time is stop year of study, no more events for person
      if (w[n] == a_end){
        continue = FALSE
      }
      else{
        # If next event is death, stop any further events
        if (n == 1){
          person$death = w[n]
          continue = FALSE
        }
        # If next event is cancer onset, generate new death waiting time
        if (n == 2){
          # Set the onset year as the current waiting time
          person$onset = w[n]
          # Set affected status to 1
          person$affected = 1
          w[1] <- min(stats::rexp(1, rate = aff_death_haz) + person$onset, a_end, w[1])
          # Remove cancer onset waiting time by setting it to the stop year
          w[n] <- a_end
        }
      }
    }
    # ind = 0 indicates that the partner has finished their events
    person$ind = 0
    # Put updated partner back in the data frame
    df[df$ID == personID, ] <- person
    return(df)
  }


  ###################### GENERATE PARTNERS AND CHILDREN #########################

  ### Generate partner

  partner <- function(personID,df){
    # Find person from the pedigree using their ID
    person = df[df$ID == personID, ]
    # Uniformly generate the birth year of the partner
    # Partner must be at least 18 and within 5 years of the person's birth year
    partner_birth <- stats::runif(1, min = max(person$birth - 5, 18), max = person$birth + 5)
    # Assign partner to have opposite sex of person
    if (person$sex == 1){
      partner_sex = 0
    }
    else{
      partner_sex = 1
    }
    # Assign new partner an ID that is the next largest of the current people
    i = max(df$ID) + 1
    new_partner <- data.frame(FamID = j, ID = i, sex = partner_sex,
                              partner = person$ID, dadID = NA, mumID = NA,
                              affected = 0, carrier = 0, birth = partner_birth,
                              onset = NA, death = NA, proband = 0,
                              event_age = NA, start_age = NA, ind = 2)
    # Add new partner to the data frame
    df = dplyr::bind_rows(df, new_partner)
    # Give assigns the person's "partner" to be the partner's ID
    df[df$ID == personID, "partner"] <- i
    return(df)
  }

  ##### Generate children

  child <- function(personID,time,df){
    # Find the person from the pedigree who is going to have a child
    person = df[df$ID == personID, ]
    # Generate a partner if the person does not have one
    if (is.na(person$partner)){
      df = partner(personID, df)
    }
    # Find the person again, now with a partner ID
    person = df[df$ID == personID, ]
    partnerID = person$partner
    # Find the partner in the pedigree using their ID
    partner = df[df$ID == partnerID, ]

    # Assign parents of the child based on person and partner sex
    if (person$sex == 1){
      dad = person
      mum = partner
    }else{
      dad = partner
      mum = person
    }
    # If one parent is a carrier, there is 50% probability that child is a carrier
    # Assign carrier status to child
    if (person$carrier == 1){
      is_carrier <- stats::rbinom(1,1,0.5)
    }else{
      is_carrier <- 0
    }
    # Randomise sex of child, 50% probability of each
    rand_sex <- stats::rbinom(1,1,0.5)
    # Give an ID to the new child which is one higher than the current largest
    i = max(df$ID) + 1
    # Create a data frame for the new child
    new_child <- data.frame(FamID = j, ID = i, sex = rand_sex, partner = NA,
                            dadID = dad$ID, mumID = mum$ID, affected = 0,
                            carrier = is_carrier, birth = time, onset = NA,
                            death = NA, proband = 0, event_age = NA,
                            start_age = NA, ind = 1)
    # Add the new child to the current pedigree
    df = dplyr::bind_rows(df, new_child)
    # Return the modified pedigree
    return(df)
  }





  ############################ GENERATE PEDIGREE ################################

  run_events <- function(dis_haz_used){
    # Uniformly generate the birth year of the founder within interval
    founder_birth <- stats::runif(1, min = f_start, max = f_end)
    # Randomise sex of founder
    rand_sex <- stats::rbinom(1,1,0.5)
    # Create a data frame of the founder, who is always a carrier
    # They have no parents and are initially unaffected (cancer free)
    df <- data.frame(FamID = j, ID = 1, sex = rand_sex, partner = NA, dadID = NA,
                     mumID = NA, affected = 0, carrier = 1, birth = founder_birth,
                     onset = NA, death = NA, proband = 0, event_age = NA,
                     start_age = NA, ind = 1)
    fin = FALSE
    # Check for people who still need to have events until all are completed
    while (fin == FALSE){
      for (i in df$ID){
        person = df[df$ID == i, ]
        # ind == 1 means this person has not had waiting times generated
        if (person$ind == 1){
          df = events(i,df,dis_haz_used)
        }
        # ind == 2 means this person is a partner needing partner waiting times
        else if (person$ind == 2){
          df = events_partner(i,df,dis_haz_used)
        }
      }
      # Once ind == 0, all people have finished their events so pedigree finished
      fin = all(df$ind == 0)
    }
    # Find carriers who had cancer onset within ascertainment interval
    asc_interval_df <- df |> dplyr::filter(.data$onset >= a_start,
        .data$onset <= a_end, .data$carrier == 1)
    # If there is at least one, assign one of these carriers as the proband
    if (nrow(asc_interval_df) > 0) {
      Proband <- asc_interval_df[sample(nrow(asc_interval_df), 1), ]
      Proband$proband = 1
      df[df$ID == Proband$ID, ] <- Proband
    }
    # Return the completed pedigree
    return(df)
  }



  ######################### ASCERTAINMENT CONDITIONS ############################

  ##### Find a pedigree that satisfies ascertainment requirements

  ascertainment <- function(aff,frailty,frail_var){
    asc = FALSE
    # Keep a record of the number of pedigrees until one ascertained
    i=0
    # Repeat until a pedigree satisfies the ascertainment requirements
    while (asc == FALSE){
      # If frailty is included, use disease hazards multiplied by the frailty term
      if (frailty == 1){
        # Frailty term
        frail <- stats::rgamma(1,1/frail_var,1/frail_var)
        dis_haz_used = dis_haz * frail
      }
      else{
        dis_haz_used = dis_haz
      }
      # Obtain a pedigree
      df = run_events(dis_haz_used)
      # Find the number of affected individuals
      num_aff = sum(df$affected == 1 & df$carrier == 1)
      # If a proband exists, the sum will be 1
      proband_exist = sum(df["proband"])
      # Ascertainment requirements are:
      # 1. Proband exists (cancer onset within ascertainment interval)
      # 2. The number of affected carriers is at least the specified amount
      if (num_aff > aff-1 & proband_exist == 1){
        # Pedigree is ascertained
        asc = TRUE
      }
      else{
        # Pedigree is not ascertained
        asc = FALSE
        # Increase number of pedigrees currently tested for ascertainment
        i=i+1
      }
    }
    return(df)
  }

  pedigree <- ascertainment(aff,frailty,frail_var)
  return(pedigree)

}
