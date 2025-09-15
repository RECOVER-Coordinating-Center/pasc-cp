# Vector of additional packages to load before executing the request
config_append('extra_packages', c('tidyr', 'stringr'))

#' Execute the request
#'
#' This function presumes the environment has been set up, and executes the
#' steps of the request.
#'
#' In addition to performing queries and analyses, the execution path in this
#' function should include periodic progress messages to the user, and logging
#' of intermediate totals and timing data through [append_sum()].
#'
#' @return The return value is dependent on the content of the request, but is
#'   typically a structure pointing to some or all of the retrieved data or
#'   analysis results.  The value is not used by the framework itself.
#' @md
.run  <- function() {

  message('Starting execution with framework version ',
          config('framework_version'))

  init_sum(cohort = 'Start', persons = 0)

  rslt <- list()
  # Developing an attrition table to show breakdown of patients meeting inclusion criteria at each step
  append_sum(cohort='Persons in the RECOVER database',
             persons=distinct_ct(cdm_tbl('person')))
  message('Determining SARS-CoV-2 Infection Status')
  # Looking for COVID-19-related indicators and flagging whether patient conclusively had COVID-19 indication
  # positive PCR test
  rslt$cohort_obs_pcr <- classify_observation(obs_ids=c(2000001530L), 
                                          val_ids=c(9191L, 2000001526L), 
                                          rule_name='positive_viral_pcr',
                                          date_bool=T,
                                          fact_concl='conclusive')
  # positive antigen test
  rslt$cohort_obs_ant <- classify_observation(obs_ids=c(2000001529L), 
                                          val_ids=c(9191L, 2000001526L), 
                                          rule_name='positive_viral_antigen',
                                          date_bool=T,
                                          fact_concl='conclusive')
  # specific COVID-19 diagnosis
  rslt$cohort_obs_specdx <- classify_observation(obs_ids=c(2000001527L), 
                                                 val_ids=c(2000001525L), 
                                                 rule_name='specific_covid19_dx',
                                                 date_bool=T,
                                                 fact_concl='conclusive')
  # COVID-19 complication code
  rslt$cohort_obs_compdx <- classify_observation(obs_ids=c(2000001527L), 
                                                 val_ids=c(2000001523L), 
                                                 rule_name='complication_covid19_dx',
                                                 date_bool=F,
                                                 fact_concl='conclusive')
  # PASC dx code
  rslt$cohort_obs_pasc <- classify_observation(obs_ids=c(2000001527L), 
                                                   val_ids=c(2000001520L), 
                                                   rule_name='pasc_misc',
                                                   date_bool=F,
                                                   fact_concl='conclusive')
  # MISC dx code
  rslt$cohort_obs_misc <- classify_observation(obs_ids=c(2000001527L), 
                                               val_ids=c(703578L), 
                                               rule_name='pasc_misc',
                                               date_bool=F,
                                               fact_concl='conclusive')
  # sequelae of other specified infectious disease (B94.8) dx code
  rslt$cohort_seq_unspec <- classify_observation(obs_ids=c(2000001527L), 
                                                 val_ids=c(2000001533L), 
                                                 rule_name='seq_unspec',
                                                 date_bool=F,
                                                 fact_concl='conclusive')
  # note that seq_unspec should only count as a qualifying event if prior to Oct 1 2021 (based on code usage for COVID-19 prior to availability of other codes)
  # if seq_unspec occurs after Oct 1 2021, another indication of COVID19 infection would be required prior
  rslt$cohort_seq_unspec_qual <- rslt$cohort_seq_unspec %>% filter(!is.na(observation_date)&observation_date<'2021-10-01')
  rslt$cohort_seq_unspec_unqual <- rslt$cohort_seq_unspec %>% filter(!is.na(observation_date)&observation_date>='2021-10-01')
  # history of COVID-19 code
  rslt$cohort_obs_histdx <- classify_observation(obs_ids=c(2000001527L), 
                                                   val_ids=c(2000001522L), 
                                                   rule_name='history_covid19_dx',
                                                   date_bool=F,
                                                 fact_concl='probable')
  # positive nucleocapsid serology test
  rslt$cohort_serology <- get_serology_positives()%>%
    mutate(rule='positive_nucleocapsid_serology',
           fact_infection='conclusive',
           date_known=FALSE)
  
  message('Build table of SARS-CoV-2 Infection observations')
  rslt$cohort_infections <- rslt$cohort_obs_pcr %>%
    dplyr::union_all(rslt$cohort_obs_ant) %>%
    dplyr::union_all(rslt$cohort_obs_specdx) %>%
    dplyr::union_all(rslt$cohort_obs_compdx) %>%
    dplyr::union_all(rslt$cohort_obs_pasc) %>%
    dplyr::union_all(rslt$cohort_obs_misc) %>%
    dplyr::union_all(rslt$cohort_seq_unspec_qual) %>%
    dplyr::union_all(rslt$cohort_obs_histdx) %>%
    dplyr::union_all(rslt$cohort_serology) %>%
    compute_new()
  
  append_sum(cohort='Patients with evidence of SARS-CoV-2 Infection',
             persons=distinct_ct(rslt$cohort_infections))
  
  # Determining index date based on either COVID-19 infection date or proxy date if patient had MISC/PASC dx with no initial COVID-19 dx
  message('Determine the qualifying event for the cohort entry')
  rslt$cohort_inf <- determine_index_date(rslt$cohort_infections)
  
  
  append_sum(cohort='Patients under 21 years of age on date of SARS-CoV-2 infection',
             persons=distinct_ct(rslt$cohort_inf))
  
  message('Restrict to patients with sufficient followup')
  rslt$cohort_visit_count <- count_visits(ce_tbl=rslt$cohort_inf,
                                          days_start=28L,
                                          days_end=179L) %>%
    filter(num_visits_window>=2L&
             num_visits_window_inperson>=1L)%>%
    compute_new()
  
  output_tbl(rslt$cohort_visit_count,
             name='cohort_visit_count')
  
  append_sum(cohort='Patients with sufficient follow-up post-SARS-CoV-2 infection',
             persons=distinct_ct(results_tbl('cohort_visit_count')))
  
  message('Determine whether patient has conclusive PASC')
  # PASC dx
  rslt$cohort_pasc_summary <- rslt$cohort_obs_pasc %>%
    group_by(person_id)%>%
    # find number of diagnoses before limiting to a certain index date
    summarise(num_pasc_dx_date=n_distinct(observation_date))%>%
    ungroup() %>%
    inner_join(rslt$cohort_obs_pasc, by = 'person_id')%>%
    mutate(ce_date=observation_date-28L)%>%
    # requiring that patient have sufficient follow up in relation to a qualifying diagnosis
    inner_join(select(results_tbl('cohort_visit_count'), c(person_id, ce_date)), by = c('person_id', 'ce_date'))%>%
    group_by(person_id, num_pasc_dx_date)%>%
    summarise(min_pasc_dx_date=min(observation_date))%>%
    ungroup() %>%
    mutate(ce_date=min_pasc_dx_date-28L)%>%
    compute_new()
  #  MISC dx
  rslt$cohort_misc_summary <- rslt$cohort_obs_misc %>%
    group_by(person_id)%>%
    # find number of diagnoses before limiting to a certain index date
    summarise(num_misc_dx_date=n_distinct(observation_date))%>%
    ungroup() %>%
    inner_join(rslt$cohort_obs_misc, by = 'person_id')%>%
    mutate(ce_date=observation_date-28L)%>%
    # requiring that patient have sufficient follow up in relation to a qualifying diagnosis
    inner_join(select(results_tbl('cohort_visit_count'), c(person_id, ce_date)), by = c('person_id', 'ce_date'))%>%
    group_by(person_id, num_misc_dx_date)%>%
    summarise(min_misc_dx_date=min(observation_date))%>%
    ungroup() %>%
    mutate(ce_date=min_misc_dx_date-28L)%>%
    compute_new()
  
  # B94.8 prior to Oct 1, 2021
  rslt$cohort_su_summary_pre<-rslt$cohort_seq_unspec_qual %>%
    group_by(person_id)%>%
    # find number of diagnoses before limiting to a certain index date
    summarise(num_su_pre_dx_date=n_distinct(observation_date))%>%
    ungroup()%>%
    inner_join(rslt$cohort_seq_unspec_qual, by = 'person_id')%>%
    mutate(ce_date=observation_date-28L)%>%
    # requiring that patient have sufficient follow up in relation to a qualifying diagnosis
    inner_join(select(results_tbl('cohort_visit_count'), c(person_id, ce_date)), by = c('person_id', 'ce_date'))%>%
    group_by(person_id, num_su_pre_dx_date)%>%
    summarise(min_su_pre_dx_date=min(observation_date))%>%
    ungroup()%>%
    mutate(ce_date=min_su_pre_dx_date-28L)%>%
    compute_new()
  # B94.8 after Oct 1, 2021
  # finding all patients who had B94.8 after Oct 1, 2021 and also had prior covid19 infection 
  rslt$cohort_su_summary_post <- rslt$cohort_seq_unspec_unqual %>% rename(b948_date=observation_date)%>%
    distinct(person_id, b948_date)%>%
    # making sure patient had COVID19 prior to B94.8 codes
    inner_join(rslt$cohort_inf %>% filter(rule!='seq_unspec'), by = 'person_id')%>%
    filter(ce_date<b948_date)%>%
    # making sure patient had sufficient follow up after the infection date
    inner_join(select(results_tbl('cohort_visit_count'), c(person_id, ce_date)), by = c('person_id', 'ce_date'))%>%
    group_by(person_id, ce_date) %>%
    summarise(num_su_dx_date=n_distinct(b948_date),
              min_su_dx_date=min(b948_date))%>%
    ungroup() %>%
    group_by(person_id)%>%
    # want to flag them as conclusive if this ever happens
    filter(num_su_dx_date==max(num_su_dx_date))%>%
    filter(ce_date==min(ce_date))%>%
    ungroup()%>%
    distinct()%>%
    compute_new()
  rslt$cohort_su_concl <-rslt$cohort_su_summary_pre %>%
    filter(num_su_pre_dx_date>=2)%>%
    dplyr::union(rslt$cohort_su_summary_post%>%filter(num_su_dx_date>=2))%>%
    group_by(person_id)%>%
    filter(ce_date==min(ce_date))%>%
    ungroup()
  rslt$cohort_su_prob <- rslt$cohort_su_summary_pre%>%filter(num_su_pre_dx_date==1)%>%
    dplyr::full_join(rslt$cohort_su_summary_post%>%filter(num_su_dx_date==1), by = c('person_id', 'ce_date'))%>%
    group_by(person_id)%>%
    filter(ce_date==min(ce_date))%>%
    ungroup()
  
  message('Build the list of conclusive PASC')
  rslt$cohort_conclusive <- rslt$cohort_pasc_summary %>% filter(num_pasc_dx_date>=2)%>%
    rename(ce_date_pasc=ce_date)%>%
    full_join(filter(rslt$cohort_misc_summary,num_misc_dx_date>=2)%>%rename(ce_date_misc=ce_date), by = 'person_id')%>%
    full_join(rslt$cohort_su_concl%>%rename(ce_date_su=ce_date), by = 'person_id')%>%
    collect()%>%
    rowwise()%>%
    # find the earliest of any of the conclusive evidence, but keep all around
    mutate(ce_date=min(ce_date_pasc, ce_date_misc, ce_date_su, na.rm = TRUE))%>%
    mutate(conclusive=TRUE)%>%
    copy_to_new(df=., name='conclusives', temporary = TRUE)
  
  message('Build the list of probable PASC based on single dx codes')
  rslt$cohort_prob <- rslt$cohort_pasc_summary %>% filter(num_pasc_dx_date==1)%>%rename(ce_date_pasc=ce_date)%>%
    full_join(filter(rslt$cohort_misc_summary,num_misc_dx_date==1)%>%rename(ce_date_misc=ce_date), by = 'person_id')%>%
    full_join(rslt$cohort_su_prob%>%rename(ce_date_su=ce_date), by = 'person_id')%>%
    anti_join(select(rslt$cohort_conclusive, person_id), by = 'person_id')%>% 
    collect()%>%
    rowwise()%>%
    # find the earliest of any of the conclusive evidence, but keep all around
    mutate(ce_date=min(ce_date_pasc, ce_date_misc, ce_date_su, na.rm = TRUE))%>%
    mutate(probable=TRUE)%>%
    copy_to_new(df=., name='probables', temporary = TRUE)
  
  message('Limit to the CED relevant to classification')
  # finding all potential CED
  rslt$cohort_id <- bind_rows((rslt$cohort_conclusive%>%rename(ce_date_conclusive=ce_date)%>%collect()), 
                              (rslt$rslt_cohort_prob%>%rename(ce_date_prob=ce_date)%>%collect()))%>%
    inner_join(select(results_tbl('cohort_visit_count'),c(person_id,site,ce_date)), by = 'person_id')%>%
    mutate(ce_date=coalesce(ce_date_conclusive, ce_date_prob, ce_date))%>% # if conclusive, use conclusive date as ce_date, otherwise probable date, otherwise other
    compute_new()
  
  message('Examining clusters')
  rslt$cohort_cluster_occs_all <- find_condition_source_occurrence(cohort=rslt$cohort_id,
                                                                   condition_tbl=cdm_tbl('condition_occurrence'),
                                                                   condition_codes=load_codeset('cluster_master_pasc', col_types = 'icccc')) %>%
    compute_new()
  
  # remove respiratory and fever concurrent with respiratory infection
  rslt$cohort_cluster_noncon <- exclude_cooccur_infx(cluster_tbl=rslt$cohort_cluster_occs_all,
                                                     condition_tbl=cdm_tbl('condition_occurrence'),
                                                     condition_codes=load_codeset('noncovid_respiratory_infections', col_types='icccccccccc'),
                                                     cluster_names=c('respiratory_signs_and_sx', 'fever_and_chills'))
  
  
  # count number and min/max dates within each phase
  rslt$cohort_cluster_summary <- characterize_clusters(ce_tbl=rslt$cohort_id,
                                                       condition_tbl=rslt$cohort_cluster_noncon) %>%
    compute_new()
  
  # limit to just post-acute and find ones with sufficient separation (note non-qualifying ones still in the table but cluster_dx_grx=0 for those)
  rslt$cohort_cluster_rolled <- rollup_cluster(cluster_summary_tbl=rslt$cohort_cluster_summary,
                                               day_diff=28L,
                                               min_n_dx=2L,
                                               phase_filt=c('post_acute'))%>%
    compute_new()
  
  rslt$cohort_conditions_towash <- find_condition_towash(ce_tbl=rslt$cohort_id,
                                                         condition_tbl=rslt$cohort_cluster_noncon,
                                                         washout_tbl=load_codeset('cluster_master_pasc_washout', col_types='icccci'))
  rslt$cohort_cluster_new <- wash_condition(conds_postacute=rslt$cohort_cluster_rolled,
                                            conds_towash=rslt$cohort_conditions_towash)
  
  # note that now only including their counts if they meet the criteria for the cluster conditions
  rslt$cohort_cluster_new_qual<-determine_cluster_qual(rslt$cohort_cluster_new)
  rslt$cohort_to_classify <- rslt$cohort_id %>%
    left_join(rslt$cohort_cluster_new_qual, by = c('person_id','ce_date'))
  message('Assign PASC certainty')
  rslt$levels_assigned<-assign_pasc_level(rslt$cohort_to_classify)

  # finding evidence of infection on index date
  # patient could have more than one, so store in separate table to keep flags table 1:patient
  rslt$cohort_index_occurrences <- rslt$cohort_inf%>%
    select(person_id, ce_date, rule, fact_infection, date_known) %>%
    inner_join(rslt$levels_assigned%>%select(person_id, ce_date),by=c('person_id', 'ce_date'))%>%
    distinct()%>%
    output_tbl('cohort_index_occurrences')
  
  output_sum()

  invisible(rslt)

}
