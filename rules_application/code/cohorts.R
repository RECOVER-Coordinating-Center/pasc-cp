#'
#' This file contains functions to identify cohorts in a study.  Its
#' contents, as well as the contents of any R files whose name begins
#' with "cohort_", will be automatically sourced when the request is
#' executed.
#'
#' For simpler requests, it will often be possible to do all cohort
#' manipulation in functions within this file.  For more complex
#' requests, it may be more readable to separate different cohorts or
#' operations on cohorts into logical groups in different files.
#'

#' Function to classify observation types in an observation table
#' @param obs_tbl table structured as the observation table
#' @param cohort table with person_ids
#' @param obs_ids list of observation_concept_ids
#' @param val_ids list of value_as_concept_ids
#' @param rule_name string to classify observation as
#' @param date_bool boolean indicating whether observation date can be relied on for the occurrence
#' @return a table with the columns in obs_tbl + the columns:
#'          `rule` name of the rule
#'          `date_known` whether or not can use the observation_date as the date of event
#'          
classify_observation <- function(obs_tbl=cdm_tbl('observation_derivation_recover'), 
                                 cohort_tbl=cdm_tbl('person'),
                                 obs_ids, val_ids, rule_name, date_bool, fact_concl){
  
cohort_tbl %>% distinct(person_id) %>%
    inner_join(obs_tbl, by = 'person_id')%>%
    filter(observation_concept_id %in% !!obs_ids &
             value_as_concept_id %in% !!val_ids) %>%
    mutate(rule=rule_name,
           date_known=date_bool,
           fact_infection=fact_concl)
}

#' Function to determine serology positives, from the treescan code
#' @param obs_der_tbl table in the format of observation_derivation_recover
get_serology_positives<-function(obs_der_tbl=cdm_tbl("observation_derivation_recover") %>% filter(observation_date<as.Date("2022-05-01"))){
  ## Positive/Negative Serology Results
  serology_tbl<-obs_der_tbl %>% 
    filter(observation_concept_id==2000001528L)%>%
    mutate(test_result=case_when(
      value_as_concept_id %in% c(9191L, 2000001526L)~"positive",
      value_as_concept_id==9189L~"negative",
      TRUE~"unknown"
    )) %>%
    filter(test_result!="unknown") %>%
    compute_new(indices=c("person_id", "visit_occurrence_id", "observation_source_concept_id"))
  
  
  measurement_tbl<-cdm_tbl('measurement_labs') %>%
    filter(measurement_date >= as.Date('2020-03-01')) %>%
    mutate(vsv_lower = tolower(value_source_value), unit_lower=tolower(unit_source_value))
  
  serology_meas<-serology_tbl %>% 
    select(-value_source_value, -unit_source_value, -unit_concept_id) %>%
    inner_join(measurement_tbl %>% 
                 select(measurement_id, measurement_date, vsv_lower, measurement_source_value, measurement_concept_id),
               by=c("observation_source_concept_id"="measurement_id")) %>%
    compute_new(indices=c("person_id", "visit_occurrence_id", "observation_source_concept_id"))
  
  serology<-serology_meas %>% left_join(select(cdm_tbl('person'),person_id,birth_date),by='person_id') %>%
    mutate(measurement_age = (measurement_date - birth_date)/365.25) %>%
    filter(measurement_age >= 0, measurement_age < 21) %>%
    mutate(age_range = case_when(measurement_age < 0.5 ~ "<6mo",
                                 measurement_age >= 0.5 & measurement_age < 5 ~ "6mo-5",
                                 measurement_age >= 5 & measurement_age < 12 ~ "5-11",
                                 measurement_age >= 12 & measurement_age < 16 ~ "12-15",
                                 measurement_age >= 16 ~ "16+"))%>% 
    mutate(sv_lower = tolower(measurement_source_value)) %>%
    mutate(measurement_concept = case_when(str_detect(sv_lower, 'nucleocapsid') ~ 2000001501,
                                           str_detect(sv_lower, 'spike') ~ 2000001502,
                                           TRUE ~ measurement_concept_id))%>% 
    mutate(serology_type = case_when(measurement_concept == 2000001501L ~ 'IgG N protein',
                                     measurement_concept %in% c(723474L, 706177L, 706181L) ~ 'IgG undifferentiated',
                                     measurement_concept %in% c(723475L, 706178L) ~ 'IgM',
                                     measurement_concept == 723473L ~ 'IgA',
                                     measurement_concept %in% c(586515L,586522L,723480L) ~ 'Ab undifferentiated',
                                     measurement_concept %in% c(36032309L, 36031956L, 36031969L) ~ 'stimulated gamma interferon release',
                                     measurement_concept %in% c(2000001502L, 1619027L, 36031734L) ~ 'S protein or RBD')) %>%
    filter(serology_type != 'stimulated gamma interferon release')
  
  serology_filtered<-serology %>%
    mutate(include_flag=case_when(
      serology_type %in% c('IgG N protein', 'IgM')~1,
      age_range=='<6mo'~1,
      age_range=='6mo-5' & measurement_date<as.Date('2022-06-18')~1,
      age_range=='16+' & measurement_date<as.Date('2020-12-12')~1,
      age_range=='12-15' & measurement_date<as.Date('2021-05-12') ~1,
      age_range=='5-11' & measurement_date<as.Date('2021-11-02')~1,
      TRUE~0
    )) %>%
    filter(include_flag==1) %>% distinct(observation_id) %>% compute_new(index="observation_id")
  
  
  serology_tbl %>% 
    filter(test_result=="positive") %>%
    select(-test_result) %>% 
    inner_join(serology_filtered, by="observation_id") %>%
    compute_new(indices=c("person_id", "visit_occurrence_id", "observation_source_concept_id")) %>%
    return()
}



#' Function to limit age on cohort entry date below a specified value
#' @param ce_tbl table with at least the columns `person_id`, `ce_date`
#' @param age_upper_years age, in years, of the maximum cohort entry years
#' @return table with the columns in the original `ce_tbl` + the columns `ce_age_days` and `ce_age_years`, limited to only patients under the age specified in the `age_upper_years` parameter
limit_ce_age <- function(ce_tbl,
                         age_upper_years){
  ce_tbl %>%
    select(person_id, ce_date) %>%
    distinct() %>%
    inner_join(select(cdm_tbl('person'), c(person_id, birth_date)), by = 'person_id') %>%
    mutate(ce_age_days=as.integer(ce_date-birth_date),
           ce_age_years=ce_age_days/365.25) %>%
    filter(ce_age_years>-28, ce_age_years<age_upper_years)%>%
    mutate(ce_age_cat=case_when(ce_age_years<1~"<1",
                                ce_age_years<5~"1-4",
                                ce_age_years<10~"5-9",
                                ce_age_years<16~"10-15",
                                TRUE~"16-21"))%>%
    inner_join(ce_tbl, by = c('person_id', 'ce_date'))
}

#' Function to count number of visits 
#'    before a cohort entry date,
#'    after a cohort entry date,
#'    in a specified time frame around a cohort entry date
#' @param ce_tbl table with the cols `person_id` and `ce_date`
#' @param visit_tbl table with visit_occurrences
#' @param days_start number of days after `ce_date` to start counting visits
#' @param days_end number of days after `ce_date` to finish counting visits
#' @return table with a count of the number of visits patient had within the `days_start` and `days_end` time frame
count_visits <- function(ce_tbl,
                         visit_tbl=cdm_tbl('visit_occurrence'),
                         days_start,
                         days_end){
  visits_distinct <- visit_tbl %>%
    distinct(person_id, site, visit_start_date, visit_concept_id)%>%
    filter(!is.na(visit_start_date))%>%# can't make a determination if there isn't a start date
    # figure out if any of the visits on the visit_start_date were in person or admitted and if so, flag as in person or admitted on that date
    mutate(in_person_flag=case_when(visit_concept_id %in%c(9201L, 9202L, 9203L, 42898160L, 44814710L,
                                                           2000000048L, 2000000088L, 581399L)~1L,
                                    TRUE~0L),
           admitted_flag=case_when(visit_concept_id%in% c(9201L, 2000000048L)~1L,
                                   TRUE~0L))%>%
    group_by(person_id, site, visit_start_date)%>%
    summarise(in_person_any=max(in_person_flag),
              admitted_any=max(admitted_flag))%>%
    ungroup()
  
  cohort_visits <- ce_tbl %>% distinct(ce_date, person_id) %>%
    inner_join(visits_distinct, by = 'person_id') %>%
    mutate(prior_visit=case_when(visit_start_date<ce_date~1L,
                                 TRUE~0L),
           post_visit=case_when(visit_start_date>=ce_date~1L,
                                TRUE~0L),
           visit_within=case_when(visit_start_date>=ce_date+days_start &
                                    visit_start_date<=ce_date+days_end~1L,
                                  TRUE~0L),
           in_person=case_when(in_person_any==1L~TRUE,
                               TRUE~FALSE),
           admitted=case_when(admitted_any==1L~TRUE,
                              TRUE~FALSE),
           in_person_flag=case_when(visit_within==1L&in_person~1L,
                                    TRUE~0L),
           admitted_flag=case_when(visit_within==1L&admitted~1L,
                                   TRUE~0L)) %>%
    group_by(person_id, site, ce_date) %>%
    summarise(num_visits_pre=sum(prior_visit),
              num_visits_post=sum(post_visit),
              num_visits_window=as.integer(sum(visit_within)),
              num_visits_window_inperson=as.integer(sum(in_person_flag)),
              num_visits_window_ip=as.integer(sum(admitted_flag)))%>%
    ungroup()
  
}
#' Function to find hits for a condition codeset in the condition_source_concept_id field
#' @param cohort table with at least a `person_id` column
#' @param condition_tbl table with condition_occurrences with a `condition_source_concept_id` field
#' @param condition_codes table with at least a `concept_id` field with the concepts of interest
#' @return a row for each occurrence in the `condition_tbl` for the `cohort` of any of the codes in `condition_codes`, 
#'              retaining the columns in the codeset (condition_codes)
find_condition_source_occurrence <- function(cohort,
                                             condition_tbl,
                                             condition_codes) {
  cohort%>% select(person_id) %>% distinct() %>%
    inner_join(condition_tbl, by = 'person_id') %>%
    inner_join(condition_codes, by = c('condition_source_concept_id'='concept_id'))
}
#' Function to classify condition occurrences grouped by cluster classifications
#' @param ce_tbl table with at least person_id and ce_date
#' @param condition_tbl table with at least person_id, condition_start_date, cluster
#' @return table with person_id, cluster, phase, and within those phase combinations:
#'           min_dx_date: first condition_start_date in the phase
#'           max_dx_date: final condition_start_date in the phase
#'           num_dx: number of condition_start_date for any of the diagnoses within the phase
#'           max_days_diff: maximum number of days between diagnoses during the phase
characterize_clusters <- function(ce_tbl,
                                  condition_tbl){
  ce_tbl %>% select(person_id, ce_date) %>% distinct() %>%
    inner_join(condition_tbl, by = 'person_id') %>%
    mutate(days_to_dx=as.integer(condition_start_date-ce_date),
           phase=case_when(days_to_dx< -28~'prior',
                           days_to_dx<= 27L~'acute',
                           days_to_dx<= 179~'post_acute',
                           days_to_dx<= 657~'chronic',
                           TRUE~'post_chronic'))%>%
    group_by(person_id, cluster, phase)%>%
    summarise(min_dx_date=min(condition_start_date, na.rm=T),
              max_dx_date=max(condition_start_date, na.rm=T),
              num_dx=as.integer(n_distinct(condition_start_date)))%>%
    ungroup()%>%
    mutate(max_days_diff=as.integer(max_dx_date-min_dx_date))
}
#' Function to summarize cluster codes per person, based on whether patient had 
#'              a certain number of diagnoses
#'              with specified days separation
#'              during specified phase
#' @param cluster_summary_tbl table with at least the cols: 
#'                            `person_id`, `cluster`, `phase`, `min_dx_date`, `num_dx`, `max_days_diff`
#' @param day_diff minimum number of days required between diagnoses
#' @param min_n_dx minimum number of diagnoses required
#' @param phase_filt vector containing strings for acceptable phase in the `phase` column
#' @return tbl with the cols:
#'           person_id
#'           num_cluster_dx: number of diagnosis dx days for any of the clusters for which patient met criteria
#'           num_clusters: number of distinct clusters for which patient met the criteria during any of the phases specified
#'           max_days_diff: maximum number of days between diagnoses within the cluster
#'           min_dx_date: minimum diagnosis date for one of the cluster diagnoses within the phase
#'           cluster_dx_grx: indicator for whether or not greater than the specified number of cluster diagnoses with specified separation
rollup_cluster <- function(cluster_summary_tbl,
                           day_diff,
                           min_n_dx,
                           phase_filt){
  cluster_summary_tbl %>%
    filter(phase %in% !!phase_filt) %>%
    mutate(cluster_dx_grx=case_when(num_dx>=min_n_dx&
             max_days_diff>=day_diff~1L,
           TRUE~0L)) %>%
    group_by(person_id) %>%
    summarise(num_cluster_dx=as.integer(sum(num_dx)),
              num_clusters=as.integer(n_distinct(cluster)),
              cluster_dx_grx=max(cluster_dx_grx),
              max_days_diff=max(max_days_diff),
              min_dx_date=min(min_dx_date))%>%
    mutate(cluster_dx_grx=case_when(cluster_dx_grx==1L~TRUE,
                                    TRUE~FALSE))%>%
    ungroup()
}


#' Function to assign conclusivity of PASC
#' @param flag_tbl table with PASC-related flags
#' @return table with the original columns in `flag_tbl` + a `pasc_level` column indicating whether PASC is
#'                                             conclusive, probable, possible, or no evidence
assign_pasc_level <- function(index_tbl){
  index_tbl %>%
    mutate(rule_num=case_when(conclusive~1L,
                              probable|cluster_dx_grx~2L,
                              TRUE~3L)) %>%
    group_by(person_id) %>%
    filter(rule_num==min(rule_num))%>%
    filter(ce_date==min(ce_date))%>%
    ungroup()%>%
    mutate(pasc_level=case_when(rule_num==1L~'conclusive',
                                rule_num==2L~'probable',
                                rule_num==3L~'no_evidence'))%>%
    distinct()
}
#' Function to determine the event that could count as the index date 
#' so that a follow up window can be established to assess for evidence of PASC,
#' taking into consideration that patient should be <21 years of age on index date
#' NOTE a patient could have more than one window of assessment
#' @param rules_tbl table with at least person_id | date_known (boolean) | observation_date
#' @return table with ce_date column containing a potentially eligible index date
determine_index_date <- function(rules_tbl){
  rules_tbl %>%
    mutate(ce_date_min=case_when(date_known~observation_date,
                                 TRUE~observation_date-28L))%>%
    inner_join(select(cdm_tbl('person'),c(person_id,birth_date)), by = 'person_id')%>%
    mutate(ce_date=case_when(ce_date_min<birth_date~birth_date,
                             TRUE~ce_date_min))%>%
    filter(ce_date>='2020-03-01')%>%
    select(-c(birth_date, ce_date_min))%>%
    limit_ce_age(., age_upper_years=21)
}
#' Exclude co-occurring infections from cluster occurrences
#'
#' @param cluster_tbl table with all cluster occurrences incl the cols person_id | cluster | condition_start_date
#' @param condition_tbl table containing all condition_occurrences
#' @param condition_codes table containing concept_ids for the infections that you want to look for
#' @param cluster_names list of the clusters that you want to remove if in the presence of an infection named in `condition_codes`
#'
#' @return new cluster occurrence table with occurrences removed
exclude_cooccur_infx <- function(cluster_tbl,
                                 condition_tbl=cdm_tbl('condition_occurrence'),
                                 condition_codes,
                                 cluster_names){
  pats <- cluster_tbl %>% distinct(person_id)
  
  # finding distinct occurrences of cluster of interest per patient
  cluster_tbl_dt <- cluster_tbl %>% 
    filter(cluster %in% !!cluster_names)%>%
    distinct(person_id, cluster, condition_start_date) %>%
    mutate(obs_start_date=condition_start_date-14L,
           obs_end_date=condition_start_date+14L)%>%
    compute_new()
  
  # find conditions and dates for removal
  coocurr_tbl_dt <- condition_tbl %>% 
    inner_join(pats, by = 'person_id')%>%
    inner_join(select(condition_codes, concept_id), by = c('condition_concept_id'='concept_id'))%>%
    distinct(person_id, condition_start_date)%>%
    rename(inf_date=condition_start_date)%>%
    compute_new()
  
  occs_to_remove <- cluster_tbl_dt %>%
    inner_join(coocurr_tbl_dt, by = 'person_id')%>%
    mutate(inf_within=case_when(inf_date>=obs_start_date&
                                  inf_date<=obs_end_date~1L,
                                TRUE~0L))%>%
    filter(inf_within==1L)
  
  # removing any occurrences that would "disqualify" a cluster occurrence
  cluster_tbl %>%
    anti_join(occs_to_remove, by = c('person_id', 'condition_start_date', 'cluster'))
}
#' Function to generate a table containing the conditions that should be washed out for the specified time frame
#' @param ce_tbl table with one row per person with the columns person_id | ce_date
#' @param condition_tbl table with cluster conditions with one row per condition_occurrence
#' @return table with conditions occurring in the washout period prior to the COVID qualifying occurrence
find_condition_towash <- function(ce_tbl,
                                  condition_tbl,
                                  washout_tbl){
  ce_tbl %>% select(person_id, ce_date) %>%
    inner_join(condition_tbl, by = 'person_id')%>%
    inner_join(select(washout_tbl, c(concept_id, washout_days)),by=c('condition_source_concept_id'='concept_id'))%>%
    filter(ce_date>condition_start_date &
             as.integer(ce_date-condition_start_date)<=washout_days)
}
#' Function to find all clusters during the post-acute period that didn't have a code occurring during the code-specific washout period
#' @param conds_post_acute table with conditions during the post-acute period, with a `cluster` column
#' @param conds_towash table with conditions occurring during the code-specific washout period, with a `cluster` column
#' @return table with conditions occurring during the post-acute period 
wash_condition <- function(conds_postacute,
                           conds_towash){
  conds_postacute %>% 
    anti_join(conds_towash, by = c('person_id', 'cluster'))
  
}
determine_cluster_qual <- function(cluster_output){
  cluster_output%>%
    group_by(person_id,ce_date) %>%
    summarise(num_cluster_dx=as.integer(sum(num_dx)),
              num_clusters=as.integer(n_distinct(cluster)),
              cluster_dx_grx=max(cluster_dx_grx),
              max_days_diff=max(max_days_diff),
              min_dx_date=min(min_dx_date))%>%
    ungroup()%>%
    mutate(cluster_dx_grx=case_when(cluster_dx_grx==1L~TRUE,
                                    TRUE~FALSE))
}