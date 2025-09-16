# pasc-cp
Computable phenotype development for Post-acute Sequelae of COVID-19 in children


# Purpose

The PASC Computable Phenotype (CP) is intended to identify pediatric patients with Long COVID and assign a level of certainty of Long COVID (conclusive, probable, no evidence). The CP was developed against the RECOVER PCORnet EHR database, a large multi-site EHR dataset. The CP incorporates both diagnosis codes for Long COVID *and* symptom clusters more common in pediatric patients.

# CP Definition

The CP is applied as follows:

1) Identify cohorts eligible for assessment:  
    1) Patients with COVID-19 infection: clinical diagnosis or PCR, antigen, or nucleocapsid serology test
    2) Age < 21 years at COVID-19 infection date
    3) Sufficient follow up: $\geq$ 2 encounters with health system between 28-179 post-COVID-19 infection
2) Classify patients as having conclusive, probable, or no evidence of Long COVID as follows:
    1) **Conclusive:** 
        - 2+ PASC or MISC diagnoses (on separate dates)
    2) **Probable:**
        - 1 PASC or MISC diagnosis OR
        - $\geq$ 2 [symptom cluster](rules_application/specs/cluster_master_pasc.csv) diagnoses separated by at least 28 days in the post acute period (28-179 post-SARS-CoV-2 infection)
            - Note: respiratory and fever diagnoses are removed if occur within 14 days of another, non-COVID-19 respiratory infection
            - Note: apply a [washout period](rules_application/specs/cluster_master_pasc_washout.csv) (30 days or 2 years, dependent on acute or chronic condition type) for conditions present prior to COVID-19 infection
    3) **No evidence:** 
        - Does not meet the "conclusive" or "probable" classifications

# Code Execution

Code is formatted using a standard R framework common across PEDSnet studies and is executable against an OMOP CDM. To execute the code, users should follow instructions in the [site](rules_application/site) directory, adjusting parameters in the [site_info](rules_application/site/site_info.R) and [run](rules_application/site/run.R) files as needed. The framework facilitates connection to a remote data source and analysis of remote or local data using functions specified in the [code](rules_application/code) directory. Study-specific functions and analyses are located in the [driver.R](rules_application/code/driver.R) and [cohorts.R](rules_application/code/cohorts.R) files. The [cohorts.R](rules_application/code/cohorts.R) file defines study-specific functions and the code/driver file executes those functions and generates results.

This code is run downstream from RECOVER pipeline processes, including the Observation Derivation RECOVER (ODR) code which generates a table called `observation_derivation_recover` which is referenced throughout the PASC CP code. The ODR code is available in the [observation_derivation_recover_ml_phenotype](https://github.com/PEDSnet/PASC/tree/main/observation_derivation_recover_ml_phenotype) repository, following the same PEDSnet standard R framework format, where executable code is in the [driver.R](https://github.com/PEDSnet/PASC/blob/main/observation_derivation_recover_ml_phenotype/code/driver.R) file


The PASC CP code is applied to the cohort in the [driver.R](rules_application/code/driver.R) file. The code contains steps to identify an eligible cohort for analyses, determine a COVID-19 infection anchor date, determine observation windows around that anchor date, and assess likelihood of Long COVID on the patient level. Documentation throughout the code explains individual steps, and all study-specific functions are in the [cohorts](rules_application/code/cohorts.R) file. 

The paper "Identifying Pediatric Long COVID: Comparing an EHR Algorithm to Manual Review" details the process of applying the CP and comparing against chart review. The paper has been accepted for publication and is currently in the production workflow (as of Sept 2025).
