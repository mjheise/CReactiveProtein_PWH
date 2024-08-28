#################################################
#                                               #  
#         ARCHES: CRP and ART Adherence         #
#                                               #
#                  MJ Heise                     #
#               August 28, 2024                 #
#                                               #
#################################################

# CODE DESCRIPTION: This code examines the relation between C-reactive protein,
# viral suppression, and antiretroviral therapy adherence in a national sample
# of sexual-minority men living with HIV.
#
# 1. READ IN DATA:
# Reads in and cleans data exported from Redcap, creates new variables for analysis.
#
# 2. TABLE 1: 
# Summarize participant demographics and health characteristics.
#
# 3. PREDICTORS OF ELEVATED CRP:
# Fit a logistic regression in which elevated CRP (>= 8) is predicted by age, 
# race/ethnicity, stimulant use, viral suppression, and ART adherence.


# Libraries
library(tidyverse) # v.2.0.0, data manipulation
library(officer) # V.0.3.15; powerpoint ggplot
library(rvg) # v.0.2.5; powerpoint ggplot
library(lubridate) # v.1.9.3; date conversion
library(emmeans) # v.1.8.9, marginal means contrasts
library(stringi) # v.1.5.0, stri_sub function
library(psych) # v. 2.3.9, describe function
library(gtsummary) # v. 1.7.2, tbl_regression function
library(zipcodeR) # v. 0.3.5, zip code mapping to categorize U.S. regions


# Functions
# Function to output odds ratios and 95% confidence intervals from model output
# -Input: Fit object (e.g., fit1) from a regression model.
# -Notes: Additional estimates will appear at the end of model output. Low 95% CI
# corresponds to the lower 95% confidence interval, Up 95% CI corresponds to the 
# upper 95% confidence interval. 

summary_or <- function(m) {
  s <- summary(m)
  
  lowCI <- (exp(s$coefficients[,1] - 1.96*s$coefficients[,2]))
  highCI <- (exp(s$coefficients[,1] + 1.96*s$coefficients[,2]))
  or <- exp(s$coefficients[,1])
  
  s$coefficients <- cbind(s$coefficients, or)
  s$coefficients <- cbind(s$coefficients, lowCI)
  s$coefficients <- cbind(s$coefficients, highCI)
  
  colnames(s$coefficients)[5] <- "Odds Ratio"
  colnames(s$coefficients)[6] <- "Low 95% CI"
  colnames(s$coefficients)[7] <- "Up 95% CI"
  return(s)
}

# Function to save ggplot object in powerpoint slide
# -Input: The ggplot object that you want to save in a powerpoint.
# -Optional inputs: Specified width and height of the outputted graph.
#  If no arguments are specified, the graph will encompass the entire
#  powerpoint slide.
# -Notes: After running the function, a window will open that 
#  allows you to select the powerpoint file. The graph will
#  save on a new slide at the end of the powerpoint.

create_pptx <- function(plt = last_plot(), path = file.choose(), width = 0, height = 0){
  if(!file.exists(path)) {
    out <- read_pptx()
  } else {
    out <- read_pptx(path)
  }
  
  if (width != 0 & height != 0) {
    out %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(value = dml(ggobj = plt), location = ph_location(left = 0, top = 0,
                                                               width = width, height = height)) %>%
      print(target = path)
  } else {
    out %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(value = dml(ggobj = plt), location = ph_location_fullsize()) %>%
      print(target = path)
    
  }
  
}

# Remove characters from the end of a string
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


#### 1. READ IN DATA ####
# Set working directory
df_fold <- ''

# List csv in wd
csvs <- list.files(path = df_fold, pattern = '\\.csv$', full.names = T)

# Read in data
rdat_original <- read.csv(csvs[1])

# Subset variables for analysis and append with event name (consent data or 
# baseline survey)
rdat_original %>%
  mutate(survey_sent = ymd(survey_sent),
         participantbirthday = ymd(participantbirthday),
         hivStatusComplete = case_when(hiv_status_submission_complete == 2 ~ 'Complete',
                                       hiv_status_submission_complete == 0 ~ 'Incomplete',
                                       hiv_status_submission_complete == 1 ~ 'Incomplete',
                                       is.na(hiv_status_submission_complete) ~ NA),
         vlComplete = case_when(blood_draw_complete == 1 ~ 1,
                                blood_draw_complete == 2 ~ 1,
                                blood_draw_complete == 3 ~ 1, 
                                .default = 0),
         subName = paste(ppt_fname, ppt_lname, collapse = ''),
         consentTemp = case_when(grepl(x = informed_consent_timestamp, 'not completed') ~ NA,
                                 (informed_consent_complete == 2) ~ stri_sub(informed_consent_timestamp, 1, 10),
                                 .default = NA),
         consentDate = ymd(consentTemp)) %>%
  rename(subNo = participant_id,
         dob_consent = participantbirthday,
         eligible = preconsent_consenter_information_complete,
         consentComplete = informed_consent_complete,
         urineComplete = urine_test_submission_complete,
         surveyComplete = baseline_survey_complete,
         yearlyCheckinComplete = yearly_checkin_complete,
         participantEmail = participant_email,
         vlOrderSent = blood_draw_order) %>%
  select(subNo, 
         participantEmail, 
         dob_consent, 
         consentDate, 
         redcap_event_name, 
         eligible, 
         consentComplete, 
         hivStatusComplete, 
         vlComplete,
         surveyComplete, 
         urine_test_sent, 
         urineComplete, 
         yearlyCheckinComplete, 
         participantzip, 
         vlOrderSent, 
         viralload, 
         survey_sent, 
         crpresult,
         urinetestresultr1, 
         urinetestresultr2, 
         hispanic, 
         raceethnicity___1, 
         raceethnicity___2, 
         raceethnicity___3, 
         raceethnicity___4, 
         raceethnicity___5, 
         raceethnicity___6, 
         raceethnicity___7, 
         raceethnicity___8, 
         income, 
         homeless_12mos, 
         cocaine_3mos, 
         meth_3mos, 
         stim_3mos, 
         assistever___4, 
         assistever___5, 
         assistever___6) %>%
  pivot_wider(id_cols = c('subNo', 'participantEmail', 'dob_consent', 'consentDate', 'participantzip', 'survey_sent'),
              names_from = 'redcap_event_name',
              values_from = c('consentComplete', 
                              'hivStatusComplete', 
                              'eligible',
                              'vlOrderSent', 
                              'vlComplete', 
                              'surveyComplete', 
                              'urine_test_sent', 
                              'urineComplete', 
                              'yearlyCheckinComplete', 
                              'viralload', 
                              'urinetestresultr1', 
                              'urinetestresultr2', 
                              'hispanic', 
                              'raceethnicity___1', 
                              'raceethnicity___2', 
                              'raceethnicity___3', 
                              'raceethnicity___4', 
                              'raceethnicity___5', 
                              'raceethnicity___6', 
                              'raceethnicity___7', 
                              'raceethnicity___8', 
                              'income', 
                              'homeless_12mos', 
                              'cocaine_3mos', 
                              'meth_3mos', 
                              'stim_3mos', 
                              'assistever___4', 
                              'assistever___5', 
                              'assistever___6', 
                              'crpresult')) -> rdat_long

# Subset variables collected during consent
rdat_long %>%
  select(subNo, 
         participantEmail, 
         dob_consent, 
         participantzip, 
         consentDate, 
         eligible_consent_process_arm_1, 
         consentComplete_consent_process_arm_1) %>%
  subset(!is.na(consentComplete_consent_process_arm_1) & participantEmail != 'kevin.sassaman@ucsf.edu') -> rdat_c

# Subset variables collected during baseline survey
rdat_long %>%
  select(subNo, 
         hivStatusComplete_baseline_arm_1, 
         vlComplete_baseline_arm_1, 
         surveyComplete_baseline_arm_1, 
         urine_test_sent_baseline_arm_1, 
         urineComplete_baseline_arm_1, 
         yearlyCheckinComplete_baseline_arm_1, 
         vlOrderSent_baseline_arm_1, 
         viralload_baseline_arm_1, 
         survey_sent,
         hispanic_baseline_arm_1, 
         raceethnicity___1_baseline_arm_1, 
         raceethnicity___2_baseline_arm_1, 
         raceethnicity___3_baseline_arm_1, 
         raceethnicity___4_baseline_arm_1, 
         raceethnicity___5_baseline_arm_1, 
         raceethnicity___6_baseline_arm_1, 
         raceethnicity___7_baseline_arm_1, 
         raceethnicity___8_baseline_arm_1, 
         income_baseline_arm_1, 
         homeless_12mos_baseline_arm_1, 
         cocaine_3mos_baseline_arm_1, 
         meth_3mos_baseline_arm_1, 
         stim_3mos_baseline_arm_1, 
         assistever___4_baseline_arm_1, 
         assistever___5_baseline_arm_1, 
         assistever___6_baseline_arm_1,
         urinetestresultr1_baseline_arm_1, 
         urinetestresultr2_baseline_arm_1, 
         crpresult_baseline_arm_1) %>%
  subset(!is.na(hivStatusComplete_baseline_arm_1) & subNo != 220000) %>%
  mutate(vlOrderSent_plus3weeks = ymd(vlOrderSent_baseline_arm_1) + days(21),
         todaysDate = today(tzone = "UTC"),
         vlOrderSent_after3weeks = as.numeric(substrRight(difftime(ymd(todaysDate), ymd(vlOrderSent_plus3weeks), units = 'days'), 4)),
         vlInProgress = case_when(vlOrderSent_after3weeks < 0 ~ 1, 
                                  .default = NA),
         surveySent_timeSince = as.numeric(substrRight(difftime(ymd(todaysDate), ymd(survey_sent), units = 'days'), 4)),
         surveyInProgress = case_when(surveySent_timeSince < 21 ~ 1, 
                                      .default = NA),
         urineSent_timeSince = as.numeric(substrRight(difftime(ymd(todaysDate), ymd(urine_test_sent_baseline_arm_1), units = 'days'), 4)),
         urineInProgress = case_when(urineSent_timeSince < 21 ~ 1, 
                                     .default = NA)) %>%
  select(subNo, 
         hivStatusComplete_baseline_arm_1, 
         vlInProgress, 
         vlComplete_baseline_arm_1, 
         surveyComplete_baseline_arm_1, 
         urineSent_timeSince, 
         urineInProgress, 
         urineComplete_baseline_arm_1, 
         yearlyCheckinComplete_baseline_arm_1, 
         viralload_baseline_arm_1, 
         survey_sent, surveyInProgress,
         hispanic_baseline_arm_1, 
         raceethnicity___1_baseline_arm_1, 
         raceethnicity___2_baseline_arm_1, 
         raceethnicity___3_baseline_arm_1, 
         raceethnicity___4_baseline_arm_1, 
         raceethnicity___5_baseline_arm_1, 
         raceethnicity___6_baseline_arm_1, 
         raceethnicity___7_baseline_arm_1, 
         raceethnicity___8_baseline_arm_1, 
         income_baseline_arm_1, 
         homeless_12mos_baseline_arm_1, 
         cocaine_3mos_baseline_arm_1, 
         meth_3mos_baseline_arm_1, 
         stim_3mos_baseline_arm_1, 
         assistever___4_baseline_arm_1, 
         assistever___5_baseline_arm_1, 
         assistever___6_baseline_arm_1,
         urinetestresultr1_baseline_arm_1, 
         urinetestresultr2_baseline_arm_1, 
         crpresult_baseline_arm_1) -> rdat_b

# Merge consent and baseline data
rdat_wide <- merge(rdat_c, rdat_b, by = 'subNo', all = T)

# Create variables for analysis (demographics, stimulant use)
rdat_wide %>%
  mutate(eligible = case_when(eligible_consent_process_arm_1 == 2 ~ 1, 
                              .default = 0),
         consentComplete = case_when(consentComplete_consent_process_arm_1 == 2 ~ 1, 
                                     .default = 0),
         vlComplete = case_when(vlComplete_baseline_arm_1 == 1 ~ 1, 
                                .default = 0),
         surveyComplete = case_when(surveyComplete_baseline_arm_1 == 2 ~ 1, 
                                    .default = 0),
         urineComplete = case_when(urineComplete_baseline_arm_1 == 2 ~ 1, 
                                   .default = 0),
         baselineYearComplete = case_when(yearlyCheckinComplete_baseline_arm_1 == 2 ~ 1, 
                                          .default = 0),
         ageYrs = time_length(difftime(as.Date(consentDate), as.Date(dob_consent)), 'years'),
         ageCategorical = case_when(ageYrs < 25 ~ '18-24 yrs',
                                    ageYrs > 24 & ageYrs < 35 ~ '25-34 yrs',
                                    ageYrs > 34 & ageYrs < 45 ~ '35-44 yrs',
                                    ageYrs > 44 & ageYrs < 55 ~ '45-54 yrs',
                                    ageYrs > 54 & ageYrs < 65 ~ '55-64 yrs',
                                    ageYrs > 65 ~ '65+ yrs',
                                    .default = NA),
         ageYrs_scale = ageYrs/10,
         race_baseline = case_when(raceethnicity___1_baseline_arm_1 == 1 ~ 'Black',
                                   raceethnicity___2_baseline_arm_1 == 1 ~ 'Black',
                                   raceethnicity___3_baseline_arm_1 == 1 ~ 'White',
                                   raceethnicity___4_baseline_arm_1 == 1 ~ 'Asian',
                                   raceethnicity___5_baseline_arm_1 == 1 ~ 'PacificIsl',
                                   raceethnicity___6_baseline_arm_1 == 1 ~ 'AmerInd',
                                   raceethnicity___7_baseline_arm_1 == 1 ~ 'MiddleEast',
                                   raceethnicity___8_baseline_arm_1 == 1 ~ 'Other',
                                   .default = NA),
         income_baseline = case_when(income_baseline_arm_1 == 1 ~ '$0',
                                     income_baseline_arm_1 == 2 ~ '$1-9,999',
                                     income_baseline_arm_1 == 3 ~ '$10,000-24,999',
                                     income_baseline_arm_1 == 4 ~ '$25,000-49,999',
                                     income_baseline_arm_1 == 5 ~ '$50,000-74,999',
                                     income_baseline_arm_1 == 6 ~ '$75,000-99,999',
                                     income_baseline_arm_1 == 7 ~ '$100,000-149,999',
                                     income_baseline_arm_1 == 8 ~ '$150,000+'),
         cocaine_3mos_baseline_arm_1 = case_when(assistever___4_baseline_arm_1 == 0 ~ 1,
                                                .default = cocaine_3mos_baseline_arm_1),
         meth_3mos_baseline_arm_1 = case_when(assistever___6_baseline_arm_1 == 0 ~ 1,
                                              .default = meth_3mos_baseline_arm_1),
         stim_3mos_baseline_arm_1 = case_when(assistever___5_baseline_arm_1 == 0 ~ 1,
                                              .default = stim_3mos_baseline_arm_1),
         stimUse3mos_baseline = case_when(cocaine_3mos_baseline_arm_1 > 1 ~ 1,
                                               meth_3mos_baseline_arm_1 > 1 ~ 1,
                                               stim_3mos_baseline_arm_1 > 1 ~ 1,
                                               .default = 0),
         stimUseEver_baseline = case_when(assistever___4_baseline_arm_1 == 1 ~ 1,
                                          assistever___5_baseline_arm_1 == 1 ~ 1,
                                          assistever___6_baseline_arm_1 == 1 ~ 1,
                                          .default = 0),
         crpresult_baseline = as.numeric(gsub('[<]', '', crpresult_baseline_arm_1)),
         crpElevated_baseline = case_when(crpresult_baseline == '<0.2' ~ 0,
                                          crpresult_baseline == '<3.0' ~ 0,
                                          crpresult_baseline == '<5.0' ~ 0,
                                          crpresult_baseline > 8 ~ 1,
                                          crpresult_baseline < 8.1 ~ 0,
                                          .default = NA),
         viralLoad_baseline = case_when(.$viralload_baseline_arm_1 == 1 ~ 'Undetectable', 
                                        .$viralload_baseline_arm_1 == 2 ~ 'Suppressed',
                                        .$viralload_baseline_arm_1 == 3 ~ 'Unsuppressed'),
         zipcode = substr(participantzip, 1, 5),
         stimUse3mos_baseline_cat = case_when(stimUse3mos_baseline == 0 ~ 'No',
                                              stimUse3mos_baseline == 1 ~ 'Yes',
                                              .default = NA),
         stimUseEver_baseline_cat = case_when(stimUseEver_baseline == 0 ~ 'No',
                                              stimUseEver_baseline == 1 ~ 'Yes',
                                              .default = NA),
         viralLoad_baseline_unsuppressed = case_when(viralLoad_baseline == 'Undetectable' ~ 0,
                                                     viralLoad_baseline == 'Suppressed' ~ 0,
                                                     viralLoad_baseline == 'Unsuppressed' ~ 1),
         viralLoad_baseline_refSuppressed = case_when(viralLoad_baseline == 'Undetectable' ~ 'Suppressed',
                                                      viralLoad_baseline == 'Suppressed' ~ 'Suppressed',
                                                      viralLoad_baseline == 'Unsuppressed' ~ 'Unsuppressed'),
         viralLoad_baseline_refSuppressed = relevel(factor(viralLoad_baseline_refSuppressed), ref = 'Suppressed'),
         viralLoad_baseline_undetectable = case_when(viralLoad_baseline == 'Undetectable' ~ 1,
                                                     viralLoad_baseline == 'Suppressed' ~ 0,
                                                     viralLoad_baseline == 'Unsuppressed' ~ 0),
         urineTest_baseline = case_when(urinetestresultr2_baseline_arm_1 == 0 ~ 'Not detected',
                                        urinetestresultr2_baseline_arm_1 == 1 ~ 'TFV detected',
                                        .default = NA),
         urineTest_baseline_refDetected = relevel(factor(urineTest_baseline), ref = 'TFV detected')) %>%
  rename(hispanic_baseline = hispanic_baseline_arm_1,
         homeless12mos_baseline = homeless_12mos_baseline_arm_1) -> rdat

# If participants haven't completed the baseline survey, replace with NA values (instead of 0s)
rdat %>%
  mutate(race_baseline = case_when(surveyComplete == 1 ~ race_baseline,
                                   .default = NA),
         hispanic_baseline = case_when(surveyComplete == 1 ~ hispanic_baseline,
                                       .default = NA),
         homeless12mos_baseline = case_when(surveyComplete == 1 ~ homeless12mos_baseline,
                                            .default = NA),
         income_baseline = case_when(surveyComplete == 1 ~ income_baseline,
                                     .default = NA),
         stimUse3mos_baseline = case_when(surveyComplete == 1 ~ stimUse3mos_baseline,
                                          .default = NA),
         stimUseEver_baseline = case_when(surveyComplete == 1 ~ stimUseEver_baseline,
                                          .default = NA),
         race_baseline_cat = case_when(race_baseline == 'White' ~ 'White',
                                       race_baseline == 'AmerInd' ~ 'Other',
                                       race_baseline == 'Asian' ~ 'Asian',
                                       race_baseline == 'Black' ~ 'Black',
                                       race_baseline == 'MiddleEast' ~ 'Other',
                                       race_baseline == 'PacificIsl' ~ 'Asian',
                                       .default = NA),
         race_baseline_refW = relevel(factor(race_baseline_cat), ref = 'White'),
         latino = case_when(is.na(hispanic_baseline) ~ 'Not latino',
                            hispanic_baseline == 0 ~ 'Not latino',
                            hispanic_baseline == 1 ~ 'Latino',
                            .default = NA),
         stimUseEver_baseline_refYes = relevel(factor(stimUseEver_baseline_cat), ref = 'Yes'),
         meth_3mo_baseline = case_when(meth_3mos_baseline_arm_1 == 1 ~ 'Never',
                                       meth_3mos_baseline_arm_1 == 2 ~ '1-2x',
                                       meth_3mos_baseline_arm_1 == 3 ~ 'Monthly',
                                       meth_3mos_baseline_arm_1 == 4 ~ 'Weekly',
                                       meth_3mos_baseline_arm_1 == 5 ~ 'Daily'),
         meth_3mo_baseline = relevel(factor(meth_3mo_baseline), ref = 'Never'),
         meth_3mo_baseline_any = case_when(meth_3mos_baseline_arm_1 == 1 ~ 'No',
                                           meth_3mos_baseline_arm_1 > 1 ~ 'Yes'),
         cocaine3mo_any = case_when(cocaine_3mos_baseline_arm_1 == 1 ~ 'Never',
                                    cocaine_3mos_baseline_arm_1 > 1 ~ 'Any use in past 3 mos'),
         meth3mo_any = case_when(meth_3mos_baseline_arm_1 == 1 ~ 'Never',
                                 meth_3mos_baseline_arm_1 > 1 ~ 'Any use in past 3 mo'),
         stimUse3mo_any = case_when(stim_3mos_baseline_arm_1 == 1 ~ 'Never',
                                    stim_3mos_baseline_arm_1 > 1 ~ 'Any use in past 3 mo')) -> dat

# Save participant zipcodes as a df
geo = reverse_zipcode(dat$zipcode)

# Merge with zip code data to map to state
dat <- merge(dat, geo, by.x = 'participantzip', by.y = 'zipcode', all.x = T)

# Create regions from states
dat %>%
  mutate(state = case_when(zipcode == 96926 ~ 'HI',
                           zipcode == 58411 ~ 'ND',
                           zipcode == 222594 ~ 'MN',
                           zipcode == 221982 ~ 'NC',
                           zipcode == 10450 ~ 'NY',
                           zipcode == 10452 ~ 'NY',
                           zipcode == 10456 ~ 'NY',
                           zipcode == 19140 ~ 'PA',
                           zipcode == 27592 ~ 'NC',
                           zipcode == 33142 ~ 'FL',
                           zipcode == 33150 ~ 'FL',
                           zipcode == 43082 ~ 'OH',
                           zipcode == 49093 ~ 'MI',
                           zipcode == 60827 ~ 'IL',
                           zipcode == 91761 ~ 'CA',
                           zipcode == 00688 ~ 'PR',
                           .default = state),
         region = case_when(
    # Northeast
    state %in% c('CT', 'ME', 'MA', 'NH', 'NY', 'RI', 'VT') ~ 'Northeast',
    
    # Midatlantic
    state %in% c('NJ', 'PA', 'DE', 'DC', 'MD', 'VA') ~ 'Midatlantic',
    
    # Midwest
    state %in% c('IL', 'IN', 'IA', 'KS', 'MI', 'MN', 'MO', 'NE', 'ND', 'OH', 'SD', 'WI') ~ 'Midwest',
    
    # South - South Atlantic
    state %in% c('FL', 'GA', 'NC', 'SC', 'WV', 'AL', 'KY', 'MS', 'TN', 'AR', 'LA', 'OK', 'TX', 'PR') ~ 'South', # Includes Puerto Rico

    # West - Mountain
    state %in% c('AZ', 'CO', 'ID', 'MT', 'NV', 'NM', 'UT', 'WY') ~ 'Mountain',
    
    # West - Pacific
    state %in% c('AK', 'CA', 'HI', 'OR', 'WA') ~ 'Pacific',
    
    # Default case (optional)
    TRUE ~ 'Unknown'
  )) -> dat

# Create df of participants who have completed their viral load and CRP
dat %>%
  filter(vlComplete == 1 & !is.na(crpElevated_baseline)) -> datv

#### 2. TABLE 1 ####
# Demographic characteristics of participants with CRP data
# Total number of participants
total = length(datv$subNo)

# Age
describe(datv$ageYrs)
quantile(datv$ageYrs, probs = c(.25, .5, .75))

# Age (categories)
datv %>%
  group_by(ageCategorical) %>%
  summarise(n = n(),
            perc = round(n/total*100))

# Ethnicity
datv %>%
  mutate(hispanic_baseline = case_when(is.na(hispanic_baseline) ~ 0,
                                       .default = hispanic_baseline)) %>%
  group_by(hispanic_baseline) %>%
  summarise(n = n(),
            perc = round(n/total*100))

# Race
datv %>%
  group_by(race_baseline) %>%
  summarise(n = n(),
            perc = round(n/total*100))

# Geographic region
datv %>%
  group_by(region) %>%
  summarise(n = n(),
            perc = round(n/total*100)) %>%
  arrange(desc(n))

# Viral suppression at blood draw
datv %>%
  group_by(viralLoad_baseline) %>%
  summarise(n = n(),
            perc = round(n/total*100))

# Stimulant use: ever
datv %>%
  group_by(stimUseEver_baseline) %>%
  summarise(n = n(),
            perc = round(n/total*100))

# Stimulant use: past 3 months
datv %>%
  group_by(stimUse3mos_baseline) %>%
  summarise(n = n(),
            perc = round(n/total*100))

# Income
datv %>%
  group_by(income_baseline) %>%
  summarise(n = n(),
            perc = round(n/total*100))

# Homelessness
datv %>%
  group_by(homeless12mos_baseline) %>%
  summarise(n = n(),
            perc = round(n/total*100))

# CRP results
datv %>%
  group_by(crpElevated_baseline) %>%
  summarise(n = n(),
            perc = round(n/total*100))

# Urine tests complete
datv %>%
  mutate(urineTestEligible = case_when(!is.na(urineSent_timeSince) ~ 1,
                                       is.na(urineSent_timeSince) ~ 0,
                                       .default = NA)) %>%
  group_by(urineComplete, urineTestEligible) %>%
  summarise(n = n(),
            perc = round(n/total*100))

# Urine tests detected
datv %>%
  group_by(urineTest_baseline) %>%
  filter(!is.na(urineTest_baseline)) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(perc = round(n/sum(n)*100))


#### 3. PREDICTORS OF ELEVATED CRP ####
# CRP predicted by age, race/ethnicity, stimulant use, TFV urine assay, viral suppression
fit_crp <- glm(crpElevated_baseline ~ ageYrs_scale + latino + race_baseline_refW + stimUse3mo_any +
                      urineTest_baseline_refDetected + viralLoad_baseline_refSuppressed, data = datv)

# Formatted regression table
fit_crp %>%
  tbl_regression(exp = T)

# Model with continuous CRP
fit_crp_cont <- glm(crpresult_baseline ~ ageYrs_scale + latino + race_baseline_refW + stimUse3mo_any +
                      urineTest_baseline_refDetected + viralLoad_baseline_refSuppressed, data = datv)

# Formatted regression table
fit_crp_cont %>%
  tbl_regression(exp = T)

# Plot continuous CRP data: Random sample from each urine assay (TFV detected/not detected)
setseed(0828)
datv_sample <- datv %>%
  filter(!is.na(urineTest_baseline_refDetected)) %>%
  group_by(urineTest_baseline_refDetected) %>%
  sample_n(69) %>%
  ungroup()

# Plot continuous sample subset (log-10 transformed CRP)
ggplot(datv_sample, aes(x = urineTest_baseline_refDetected, y = log10(crpresult_baseline), color = viralLoad_baseline_refSuppressed)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 3) +
  theme_minimal() +
  labs(x = "Tenofovir Detected in Urine Test", 
       y = "CRP Result (log)", 
       color = "Viral Load") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot continous CRP data: All data
datv_plot <- datv %>%
  filter(!is.na(urineTest_baseline_refDetected))

# Plot continuous sample (log-10 transformed CRP)
fig1 <- ggplot(datv_plot, aes(x = urineTest_baseline_refDetected, y = log10(crpresult_baseline), color = viralLoad_baseline_refSuppressed)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 3) +
  theme_minimal() +
  labs(x = "Tenofovir Detected in Urine Test", 
       y = "CRP Result (log)", 
       color = "Viral Load") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export to powerpoint slide
create_pptx(fig1)
