# gig_rate_calculation
Calculation of Chl-a concentration and phytoplankton growth and grazing rates for the Growing In Gansett (GIG) project.

All the raw fluorescence data for each experiment are stored in the GIG-raw-data-file.csv

1) Calculation of Chl-a concentrations using the **chl_grazing_experiment_chl_calc.m** Matlab script to calculate the Chl-a and Phaoe pigments concentration, after calculation of Fo-blank, Fa-blank and Fo/Fa ratios.\
*Inputs: GIG-raw-data-file.csv file with Fo, Fa, blank values and calibration coefficients Fs and r*\
*Outputs: GIG-chl-calc.csv file*

2) Special calculation of Chl-a concentrations for few experiments where some problem occured using the **chl_grazing_experiment_chl_calc_special.m** Matlab script. All exceptions are detailed in the header of the Matlab script.\
*Inputs: GIG-chl-calc.csv file*\
*Outputs: GIG-chl-calc-special.csv file*

3) First step of the Chl-a data quality check (aka cleaning) using the **chl_grazing_experiment_fofa_cleaning** Matlab script. This first step is based on 2 criteria, for each experiment and filter type:
  - Fo/Fa within 1-3 range,
  - Fo/Fa witihn +/- 2 StdDev confidence interval for a given type of filter
All the values that don't fit these criteria are flagged as questionable with a iode_quality_flag = 3
All the other values are flagged as good with a iode_quality_flag = 1\
*Input: GIG-chl-calc-special.csv file*\
*Outputs:GIG-fofa-clean.csv file*

4) Second step of the Chl-a data quality check (aka cleaning) using the **chl_grazing_experiment_chl_conc_cleaning.m** Matlab script. This second step is based on 2 criteria:
  - All negative Chl-a conc are flagged as questionnable/suspect with a iode_quality_flag = 3
  - For each station/depth/treatment/triplicate values:
Each triplicate value should stand in the +/- 2 x %CV of the mean values of the triplicate with a QC flag=1. %CV is considered as the mean %CV obtained on a given type of filter (GFF/10um) at T0 and at TF.
All the values that don't fit these criteria are flagged as questionable with a iode_quality_flag = 3
All the other values are flagged as good with a iode_quality_flag = 1\
*Input: GIG-fofa-clean.csv file*\
*Outputs: GIG-chl-conc-clean.csv file.*

5) Calculation of the apparent growth rates k for each treatment and filter types using the **chl_grazing_experiments_k_values.m** Matlab script.
It creates a new table gathering all the caluclated k-values (apparent groth rates) for each tretment and flter type.\
The duration of each incuabtion is calculated from the date_time_utc_end and the date_time_utc_start for each experiment.\
All the mean T0 Chl-a conc are calculated from the value with iode_quality_flag = 1\
T0 WSW mean Chla (Total = Chla, >10um = Chlau10, <10um = Chlad10) are reported in the new table, in addition to the % of Chl-a </>10um (Chlad10per and Chlau10per). Chlad10 = Chla - Chlau10. If Chlad10 < 0, Chlad10 = 0, Chlad10per = 0% and Chlau10per = 100%. \
k are caluclated from each TF values. Up to 6 k values are then obatined for each treatment (dilution/nutrient/light) and filter type (>0&<200, >10&<200, >0 and >0&<10 ).\
Based on >0&200 (GFF) and >10&<200 (10um filters), >0&<10 (from difference between Chl-a conc. in >0&<200 and >10&<200 size fractions) Chl-a values are calculated and then corresponding k values. For each triplicate, Chl-a d10 triplcate values are calculated as the diffrenece between the mean Chl-a value on >0&<200 and individual triplicate Chl-a values of >10&<200.\ 
Only data with QC = 1 are considered. If Chl-a d10 <0, then Chl-a conc and k = NaN.\
*Input: CRUSIE-chl-grazing-experiments-clean.csv files*\
*Outputs: CRUISE-chla-grazing-experiments-k-values.csv files.*\

6) Calculation of the phytoplankton growth rates (mu0), microzooplankton grazing rates (g), apparent growth rates in nonamended nutrient treatment (wsw NoN, kNoN) and phytoplankton growth rates in nutrient amended treatments (muN = g + kN) using the **chl_grazing_experiments_rates.m** Matlab script.\
Associted errors (std) were estimated for all these rates. Note that kNoN and muN were estimated only when there was apparent nutrient limitation.\
Create a new table gathering all the caluclated rates for each flter type from the k (apparent growth rates) values and the associated data and metadata.\
Phytoplankton growth rates and protist grazing rates were estimated from 24 h changes in Chl-a concentration. For each incubation bottle, the apparent growth rates (k, d^-1) were calculated as:\ 
k=1⁄t×ln(C_t⁄C_0) \
where t is the incubation time (d) and C_t and C_0 the final and initial Chl-a concentration (µg L^-1), respectively.\
Protist grazing rates (g, d^-1) were estimated with the equation:\
g=((k_d-k_N))⁄((1-x))\
where k_d and k_N are the apparent growth rates in 20WSW and WSW nutrient amended treatments, respectively, and x is the achieved fraction of WSW in the diluted treatment calculated from T0 Chl-a in 20WSW and WSW. \
Accordingly, the instantaneous, or in situ, growth rate (mu_0, d^-1) was estimated as:\
mu_0=g+k_NoN\
where k_NoN is apparent phytoplankton growth rate k without nutrient addition.\
The potential for nutrient limitation was assessed by comparing apparent phytoplankton growth rates k in nutrient amended (k_N) and nonamended (k_NoN) replicates using a paired t-test. If a significant difference was found (p below 0.05) between k_N and k_NoN, nutrient-amended growth rates (mu_N, d^-1) were also calculated as:\
mu_N = g + k_N. \
Otherwise, all k_N and k_NoN triplicate of replicate values were used to calculate both g and mu_0.\
The uncertainty of g estimates was quantified using the standard error of the slope fit from a linear regression between replicate k values and dilution levels. When the slope obtained was not significantly different from zero (p higher than 0.05), g was set to 0. \
Thus, the average k_N represented mu_N and the average k_NoN represented mu_0. \
A significant positive slope (i.e. higher growth in the WSW treatment than in the diluted) represents a violation of the method’s assumption. \
In such cases, g was reported as undetermined, and k in the undiluted bottles represented mu_N and mu_0. Uncertainties relative to mu_N and mu_0 were estimated from the standard deviations observed on k_N and k_NoN triplicate values.\
*Input: GIG-k-values.csv files*\
*Outputs: GIG-rates.csv files.*\

7) Quality check of the rates and renaming of some values of the rate data based using the **chl_grazing_experiments_rates_qc.m** Matlab script following these criteria:
  - Get only 2 decimal digits for the numeric values in the table
  - Change grazing rates < 0 (and g_std) to n/d and change grazing rates = NaN to n/d
  - Change muN = NaN (and mu_N_std) to n/n
  - if dilution (dilution level for the dilution experiment) > 0.4 (40%), iode_quality_flag (QC flag) = 3 (questionable). The optimal dilution level for the 2-points method is <40% (Morison and Menden-Deuer, 2017)
  - if Chlad10 (<10um) or Chlau10 (>10um) concentrations are < 0.02 mg m-3, the rates for these size fractions are considered questionable (iode_quality_flag (QC flag) = 3)
  - if Chlad10per (<10um) or Chlau10per (>10um) relative contribution to total Chl-a are < 0.02 (2%), the rates for these size fractions are considered questionable (iode_quality_flag (QC flag) = 3)
  - Change mu0 = NaN (and mu_N_std) to n/d.
*Input: GIG-rates.csv files*\
*Outputs: GIG-rates-qc.csv files.*\
