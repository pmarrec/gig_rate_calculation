%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for the Quality Check (QC) and renaming  of some values 
% of the rate data based on the following criteria:
%
% 0) Get only 2 decimal digits for the numeric values in the table
%
% 1) Change grazing rates < 0 (and g_std) to n/d and change grazing rates =
% NaN to n/d (and g_std)
%
% 2) Change muN = NaN (and mu_N_std) to n/n
%
% 3) Change mu0 = NaN (and mu_N_std) to n/d.
%
% 4) if dilution (dilution level for the dilution experiment) > 0.4 (40%), 
% iode_quality_flag (QC flag) = 3 (questionable). The optimal dilution
% level for the 2-points method is <40% (Morison and Menden-Deuer, 2017)
%
% 5) if Chlad10 (<10um) or Chlau10 (>10um) concentrations are < 0.02 mg m-3, the rates for
% these size fractions are considered questionable (iode_quality_flag (QC
% flag) = 3)
%
% 6) if Chlad10per (<10um) or Chlau10per (>10um) relative contribution to total Chl-a 
% are < 0.02 (2%), the rates for these size fractions are considered questionable (iode_quality_flag (QC
% flag) = 3)
%
% Input: GIG-rates.csv file
%
% Outputs: GIG-rates-qc.csv file.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 5/9/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\Andria\Data Validation\';

%Open the raw csv file
T1=readtable("GIG-rates.csv");

%Create a new column "iode_quality_flag" with flag = 1
T1.iode_quality_flag=ones(height(T1),1);

% 0) Ge only 2 decimal digits for the numeric values in the table
T1.temperature_sampling=round(T1.temperature_sampling,2);
T1.temperature_incubation_avg=round(T1.temperature_incubation_avg,2);
T1.Chla_avg=round(T1.Chla_avg,2);
T1.Chla_sd=round(T1.Chla_sd,2);
T1.Chlad10=round(T1.Chlad10,2);
T1.Chlau10=round(T1.Chlau10,2);
T1.Chlad10per=round(T1.Chlad10per,2);
T1.Chlau10per=round(T1.Chlau10per,2);
T1.duration_incubation=round(T1.duration_incubation,2);
T1.dilution=round(T1.dilution,2);
T1.mu_0=round(T1.mu_0,2);
T1.mu_0_std=round(T1.mu_0_std,2);
T1.grazing=round(T1.grazing,2);
T1.grazing_std=round(T1.grazing_std,2);
T1.mu_N=round(T1.mu_N,2);
T1.mu_N_std=round(T1.mu_N_std,2);


% 1) Change grazing rates < 0 (and g_std) to n/d and grazing rates = NaN to
% n/d (and g_std)
a1=find(T1.grazing<0|isnan(T1.grazing));
T1.grazing=num2cell(T1.grazing);
T1.grazing_std=num2cell(T1.grazing_std);
for n1=1:length(a1)
T1.grazing{a1(n1)}={'n/d'};
T1.grazing_std{a1(n1)}={'n/d'};
end

% 2) Change mu_N = NaN (an mu_N_std) to n/n
a2=find(isnan(T1.mu_N));
T1.mu_N=num2cell(T1.mu_N);
T1.mu_N_std=num2cell(T1.mu_N_std);
for n2=1:length(a2)
T1.mu_N{a2(n2)}={'n/n'};
T1.mu_N_std{a2(n2)}={'n/n'};
end

% 3) Change mu0 = NaN (and mu_N_std) to n/d. \
a3=find(isnan(T1.mu_0));
T1.mu_0=num2cell(T1.mu_0);
T1.mu_0_std=num2cell(T1.mu_0_std);
for n3=1:length(a3)
T1.mu_0{a3(n3)}={'n/d'};
T1.mu_0_std{a3(n3)}={'n/d'};
end

% 4) if dilution (dilution level for the dilution experiment) > 0.4 (40%), 
% iode_quality_flag (QC flag) = 3 (questionable). The optimal dilution
% level for the 2-points method is <40% (Morison and Menden-Deuer, 2017)
a4=(T1.dilution>0.4);
T1.iode_quality_flag(a4)=3;

% 5) if Chlad10 (<10um) or Chlau10 (>10um) concentrations are < 0.02 mg m-3, the rates for
% these size fractions are considered questionable (iode_quality_flag (QC
% flag) = 3)
a51=(strcmp(T1.size_fraction,'>0<10'))&(T1.Chlad10<0.02);
T1.iode_quality_flag(a51)=3;
a52=(strcmp(T1.size_fraction,'>10&<200'))&(T1.Chlau10<0.02);
T1.iode_quality_flag(a52)=3;

% 6) if Chlad10per (<10um) or Chlau10per (>10um) relative contribution to total Chl-a 
% are < 0.02 (2%), the rates for these size fractions are considered questionable (iode_quality_flag (QC
% flag) = 3)
a61=(strcmp(T1.size_fraction,'>0<10'))&(T1.Chlad10per<0.02);
T1.iode_quality_flag(a61)=3;
a62=(strcmp(T1.size_fraction,'>10&<200'))&(T1.Chlau10per<0.02);
T1.iode_quality_flag(a62)=3;



%Save the new table csv file
newtablename=strcat(rep,'GIG-rates-qc.csv');
writetable(T1,newtablename);