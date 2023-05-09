%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for Chl-a conc cleaning of the raw Chl-a
% obtained during dilution (grazing) experiments for the GIG project.
%
% Chl-a conc cleaning:
% All negative Chl-a conc are flagged as questionnable/suspect with a
% iode_quality_flag = 3
% Criteria, for each experiment/treatment/triplicate values:
% Each triplicate value standing in the +/- 2 x %CV of
% the mean values of the triplicate assigned with a QC flag=1
% %CV is considered as the mean %CV obtained on a given type of
% filter (GFF/10um) at T0 AND at TF.
% All the values that don't fit these criteria are flagged as questionable
% with a iode_quality_flag = 3
% All the other values are flagged as good with a iode_quality_flag = 1
%
% Input: GIG-chl-fofa-clean.csv files
%
% Outputs: GIG-chl-conc-clean.csv files.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 5/2/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

%Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\Andria\Data Validation\';

%Open the raw csv file
T1=readtable("GIG-chl-fofa-clean.csv");

%Flag as questionable/suspect all <0 chl conc value
for n=1:height(T1)
    if T1.chl(n)<0
        T1.iode_quality_flag(n)=3;
    end
end

%identify each unique experiment
a1=unique(T1.date);

%find the rows corresponding to the corresponding experiment
for n1=1:length(a1)
    b1=(T1.date==a1(n1));

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %Chl-a conc cleaning
    % Criteria, for each experiment/triplicate values:
    % Each triplicate value should stand in the +/- !.5 x %CV of
    % the mean values of the triplicate with a QC flag=1
    % %CV is considered as the mean %CV obtained on a given type of
    % filter (GFF/10um) at T0 AND at TF.
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %>0&<200 filters (GFF)
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %Identify all values obtained with >0&<200 filters at T0
    c1=b1 & strcmp(T1.filter_size,'>0&<200');

    %Identify as group the dil and wsw values
    [d1,D1,D2,D3]=findgroups(T1.dilution(c1),T1.nutrient_treatment(c1),T1.replicate_bottle(c1));

    %Create a new vector to store the mean values
    chl_avg=nan(max(d1),1);
    %Create a new vector to store the stdev values
    chl_std=nan(max(d1),1);
    %Create a new vector to store the %CV values
    chl_cv=nan(max(d1),1);

    for m1=1:max(d1)
        %Identify the values of the given group with a
        %iode_quality_flag of 1
        c2=c1 & strcmp(T1.dilution,D1(m1)) & strcmp(T1.nutrient_treatment,D2(m1)) ...
            & strcmp(T1.replicate_bottle,D3(m1)) & T1.iode_quality_flag==1;
        chl_avg(m1)=mean(T1.chl(c2));
        chl_std(m1)=std(T1.chl(c2));
        chl_cv(m1)=chl_std(m1)/chl_avg(m1);
    end

    %Define the lower limit values
    chl_llim=chl_avg-2*mean(chl_cv)*chl_avg;
    %Define the upper limit values
    chl_ulim=chl_avg+2*mean(chl_cv)*chl_avg;

    %Check if the values for each type of treatment are in the
    %confidence interval llim<x<ulim. If they are not, values are
    %flag as questionable/suspect
    for m1=1:max(d1)
        c2=c1 & strcmp(T1.dilution,D1(m1)) & strcmp(T1.nutrient_treatment,D2(m1)) ...
            & strcmp(T1.replicate_bottle,D3(m1)) & T1.iode_quality_flag==1;
        for M1=1:length(c2)
            if c2(M1)==1
                if (T1.chl(M1)<chl_llim(m1)) || (T1.chl(M1)>chl_ulim(m1))
                    T1.iode_quality_flag(M1)=3;
                end
            end
        end
    end

    clear chl_avg chl_std chl_llim chl_ulim chl_cv

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %T0 and >10&<200 filters (10um)
    %%%%%%%%%%%%%%%%%%%%%%%%%%

        %Identify all values obtained with >0&<200 filters at T0
    c1=b1 & strcmp(T1.filter_size,'>10&<200');

    %Identify as group the dil and wsw values
    [d1,D1,D2,D3]=findgroups(T1.dilution(c1),T1.nutrient_treatment(c1),T1.replicate_bottle(c1));

    %Create a new vector to store the mean values
    chl_avg=nan(max(d1),1);
    %Create a new vector to store the stdev values
    chl_std=nan(max(d1),1);
    %Create a new vector to store the %CV values
    chl_cv=nan(max(d1),1);

    for m1=1:max(d1)
        %Identify the values of the given group with a
        %iode_quality_flag of 1
        c2=c1 & strcmp(T1.dilution,D1(m1)) & strcmp(T1.nutrient_treatment,D2(m1)) ...
            & strcmp(T1.replicate_bottle,D3(m1)) & T1.iode_quality_flag==1;
        chl_avg(m1)=mean(T1.chl(c2));
        chl_std(m1)=std(T1.chl(c2));
        chl_cv(m1)=chl_std(m1)/chl_avg(m1);
    end

    %Define the lower limit values
    chl_llim=chl_avg-2*mean(chl_cv)*chl_avg;
    %Define the upper limit values
    chl_ulim=chl_avg+2*mean(chl_cv)*chl_avg;

    %Check if the values for each type of treatment are in the
    %confidence interval llim<x<ulim. If they are not, values are
    %flag as questionable/suspect
    for m1=1:max(d1)
        c2=c1 & strcmp(T1.dilution,D1(m1)) & strcmp(T1.nutrient_treatment,D2(m1)) ...
            & strcmp(T1.replicate_bottle,D3(m1)) & T1.iode_quality_flag==1;
        for M1=1:length(c2)
            if c2(M1)==1
                if (T1.chl(M1)<chl_llim(m1)) || (T1.chl(M1)>chl_ulim(m1))
                    T1.iode_quality_flag(M1)=3;
                end
            end
        end
    end

    clear chl_avg chl_std chl_llim chl_ulim chl_cv

    
end


%Save the new GIG-chl-fofa-clean.csv
writetable(T1,'GIG-chl-conc-clean.csv')





