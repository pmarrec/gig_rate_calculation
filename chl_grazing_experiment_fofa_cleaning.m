%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for Chl-a conc cleaning of the raw Chl-a
% obtained during dilution (grazing) experiments based on Fo/Fa ratios
% for the GIG project.
%
% Fo/Fa cleaning:
% 2 criteria, for each cast/depth:
% 1) within 1-3 range,
% 2) witihn +/- 1.5 StdDev confidence interval for a given type of filter
% GFF after screening with 200um mesh = >0&<200um
% 10um after screening with 200um mesh = >10&<200um
% All the values that don't fit these criteria are flagged as questionable
% with a iode_quality_flag = 3
% All the other values are flagged as good with a iode_quality_flag = 1
%
% Input: GIG-chl-calc.csv files with Chl-a
%        calculated from the chl_grazing_experiment_chl_calc.m script
%
% Outputs: GIG-chl-fofa-clean.csv files.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 5/2/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

%Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\Andria\Data Validation\';

%Open the raw csv file
T1=readtable("GIG-chl-calc-special.csv");

%identify each unique experiment
a1=unique(T1.date);

%find the rows corresponding to the corresponding experiment
for n1=1:length(a1)
    b1=(T1.date==a1(n1));

    %for each sampling experiment, indetify the unique filter sizes
    a2=unique(T1.filter_size(b1));

    for n2=1:length(a2)
        %Get the index of all the values from a given size
        %fraction
        b2=b1 & strcmp(T1.filter_size,a2(n2));

        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Fo/Fa cleaning
        % 2 criteria, for each cast/depth:
        % 1) within 1-3 range,
        % 2) witihn +/- 1.5 StdDev confidence interval for a given type of filter
        % GFF after screening with 200um mesh = >0&<200um
        % 10um after screening with 200um mesh = >10&<200um
        % GFF without screening with 200um mesh = >0 (for EN627 L11-B)
        % GFF with screening with 10um mesh = >0&<10um (for EN668)
        %%%%%%%%%%%%%%%%%%%%%%%%%

        %1st step, QC based on FoFa ratios
        FoFa=T1.fo_fa(b2);
        %Discard all values 1<FoFa<3
        FoFa(FoFa>3)=[];
        FoFa(FoFa<1)=[];
        FoFa(isnan(FoFa))=[];
        %Get the mean/stddev and the upper and lower limits of the
        %confidence interval
        FoFa_avg=mean(FoFa);
        FoFa_std=std(FoFa);
        FoFa_ulim=FoFa_avg+1.5*FoFa_std;
        FoFa_llim=FoFa_avg-1.5*FoFa_std;

        %Assigned a iode_quality flag (1=good, 3=questionable/suspect)
        %to the data based on the nan values, the threshold (1<x<3)
        % and the upper/lower limits defined (ulim and llim)
        for n3=1:length(b2)
            if b2(n3)==1
                if (isnan(T1.fo_fa(n3))) || (T1.fo_fa(n3)>3) ...
                        || (T1.fo_fa(n3)<1) || (T1.fo_fa(n3)>FoFa_ulim) ...
                        || (T1.fo_fa(n3)<FoFa_llim)
                    T1.iode_quality_flag(n3)=3;
                else
                    T1.iode_quality_flag(n3)=1;
                end
            end
        end


        clear FoFa FoFa_avg FoFa_std FoFa_ulim FoFa_llim


    end
end

%Save the new GIG-chl-fofa-clean.csv
writetable(T1,'GIG-chl-fofa-clean.csv')




