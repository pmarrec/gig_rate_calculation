%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for calculation of the rates obtained during dilution (grazing)
% experiments for GIG project.
%
% Create a new table gathering all the caluclated rates
% for each tretment and flter type from the k (apparent growth rates) values
% and the associated data and metadata.
% Phytoplankton growth rates and protist grazing rates were estimated from
% 24 h changes in Chl-a concentration. For each incubation bottle,
% the apparent growth rates (k, d^-1) were calculated as: k=1⁄t×ln(C_t⁄C_0)
% where t is the incubation time (d) and C_t and C_0 the final and
% initial Chl-a concentration (µg L^-1), respectively.
% Protist grazing rates (g, d^-1) were estimated with the equation:
% g=((k_d-k_N))⁄((1-x))
% where k_d and k_N are the apparent growth rates in 20WSW and WSW nutrient
% amended treatments, respectively, and x is the achieved fraction of WSW
% in the diluted treatment calculated from T0 Chl-a in 20WSW and WSW. Accordingly,
% the instantaneous, or in situ, growth rate (mu_0, d^-1) was estimated as:
% mu_0=g+k_NoN
% where k_NoN is apparent phytoplankton growth rate k without nutrient addition.
% The potential for nutrient limitation was assessed by comparing apparent
% phytoplankton growth rates k in nutrient amended (k_N) and nonamended
% (k_NoN) replicates using a paired t-test. If a significant difference was
% found (p below 0.05) between k_N and k_NoN, nutrient-amended growth rates
% (mu_N, d^-1) were also calculated as mu_N = g + k_N. Otherwise, all k_N
% and k_NoN triplicate of replicate values were used to calculate both g and mu_0.
% When size fractionation at 10 µm was performed only on nutrient amended samples,
% growth rates reported on greater than 10 µm and less than 10 µm fractions
% in this study were nutrient-amended growth rates (mu_N) when nutrient
% limitation was observed. If no nutrient limitation was observed,
% mu_N obtained is equivalent to mu_0.
% The uncertainty of g estimates was quantified using the standard error
% of the slope fit from a linear regression between replicate k values
% and dilution levels. When the slope obtained was not significantly
% different from zero (p higher than 0.05), g was set to 0. Thus,
% the average k_N represented mu_N and the average k_NoN represented mu_0.
% A significant positive slope (i.e. higher growth in the WSW treatment
% than in the diluted) represents a violation of the method’s assumption.
% In such cases, g was reported as undetermined, and k in the undiluted
% bottles represented mu_N and mu_0. Uncertainties relative to mu_N and
% mu_0 were estimated from the standard deviations observed on k_N and
% k_NoN triplicate values.
%
%
% Input: GIG-k-values.csv files
%
% Outputs: GIG-rates.csv files.
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
T1=readtable("GIG-k-values.csv");

%identify each unique experiment
a1=unique(T1.date);


% Create a table for each cruise to store the rates
% for each cast/niskin
% The number of colums depends on the number of treatments (nutrient, light)
% and on the filter types (>0&<200, >0&<200, but also >0 for EN627-L11-B and
% >0&<10 for EN668 and will be implemented for each cruise

% Table T2 for with the good nb of rows and columns for rates
T2=table('Size',[length(a1)*3 21],'VariableTypes',...
    {'double','string','string','string',...
    'double','double','string',...
    'double','double','double','double','double',...
    'double','double','double',...
    'double','double','double','double','double','double'},...
    'VariableNames',{'date','date_time_utc_sampling','date_time_utc_start','date_time_utc_end',...
    'temperature_sampling','temperature_incubation_avg',...
    'size_fraction','Chla_avg','Chla_sd','Chlad10','Chlau10','Chlad10per','Chlau10per',...
    'duration_incubation','dilution',...
    'mu_0','mu_0_std','grazing','grazing_std','mu_N','mu_N_std'});

% Set up a loop counter for indexing the unique experiment
cnt1 = 0;
%find the rows corresponding to the corresponding experiment
for n1=1:length(a1)

    b1=(T1.date==a1(n1));

    %for each sampling depth, indetify the unique filter sizes
    a2=unique(T1.size_fraction(b1));

    for n2=1:length(a2)

        cnt1 = cnt1 +1;
        %Get the index of all the values from a given size
        %fraction
        b2=b1 & strcmp(T1.size_fraction,a2(n2));
        % Duration of the incubation
        B2 = find(b2, 1, 'first');%find the first occurence of b4=1

        % Correpsonding cruise/cast/niskin/date_time...
        T2.date(cnt1)=T1.date(B2);
        T2.date_time_utc_sampling(cnt1)=T1.date_time_utc_sampling(B2);
        T2.date_time_utc_start(cnt1)=T1.date_time_utc_start(B2);
        T2.date_time_utc_end(cnt1)=T1.date_time_utc_end(B2);
        T2.temperature_sampling(cnt1)=T1.temperature_sampling(B2);
        T2.temperature_incubation_avg(cnt1)=T1.temperature_incubation_avg(B2);
        T2.duration_incubation(cnt1)=T1.duration_incubation(B2);
        T2.dilution(cnt1)=T1.dilution(B2);
        T2.size_fraction(cnt1)=T1.size_fraction(B2);
        T2.Chla_avg(cnt1)=T1.Chla_avg(B2);
        T2.Chla_sd(cnt1)=T1.Chla_sd(B2);
        T2.Chlau10(cnt1)=T1.Chlau10(B2);
        T2.Chlad10(cnt1)=T1.Chlad10(B2);
        T2.Chlau10per(cnt1)=T1.Chlau10per(B2);
        T2.Chlad10per(cnt1)=T1.Chlad10per(B2);

        d1=[T1.dilution(b2);ones(12,1)];%When no nutrient limitation
        d2=[T1.dilution(b2);ones(6,1)];%When nutrient limitation

        k_dil=T1.k_dil(b2);
        k_dil_nan=isnan(k_dil);
        k_wsw_NoN=T1.k_wsw_NoN(b2);
        k_wsw_NoN_nan=isnan(k_wsw_NoN);
        k_wsw_N=T1.k_wsw_N(b2);
        k_wsw_N_nan=isnan(k_wsw_N);

        if sum(k_dil_nan)<length(k_dil) && ...
                (sum(k_wsw_NoN_nan)<length(k_wsw_NoN) || sum(k_wsw_N_nan)<length(k_wsw_N)) %test if k values for NoN and N treatments are available

            %Test if Nutrient limitation
            [h,p]=ttest2(k_wsw_N,k_wsw_NoN,'Tail','right','Vartype','equal');%Need to define the different argument to consider in the paired t-test


            if h==0%No nutrient Limitation

                k=[k_dil;k_wsw_NoN;k_wsw_N];
                g=(mean(k_dil,'omitnan')-mean([k_wsw_NoN;k_wsw_N],'omitnan'))/(1-T1.dilution(B2));

                if g<0
                    %mu is then equal to k_100
                    mu=mean([k_wsw_NoN;k_wsw_N],'omitnan');

                    %test if g is significantly different from
                    %0 with a ttest between k_dil and k_100
                    [h1,p1]=ttest2(k_dil,[k_wsw_NoN;k_wsw_N],'Tail','both','Vartype','equal');

                    if h1==0% No signiicant diference between k_dil and k_100
                        T2.grazing(cnt1)=0;%set g = 0
                    else% Significant difference between k_dil and k_100
                        T2.grazing(cnt1)=g;%keep the <0 values, will be considered as non determined in the rate QC
                    end

                else%g is positive
                    T2.grazing(cnt1)=g;
                    %mu is then equal to k_100 + g
                    mu=mean([k_wsw_NoN;k_wsw_N],'omitnan')+g;
                end

                T2.mu_0(cnt1)=mu;
                mu_stdev=max(std(k_dil,'omitnan'),std([k_wsw_NoN;k_wsw_N],'omitnan'));
                %mu_stdev=mdl.Coefficients{1,2};
                T2.mu_0_std(cnt1)=mu_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                g_stdev=mean([std(k_dil,'omitnan'),std([k_wsw_NoN;k_wsw_N],'omitnan')],'omitnan');
                T2.grazing_std(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                %No-nutrietn limitation, mu_N = NaN
                T2.mu_N(cnt1)=nan;
                T2.mu_N_std(cnt1)=nan;

                clear k mdl mu g mu_stdev g_stdev

            else %Nutirent Limited


                %N ammended samples are used to compute muN and g
                k=[k_dil;k_wsw_N];
                g=(mean(k_dil,'omitnan')-mean(k_wsw_N,'omitnan'))/(1-T1.dilution(B2));

                if g<0
                    %mu_N = k_wsw_N
                    muN=mean(k_wsw_N,'omitnan');
                    %test if g is significantly different from
                    %0 with a ttest between k_dil and k_100
                    [h1,p1]=ttest2(k_dil,k_wsw_N,'Tail','both','Vartype','equal');

                    if h1==0% No signiicant diference between k_dil and k_100
                        T2.grazing(cnt1)=0;%set g = 0
                    else% Significant difference between k_dil and k_100
                        T2.grazing(cnt1)=g;%keep the <0 values, will be considered as non determined in the rate QC
                    end

                else%g is positive
                    T2.grazing(cnt1)=g;
                    muN=mean(k_wsw_N,'omitnan')+g;
                end
                T2.mu_N(cnt1)=muN;% mu = y-intercept = g + average k(1)
                muN_stdev=max(std(k_dil,'omitnan'),std(k_wsw_N,'omitnan'));
                T2.mu_N_std(cnt1)=muN_stdev;% mu StdDev = Standard Error on y-intercept from the linear model
                g_stdev=mean([std(k_dil,'omitnan'),std(k_wsw_N,'omitnan')],'omitnan');
                T2.grazing_std(cnt1)=g_stdev;% mu StdDev = Standard Error on the slope from the linear model

                clear d mdl muN muN_stdev g_stdev %We keep g for muNoN computation

                %Computation of muNoN (in-situ growth rate) from k(1)NoN
                %and g calculated from N amended samples

                kNoN=mean(k_wsw_NoN,'omitnan');%k(1)NoN = mean of k(1)NoN (k(d) not included)

                if g<0
                    muNoN=kNoN;
                else
                    muNoN=kNoN+g;
                end

                T2.mu_0(cnt1)=muNoN;
                T2.mu_0_std(cnt1)=max(std(k_dil,'omitnan'),std(k_wsw_NoN,'omitnan'));%Find a better way to estimate StdDev on muNoN
            end

            clear h p k_dil k_dil_nan k_wsw_NoN k_wsw_NoN_nan k_wsw_N k_wsw_N_nan g kNoN kNoN_stdev muNoN muNoN_stdev

        else%Nan Values if only NaN values for k
            T2.mu_0(cnt1)=nan;
            T2.mu_0_std(cnt1)=nan;
            T2.grazing(cnt1)=nan;
            T2.grazing_std(cnt1)=nan;
            T2.mu_N(cnt1)=nan;
            T2.mu_N_std(cnt1)=nan;

        end

    end

end

writetable(T2,'GIG-rates.csv')