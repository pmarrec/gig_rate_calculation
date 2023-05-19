%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for calculation of the k values (apparent growth rate)
% obtained during dilution (grazing) experiments for GIG project.
%
% Create a new table gathering all the caluclated k-values
% (apparent groth rates) for each tretment and flter type.
% The duration of each incuabtion is calculated from the date_time_utc_end
% and the date_time_utc_start for each experiment.
% Only values with iode_quality_flag = 1 are considered
% T0 WSW mean Chla (Total = Chla, >10um = Chlau10, <10um = Chlad10) are
% reported in the new table, in addition to the % of Chl-a </>10um
% (Chlad10per and Chlau10per).
% Chlad10 = Chla - Chlau10. If Chlad10 < 0, Chlad10 = 0, Chlad10per = 0%
% and Chlau10per = 100%.
% k are caluclated from each TF values.
% Up to 6 k values are then obatined
% for each treatment (dilution/nutrient/light) and filter type (>0&<200,
% >10&<200).
% Based on >0&200 (GFF) and >10&<200 (10um filters, u10 = up 10),
% >0&<10 (d10 = down 10) Chl-a values are calculated and then corresponding
% k values. For each triplicate, Chl-a d10 triplcate values are calculated
% as the diffrenece between the mean Chl-a value on >0&<200 and individual
% triplicate Chl-a values of >10&<200. Only dat with QC = 1 are considered.
% If Chl-a d10 <0, then Chl-a conc and k = NaN.
%
% Input: GIG-chl-conc-clean.csv files
%
% Outputs: GIG-k-values.csv files.
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
T1=readtable("GIG-chl-conc-clean.csv");

% Only consider the TF values
b0=strcmp(T1.T0_TF,'TF');

% Get the number of station/depth (cast/niskin) sampled during the
% cruise and the number of different treatment (light/filter size)
[c1,C1,C2]=findgroups(T1.date(b0),T1.filter_size(b0));

nrow=6;

%find the number of >10&<200 "groups" in C3 to then add the requested
%amount of rows for >0&<10 k values
f0=strcmp(C2,'>10&<200');
F0=sum(f0);

% Start with a table T2 with the good nb of rows
% Nb of rows = number of of different station/depth/filter size/light
% treatment + F0 nb of >0&<10  multiply by nrows (6 or 9 depending on
% cruises)
T2=table('Size',[(length(C1)+F0)*nrow 20],'VariableTypes',...
    {'double','string','string','string',...
    'double','double','string',...
    'string','string','double','double','double','double','double',...
    'double','double','double','double','double','double'},...
    'VariableNames',{'date',...
    'date_time_utc_sampling','date_time_utc_start','date_time_utc_end',...
    'temperature_sampling','temperature_incubation_avg', ...
    'size_fraction','replicate_bottle','replicate_chl',...
    'Chla_avg','Chla_sd','Chlad10','Chlau10','Chlad10per','Chlau10per',...
    'duration_incubation','dilution','k_dil','k_wsw_NoN','k_wsw_N'});

%identify each unique experiment
a1=unique(T1.date);

% Set up a loop counter for indexing the unique experiments
cnt1 = 0;

%find the rows corresponding to the corresponding experiment
for n1=1:length(a1)
    b1=(T1.date==a1(n1));
    %Get T0 WSW Chl-a Total, Chl-a<10um and Chl-a>10um, %<10um and
    % %>10um
    e1=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
        & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
    CHLA=mean(T1.chl(e1));
    CHLA_SD=std(T1.chl(e1));

    e2=b1 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') ...
        & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
    CHLAu10=mean(T1.chl(e2));

    a2=unique(T1.filter_size(b1));

    for n2=1:length(a2)

        cnt1 = cnt1 +1;
        %Get the index of all the values from a given size
        %fraction
        b2=b1 & strcmp(T1.filter_size,a2(n2));
        % Duration of the incubation
        B2 = find(b2, 1, 'first');%find the first occurence of b2=1
        Tinc=datenum(T1.date_time_utc_end(B2))-datenum(T1.date_time_utc_start(B2));
        %Define where to store the first value for this given cast/depth
        T2_start=(cnt1*nrow-(nrow-1));
        %Define where to store the last value for this given cast/depth
        T2_end=(cnt1*nrow);
        % Correpsonding cruise/cast/niskin/date_time...
        T2.date(T2_start:T2_end)=T1.date(B2);
        T2.date_time_utc_sampling(T2_start:T2_end)=T1.date_time_utc_sampling(B2);
        T2.date_time_utc_start(T2_start:T2_end)=T1.date_time_utc_start(B2);
        T2.date_time_utc_end(T2_start:T2_end)=T1.date_time_utc_end(B2);
        T2.temperature_sampling(T2_start:T2_end)=T1.temperature_sampling(B2);
        T2.temperature_incubation_avg(T2_start:T2_end)=T1.temperature_incubation_avg(B2);
        T2.duration_incubation(T2_start:T2_end)=Tinc;
        T2.replicate_bottle(T2_start:T2_end)={'a';'a';'a';'b';'b';'b'};
        T2.replicate_chl(T2_start:T2_end)={'a';'b';'c';'a';'b';'c'};
        T2.size_fraction(T2_start:T2_end)=T1.filter_size(B2);
        T2.Chla_avg(T2_start:T2_end)=CHLA;
        T2.Chla_sd(T2_start:T2_end)=CHLA_SD;
        T2.Chlau10(T2_start:T2_end)=CHLAu10;
        CHLAd10=CHLA-CHLAu10;
        %correct if Chla '>0&<200' by difference < 0 and then
        %compute the % of Chl-a in each size fraction
        if CHLAd10<0
            T2.Chlad10(T2_start:T2_end)=0;
            T2.Chlau10per(T2_start:T2_end)=1;
            T2.Chlad10per(T2_start:T2_end)=0;
        else
            T2.Chlad10(T2_start:T2_end)=CHLAd10;
            T2.Chlau10per(T2_start:T2_end)=CHLAu10/CHLA;
            T2.Chlad10per(T2_start:T2_end)=CHLAd10/CHLA;
        end
        %if Chla '>10&<200' > Chla '>0&<200'
        if CHLAu10>CHLA
            CHLAu10=CHLA;
            T2.Chlau10(T2_start:T2_end)=CHLAu10;
        end

        %Get the dilution level from >0&<200 filters at T0
        %This dilution level also used for size fractions
        % Identify all values obtained with >0&<200 filters at T0 dil
        % with a iode_quality_flag = 1
        c1=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
            & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
        CHL_T0_dil=mean(T1.chl(c1));
        % Mean Chl-a T0 wsw >0&<200
        c2=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
            & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
        CHL_T0_wsw=mean(T1.chl(c2));
        % Dilution level
        T2.dilution(T2_start:T2_end)=CHL_T0_dil/CHL_T0_wsw;

        clear c1 c2


        % Define the T0 Chl-a concentrations to use depending on
        % the filter type
        % Identify all values obtained with these filters at T0 dil
        % with a iode_quality_flag = 1
        c1=b1 & strcmp(T1.filter_size,a2(n2)) & strcmp(T1.T0_TF,'T0') ...
            & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
        chl_T0_dil=mean(T1.chl(c1));
        % Mean Chl-a T0 wsw >0&<200
        c2=b1 & strcmp(T1.filter_size,a2(n2)) & strcmp(T1.T0_TF,'T0') ...
            & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
        chl_T0_wsw=mean(T1.chl(c2));

        % k-values dil
        c3= b0 & b2 & strcmp(T1.dilution,'dil');
        C3=find(c3==1);
        if ~isempty(C3)
            T2.k_dil(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_dil);
            C4=find(T1.iode_quality_flag(C3)==3);
            T2.k_dil(T2_start+C4-1)=nan;
        else
            T2.k_dil(T2_start:T2_end)=nan;
        end
        clear c3 C3 C4
        % k-values wsw
        c3=b0 & b2 & strcmp(T1.dilution,'wsw') & strcmp(T1.nutrient_treatment,'NoN');
        C3=find(c3==1);
        if ~isempty(C3)
            T2.k_wsw_NoN(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
            C4=find(T1.iode_quality_flag(C3)==3);
            T2.k_wsw_NoN(T2_start+C4-1)=nan;
        else
            T2.k_wsw_NoN(T2_start:T2_end)=nan;
        end
        clear c3 C3 C4
        % k-values wsw + N High Light (65% or 100% for EN644)
        c3=b0 & b2 & strcmp(T1.dilution,'wsw') & strcmp(T1.nutrient_treatment,'N');
        C3=find(c3==1);
        if ~isempty(C3)
            T2.k_wsw_N(T2_start:T2_end)=1/Tinc*log(T1.chl(C3)./chl_T0_wsw);
            C4=find(T1.iode_quality_flag(C3)==3);
            T2.k_wsw_N(T2_start+C4-1)=nan;
        else
            T2.k_wsw_N(T2_start:T2_end)=nan;
        end
        clear c3 C3 C4

    end

    %>0&<10 data is obtained by the difference between >0&<200 and
    %>10&<200.
    if sum(ismember(a2,'>10&<200'))==1
        cnt1 = cnt1 +1;

        %Get the index of all the values from the >0&<200
        %filter_size
        b21=b1 & strcmp(T1.filter_size,'>0&<200');
        % Duration of the incubation
        B21 = find(b21, 1, 'first');%find the first occurence of b4=1
        Tinc=datenum(T1.date_time_utc_end(B21))-datenum(T1.date_time_utc_start(B21));
        %Define where to store the first value for this given cast/depth
        T2_start=(cnt1*nrow-(nrow-1));
        %Define where to store the last value for this given cast/depth
        T2_end=(cnt1*nrow);
        % Correpsonding cruise/cast/niskin/date_time...
        T2.date(T2_start:T2_end)=T1.date(B21);
        T2.date_time_utc_sampling(T2_start:T2_end)=T1.date_time_utc_sampling(B21);
        T2.date_time_utc_start(T2_start:T2_end)=T1.date_time_utc_start(B21);
        T2.date_time_utc_end(T2_start:T2_end)=T1.date_time_utc_end(B21);
        T2.temperature_sampling(T2_start:T2_end)=T1.temperature_sampling(B21);
        T2.temperature_incubation_avg(T2_start:T2_end)=T1.temperature_incubation_avg(B21);
        T2.duration_incubation(T2_start:T2_end)=Tinc;
        T2.size_fraction(T2_start:T2_end)='>0&<10';%Name '<10' for the moment, will be rename at the end of the script
        T2.replicate_bottle(T2_start:T2_end)={'a';'a';'a';'b';'b';'b'};
        T2.replicate_chl(T2_start:T2_end)={'a';'b';'c';'a';'b';'c'};
        T2.Chla_avg(T2_start:T2_end)=CHLA;
        T2.Chla_sd(T2_start:T2_end)=CHLA_SD;
        T2.Chlau10(T2_start:T2_end)=CHLAu10;
        CHLAd10=CHLA-CHLAu10;
        %correct if Chla '>0&<200' by difference < 0 and then
        %compute the % of Chl-a in each size fraction
        if CHLAd10<0
            T2.Chlad10(T2_start:T2_end)=0;
            T2.Chlau10per(T2_start:T2_end)=1;
            T2.Chlad10per(T2_start:T2_end)=0;
        else
            T2.Chlad10(T2_start:T2_end)=CHLAd10;
            T2.Chlau10per(T2_start:T2_end)=CHLAu10/CHLA;
            T2.Chlad10per(T2_start:T2_end)=CHLAd10/CHLA;
        end
        %if Chla '>10&<200' > Chla '>0&<200'
        if CHLAu10>CHLA
            CHLAu10=CHLA;
            T2.Chlau10(T2_start:T2_end)=CHLAu10;
        end

        %Get the dilution level from >0&<200 filters at T0
        %This dilution level also used for size fractions
        % Identify all values obtained with >0&<200 filters at T0 dil
        % with a iode_quality_flag = 1
        c1=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
            & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
        CHL_T0_dil=mean(T1.chl(c1));
        % Mean Chl-a T0 wsw >0&<200
        c2=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
            & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
        CHL_T0_wsw=mean(T1.chl(c2));
        % Dilution level
        T2.dilution(T2_start:T2_end)=CHL_T0_dil/CHL_T0_wsw;

        clear c1 c2


        % Identify all values obtained with >0&<200 filters at T0 dil
        % with a iode_quality_flag = 1
        c1=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
            & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
        chl_T0_dil=mean(T1.chl(c1));
        % Mean Chl-a T0 wsw >0&<200
        c2=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') ...
            & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
        chl_T0_wsw=mean(T1.chl(c2));

        % Identify all values obtained with >10&<200 filters at T0 dil
        % with a iode_quality_flag = 1
        c1=b1 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') ...
            & strcmp(T1.dilution,'dil') & T1.iode_quality_flag==1;
        chl_T0_dil_u10=mean(T1.chl(c1));
        % Mean Chl-a T0 wsw >0&<200
        c2=b1 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') ...
            & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
        chl_T0_wsw_u10=mean(T1.chl(c2));

        % Mean Chl-a T0 dil <10um (d10)
        chl_T0_dil_d10=chl_T0_dil-chl_T0_dil_u10;
        if chl_T0_dil_d10<0
            chl_T0_dil_d10=nan;
        end
        % Mean Chl-a T0 dil <10um (d10)
        chl_T0_wsw_d10=chl_T0_wsw-chl_T0_wsw_u10;
        if chl_T0_wsw_d10<0
            chl_T0_wsw_d10=nan;
        end

        % k-values dil for bottle a
        c31= b0 & b1 & strcmp(T1.filter_size,'>10&<200') ...
            & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_bottle,'a');
        c32= b0 & b1 & strcmp(T1.filter_size,'>0&<200') ...
            & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_bottle,'a');
        chl_TF=mean(T1.chl(c32));
        C3=find(c31==1);
        if ~isempty(C3)
            chl_TF_d10=chl_TF-T1.chl(C3);
            c4=chl_TF_d10<0;
            chl_TF_d10(c4)=nan;
            T2.k_dil(T2_start:T2_start+2)=1/Tinc*log(chl_TF_d10./chl_T0_dil_d10);
            C4=find(T1.iode_quality_flag(C3)==3);
            T2.k_dil(T2_start+C4-1)=nan;
        else
            T2.k_dil(T2_start:T2_start+2)=nan;
        end
        clear c31 c32 chl_TF C3 chl_TF_d10

        % k-values dil for bottle b
        c31= b0 & b1 & strcmp(T1.filter_size,'>10&<200') ...
            & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_bottle,'b');
        c32= b0 & b1 & strcmp(T1.filter_size,'>0&<200') ...
            & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_bottle,'b');
        chl_TF=mean(T1.chl(c32));
        C3=find(c31==1);
        if ~isempty(C3)
            chl_TF_d10=chl_TF-T1.chl(C3);
            c4=chl_TF_d10<0;
            chl_TF_d10(c4)=nan;
            T2.k_dil(T2_start+3:T2_start+5)=1/Tinc*log(chl_TF_d10./chl_T0_dil_d10);
            C4=find(T1.iode_quality_flag(C3)==3);
            T2.k_dil(T2_start+C4-1)=nan;
        else
            T2.k_dil(T2_start+3:T2_start+5)=nan;
        end
        clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10


        % k-values wsw NoN for bottle a
        c31= b0 & b1 & strcmp(T1.filter_size,'>10&<200') ...
            & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'a')...
            & strcmp(T1.nutrient_treatment,'NoN');
        c32= b0 & b1 & strcmp(T1.filter_size,'>0&<200') ...
            & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'a')...
            & strcmp(T1.nutrient_treatment,'NoN');
        chl_TF=mean(T1.chl(c32));
        C3=find(c31==1);
        if ~isempty(C3)
            chl_TF_d10=chl_TF-T1.chl(C3);
            c4=chl_TF_d10<0;
            chl_TF_d10(c4)=nan;
                T2.k_wsw_NoN(T2_start:T2_start+2)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_NoN(T2_start+C4-1)=nan;
        else
            T2.k_wsw_NoN(T2_start:T2_end)=nan;
        end
        clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10

        % k-values wsw NoN for bottle b
        c31= b0 & b1 & strcmp(T1.filter_size,'>10&<200') ...
            & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'b')...
            & strcmp(T1.nutrient_treatment,'NoN');
        c32= b0 & b1 & strcmp(T1.filter_size,'>0&<200') ...
            & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'b')...
            & strcmp(T1.nutrient_treatment,'NoN');
        chl_TF=mean(T1.chl(c32));
        C3=find(c31==1);
        if ~isempty(C3)
            chl_TF_d10=chl_TF-T1.chl(C3);
            c4=chl_TF_d10<0;
            chl_TF_d10(c4)=nan;        
                T2.k_wsw_NoN(T2_start+3:T2_end)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
                C4=find(T1.iode_quality_flag(C3)==3);
                T2.k_wsw_NoN(T2_start+C4-1)=nan;
            
        else
            T2.k_wsw_NoN(T2_start+3:T2_end)=nan;
        end
        clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10

        % k-values wsw + N for bottle a
        c31= b0 & b1 & strcmp(T1.filter_size,'>10&<200') ...
            & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'a')...
            & strcmp(T1.nutrient_treatment,'N');
        c32= b0 & b1 & strcmp(T1.filter_size,'>0&<200') ...
            & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'a')...
            & strcmp(T1.nutrient_treatment,'N');
        chl_TF=mean(T1.chl(c32));
        C3=find(c31==1);
        if ~isempty(C3)
            chl_TF_d10=chl_TF-T1.chl(C3);
            c4=chl_TF_d10<0;
            chl_TF_d10(c4)=nan;
            T2.k_wsw_N(T2_start:T2_start+2)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
            C4=find(T1.iode_quality_flag(C3)==3);
            T2.k_wsw_N(T2_start+C4-1)=nan;
        else
            T2.k_wsw_N(T2_start:T2_start+2)=nan;
        end
        clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10

        % k-values wsw + N for bottle b
        c31= b0 & b1 & strcmp(T1.filter_size,'>10&<200') ...
            & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'b')...
            & strcmp(T1.nutrient_treatment,'N');
        c32= b0 & b1 & strcmp(T1.filter_size,'>0&<200') ...
            & strcmp(T1.dilution,'wsw') & strcmp(T1.replicate_bottle,'b')...
            & strcmp(T1.nutrient_treatment,'N');
        chl_TF=mean(T1.chl(c32));
        C3=find(c31==1);
        if ~isempty(C3)
            chl_TF_d10=chl_TF-T1.chl(C3);
            c4=chl_TF_d10<0;
            chl_TF_d10(c4)=nan;
            T2.k_wsw_N(T2_start+3:T2_start+5)=1/Tinc*log(chl_TF_d10./chl_T0_wsw_d10);
            C4=find(T1.iode_quality_flag(C3)==3);
            T2.k_wsw_N(T2_start+C4-1)=nan;
        else
            T2.k_wsw_N(T2_start+3:T2_start+5)=nan;
        end
        clear c31 c32 c4 chl_TF C3 C4 chl_TF_d10

    end

end



%Save the new GIG-chl-fofa-clean.csv
writetable(T2,'GIG-k-values.csv')



