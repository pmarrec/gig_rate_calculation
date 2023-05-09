%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for Chl-a conc calculation of given experiments where some
% problems occured
%
% Input: GIG-chl-chl-calc.csv files
%
% Outputs: GIG-chl-chl-calc-clean.csv files.
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
T1=readtable("GIG-chl-calc.csv");

%identify each unique experiment
a1=unique(T1.date);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%20220314
%Bad f0/fa at T0, fo looks good, but fa values are bad
%Use fo/fa TF (for each filter type) to compute fa-blank values at T0
%Identify all values obtained on 20220314
n=4; b1=(T1.date==a1(n));

%Identify all values obtained with at T0 with >0&<200 filters
c11=b1 & strcmp(T1.T0_TF,'T0') & strcmp(T1.filter_size,'>0&<200');
%Identify all values obtained at TF with >0&<200 filters
c21=b1 & strcmp(T1.T0_TF,'TF') & strcmp(T1.filter_size,'>0&<200');
%Identify all values obtained with at T0 with >10&<200 filters
c12=b1 & strcmp(T1.T0_TF,'T0') & strcmp(T1.filter_size,'>10&<200');
%Identify all values obtained at TF with >10&<200 filters
c22=b1 & strcmp(T1.T0_TF,'TF') & strcmp(T1.filter_size,'>10&<200');

%Get Fo/Fa TF for >0&<200 filters
fofa21=mean(T1.fo_fa(c21));
%Replace Fo/Fa T0 by the mean TF value
T1.fo_fa(c11)=fofa21;
%Get Fo/Fa TF for >10&<200 filters
fofa22=mean(T1.fo_fa(c22));
%Replace Fo/Fa T0 by the mean TF value
T1.fo_fa(c12)=fofa21;

%Get new fa_blank values at T0 for >0&<200 filters
T1.fa_blank(c11)=T1.fo_blank(c11)./fofa21;
%Get new fa_blank values at T0 for >10&<200 filters
T1.fa_blank(c12)=T1.fo_blank(c12)./fofa22;

%Calculation of chl
T1.chl=T1.Fs.*(T1.r./(T1.r-1)).*(T1.fo_blank-T1.fa_blank)...
    .*T1.vol_extracted.*T1.dilution_factor./T1.vol_filtered;
%Calculation of phaeo
T1.phaeo=T1.Fs.*(T1.r./(T1.r-1)).*(T1.r.*T1.fa_blank-T1.fo_blank)...
    .*T1.vol_extracted.*T1.dilution_factor./T1.vol_filtered;

clear n b1 c11 c12 c21 c22 fofa21 fofa22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%20221010
%Bad f0/fa at T0 for WSW, fa WSW values are kind of the same as f0 dil values
%Use fo/fa TF (for each filter type) to compute fa-blank values from fo/fa 
%at T0 WSW
%Identify all values obtained on 20220314
n=31; b1=(T1.date==a1(n));

%Identify all values obtained with at T0 with >0&<200 filters for WSW
c11=b1 & strcmp(T1.T0_TF,'T0') & strcmp(T1.filter_size,'>0&<200')...
    & strcmp(T1.dilution,'wsw');
%Identify all values obtained at TF with >0&<200 filters
c21=b1 & strcmp(T1.T0_TF,'TF') & strcmp(T1.filter_size,'>0&<200');
%Identify all values obtained with at T0 with >10&<200 filters for WSW
c12=b1 & strcmp(T1.T0_TF,'T0') & strcmp(T1.filter_size,'>10&<200')...
    & strcmp(T1.dilution,'wsw');
%Identify all values obtained at TF with >10&<200 filters
c22=b1 & strcmp(T1.T0_TF,'TF') & strcmp(T1.filter_size,'>10&<200');

%Get Fo/Fa TF for >0&<200 filters
fofa21=mean(T1.fo_fa(c21));
%Replace Fo/Fa T0 by the mean TF value
T1.fo_fa(c11)=fofa21;
%Get Fo/Fa TF for >10&<200 filters
fofa22=mean(T1.fo_fa(c22));
%Replace Fo/Fa T0 by the mean TF value
T1.fo_fa(c12)=fofa21;

%Get new fa_blank values at T0 for >0&<200 filters
T1.fa_blank(c11)=T1.fo_blank(c11)./fofa21;
%Get new fa_blank values at T0 for >10&<200 filters
T1.fa_blank(c12)=T1.fo_blank(c12)./fofa22;

%Calculation of chl
T1.chl=T1.Fs.*(T1.r./(T1.r-1)).*(T1.fo_blank-T1.fa_blank)...
    .*T1.vol_extracted.*T1.dilution_factor./T1.vol_filtered;
%Calculation of phaeo
T1.phaeo=T1.Fs.*(T1.r./(T1.r-1)).*(T1.r.*T1.fa_blank-T1.fo_blank)...
    .*T1.vol_extracted.*T1.dilution_factor./T1.vol_filtered;

clear n b1 c11 c12 c21 c22 fofa21 fofa22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%20221128
%Bad f0/fa at T0 for WSW, fa WSW values are kind of the same as f0 dil values
%Use fo/fa TF (for each filter type) to compute fa-blank values from fo/fa 
%at T0 WSW
%Identify all values obtained on 20220314
n=38; b1=(T1.date==a1(n));

%Identify all values obtained with at T0 with >0&<200 filters for WSW
c11=b1 & strcmp(T1.T0_TF,'T0') & strcmp(T1.filter_size,'>0&<200')...
    & strcmp(T1.dilution,'dil');
%Identify all values obtained at TF with >0&<200 filters
c21=b1 & strcmp(T1.T0_TF,'TF') & strcmp(T1.filter_size,'>0&<200');

%Get Fo/Fa TF for >0&<200 filters
fofa21=mean(T1.fo_fa(c21));
%Replace Fo/Fa T0 by the mean TF value
T1.fo_fa(c11)=fofa21;

%Get new fa_blank values at T0 for >0&<200 filters
T1.fa_blank(c11)=T1.fo_blank(c11)./fofa21;

%Calculation of chl
T1.chl=T1.Fs.*(T1.r./(T1.r-1)).*(T1.fo_blank-T1.fa_blank)...
    .*T1.vol_extracted.*T1.dilution_factor./T1.vol_filtered;
%Calculation of phaeo
T1.phaeo=T1.Fs.*(T1.r./(T1.r-1)).*(T1.r.*T1.fa_blank-T1.fo_blank)...
    .*T1.vol_extracted.*T1.dilution_factor./T1.vol_filtered;

clear n b1 c11 c21 fofa21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%20221226
%Bad T0 WSW >10um values. Use of the %>10um 20% to compute T0 WSW >10um.
%Identify all values obtained on 20220314
n=41; b1=(T1.date==a1(n));

%Identify all values obtained with at T0 with >0&<200 filters for WSW
c11=b1 & strcmp(T1.T0_TF,'T0') & strcmp(T1.filter_size,'>0&<200')...
    & strcmp(T1.dilution,'wsw');
%Identify all values obtained with at T0 with >10&<200 filters for WSW
c1=b1 & strcmp(T1.T0_TF,'T0') & strcmp(T1.filter_size,'>10&<200')...
    & strcmp(T1.dilution,'wsw');
%Identify all values obtained with at T0 with >0&<200 filters for 20%
c12=b1 & strcmp(T1.T0_TF,'T0') & strcmp(T1.filter_size,'>0&<200')...
    & strcmp(T1.dilution,'dil') & strcmp(T1.replicate_chl,'c');
%Identify all values obtained with at T0 with >10&<200 filters for 20%
c13=b1 & strcmp(T1.T0_TF,'T0') & strcmp(T1.filter_size,'>10&<200')...
    & strcmp(T1.dilution,'dil');
%Calculate the >10%/<10% from the 20% dilution
u10per=mean(T1.chl(c13),'omitnan')./mean(T1.chl(c12),'omitnan');

%Calculation of chl
T1.chl(c1)=mean(T1.chl(c11),'omitnan').*u10per;
%Calculation of phaeo
T1.phaeo(c1)=mean(T1.phaeo(c11),'omitnan').*u10per;

clear n b1 c11 c12 c21 c22 fofa21 fofa22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save the new GIG_chl-calc,csv
writetable(T1,'GIG-chl-calc-special.csv')

