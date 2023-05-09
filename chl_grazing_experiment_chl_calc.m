%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for Chl-a conc computation from the raw fluorescence Fo and
% Fa obtained during dilution (grazing) experiments of the GIG project.
%
% Chl-a computation:
% Chl-a = [Fs ([r/(r-1)] (Fo-Fa))] * vol_extracted / vol_filtered
% Phaeo = [Fs([r/(r-1)](r*Fa-Fo))] * vol_extracted / vol_filtered
%
% Input: GIG_raw_data_file.csv files with Fo, Fa,
%        blank values and calibration coefficients Fs and r
%
% Outputs: GIG-chl-calc.csv files.
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
T=readtable("GIG_raw_data_file.csv");

%Calculation of fo_blank
T.fo_blank=T.fo-T.blank;
%Calculation of fa_blank
T.fa_blank=T.fa-T.blank_acid;
%Calculation of fo_fa
T.fo_fa=T.fo_blank./T.fa_blank;
%Calculation of chl
T.chl=T.Fs.*(T.r./(T.r-1)).*(T.fo_blank-T.fa_blank)...
    .*T.vol_extracted.*T.dilution_factor./T.vol_filtered;
%Calculation of phaeo
T.phaeo=T.Fs.*(T.r./(T.r-1)).*(T.r.*T.fa_blank-T.fo_blank)...
    .*T.vol_extracted.*T.dilution_factor./T.vol_filtered;

%Save the new GIG_chl-calc,csv
writetable(T,'GIG-chl-calc.csv')
