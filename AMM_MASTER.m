%% AMM MASTER SCRIPT
% This script will process the data for the auditory motion metaphor study
%
% Written by Anja Jamrozik and Andrew S Bock Aug 2016


%% Set up initial variables
dataDir = '/data/jag/ajamrozik/Data/AMMet/';
subjDirs = listdir(dataDir,'dirs');


%% Sort the dicoms and convert to nifti

%%%% add these scripts %%%%

%% Run Freesurfer's recon-all


%% Anything else we did 

%%% e.g. slice timing correction? %%%

%%% e.g. Joe's scripts %%% this could just be a description, or a link to
%%% the python scripts

%%% feat stuff %%%

%% Register the functional data to the freesurfer anatomical data
% Use Freesurfer's 'bbregister' function

for i = 1:length(subjDirs)
    subDir = fullfile(dataDir,subjDirs{i});
    sessDir = listdir(subDir,'dirs');
    
    bbregister(subject_name,filefor_reg,bbreg_out_file,acq_type);
end


