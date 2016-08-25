%% AMM MASTER SCRIPT
% This script will process the data for the auditory motion metaphor study
%
% Written by Anja Jamrozik and Andrew S Bock Aug 2016


%% Set up initial variables
dataDir = '/data/jag/ajamrozik/Data/AMMet/';
subjDirs = listdir(dataDir,'dirs');
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
shellScriptDir = '/data/jag/ajamrozik/clusterscripts/AMMet/reconallScripts';
if ~exist(shellScriptDir,'dir')
    mkdir(shellScriptDir);
end
logDir = '/data/jag/ajamrozik/LOGS';
mem = 20; % memory for reconall
%% Sort the dicoms and convert to nifti

%%%% add these scripts %%%%

%% Run Freesurfer's recon-all
% We ran some of the subjects already, so we just want to run the remained
fsSubs = listdir(SUBJECTS_DIR,'dirs');
for i = 1:length(subjDirs)
    sessDir = listdir(fullfile(dataDir,subjDirs{i}),'dirs');
    if ~ismember(sessDir,fsSubs)
        % make a recon-all shell script
        inAnat = fullfile(dataDir,subjDirs{i},sessDir{1},'MPRAGE/001/ACPC/MPRAGE.ACPC.nii.gz');
        reconallText = ['recon-all -i ' inAnat ' -s ' sessDir{1} ' -all'];
        fname = fullfile(shellScriptDir,[sessDir{1} '_reconall.sh']);
        fid = fopen(fname,'w');
        fprintf(fid,reconallText);
        fclose(fid);
    end
end
% make submit script
allScripts = listdir(shellScriptDir,'files');
fname = fullfile(shellScriptDir,'submit_reconall.sh');
fid = fopen(fname,'w');
for i = 1:length(allScripts);
fprintf(fid,['qsub -l h_vmem=' num2str(mem) ...
    '.2G,s_vmem=' num2str(mem) 'G -e ' logDir ' -o ' logDir ' ' ...
    fullfile(shellScriptDir,allScripts{i}) '\n']);
end
fclose(fid);

%%% The above will make shell scripts and a 'submit' script. Run the 'submit' script from chead %%%

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