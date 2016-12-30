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

%% Register functional runs to Freesurfer anatomical
for i = 1:length(subjDirs)
    tDir = listdir(fullfile(dataDir,subjDirs{i}),'dirs');
    sessionDir = fullfile(dataDir,subjDirs{i},tDir{1});
    [~,subject_name] = fileparts(sessionDir);
    boldDirs = find_bold(sessionDir);
    for j = 1:length(boldDirs)
        featDirs = listdir(fullfile(sessionDir,boldDirs{j},'*.feat'),'dirs');
        regFile = fullfile(sessionDir,boldDirs{j},featDirs{end},'mean_func.nii.gz');
        bbreg_out_file = fullfile(sessionDir,boldDirs{j},'func_bbreg.dat');
        bbregister(subject_name,regFile,bbreg_out_file,'t2');
    end
end
%% Check the results of bbregister
for i = 1:length(subjDirs)
    tDir = listdir(fullfile(dataDir,subjDirs{i}),'dirs');
    sessionDir = fullfile(dataDir,subjDirs{i},tDir{1});
    boldDirs = find_bold(sessionDir);
    for j = 1:length(boldDirs)
        minFile = fullfile(sessionDir,boldDirs{j},'func_bbreg.dat.mincost');
        if exist(minFile,'file')
            tmp = load(minFile);
            if tmp(1) > 0.7
                disp([fullfile(sessionDir,boldDirs{j}) ' > 0.7']);
            end
        else
            disp([minFile ' doesn''t exist!']);
        end
    end
end
%% Project stats to surface
for i = 1:length(subjDirs)
    tDir = listdir(fullfile(dataDir,subjDirs{i}),'dirs');
    sessionDir = fullfile(dataDir,subjDirs{i},tDir{1});
    boldDirs = find_bold(sessionDir);
    for j = 1:length(boldDirs)
        featDirs = listdir(fullfile(sessionDir,boldDirs{j},'*.feat'),'dirs');
        bbreg_out_file = fullfile(sessionDir,boldDirs{j},'func_bbreg.dat');
        if exist(bbreg_out_file,'file')
            for k = 1:length(featDirs)
                statDir = fullfile(sessionDir,boldDirs{j},featDirs{k},'stats');
                statFiles = listdir(fullfile(statDir,'zstat*.nii.gz'),'files');
                for l = 1:length(statFiles)
                    tmp = strfind(statFiles{l},'.nii.gz');
                    baseName = statFiles{l}(1:(tmp-1));
                    % left hemisphere
                    system(['mri_vol2surf --mov ' ...
                        fullfile(statDir,statFiles{l}) ' --reg ' bbreg_out_file ...
                        ' --hemi lh --o ' fullfile(statDir,[baseName '.surf.lh.nii.gz'])]);
                    % right hemisphere
                    system(['mri_vol2surf --mov ' ...
                        fullfile(statDir,statFiles{l}) ' --reg ' bbreg_out_file ...
                        ' --hemi rh --o ' fullfile(statDir,[baseName '.surf.rh.nii.gz'])]);
                end
            end
        end
    end
end