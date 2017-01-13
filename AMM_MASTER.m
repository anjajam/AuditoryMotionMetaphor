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
                system(['rm ' fullfile(statDir,'*surf*.nii.gz')]);
                % z-stats
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
                % f-stat
                fstatFiles = listdir(fullfile(statDir,'zfstat*.nii.gz'),'files');
                if ~isempty(fstatFiles)
                    for l = 1:length(fstatFiles)
                        tmp = strfind(fstatFiles{l},'.nii.gz');
                        baseName = fstatFiles{l}(1:(tmp-1));
                        % left hemisphere
                        system(['mri_vol2surf --mov ' ...
                            fullfile(statDir,fstatFiles{l}) ' --reg ' bbreg_out_file ...
                            ' --hemi lh --o ' fullfile(statDir,[baseName '.surf.lh.nii.gz'])]);
                        % right hemisphere
                        system(['mri_vol2surf --mov ' ...
                            fullfile(statDir,fstatFiles{l}) ' --reg ' bbreg_out_file ...
                            ' --hemi rh --o ' fullfile(statDir,[baseName '.surf.rh.nii.gz'])]);
                    end
                end
            end
        end
    end
end
%% Project Heschl's Gyrus ROIs (FSL) to subject space
for i = 1:length(subjDirs)
    tDir = listdir(fullfile(dataDir,subjDirs{i}),'dirs');
    sessionDir = fullfile(dataDir,subjDirs{i},tDir{1});
    [~,subject_name] = fileparts(sessionDir);
    subjDir = fullfile(SUBJECTS_DIR,subject_name);
    if exist(subjDir,'dir')
        % left hemisphere
        sval = '/data/jag/ajamrozik/Data/fsaverageROIs/HeschlsGyrus.surf.lh.nii.gz';
        tval = fullfile(SUBJECTS_DIR,subject_name,'label','HeschlsGyrus.surf.lh.nii.gz');
        mri_surf2surf('fsaverage',subject_name,sval,tval,'lh');
        % right hemisphere
        sval = '/data/jag/ajamrozik/Data/fsaverageROIs/HeschlsGyrus.surf.rh.nii.gz';
        tval = fullfile(SUBJECTS_DIR,subject_name,'label','HeschlsGyrus.surf.rh.nii.gz');
        mri_surf2surf('fsaverage',subject_name,sval,tval,'rh');
    end
end
%% Make the cross-correlation matrices
for i = 1:length(subjDirs)
    ct = 0;
    % Get the ROI vertices
    tDir = listdir(fullfile(dataDir,subjDirs{i}),'dirs');
    sessionDir = fullfile(dataDir,subjDirs{i},tDir{1});
    [~,subject_name] = fileparts(sessionDir);
    subjDir = fullfile(SUBJECTS_DIR,subject_name);
    if exist(subjDir,'dir')
        % Primary auditory
        locAudDir = listdir(fullfile(sessionDir,'*BOLD_LOCAUD'),'dirs');
        
        lhH = load_nifti(fullfile(SUBJECTS_DIR,subject_name,'label',...
            'HeschlsGyrus.surf.lh.nii.gz'));
        lhAud = load_nifti(fullfile(sessionDir,locAudDir{1},'OUTPUT.feat','stats',...
            'zstat1.surf.lh.nii.gz'));
        lhAudQ = quantile(lhAud.vol,0.9);
        lhAudVerts = find(lhH.vol > 5 & lhAud.vol > lhAudQ); % > 5% probability
        
        rhH = load_nifti(fullfile(SUBJECTS_DIR,subject_name,'label',...
            'HeschlsGyrus.surf.rh.nii.gz'));
        rhAud = load_nifti(fullfile(sessionDir,locAudDir{1},'OUTPUT.feat','stats',...
            'zstat1.surf.rh.nii.gz'));
        rhAudQ = quantile(rhAud.vol,0.9);
        rhAudVerts = find(rhH.vol > 5 & rhAud.vol > rhAudQ); % > 5% probability
        
        % Primary motor
        locMotDir = listdir(fullfile(sessionDir,'*BOLD_LOCMOTOR'),'dirs');
        
        lhM1 = read_label(subject_name,'lh.BA4a');
        lhM2 = read_label(subject_name,'lh.BA4p');
        lhM = [lhM1(:,1) + 1;lhM2(:,1) + 1]; % add 1 to 0-based index
        lhMot = load_nifti(fullfile(sessionDir,locMotDir{1},'OUTPUT.feat','stats',...
            'zfstat1.surf.lh.nii.gz'));
        lhMotQ = quantile(lhMot.vol,0.9);
        lhMotVerts = find(lhMot.vol > lhMotQ & ismember(1:length(lhMot.vol),lhM));
        
        rhM1 = read_label(subject_name,'rh.BA4a');
        rhM2 = read_label(subject_name,'rh.BA4p');
        rhM = [rhM1(:,1) + 1;rhM2(:,1) + 1]; % add 1 to 0-based index
        rhMot = load_nifti(fullfile(sessionDir,locMotDir{1},'OUTPUT.feat','stats',...
            'zfstat1.surf.rh.nii.gz'));
        rhMotQ = quantile(rhMot.vol,0.9);
        rhMotVerts = find(rhMot.vol > rhMotQ & ismember(1:length(rhMot.vol),rhM));
        
        % MT
        locVisMotDir = listdir(fullfile(sessionDir,'*BOLD_LOCVISMOT'),'dirs');
        
        lhMT1 = read_label(subject_name,'lh.MT');
        lhMTv = lhMT1(:,1) + 1; % add 1 to 0-based index
        lhMT = load_nifti(fullfile(sessionDir,locVisMotDir{1},'OUTPUT.feat','stats',...
            'zstat3.surf.lh.nii.gz'));
        lhMTQ = quantile(lhMT.vol,0.9);
        lhMTVerts = find(lhMT.vol > lhMTQ & ismember(1:length(lhMT.vol),lhMTv));
        
        rhMT1 = read_label(subject_name,'rh.MT');
        rhMTv = rhMT1(:,1) + 1; % add 1 to 0-based index
        rhMT = load_nifti(fullfile(sessionDir,locVisMotDir{1},'OUTPUT.feat','stats',...
            'zstat3.surf.rh.nii.gz'));
        rhMTQ = quantile(rhMT.vol,0.9);
        rhMTVerts = find(rhMT.vol > rhMTQ & ismember(1:length(rhMT.vol),rhMTv));
        
        % Correlate each sentence run with the auditory and video runs
        metDirs = listdir(fullfile(sessionDir,'*BOLD_MET*'),'dirs');
        
        soundDir = listdir(fullfile(sessionDir,'*BOLD_SOUND'),'dirs');
        lhSoundStats = listdir(fullfile(sessionDir,soundDir{1},'OUTPUT.feat','stats',...
            'zstat*.surf.lh.nii.gz'),'files');
        rhSoundStats = listdir(fullfile(sessionDir,soundDir{1},'OUTPUT.feat','stats',...
            'zstat*.surf.rh.nii.gz'),'files');
        
        videoDir = listdir(fullfile(sessionDir,'*BOLD_VIDEO'),'dirs');
        lhVideoStats = listdir(fullfile(sessionDir,videoDir{1},'OUTPUT.feat','stats',...
            'zstat*.surf.lh.nii.gz'),'files');
        rhVideoStats = listdir(fullfile(sessionDir,videoDir{1},'OUTPUT.feat','stats',...
            'zstat*.surf.rh.nii.gz'),'files');
        
        for j = 1:length(metDirs)
            indivDir = listdir(fullfile(sessionDir,metDirs{j},'*INDIV.feat'),'dirs');
            lhzstats = listdir(fullfile(sessionDir,metDirs{j},indivDir{1},...
                'zstat*.surf.lh.nii.gz'),'files');
            rhzstats = listdir(fullfile(sessionDir,metDirs{j},indivDir{1},...
                'zstat*.surf.rh.nii.gz'),'files');
            for k = 1:length(lhzstats)
                lhMet = load_nifti(fullfile(sessionDir,metDirs{j},indivDir{1},...
                    lhzstats{k}));
                rhMet = load_nifti(fullfile(sessionDir,metDirs{j},indivDir{1},...
                    rhzstats{k}));
                % Sound
                for l = 1:length(lhSoundStats)
                    lhSound = load_nifti(fullfile(sessionDir,soundDir{1},...
                        'OUTPUT.feat','stats',lhSoundStats{l}));
                    rhSound = load_nifti(fullfile(sessionDir,soundDir{1},...
                        'OUTPUT.feat','stats',rhSoundStats{l}));
                    % Primary Auditory
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        lhMet.vol(lhAudVerts),lhSound.vol(lhAudVerts));
                    labels = '';
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_sound%02d_lhAud_',j,k,l)];
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        rhMet.vol(rhAudVerts),rhSound.vol(rhAudVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_sound%02d_rhAud',j,k,l)];
                    % Primary Motor
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        lhMet.vol(lhMotVerts),lhSound.vol(lhMotVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_sound%02d_lhMot',j,k,l)];
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        rhMet.vol(rhMotVerts),rhSound.vol(rhMotVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_sound%02d_rhMot',j,k,l)];
                    % MT
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        lhMet.vol(lhMTVerts),lhSound.vol(lhMTVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_sound%02d_lhMT',j,k,l)];
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        rhMet.vol(rhMTVerts),rhSound.vol(rhMTVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_sound%02d_rhMT',j,k,l)];
                end
                % Video
                for l = 1:length(lhVideoStats)
                    lhVideo = load_nifti(fullfile(sessionDir,videoDir{1},...
                        'OUTPUT.feat','stats',lhVideoStats{l}));
                    rhVideo = load_nifti(fullfile(sessionDir,videoDir{1},...
                        'OUTPUT.feat','stats',rhVideoStats{l}));
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        lhMet.vol(lhAudVerts),lhVideo.vol(lhAudVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_video%02d_lhAud',j,k,l)];
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        rhMet.vol(rhAudVerts),rhVideo.vol(rhAudVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_video%02d_rhAud',j,k,l)];
                    % Primary Motor
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        lhMet.vol(lhMotVerts),lhVideo.vol(lhMotVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_video%02d_lhMot',j,k,l)];
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        rhMet.vol(rhMotVerts),rhVideo.vol(rhMotVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_video%02d_rhMot',j,k,l)];
                    % MT
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        lhMet.vol(lhMTVerts),lhVideo.vol(lhMTVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_video%02d_lhMT',j,k,l)];
                    ct = ct + 1;
                    corrVals(i,ct) = corr(...
                        rhMet.vol(rhMTVerts),rhVideo.vol(rhMTVerts));
                    labels = [labels ' ' ...
                        sprintf('run%02d_sentence%02d_video%02d_rhMT',j,k,l)];
                end
            end
        end
    end
end