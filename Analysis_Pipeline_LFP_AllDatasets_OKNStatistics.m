clear all
close all
clc

% This script does all the LFP analysis from the structures generated using
% Pipeline_Analysis

%% Enumerate the datasets

nDatasets = 6;
directories.taskdirPFC{1} = 'E:\Data\H07\12-06-2016\PFC\Bfsgrad1';
recDate{1} = '12062016';
fileID{1} = '12-06-2016';
subjID{1} = 'Hayo';
directories.taskdirPFC{2} = 'E:\Data\H07\13-07-2016\PFC\Bfsgrad1';
recDate{2} = '13072016';
fileID{2} = '13-07-2016';
subjID{2} = 'Hayo';
directories.taskdirPFC{3} = 'E:\Data\H07\20161019\PFC\Bfsgrad1';
recDate{3} = '19102016';
fileID{3} = '20161019';
subjID{3} = 'Hayo';
directories.taskdirPFC{4} = 'E:\Data\H07\20161025\PFC\Bfsgrad1';
recDate{4} = '25102016';
fileID{4} = '20161025';
subjID{4} = 'Hayo';
directories.taskdirPFC{5} = 'E:\Data\A11\20170305\PFC\Bfsgrad1';
recDate{5} = '05032017';
fileID{5} = '20170305';
subjID{5} = 'Anton';
directories.taskdirPFC{6} = 'E:\Data\A11\20170302\PFC\Bfsgrad1';
recDate{6} = '02032017';
fileID{6} = '20170302';
subjID{6} = 'Anton';

%% Compute Spectrograms, LFP event statistics and event correlations for all these datasets in one bigass loop

%params
domdur = [500 1000 1500 2000]; % in ms. Take half for samples in LFP;
bands = {'delta','theta','low','beta'};
durs.switch = 250;
params.conditions = 8;
params.elecs = 96;

for iDataset = 1:nDatasets
    
    % Call relevant files
    eventfile = 'finalevents_audio.mat';
    spikeFile = 'jMUSpikesByTime.mat';
    qnxfile = [subjID{iDataset} '_' recDate{iDataset} '_' 'Bfsgrad1' '.dgz'];
    
    % Load relevant LFPs and Spikes
    cd(directories.taskdirPFC{iDataset})
    load('lfpByTrial.mat') % This is massive. We should find a better solution
    load('jMUSpikesByTime.mat')
    
    fprintf('Processing Dataset %d of : %d \n',iDataset,nDatasets)
    
    % Compute Statistics and Spectrograms
    
    sav_dir = [directories.taskdirPFC{iDataset} '\OKNStatistics'];
    mkdir(sav_dir)
    if strcmp(subjID{iDataset},'Hayo')==1
        monkID = 'H07';
    else
        monkID = 'A11';
    end
    
    % Compute Spectrograms
    
    for iDomDur = 1:length(domdur)
        cd(directories.taskdirPFC{iDataset})
        durs.domBehind = domdur(iDomDur);
        durs.domForward = domdur(iDomDur);
        %Spectrograms
        trialInformation = collectTrialInformationHayo(params,qnxfile,eventfile,'pfc','lfp');
        [lfpActivity,endDomTimes] = collectCleanDominancesLFPHayo(params,LFP,jMUspikes,durs);
        tag = [num2str(durs.domBehind/1000) 's_back_' num2str(durs.domForward/1000) 's'];
        [oknStatistics] = computeOKNStatistics(lfpActivity,endDomTimes,durs,sav_dir,tag);
    end
    
end
