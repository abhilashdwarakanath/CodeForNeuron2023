clear all
close all
clc

% This script does all the LFP analysis from the structures generated using
% Pipeline_Analysis

%% Enumerate the datasets

nDatasets = 6;
directories.taskdirPFC{1} = 'B:\H07\12-06-2016\PFC\Bfsgrad1';
recDate{1} = '12062016';
fileID{1} = '12-06-2016';
subjID{1} = 'Hayo';
directories.taskdirPFC{2} = 'B:\H07\13-07-2016\PFC\Bfsgrad1';
recDate{2} = '13072016';
fileID{2} = '13-07-2016';
subjID{2} = 'Hayo';
directories.taskdirPFC{3} = 'B:\H07\20161019\PFC\Bfsgrad1';
recDate{3} = '19102016';
fileID{3} = '20161019';
subjID{3} = 'Hayo';
directories.taskdirPFC{4} = 'B:\H07\20161025\PFC\Bfsgrad1';
recDate{4} = '25102016';
fileID{4} = '20161025';
subjID{4} = 'Hayo';
directories.taskdirPFC{5} = 'B:\A11\20170305\PFC\Bfsgrad1';
recDate{5} = '05032017';
fileID{5} = '20170305';
subjID{5} = 'Anton';
directories.taskdirPFC{6} = 'B:\A11\20170302\PFC\Bfsgrad1';
recDate{6} = '02032017';
fileID{6} = '20170302';
subjID{6} = 'Anton';

%% Compute Spectrograms, LFP event statistics and event correlations for all these datasets in one bigass loop

%params
domdur = [500]; % in ms. Take half for samples in LFP;
bands = {'beta'};
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
    
    sav_dir = [directories.taskdirPFC{iDataset} '\LFPPiecemeals'];
    mkdir(sav_dir)
    if strcmp(subjID{iDataset},'Hayo')==1
        monkID = 'H07';
    else
        monkID = 'A11';
    end
    
    % Compute Spectrograms
    
%     for iDomDur = 1:length(domdur)
%         cd(directories.taskdirPFC{iDataset})
%         durs.domBehind = domdur(iDomDur);
%         durs.domForward = domdur(iDomDur);
%         durs.switch = domdur(iDomDur);
%         %Spectrograms
%         trialInformation = collectTrialInformationHayo(params,qnxfile,eventfile,'pfc','lfp');
%         [lfpActivity] = collectPiecemealLFP(params,LFP,jMUspikes,durs);
%         nTransitions_BR_90 = 0;
%         nTransitions_BR_270 = 0;
%         for iConds = 1:length(lfpActivity.validOKN.BR.data.dom90)
%             nTransitions_BR_90 = nTransitions_BR_90+length(lfpActivity.validOKN.BR.data.dom90{iConds});
%         end
%         for iConds = 1:length(lfpActivity.validOKN.BR.data.dom270)
%             nTransitions_BR_270 = nTransitions_BR_270+length(lfpActivity.validOKN.BR.data.dom270{iConds});
%         end
%         save(['nTransitionsPiecemeals_' num2str(domdur(iDomDur)/1000) 's.mat'],'nTransitions_BR_90','nTransitions_BR_270','-v7.3');
%         tag = [num2str(durs.domBehind/1000) 's_back_' num2str(durs.domForward/1000) 's'];
%         [~] = computePiecemealSpecgrams(lfpActivity,sav_dir,durs,tag);
%         %[~] = computePiecemealSpecgramsTbyT(lfpActivity,sav_dir,durs,tag);
%         %lfpEventTriggeredTracesPiecemeal(lfpActivity, durs, sav_dir);
%         %lfpEventTriggeredTracesPiecemealPLV(lfpActivity, durs, sav_dir)
%     end
    
    % Statistics
    
    sav_dir = [directories.taskdirPFC{iDataset} '\LFPStatisticsPiecemeals'];
    mkdir(sav_dir)
    if strcmp(subjID{iDataset},'Hayo')==1
        monkID = 'H07';
    else
        monkID = 'A11';
    end
    
    for iBand = 1:length(bands)
        fprintf('Computing LFP Statistics for band %d of : %d \n',iBand,length(bands))
        fprintf([bands{iBand} '\n'])
        for iDomDur = 1:length(domdur)
            fprintf('Computing LFP Statistics for CleanDom Duration %d of : %d \n',iDomDur,length(domdur)-1)
            fprintf([num2str(domdur(iDomDur)) 'ms\n'])
            % Collect clean LFPs
            cd(directories.taskdirPFC{iDataset})
            durs.domForward = 500;
            durs.domBehind = 500;
            durs.switch = domdur(iDomDur);
            trialInformation = collectTrialInformationMM(params,qnxfile,eventfile,'pfc','lfp');
            [lfpActivity] = collectPiecemealLFP(params,LFP,jMUspikes,durs);
            % Characterisation
            lfpEventStatisticsPiecemeal(lfpActivity,durs, sav_dir,'bfsgrad',bands{iBand});
            lfpEventRateInTimePiecemealsAllChans(lfpActivity, durs, sav_dir,'bfsgrad',bands{iBand})
        end
    end
    
    
end
