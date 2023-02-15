clear all
close all
clc

% This script does all the LFP analysis from the structures generated using
% Pipeline_Analysis

%% Enumerate the datasets

nDatasets = 6;
directories.taskdirPFC{1} = 'B:\H07\12-06-2016\PFC\Bfsgrad1';
recDate{1} = '12062016';
fileID{1} = 'H07\12-06-2016';
subjID{1} = 'Hayo';
directories.taskdirPFC{2} = 'B:\H07\13-07-2016\PFC\Bfsgrad1';
recDate{2} = '13072016';
fileID{2} = 'H07\13-07-2016';
subjID{2} = 'Hayo';
directories.taskdirPFC{3} = 'B:\H07\20161019\PFC\Bfsgrad1';
recDate{3} = '19102016';
fileID{3} = 'H07\20161019';
subjID{3} = 'Hayo';
directories.taskdirPFC{4} = 'B:\H07\20161025\PFC\Bfsgrad1';
recDate{4} = '25102016';
fileID{4} = 'H07\20161025';
subjID{4} = 'Hayo';
directories.taskdirPFC{5} = 'B:\A11\20170305\PFC\Bfsgrad1';
recDate{5} = '05032017';
fileID{5} = 'A11\20170305';
subjID{5} = 'Anton';
directories.taskdirPFC{6} = 'B:\A11\20170302\PFC\Bfsgrad1';
recDate{6} = '02032017';
fileID{6} = 'A11\20170302';
subjID{6} = 'Anton';

%% Compute SFC

%params
domdur = [1000]; % in ms. Take half for samples in LFP;
durs.switch = 250;
bands = {'low','beta'};
params.conditions = 8;
params.elecs = 96;

for iDataset = 2:nDatasets
    
    % Call relevant files
    eventfile = 'finalevents_audio.mat';
    spikeFile = 'jMUDominancesByTime.mat';
    qnxfile = [subjID{iDataset} '_' recDate{iDataset} '_' 'Bfsgrad1' '.dgz'];
    
    % Load relevant LFPs and Spikes
    cd(directories.taskdirPFC{iDataset})
    load('lfpByTrial.mat') % This is massive. We should find a better solution
    load('jMUSpikesByTime.mat')
    load('SUMUSpikesByTime.mat')
    
    num_su = 0;
    SUspikes_t = struct;
    SUspikes_t = [];
    spk_flashDom_su = [];

    % collect all the data
    
    for iChan = 1:length(SUMUspikes.data)
        
        if isstruct(SUMUspikes.data{iChan})
            
            for iSu = 1:size(SUMUspikes.data{iChan}.spikesUnaligned,2)
                
                num_su = num_su + 1;
                
                SUspikes_t.data{num_su}.spikesUnaligned = SUMUspikes.data{iChan}.spikesUnaligned{iSu};
                spk_flashDom_su{num_su} = SUMUspikes.data{iChan}.spikesUnaligned{iSu};
                
                %     keep information about channel number and SU_idx
                SUspikes_t.Info(num_su,:) = [iChan iSu];
                SUspikes_t.data{num_su}.okn.trace = SUMUspikes.data{iChan}.okn.trace;
            end
            
        end
        
    end
    
    nSU = num_su;
    
    params.nSU = nSU;
    
    fprintf('Processing Dataset %d of : %d \n',iDataset,nDatasets)
    
    % Compute Statistics and Spectrograms
    
    sav_dir = [directories.taskdirPFC{iDataset} '\SFC'];
    mkdir(sav_dir)
    if strcmp(subjID{iDataset},'Hayo')==1
        monkID = 'H07';
    else
        monkID = 'A11';
    end
    
    
    % Compute Spectrograms
    
    for iDomDur = 1:length(domdur)
        cd(directories.taskdirPFC{iDataset})
        durs.domForward = domdur(iDomDur);
        durs.domBehind = domdur(iDomDur);
        trialInformation = collectTrialInformationMM(params,qnxfile,eventfile,'pfc','lfp');
        [lfpActivityMM] = collectCleanDominancesLFPMM(params,LFP,jMUspikes,durs);
        [spikingActivityMM] = collectCleanDominancesSpikesSUMM(params,SUspikes_t,durs);
        %Spectrograms
        tag = [num2str(durs.domBehind/1000) 's_back_' num2str(durs.domForward/1000) 's'];
        %lfpTransitionsTbyT_SFC(lfpActivityMM,spikingActivityMM,sav_dir,durs,fileID{iDataset});
        lfpTransitionsTbyT_SFCohgram(lfpActivityMM,spikingActivityMM,sav_dir,durs,fileID{iDataset});
        %lfpTransitionsTbyT_SFC_Individual(lfpActivity,spikingActivity,sav_dir,durs,fileID{iDataset});
        %eventTriggered_SFC(lfpActivity,spikingActivity,sav_dir,durs,fileID{iDataset});
        %lfpTransitionsTbyT_SFC_globalLFP_ensembleSpiking(lfpActivityMM,spikingActivityMM,sav_dir,durs,fileID{iDataset});
        %lfpTransitionsTbyT_SFCohgram_globalLFP_ensembleSpiking(lfpActivityMM,spikingActivityMM,sav_dir,durs,fileID{iDataset});
    end
    
end
