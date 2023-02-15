clear all
clc
close all

% This is the pipeline that post-processes the data collected by the AG
% Fanis from 2 Utah Arrays.

% Abhilash Dwarakanath, MPI for biological cybernetics, April 2017.

%% Specify directories

subjectName = 'Anton';% Only change monkey name here

if strcmp(subjectName,'Anton')==1
    
    directories.recording = 'B:\A11\20170302';
    recDate = '02032017';
    subjID = 'Anton';
    
    cd(directories.recording) % Enter the data directory here % On Node 2
    if exist('PFC') ~=7 && exist('QNX') ~=7
        mkdir PFC
        mkdir QNX
    end
    
    directories.PFC = [directories.recording '\PFC'];
    directories.QNX = [directories.recording '\QNX'];
else
    directories.recording = 'B:\H07\20161022';
    recDate = '22102016';
    subjID = 'Hayo';
    cd(directories.recording) % Enter the data directory here % On Node 2
    if exist('PFC') ~=7 && exist('PPC') ~=7 && exist('QNX') ~=7
        mkdir PFC
        mkdir PPC
        mkdir QNX
    end
    
    directories.PFC = [directories.recording '\PFC'];
    directories.PPC = [directories.recording '\PPC'];
    directories.QNX = [directories.recording '\QNX'];
    
end

%% Get the event file and parse it

% NEV files
nevfiles = dir( [directories.recording '/*.nev']); %get files matching pattern
nevfiles = {nevfiles.name};

if strcmp(subjectName,'Hayo')==1
    
    cd(directories.recording)
    PFCevents = openNEV(['./' nevfiles{1}]);
    PPCevents = openNEV(['./' nevfiles{2}]);
    
else
    
    cd(directories.recording)
    PFCevents = openNEV(['./' nevfiles{1}]);
    
end


% %% Get the event file and parse it
%
% % NEV files
% nevfiles = dir( [directories.recording '/*.nev']); %get files matching pattern
% nevfiles = {nevfiles.name};
%
% if strcmp(subjectName,'Hayo')==1
%
%     cd(directories.recording)
%     PFCevents = openNEV(['./' nevfiles{1}]);
%     PPCevents = openNEV(['./' nevfiles{2}]);
%
% else
%     cd(directories.recording)
%     if ~exist([nevfiles{1}(1:end-4) '.mat'],'file')
%
%         PFCevents = openNEV(['./' nevfiles{1}]);
%
%     else
%         load([nevfiles{1}(1:end-4) '.mat']);
%         PFCevents = NEV;
%     end
%
% end

%% Get the events for a particular experiment

%tagsNeeded = [{'Bfsgrad1'}, {'Bfsgrad1 end'}]; % INPUT THE REQUIRED TAG FROM THE ABOVE LIST HERE!! CHECK THE EXCEL LOGFILE BEFORE YOU INPUT THE TAGS
%tagsNeeded = [{'Bfsnatur1'}, {'Bfsnatur1 end'}];
%tagsNeeded = [{'Spontaneous'}, {'Spontaneous end'}]; % INPUT THE REQUIRED TAG FROM THE ABOVE LIST HERE!! CHECK THE EXCEL LOGFILE BEFORE YOU INPUT THE TAGS
%tagsNeeded = [{'OriDisk1'}, {'OriDisk1 end'}];
tagsNeeded = [{'OriDiskf'}, {'OriDiskfixoff1 end'}];
%tagsNeeded = [{'AuSeq1_4'}, {'AuSeq1_4 end'}];
%tagsNeeded = [{'ImRivRE75_1'}, {'ImRivRE75_1 end'}];
%tagsNeeded = [{'BfsgradFs2FixOn_1'}, {'BfsgradFs2FixOn_1 end'}];

if strcmp(subjectName,'Hayo')==1
    
    % Create a specific task directory
    directories.taskdirPFC = [directories.PFC '\' tagsNeeded{1}];
    directories.taskdirPPC = [directories.PPC '\' tagsNeeded{1}];
    cd(directories.PFC)
    if ~exist(directories.taskdirPFC)
        mkdir(directories.taskdirPFC)
    end
    cd(directories.PPC)
    if ~exist(directories.taskdirPPC)
        mkdir(directories.taskdirPPC)
    end
else
    % Create a specific task directory
    directories.taskdirPFC = [directories.PFC '\' tagsNeeded{1}];
    
    cd(directories.PFC)
    if ~exist(directories.taskdirPFC)
        mkdir(directories.taskdirPFC)
    end
    
end

%Use ismember to get indices of the epochs that did happen in the experiment
expEpochTags = cellstr(PFCevents.Data.Comments.Text);

[ispresent,epochLogical] = ismember(tagsNeeded,expEpochTags);

if epochLogical(2)==0
    epochLogical(2) = epochLogical(1)+1;
end

% DUE TO A RECORDING BUG; WE WILL HAVE TO CHECK THIS AND CHANGE THE INDICES
% MANUALLY HERE. CONSULT THE LAB RECORD IF UNSURE!!!!!

epochinds = PFCevents.Data.Comments.TimeStamp(epochLogical);

if strcmp(tagsNeeded{1},'OriDiskf')
    tagsNeeded = [{'OriDiskfixoff1'}, {'OriDiskfixoff1 end'}];
    %tagsNeeded = [{'OriDiskFixOn1'}, {'OriDiskFixOn1 end'}];
end

%% Define parameters

global params % use global variable params so that PARFOR doesn't get cranky

% General params
params.fs = 30000; % in Hz
params.chans = 128;
params.elecs = 96;
params.nAnalogIns = 16;
params.clean = 0; % In samples. Usually for NSP1 this is the "Timestamp" value.

% For LFP stuff

params.targetFreq = 500;
params.trialLength = 10; % in seconds!! This is either the length of your trial from the DGZ file or the pieces you use to compute the coherence
params.tapers = 5/2;
params.decimCoeff = params.fs/params.targetFreq; % For reducing sampling rate
params.decimFacs = computeDecimationFactors(params.decimCoeff);

% For spike stuff

params.passband = [600 3000]/(params.fs/2);
%params.passband = [600 5800]/(params.fs/2);
params.stopband = [400 6000]/(params.fs/2);
params.passripple = 0.002; % Maximum allowed dB of the pass-band ripple
params.stopatten = 60;
params.minSpkNum = 4; % Currently this is set to number of spikes. We can later threshold it in Hz
params.numFeatures = 3; % THIS HAS TO BE LESS THAN OR EQUAL TO THE NUMBER OF SPIKES!!!!
params.stdmin = -5; % First number is for negative thresholded channels
params.stdmax = -50;                     % maximum threshold for detection
params.refractory = 0.5; %ms
params.beforeSpike = 0.5;
params.afterSpike = 1;
params.min_spk = 1;
params.filterorder = 2; % We will use a 2nd order Butterworth filter.
params.offs = 4;
params.thrsh = 60;

% For PSTH stuff
params.kerneltype = 'alpha'; %
params.psthbinsize = 0.025; % in SECONDS please
params.conditions = 8;

% for correlations stuff
params.binsize = 0.015; % in SECONDS please
params.maxlag = 0.5; % in SECONDS

max_memo_GB = 8;
notchfilter = 0;

%% Load electrode maps

if strcmp(subjectName,'Hayo')==1
    
    HayoPFCmap=[NaN	88	78	68	58	48	38	28	18	NaN
        96	87	77	67	57	47	37	27	17	8
        95	86	76	66	56	46	36	26	16	7
        94	85	75	65	55	45	35	25	15	6
        93	84	74	64	54	44	34	24	14	5
        92	83	73	63	53	43	33	23	13	4
        91	82	72	62	52	42	32	22	12	3
        90	81	71	61	51	41	31	21	11	2
        89	80	70	60	50	40	30	20	10	1
        NaN	79	69	59	49	39	29	19	9	NaN];
    
    HayoPPCmap = [NaN    88    78    68    58    48    38    28    18   8
        96    87    77    67    57    47    37    27    17     NaN
        95    86    76    66    56    46    36    26    16     7
        94    85    75    65    55    45    35    25    15     6
        93    84    74    64    54    44    34    24    14     5
        92    83    73    63    53    43    33    23    13     4
        91    82    72    62    52    42    32    22    12     3
        90    81    71    61    51    41    31    21    11     2
        89    80    70    60    50    40    30    20    10     1
        NaN    79    69    59    49    39    29    19     9   NaN];
    
else
    
    AntonPFCmap = [NaN    88    78    68    58    48    38    28    18   8
        96    87    77    67    57    47    37    27    17     7
        95    86    76    66    56    46    36    26    16     6
        94    85    75    65    55    45    35    25    15     5
        93    84    74    64    54    44    34    24    14     NaN
        92    83    73    63    53    43    33    23    13     4
        91    82    72    62    52    42    32    22    12     3
        90    81    71    61    51    41    31    21    11     2
        89    80    70    60    50    40    30    20    10     1
        NaN    79    69    59    49    39    29    19     9   NaN];
    
end

%% Load QNX files and setup analysis

% Get the starting and ending timestamps of the particular chosen
% experiment

ss = epochinds(1);
es = epochinds(2);

% Create folders and copy QNX files

if strcmp(subjectName,'Hayo')==1
    
    directories.parentQnx = 'L:\projects\AbhilashD\Recordings\QNX\Hayo';
    qnxfile = [subjID '_' recDate '_' tagsNeeded{1} '.dgz'];
    cd(directories.QNX)
    copyfile([directories.parentQnx '\' qnxfile], qnxfile);
    cd(directories.taskdirPFC)
    copyfile([directories.parentQnx '\' qnxfile], qnxfile);
    cd(directories.taskdirPPC)
    copyfile([directories.parentQnx '\' qnxfile], qnxfile);
    
    % mapping
    
    cd('B:\H07')
    load('H07ElectrodeInfo.mat')
    cd(directories.recording)
    
else
    
    directories.parentQnx = 'L:\projects\AbhilashD\Recordings\QNX\Anton';
    qnxfile = [subjID '_' recDate '_' tagsNeeded{1} '.dgz'];
    cd(directories.QNX)
    copyfile([directories.parentQnx '\' qnxfile], qnxfile);
    cd(directories.taskdirPFC)
    copyfile([directories.parentQnx '\' qnxfile], qnxfile);
    
    % mapping
    
    cd('B:\A11')
    load('A11ElectrodeInfo.mat')
    cd(directories.recording)
    
end

%% Read and decimate the data for Travelling-wave analysis if Spontaneous. Make a different Pipeline for this!!!

if tagsNeeded{1} == 'Spontaneous' || tagsNeeded{1} == 'spontaneous'
    
    if strcmp(subjectName,'Hayo')==1
        
        if exist('epochinds','var')
            
            espfc = endSampPFC;
            startSamp = epochinds(1);
            endSampPFC = epochinds(2);
            
        end
        
        
        cd(directories.PPC)
        for chans = 1:params.elecs
            data(chans,:)=readDecimate(params,chans,startSamp,endSampPPC);
        end
        
        cd(directories.PPC)
        LFP.data = data; clear data;
        LFP.map = HayoPFCmap;
        LFP.dx = 1/(params.fs/60);
        LFP.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        LFP.date = directories.recording(8:end);
        LFP.subject = 'H07';
        LFP.site = 'PPC';
        LFP.expttype = 'spontaneous';
        LFP.rawdatadir = directories.recording;
        save('LFPPFC.mat','LFP','-v7.3');
        cd(directories.recording)
        
    else
        
        cd(directories.PFC)
        parfor chans = 1:params.elecs
            data(chans,:)=readDecimate(params,chans,startSamp,endSampPFC);
        end
        
        cd(directories.PFC)
        LFP.data = data; clear data;
        LFP.map = AntonPFCmap;
        LFP.dx = 1/(params.fs/4);
        LFP.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        LFP.date = directories.recording(8:end);
        LFP.subject = 'A11';
        LFP.site = 'PFC';
        LFP.expttype = 'spontaneous';
        LFP.rawdatadir = directories.recording;
        save('LFPPFC.mat','LFP','-v7.3');
        cd(directories.recording)
        
    end
    
else
    
    fprintf('Not a spontaneous recording, moving on....\n');
    
end

%% Get SU Spikes Spontaneous

% ONLY EXTRACTION OF SORTED SU SPIKES HAS BEEN IMPLEMENTED. Copy the
% contents of extractRSRecSUKK.m to a new function-m file and just remove
% the looping over units. Look below in this script where we collect
% jMUSpikes and use that function to adapt the looping. E.g. -
% extractTrialsRecjMUHayo.m

if strcmp(subjectName,'Hayo')==1
    
    cd(directories.taskdirPFC)
    if exist('SUSpikesRS.mat')==0
        
        if exist('spikesPFC') == 0
            
            cd(directories.PFC)
            load('PFCspikes.mat')
            
        end
        
        for chans = 1:params.elecs
            cd([directories.PFC])
            fprintf('Extracting trials in Channel : %d \n',chans);
            [SUspikesByTime{chans}] = extractRSRecSUKK(spikesPFC{chans},chans);
        end
        
        cd(directories.taskdirPFC)
        SUspikes.data = SUspikesByTime; clear data;
        SUspikes.map = HayoPFCmap; % We need to change this
        SUspikes.dx = 1/(params.fs);
        SUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        SUspikes.date = directories.recording(8:end);
        SUspikes.subject = 'H07';
        SUspikes.site = 'PFC';
        SUspikes.chanElecs = H07.PFC;
        SUspikes.expttype = tagsNeeded{1}(1:end-1);
        SUspikes.rawdatadir = directories.recording;
        SUspikes.params = params;
        save('SUSpikesRS.mat','SUspikes','-v7.3');
        
    else
        fprintf('SUSpikesByTime has already been created. Moving on...\n')
        
    end
    clear SUspikes;
    cd(directories.recording)
    
    cd(directories.taskdirPPC)
    if exist('SUSpikesByTime.mat')==0
        
        if exist('spikesPPC') == 0
            
            cd(directories.PPC)
            load('PPCspikes.mat')
            
        end
        
        for chans = 1:params.elecs
            cd([directories.PPC])
            fprintf('Extracting trials in Channel : %d \n',chans);
            [SUspikesByTime{chans}] = extractRSRecSUKK(spikesPPC{chans},chans);
        end
        
        cd(directories.taskdirPPC)
        SUspikes.data = SUspikesByTime; clear data;
        SUspikes.map = HayoPPCmap; % We need to change this
        SUspikes.dx = 1/(params.fs);
        SUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        SUspikes.date = directories.recording(8:end);
        SUspikes.subject = 'H07';
        SUspikes.site = 'PPC';
        SUspikes.chanElecs = H07.PPC;
        SUspikes.expttype = tagsNeeded{1}(1:end-1);
        SUspikes.rawdatadir = directories.recording;
        SUspikes.params = params;
        save('SUSpikesByTime.mat','SUspikes','-v7.3');
    else
        fprintf('SUSpikesByTime has already been created. Moving on...\n')
        
    end
    clear SUSpikes;
    cd(directories.recording)
    
else
    
    cd(directories.taskdirPFC)
    if exist('SUSpikesByTime.mat')==0
        
        if exist('spikesPFC') == 0
            
            cd(directories.PFC)
            load('PFCspikes.mat')
            
        end
        
        for chans = 1:params.elecs
            cd([directories.PFC])
            fprintf('Extracting trials in Channel : %d \n',chans);
            [SUspikesByTime{chans}] = extractRSRecSUKK(spikesPFC{chans},chans);
        end
        
        cd(directories.taskdirPFC)
        SUspikes.data = SUspikesByTime; clear data;
        SUspikes.map = AntonPFCmap; % We need to change this
        SUspikes.dx = 1/(params.fs);
        SUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        SUspikes.date = directories.recording(8:end);
        SUspikes.subject = 'A11';
        SUspikes.site = 'PFC';
        SUspikes.chanElecs = A11.PFC;
        SUspikes.expttype = tagsNeeded{1}(1:end-1);
        SUspikes.rawdatadir = directories.recording;
        SUspikes.params = params;
        save('SUSpikesRS.mat','SUspikes','-v7.3');
        
    else
        fprintf('SUSpikesByTime has already been created. Moving on...\n')
        
    end
    clear SUspikes;
    cd(directories.recording)
    
end

%% Check synchronisation and get delays

if strcmp(subjectName,'Hayo')==1
    
    cd(directories.taskdirPFC);
    [trial_lengths,delay] = checkTrialSync(PFCevents,PPCevents,qnxfile,ss,es);
    cd(directories.taskdirPPC)
    save('delay.mat','delay','-v7.3')
    cd(directories.taskdirPFC)
else
    delay=0;
    fprintf('Only 1 active array; no need to check synchornisation...\n')
    
end

%% Get jMU spikes

if strcmp(subjectName,'Hayo')==1
    
    cd(directories.taskdirPFC)
    if exist('jMUSpikesByTime.mat')==0
        cd(directories.taskdirPFC)
        copyfile([directories.PFC '\NSX130.NC5'],'NSX130.NC5'); %Y coord EM file
        copyfile([directories.PFC '\NSX129.NC5'],'NSX129.NC5'); % X coord EM file
        copyfile([directories.PFC '\Bfsgrad1\NSX_TimeStamps.mat'],'NSX_TimeStamps.mat');
        copyfile('L:\projects\AbhilashD\BlackrockPreProcessingCodes\switchTypes.txt','switchTypes.txt')
        
        if exist('spikesPFC') == 0
            
            cd(directories.PFC)
            load('PFCspikes.mat')
            
        end
        
        cd(directories.taskdirPFC)
        parfor chans = 1:params.elecs %This can be changed to parfar
            fprintf('Extracting trials in Channel : %d \n',chans);
            cd(directories.PFC)
            [jMUspikesByTime{chans},times{chans}] = extractTrialsRecjMUHayo(PFCevents,spikesPFC{chans},qnxfile,ss,es,chans,directories,'jMU','pfc',delay);
        end
        
        cd(directories.taskdirPFC)
        jMUspikes.data = jMUspikesByTime; clear data;
        jMUspikes.map = HayoPFCmap;
        jMUspikes.dx = 1; %(1ms)
        jMUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude.
        jMUspikes.date = directories.recording(8:end);
        jMUspikes.subject = 'H07';
        jMUspikes.site = 'PFC';
        jMUspikes.chanElecs = H07.PFC;
        jMUspikes.expttype = tagsNeeded{1}(1:end-1);
        jMUspikes.rawdatadir = directories.recording;
        jMUspikes.params = params;
        jMUspikes.events = times{1};
        save('jMUSpikesByTime.mat','jMUspikes','-v7.3');
        cd(directories.recording)
        
        clear jMUspikes; clear jMUspikesByTime;
        
    else
        fprintf('jMUSpikesByTime has already been created. Moving on...\n')
        
    end
    % For PPC
    
    cd(directories.taskdirPPC)
    
    if exist('jMUSpikesByTime.mat')==0
        cd(directories.taskdirPPC)
        copyfile([directories.PFC '\Bfsgrad1\finalevents_audio.mat'],'finalevents_audio.mat');
        copyfile([directories.PFC '\NSX130.NC5'],'NSX130.NC5');
        copyfile([directories.PFC '\NSX129.NC5'],'NSX129.NC5');
        copyfile([directories.PFC '\Bfsgrad1\NSX_TimeStamps.mat'],'NSX_TimeStamps.mat');
        copyfile('L:\projects\AbhilashD\BlackrockPreProcessingCodes\switchTypes.txt','switchTypes.txt')
        
        if exist('spikesPPC') == 0
            
            cd(directories.PPC)
            load('PPCspikes.mat')
            
        end
        
        cd(directories.taskdirPPC)
        parfor chans = 1:params.elecs
            fprintf('Extracting trials in Channel : %d \n',chans);
            cd(directories.PPC)
            [jMUspikesByTime{chans}] = extractTrialsRecjMUHayo(PFCevents,spikesPPC{chans},qnxfile,ss,es,chans,directories,'jMU','ppc',delay);
        end
        
        cd(directories.taskdirPPC)
        jMUspikes.data = jMUspikesByTime; clear data;
        jMUspikes.map = HayoPPCmap;
        jMUspikes.dx = 1;
        jMUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        jMUspikes.date = directories.recording(8:end);
        jMUspikes.subject = 'H07';
        jMUspikes.site = 'PPC';
        jMUspikes.chanElecs = H07.PPC;
        jMUspikes.expttype = tagsNeeded{1}(1:end-1);
        jMUspikes.rawdatadir = directories.recording;
        jMUspikes.params = params;
        save('jMUSpikesByTime.mat','jMUspikes','-v7.3');
        cd(directories.recording)
        
        clear jMUspikes;clear jMUspikesByTime;
        
    else
        fprintf('jMUSpikesByTime has already been created. Moving on...\n')
        
    end
    
else
    
    cd(directories.taskdirPFC)
    
    if exist('jMUSpikesByTime.mat')==0
        cd(directories.taskdirPFC)
        copyfile([directories.PFC '\NSX130.NC5'],'NSX130.NC5');
        copyfile([directories.PFC '\NSX129.NC5'],'NSX129.NC5');
        %copyfile([directories.PFC '\NSX_TimeStamps.mat'],'NSX_TimeStamps.mat');
        copyfile([directories.PFC '\Bfsgrad1\NSX_TimeStamps.mat'],'NSX_TimeStamps.mat');
        copyfile('L:\projects\AbhilashD\BlackrockPreProcessingCodes\switchTypes.txt','switchTypes.txt')
        
        if exist('spikesPFC') == 0
            
            cd(directories.PFC)
            load('PFCspikes.mat')
            
        end
        
        cd(directories.taskdirPFC)
        delay = 0;
        parfor chans = 1:params.elecs
            fprintf('Extracting trials in Channel : %d \n',chans);
            cd(directories.PFC)
            [jMUspikesByTime{chans},times{chans}] = extractTrialsRecjMUHayo(PFCevents,spikesPFC{chans},qnxfile,ss,es,chans,directories,'jMU','pfc',delay);
        end
        
        cd(directories.taskdirPFC)
        jMUspikes.data = jMUspikesByTime; clear data;
        jMUspikes.map = AntonPFCmap;
        jMUspikes.dx = 1;
        jMUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        jMUspikes.date = directories.recording(8:end);
        jMUspikes.subject = 'A11';
        jMUspikes.site = 'PFC';
        jMUspikes.chanElecs = A11.PFC;
        jMUspikes.expttype = tagsNeeded{1}(1:end-1);
        jMUspikes.rawdatadir = directories.recording;
        jMUspikes.params = params;
        jMUspikes.events = times{1};
        save('jMUSpikesByTime.mat','jMUspikes','-v7.3');
        cd(directories.recording)
        
        clear jMUspikes;
        
    else
        fprintf('jMUSpikesByTime has already been created. Moving on...\n')
        
    end
    
end

%% Get sorted SU spikes

if strcmp(subjectName,'Hayo')==1
    
    cd(directories.taskdirPFC)
    if exist('SUSpikesByTime.mat')==0
        
        if exist('spikesPFC') == 0
            
            cd(directories.PFC)
            load('PFCspikes.mat')
            
        end
        
        parfor chans = 1:params.elecs
            cd([directories.PFC])
            fprintf('Extracting trials in Channel : %d \n',chans);
            [SUspikesByTime{chans}] = extractTrialsRecSUHayo(PFCevents,spikesPFC{chans},qnxfile,ss,es,chans,directories,'pfc',delay);
        end
        
        cd(directories.taskdirPFC)
        SUspikes.data = SUspikesByTime; clear data;
        SUspikes.map = HayoPFCmap; % We need to change this
        SUspikes.dx = 1/(params.fs);
        SUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        SUspikes.date = directories.recording(8:end);
        SUspikes.subject = 'H07';
        SUspikes.site = 'PFC';
        SUspikes.chanElecs = H07.PFC;
        SUspikes.expttype = tagsNeeded{1}(1:end-1);
        SUspikes.rawdatadir = directories.recording;
        SUspikes.params = params;
        save('SUSpikesByTime.mat','SUspikes','-v7.3');
        
    else
        fprintf('SUSpikesByTime has already been created. Moving on...\n')
        
    end
    clear SUspikes;
    cd(directories.recording)
    
    cd(directories.taskdirPPC)
    if exist('SUSpikesByTime.mat')==0
        
        if exist('spikesPPC') == 0
            
            cd(directories.PPC)
            load('PPCspikes.mat')
            
        end
        
        for chans = 1:params.elecs
            cd([directories.PPC])
            fprintf('Extracting trials in Channel : %d \n',chans);
            [SUspikesByTime{chans}] = extractTrialsRecSUHayo(PFCevents,spikesPPC{chans},qnxfile,ss,es,chans,directories,'ppc',delay);
        end
        
        cd(directories.taskdirPPC)
        SUspikes.data = SUspikesByTime; clear data;
        SUspikes.map = HayoPPCmap; % We need to change this
        SUspikes.dx = 1/(params.fs);
        SUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        SUspikes.date = directories.recording(8:end);
        SUspikes.subject = 'H07';
        SUspikes.site = 'PPC';
        SUspikes.chanElecs = H07.PPC;
        SUspikes.expttype = tagsNeeded{1}(1:end-1);
        SUspikes.rawdatadir = directories.recording;
        SUspikes.params = params;
        save('SUSpikesByTime.mat','SUspikes','-v7.3');
    else
        fprintf('SUSpikesByTime has already been created. Moving on...\n')
        
    end
    clear SUSpikes;
    cd(directories.recording)
    
else
    
    cd(directories.taskdirPFC)
    if exist('SUSpikesByTime.mat')==0
        
        if exist('spikesPFC') == 0
            
            cd(directories.PFC)
            load('PFCspikes.mat')
            
        end
        
        for chans = 1:params.elecs
            cd([directories.PFC])
            fprintf('Extracting trials in Channel : %d \n',chans);
            [SUspikesByTime{chans}] = extractTrialsRecSUHayo(PFCevents,spikesPFC{chans},qnxfile,ss,es,chans,directories,'pfc',delay);
        end
        
        cd(directories.taskdirPFC)
        SUspikes.data = SUspikesByTime; clear data;
        SUspikes.map = AntonPFCmap; % We need to change this
        SUspikes.dx = 1;
        SUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        SUspikes.date = directories.recording(8:end);
        SUspikes.subject = 'A11';
        SUspikes.site = 'PFC';
        SUspikes.chanElecs = A11.PFC;
        SUspikes.expttype = tagsNeeded{1}(1:end-1);
        SUspikes.rawdatadir = directories.recording;
        SUspikes.params = params;
        save('SUSpikesByTime.mat','SUspikes','-v7.3');
        
    else
        fprintf('SUSpikesByTime has already been created. Moving on...\n')
        
    end
    clear SUspikes;
    cd(directories.recording)
    
end

%% Get SU and MU spikes separate but in the same structure

if strcmp(subjectName,'Hayo')==1
    
    cd(directories.taskdirPFC)
    if exist('SUMUSpikesByTime.mat')==0
        
        if exist('spikesPFC') == 0
            
            cd(directories.PFC)
            load('PFCspikes.mat')
            
        end
        
        for chans = 1:params.elecs
            cd([directories.PFC])
            fprintf('Extracting trials in Channel : %d \n',chans);
            [SUspikesByTime{chans}] = extractTrialsRecSUMUHayo(PFCevents,spikesPFC{chans},qnxfile,ss,es,chans,directories,'pfc',delay);
        end
        
        cd(directories.taskdirPFC)
        SUMUspikes.data = SUspikesByTime; clear data;
        SUMUspikes.map = HayoPFCmap; % We need to change this
        SUMUspikes.dx = 1/(params.fs);
        SUMUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        SUMUspikes.date = directories.recording(8:end);
        SUMUspikes.subject = 'H07';
        SUMUspikes.site = 'PFC';
        SUMUspikes.chanElecs = H07.PFC;
        SUMUspikes.expttype = tagsNeeded{1}(1:end-1);
        SUMUspikes.rawdatadir = directories.recording;
        SUMUspikes.params = params;
        save('SUMUSpikesByTime.mat','SUMUspikes','-v7.3');
        
    else
        fprintf('SUMUSpikesByTime has already been created. Moving on...\n')
        
    end
    clear SUMUspikes;
    cd(directories.recording)
    
%     cd(directories.taskdirPPC)
%     if exist('SUMUSpikesByTime.mat')==0
%         
%         if exist('spikesPPC') == 0
%             
%             cd(directories.PPC)
%             load('PPCspikes.mat')
%             
%         end
%         
%         for chans = 1:params.elecs
%             cd([directories.PPC])
%             fprintf('Extracting trials in Channel : %d \n',chans);
%             [SUMUspikesByTime{chans}] = extractTrialsRecSUMUHayo(PFCevents,spikesPPC{chans},qnxfile,ss,es,chans,directories,'ppc',delay);
%         end
%         
%         cd(directories.taskdirPPC)
%         SUMUspikes.data = SUMUspikesByTime; clear data;
%         SUMUspikes.map = HayoPPCmap; % We need to change this
%         SUMUspikes.dx = 1/(params.fs);
%         SUMUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
%         SUMUspikes.date = directories.recording(8:end);
%         SUMUspikes.subject = 'H07';
%         SUMUspikes.site = 'PPC';
%         SUMUspikes.chanElecs = H07.PPC;
%         SUMUspikes.expttype = tagsNeeded{1}(1:end-1);
%         SUMUspikes.rawdatadir = directories.recording;
%         SUMUspikes.params = params;
%         save('SUMUSpikesByTime.mat','SUMUspikes','-v7.3');
%     else
%         fprintf('SUMUSpikesByTime has already been created. Moving on...\n')
%         
%     end
%     clear SUMUSpikes;
%     cd(directories.recording)
    
else
    
    cd(directories.taskdirPFC)
    if exist('SUMUSpikesByTime.mat')==0
        
        if exist('spikesPFC') == 0
            
            cd(directories.PFC)
            load('PFCspikes.mat')
            
        end
        
        parfor chans = 1:params.elecs
            cd([directories.PFC])
            fprintf('Extracting trials in Channel : %d \n',chans);
            [SUMUspikesByTime{chans}] = extractTrialsRecSUMUHayo(PFCevents,spikesPFC{chans},qnxfile,ss,es,chans,directories,'pfc',delay);
        end
        
        cd(directories.taskdirPFC)
        SUMUspikes.data = SUMUspikesByTime; clear data;
        SUMUspikes.map = AntonPFCmap; % We need to change this
        SUMUspikes.dx = 1;
        SUMUspikes.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        SUMUspikes.date = directories.recording(8:end);
        SUMUspikes.subject = 'A11';
        SUMUspikes.site = 'PFC';
        SUMUspikes.chanElecs = A11.PFC;
        SUMUspikes.expttype = tagsNeeded{1}(1:end-1);
        SUMUspikes.rawdatadir = directories.recording;
        SUMUspikes.params = params;
        save('SUMUSpikesByTime.mat','SUMUspikes','-v7.3');
        
    else
        fprintf('SUMUSpikesByTime has already been created. Moving on...\n')
        
    end
    clear SUMUspikes;
    cd(directories.recording)
    
end

%% Extract LFP

if strcmp(subjectName,'Hayo')==1
    
    % PFC
    cd(directories.taskdirPFC)
    if exist('lfpByTrial.mat')==0
        
        cd(directories.PFC)
        parfor chans = 1:params.elecs
            filenameraw = ['NSX' num2str(chans) '.NC5'];
            fprintf('Extracting trials in Channel : %d \n',chans);
            [lfp{chans},LFPtimes{chans}] = extractTrialsLFPHayo(params,filenameraw,PFCevents,qnxfile,ss,es,directories,'pfc',delay);
        end
        
        cd(directories.taskdirPFC)
        LFP.data = lfp; clear lfp;
        LFP.map = HayoPFCmap;
        LFP.dx = 1/500;
        LFP.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        LFP.date = directories.recording(8:end);
        LFP.subject = 'H07';
        LFP.site = 'PFC';
        LFP.chanElecs = H07.PFC;
        LFP.expttype = tagsNeeded{1}(1:end-1);
        LFP.events = LFPtimes{1};
        LFP.epochStart = ss;
        LFP.rawdatadir = directories.recording;
        save('lfpByTrial.mat','LFP','-v7.3');
        clear LFP
    else
        fprintf('LFPs by Trial have been extracted. Moving on...\n')
        clear lfp;
        cd(directories.recording)
    end
    
    % PPC
    
    cd(directories.taskdirPPC)
    if exist('lfpByTrial.mat')==0
        
        cd(directories.PPC)
        for chans = 1:params.elecs
            filenameraw = ['NSX' num2str(chans) '.NC5'];
            fprintf('Extracting trials in Channel : %d \n',chans);
            [lfp{chans},LFPtimes{1}] = extractTrialsLFPHayo(params,filenameraw,PFCevents,qnxfile,ss,es,directories,'ppc',delay);
        end
        
        cd(directories.taskdirPPC)
        LFP.data = lfp; clear lfp;
        LFP.map = HayoPPCmap;
        LFP.dx = 1/500;
        LFP.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        LFP.date = directories.recording(8:end);
        LFP.subject = 'H07';
        LFP.site = 'PPC';
        LFP.chanElecs = H07.PPC;
        LFP.events = LFPtimes{1};
        LFP.epochStart = ss;
        LFP.expttype = tagsNeeded{1}(1:end-1);
        LFP.rawdatadir = directories.recording;
        save('lfpByTrial.mat','LFP','-v7.3');
        clear LFP;
    else
        fprintf('LFPs by Trial have been extracted. Moving on...\n')
        clear lfp;
        cd(directories.recording)
    end
    
else
    
    % PFC
    
    cd(directories.taskdirPFC)
    if exist('lfpByTrial.mat')==0
        
        cd(directories.PFC)
        parfor chans = 1:params.elecs
            filenameraw = ['NSX' num2str(chans) '.NC5'];
            fprintf('Extracting trials in Channel : %d \n',chans);
            [lfp{chans},LFPtimes{chans}] = extractTrialsLFPHayo(params,filenameraw,PFCevents,qnxfile,ss,es,directories,'pfc',0);
        end
        
        cd(directories.taskdirPFC)
        LFP.data = lfp; clear lfp;
        LFP.map = AntonPFCmap;
        LFP.dx = 1/500;
        LFP.exclchan = zeros(1,params.elecs); % A 0 means don't exclude
        LFP.date = directories.recording(8:end);
        LFP.subject = 'A11';
        LFP.site = 'PFC';
        LFP.chanElecs = A11.PFC;
        LFP.expttype = tagsNeeded{1}(1:end-1);
        LFP.events = LFPtimes{1};
        LFP.epochStart = ss;
        LFP.rawdatadir = directories.recording;
        save('lfpByTrial.mat','LFP','-v7.3');
        clear LFP;
    else
        fprintf('LFPs by Trial have been extracted. Moving on...\n')
        clear lfp;
        cd(directories.recording)
    end
    
end

%% Extract Photodiode Signals


% PFC
cd(directories.taskdirPFC)
if exist('pdByTrial.mat')==0
    
    cd(directories.PFC)
    recfile{1} = 'NSX133.NC5';
    recfile{2} = 'NSX134.NC5';
    [pdByTrials,times] = extractTrialsPDSignal(params,recfile,PFCevents,qnxfile,ss,es,directories,'pfc',0);
    save('pdByTrials.mat','pdByTrials','-v7.3')
else
    
    disp('PD Signal already extracted!')
    
end


%% Prepare for Marking

if exist('jMUspikes') == 0
    
    cd(directories.taskdirPFC)
    load('jMUSpikesByTime.mat')
    
end

cd(directories.taskdirPFC)
copyfile([directories.PFC '\NSX130.NC5'],'NSX130.NC5');
copyfile('L:\projects\AbhilashD\BlackrockPreProcessingCodes\switchTypes.txt','switchTypes.txt')
events_samples = jMUspikes.events;
if iscell(events_samples)
    events_samples = events_samples{1};
end
if ~exist('finalevents_audio.mat')==1
    markEventsGUIFullRec(params,events_samples,tagsNeeded,directories);
else
    sprintf('Marking is done? Check!\n')
end
