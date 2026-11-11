%% Concatenate LFP of social interaction test phases
% Project: Adriana PhD, Felipe Gomes
% 
% Pipeline:
%   Set paths and general parameters
%   Define data ID
%   Find TTLs and define SIT phases
%   Load, check, and concatenate LFP
%   Load, check, and concatenate ACC
%   Save SIT phases indices, ACC, and LFP
%  
% Instructions:
%   - Set once the correct paths and folders to your data and scripts
%   - For every recording, define the data ID (subject_group) 
%   - Run each section of the whole script per recording from the start
% 
% Authors: Danilo Benette Marques, Matheus Teixeira Rossignoli, Adriana
%   Souza
% Last update: 26-06-2025 15h47

% %% Set paths
clear all

%Set main folder
sitephysdir = ['C:\Users\User\Desktop\MESTRADO\ELETROFISIOLOGIA\ANALISE\MATLAB\ROTINAS\posdoutorado' ...
    '\parcerias\Felipe Gomes\Adriana'];
cd(sitephysdir)

%Add scripts to path
addpath(genpath(fullfile(sitephysdir, 'scripts')))

%Define saved data path
savedatadir = fullfile(sitephysdir,'data');

% %% General parameters
Fs = 500;
gain = 1/1000; %uV -> mV

%SIT phases durations
Nsitphases = 4; %number of phases
% Phase 1: habituation
sitphases_duration(1) = 30*60*Fs; 
% Phase 2: empty cage
sitphases_duration(2) = 10*60*Fs; 
% Phase 3: intruder
sitphases_duration(3) = 20*60*Fs; 
% Phase 4: post
sitphases_duration(4) = 30*60*Fs;


%% Define subject data ID and path
subjID = 'M1'; %subject ID: M#
sexID = 'M'; %sex ID: M | F
groupID = 'VEH'; %group ID: VAL | VEH

dataID = strcat(subjID,'_',sexID,'_',groupID)

%Select subject ephys data folder
filepath = uigetdir(sitephysdir,['Select ' dataID ' folder']);


%% Load and find TTL

labels = {'ADC1','_36.continuous'}; %labels to load

[ttl,t,info] = load_open_ephys_data_dbm(filepath,[],labels,gain,Fs);

Ns = length(ttl); %total number of samples


% %% Find TTL
idx_ttl = findpulses(ttl,4.9*Fs); %automatically find pulses >4.9s
Nttl = numel(idx_ttl);
fprintf('Found %d TTLs \n',Nttl)

%Check TTL
figure,
plot(t,ttl,'k')
    hold on,scatter(t(idx_ttl),ttl(idx_ttl),'sr','filled')
    xlabel('Tempo (s)')
    ylabel('mV')
    pause(0)

%Manual fixes of TTL locations
manualttl = input('Need to locate TTLs manually? Y/N or Enter: ','s');
if strcmp(manualttl,'Y')
    close

    sigapp = signaltsmarker(ttl,Fs); %App for manual locating and clearing of points
    sigapp.Markers{2,1} = idx_ttl;
    sigapp.Markers{3,1} = idx_ttl/Fs;
    sigapp.plotSignal; 

    uiwait(sigapp.SignalTimestampMarkerappUIFigure) %waits until close app

    idx_ttl = signalmarkers{2,1};
end

%% Get SIT phases
%Get start-stop indices of SIT phases
% start: TTLs
% stop: phase duration
clear idx_sitphases
for sitphase = 1:Nsitphases
    idx_sitphases(sitphase,1) = idx_ttl(sitphase);
    idx_sitphases(sitphase,2) = idx_ttl(sitphase) + sitphases_duration(sitphase)-1;
end
    
%OBS: teve um rato com 15min de empty e 15min social. preencher com Nan (se
%   passar do next ttl)

%% Load LFP

% %% Option 1) Load sequentially or keep missing (if missing channels)
%Region labels
labels_pfc = {'PFC'};
labels_vhpc = {'HipV'};
labels_nac = {'Nac'};
labels_vta = {'VTA'};
Nlfp = 4;

LFP = nan(Ns,Nlfp); %initialize with missing

try [lfp,t,info] = load_open_ephys_data_dbm(filepath,[],labels_pfc,gain,Fs); LFP(:,1) = lfp; 
catch warning('Did not load PFC'); end

try [lfp,t,info] = load_open_ephys_data_dbm(filepath,[],labels_vhpc,gain,Fs); LFP(:,2) = lfp; 
catch warning('Did not load vHPC'); end

try [lfp,t,info] = load_open_ephys_data_dbm(filepath,[],labels_nac,gain,Fs); LFP(:,3) = lfp; 
catch warning('Did not load NAc'); end

try [lfp,t,info] = load_open_ephys_data_dbm(filepath,[],labels_vta,gain,Fs); LFP(:,4) = lfp; 
catch warning('Did not load VTA'); end

% %% Option 2) Load all at once (if all folders have all recordings)
% labels = {'PFC',... %PFC
%     'dHPC',... %dorsal hippocampus
%     'vHPC',... %ventral hippocampus
%     'NAc',... %nucleus accumbens
%     'VTA'}; %VTA
% 
% [LFP,t,info] = load_open_ephys_data_dbm(filepath,[],labels,gain,Fs);

%Check LFP
figure
h=plotmultisignals(t,LFP'); set(h,'color','k')
hold on,xline(t(idx_sitphases(:,1)),'b--')
hold on,xline(t(idx_sitphases(:,2)),'r--')
title('LFP')
    xlabel('Tempo (s)')
    ylabel('mV')
    pause(0)

%% Concatenate LFP and time per SIT phase
% LFP_hab = LFP(idx_sitphases(1,1):idx_sitphases(1,2),:);
% LFP_cage = LFP(idx_sitphases(2,1):idx_sitphases(2,2),:);
% LFP_social = LFP(idx_sitphases(3,1):idx_sitphases(3,2),:);
% LFP_post = LFP(idx_sitphases(4,1):idx_sitphases(4,2),:);
% 
% t_hab = t(idx_sitphases(1,1):idx_sitphases(1,2),1);
% t_cage = t(idx_sitphases(2,1):idx_sitphases(2,2),1);
% t_social = t(idx_sitphases(3,1):idx_sitphases(3,2),1);
% t_post = t(idx_sitphases(4,1):idx_sitphases(4,2),1);

%% Concatenate LFP and time per SIT phase (with padding)
% OBS: actually does not need padding before
npad = 10*60*Fs;
idx_sitphases_pad = idx_sitphases + npad;

%fill sides with NaN to allow standardized concatenation
LFP = [nan(npad,Nlfp) ; LFP; nan(npad,Nlfp)]; 
t = [nan(npad,1) ; shiftdim(t); nan(npad,1)]; %some ephys load in column vectors (?! weird)

LFP_hab = LFP(idx_sitphases_pad(1,1):idx_sitphases_pad(1,2),:);
LFP_cage = LFP(idx_sitphases_pad(2,1):idx_sitphases_pad(2,2),:);
LFP_social = LFP(idx_sitphases_pad(3,1):idx_sitphases_pad(3,2),:);
LFP_post = LFP(idx_sitphases_pad(4,1):idx_sitphases_pad(4,2),:);

t_hab = t(idx_sitphases(1,1):idx_sitphases(1,2),1);
t_cage = t(idx_sitphases(2,1):idx_sitphases(2,2),1);
t_social = t(idx_sitphases(3,1):idx_sitphases(3,2),1);
t_post = t(idx_sitphases(4,1):idx_sitphases(4,2),1);

%% Ignore overlapping with next phase
%OBS: we can do it on the fly
% % Turn end of phase overlapping with start of next phase to NaN
% 
% noverlapphases = cell(3,1); %empty
% for idx_sitphase = 1:Nsitphases-1 %hab to social
% 
%     %If end of phase N is after start of phase N+1
%     if idx_sitphases(idx_sitphase,2)>idx_sitphases(idx_sitphase+1,1)
%         noverlapphases{idx_sitphase} = idx_sitphases(idx_sitphase,2) - idx_sitphases(idx_sitphase+1,1);
%     end
% end
% 
% %Turn overlapping end of each phase to missing
% LFP_hab(end-noverlapphases{1}+1:end,:) = NaN;
% LFP_cage(end-noverlapphases{2}+1:end,:) = NaN;
% LFP_social(end-noverlapphases{3}+1:end,:) = NaN;
% 
% t_hab(end-noverlapphases{1}+1:end) = NaN;
% t_cage(end-noverlapphases{2}+1:end) = NaN;
% t_social(end-noverlapphases{3}+1:end) = NaN;

%% Load accelerometer

labels = {'AUX1','_33.continuous',... %AUX1
    'AUX2','_34.continuous',... %AUX2
    'AUX3','_35.continuous',}; %AUX3

[ACC,t,info] = load_open_ephys_data_dbm(filepath,[],labels,[],Fs);
Nacc = 3;

figure
plotmultisignals(t,ACC');
hold on,xline(idx_sitphases(:),'k--')
    title('Accelerometer')
    xlabel('Tempo (s)')
    ylabel('g')
    pause(0)

%% Concatenate ACC over SIT phases
% ACC_hab = ACC(idx_sitphases(1,1):idx_sitphases(1,2),:);
% ACC_cage = ACC(idx_sitphases(2,1):idx_sitphases(2,2),:);
% ACC_social = ACC(idx_sitphases(3,1):idx_sitphases(3,2),:);
% ACC_post = ACC(idx_sitphases(4,1):idx_sitphases(4,2),:);

%% Concatenate ACC over SIT phases (with padding)
npad = 10*60*Fs;
idx_sitphases_pad = idx_sitphases + npad;

%fill sides with NaN to allow standardized concatenation
ACC = [nan(npad,Nacc) ; ACC; nan(npad,Nacc)]; 

ACC_hab = ACC(idx_sitphases_pad(1,1):idx_sitphases_pad(1,2),:);
ACC_cage = ACC(idx_sitphases_pad(2,1):idx_sitphases_pad(2,2),:);
ACC_social = ACC(idx_sitphases_pad(3,1):idx_sitphases_pad(3,2),:);
ACC_post = ACC(idx_sitphases_pad(4,1):idx_sitphases_pad(4,2),:);

% Ignore overlapping with next phase
ACC_hab(end-noverlapphases{1}+1:end,:) = NaN;
ACC_cage(end-noverlapphases{2}+1:end,:) = NaN;
ACC_social(end-noverlapphases{3}+1:end,:) = NaN;

%% Save data
eval([dataID '_idx_sitphases = idx_sitphases;']);

eval([dataID '_LFP_hab = LFP_hab;']);
eval([dataID '_LFP_cage = LFP_cage;']);
eval([dataID '_LFP_social = LFP_social;']);
eval([dataID '_LFP_post = LFP_post;']);

eval([dataID '_t_hab = t_hab;']);
eval([dataID '_t_cage = t_cage;']);
eval([dataID '_t_social = t_social;']);
eval([dataID '_t_post = t_post;']);

try %if there is acc data
eval([dataID '_ACC_hab = ACC_hab;']);
eval([dataID '_ACC_cage = ACC_cage;']);
eval([dataID '_ACC_social = ACC_social;']);
eval([dataID '_ACC_post = ACC_post;']);
end

%Save data to specified folder
disp(['Saving LFP of ' dataID])
save([savedatadir '\' dataID '_lfp.mat'],'-regexp',dataID);
disp(['Saved LFP of ' dataID])
