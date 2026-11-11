% ANALYSIS OF LFP VEH x VPA
% Project: PhD's Adriana Souza
% 
% Pipeline:
% 
% Authors: Danilo Benette Marques, Adriana Souza. 2025
% 
% Last update: 08-10-2025

cd('D:\FGomesLab\Adriana Souza')

%% Load and concatenate
% Load data to workspace

% clear idx_train tf_train
clearvars -except *_LFP_* *_ACC_* *_t_* *_idx_sitphases

varNames = who; %keep true names for eval

%Get subejects ID
[subjNames,~,~,] = cellfun(@getsubjid,varNames,'un',0);

subjNames = unique(subjNames);
[subjNames,groupNames,sexNames] = cellfun(@getsubjid,strcat(subjNames,'_'),'un',0);

%Run concat 
for idx_subj = 1:numel(subjNames)
subjName = subjNames{idx_subj};

idx_vars_subj = find(contains(varNames,subjName));

    %---LFP---
    %hab
    idx_var = find(contains(varNames,[subjName '_LFP_hab']));
    if idx_var
        LFP{idx_subj,1} = eval(varNames{idx_var});
    else
        LFP{idx_subj,1} = [];
    end

    %cage
    idx_var = find(contains(varNames,[subjName '_LFP_cage']));
    if idx_var
        LFP{idx_subj,2} = eval(varNames{idx_var});
    else
        LFP{idx_subj,2} = [];
    end

    %social
    idx_var = find(contains(varNames,[subjName '_LFP_social']));
    if idx_var
        LFP{idx_subj,3} = eval(varNames{idx_var});
    else
        LFP{idx_subj,3} = [];
    end

    %post
    idx_var = find(contains(varNames,[subjName '_LFP_post']));
    if idx_var
        LFP{idx_subj,4} = eval(varNames{idx_var});
    else
        LFP{idx_subj,4} = [];
    end

    %---ACC---
    %hab
    idx_var = find(contains(varNames,[subjName '_ACC_hab']));
    if idx_var
        ACC{idx_subj,1} = eval(varNames{idx_var});
    else
        ACC{idx_subj,1} = [];
    end

    %cage
    idx_var = find(contains(varNames,[subjName '_ACC_cage']));
    if idx_var
        ACC{idx_subj,2} = eval(varNames{idx_var});
    else
        ACC{idx_subj,2} = [];
    end

    %social
    idx_var = find(contains(varNames,[subjName '_ACC_social']));
    if idx_var
        ACC{idx_subj,3} = eval(varNames{idx_var});
    else
        ACC{idx_subj,3} = [];
    end

    %post
    idx_var = find(contains(varNames,[subjName '_ACC_post']));
    if idx_var
        ACC{idx_subj,4} = eval(varNames{idx_var});
    else
        ACC{idx_subj,4} = [];
    end
        
    %---Idx_sitphases---
    idx_var = find(contains(varNames,[subjName '_idx_sitphases']));
    if idx_var
        Idx_sitphases{idx_subj,1} = eval(varNames{idx_var});
    else
        Idx_sitphases{idx_subj,1} = [];
    end         
        
end %subject

%% Identifications and general parameters
idx_incl = true(size(LFP,1),1);

%Identify groups
groupNames = categorical(groupNames);
% GroupNames = unique(groupNames)';
GroupNames = categorical(["SAL","VPA"]);

groupcolor = zeros(numel(GroupNames),3);
try groupcolor(contains(string(GroupNames),'SAL'),:) = [.2 .2 .2]; end
try groupcolor(contains(string(GroupNames),'VPA'),:) = [1 .25 .25]; end

%Identify sex
sexNames = categorical(sexNames);
SexNames = unique(sexNames)';

% sexcolor = zeros(numel(SexNames),3);
% try sexcolor(contains(string(SexNames),'_M_'),:) = [.25 .25 1]; end
% try sexcolor(contains(string(SexNames),'_F_'),:) = [1 .25 .25]; end
sexcolor(1) = 1;
sexcolor(2) = .5;

groupNames2 = strcat(string(groupNames),'_',string(sexNames));
groupNames2 = categorical(groupNames2);
GroupNames2 = unique(groupNames2);
groupcolor2 = zeros(2,3,2);
try groupcolor2(1,:,1) = [0 0 0]; end %VEH_F
try groupcolor2(1,:,2) = [.25 .25 .25]; end %VEH_M
try groupcolor2(2,:,1) = [1 .25 .25]; end %VPA_F
try groupcolor2(2,:,2) = [.25 .75 .25]; end %VPA_M

%Region labels
regionlabels{1} = 'PFC'; %PFC
regionlabels{2} = 'vHip'; %HipV
regionlabels{3} = 'NAc'; %NAc
regionlabels{4} = 'VTA'; %VTA

%Identify noise or missing epochs
% Idx_noise = cellfun(@isnan,Lfp,'un',0); %get already NaNs

%General parameters
Fs = 500;

%SIT phases labels
sitphaselabels{1} = 'HAB'; %HAB
sitphaselabels{2} = 'EC'; %JAULA
sitphaselabels{3} = 'SOCIAL'; %SOCIAL
sitphaselabels{4} = 'POST'; %PÓS

%OBS: SIT phases durations
% Phase 1: hab: 30min
% Phase 2: cage: 10min
% Phase 3: social: 20min
% Phase 4: post: 30min

%Time of interest of each phase
Nsitphases = 4; %number of phases
toi_sitphases(1,1) = 10*60*Fs; %end
toi_sitphases(2,1) = 10*60*Fs; %beginning
toi_sitphases(3,1) = 20*60*Fs; %beginning
toi_sitphases(4,1) = 10*60*Fs; %beginning

%Segmentation
L = 5*Fs; 
noverlap = 0*Fs;

% %% Ignore overlapping with next phase
% Turn end of phase overlapping with start of next phase to NaN
for idx_subj = 1:size(LFP,1)
    
    idx_sitphases = Idx_sitphases{idx_subj};

    noverlapphases = cell(3,1); %empty
    for idx_sitphase = 1:Nsitphases-1 %hab to social
    
        %If end of phase N is after start of phase N+1
        if idx_sitphases(idx_sitphase,2)>idx_sitphases(idx_sitphase+1,1)
            noverlapphases{idx_sitphase} = idx_sitphases(idx_sitphase,2) - idx_sitphases(idx_sitphase+1,1);
        end
    end
    
    %Turn overlapping end of each phase to missing
    LFP{idx_subj,1}(end-noverlapphases{1}+1:end,:) = NaN;
    LFP{idx_subj,2}(end-noverlapphases{2}+1:end,:) = NaN;
    LFP{idx_subj,3}(end-noverlapphases{3}+1:end,:) = NaN;

    % ACC{idx_subj,1}(end-noverlapphases{1}+1:end,:) = NaN;
    % ACC{idx_subj,2}(end-noverlapphases{2}+1:end,:) = NaN;
    % ACC{idx_subj,3}(end-noverlapphases{3}+1:end,:) = NaN;
end

% %% Cut LFP within time of interest per SIT phase
for idx_subj = 1:size(LFP,1)

    LFP{idx_subj,1} = LFP{idx_subj,1}(end-toi_sitphases(1)+1:end,:); %hab = 1:toi (end)
    LFP{idx_subj,2} = LFP{idx_subj,2}(1:toi_sitphases(2),:); %cage = end-toi+1 (beginning)
    LFP{idx_subj,3} = LFP{idx_subj,3}(1:toi_sitphases(3),:); %cage = end-toi+1 (beginning)
    LFP{idx_subj,4} = LFP{idx_subj,4}(1:toi_sitphases(4),:); %cage = end-toi+1 (beginning)

end

%Define indices of starts of each phase
clear idx_sitphases
idx_sitphases(1,1) = 1;
idx_sitphases(2,1) = toi_sitphases(1)+1;
idx_sitphases(3,1) = sum(toi_sitphases(1:2))+1;
idx_sitphases(4,1) = sum(toi_sitphases(1:3))+1;

%% Filter signals per SIT phase
locutof = 1;
hicutof = 50;

% %Adjust gain to mV
% LFP = cellfun(@(x) x*1000,LFP,'un',0);

for idx_subj = 1:size(LFP,1)
    for j = 1:size(LFP,2)
        fprintf('Filtering subject %d/%d phase %d \n',idx_subj,size(LFP,1),j)

        lfp = LFP{idx_subj,j}';

        %Find missing
        idx_missing_ch = all(isnan(lfp),2);
        idx_missing_samples = all(isnan(lfp),1);

        %Fill missing with zeros to apply filter
        lfp(idx_missing_ch,:) = 0;
        lfp(:,idx_missing_samples) = 0;

        %Filter
        lfp = eegfilt(lfp,Fs,locutof,0);
        lfp = eegfilt(lfp,Fs,0,hicutof);
    
        %Return missing
        lfp(idx_missing_ch,:) = NaN;
        lfp(:,idx_missing_samples) = NaN;
        
        LFP{idx_subj,j} = lfp';
    end
end

%% Concatenate continuous LFP (channels,samples,subjects)
if ~iscell(LFP)
    error('LFP is already continuous')
end

disp('Getting continuous LFP')
clear Lfp
for idx_subj = 1:size(LFP,1)

    Lfp(:,:,idx_subj) = cell2mat(cellfun(@transpose,LFP(idx_subj,:),'un',false));

end
LFP = Lfp; clear Lfp;

%Create time vector
dt = 1/Fs;
t = 0:dt:size(LFP,2)*dt-dt;
t = t/60; %min

t0 = t(idx_sitphases(3)); %reference to intruder
t = t-t0; 

%% Plot example LFP
idx_subj = 11; %subj legais: 11
lfp = LFP(:,:,idx_subj);

figure,
h=plotmultisignals(t*60,lfp,-.5);
    set(h,'color','k','linewidth',1)
    figurapadraoas
xline(t(idx_sitphases)*60,'--k')
return

xlim([0 10]) %subj 11 0:10 vs -10:0
scalebar(1,'1 seg',.5,'0,5 mV'); 
set(gca,'Visible','off')

%% Segmentation
clear LFPseg
disp('Segmenting LFP')
[LFPseg,tseg] = segmentation(LFP,L,noverlap);
disp('Segmented')

tseg = tseg(1,:)';
tseg = tseg/Fs; %seg
tseg = tseg/60; %min
tseg = tseg-t0; 

%Get segment phase ID
clear Idx_sitphases
Idx_sitphases = zeros(1,length(tseg));
Idx_sitphases(tseg>=t(idx_sitphases(1)))=1;
Idx_sitphases(tseg>=t(idx_sitphases(2)))=2;
Idx_sitphases(tseg>=t(idx_sitphases(3)))=3;
Idx_sitphases(tseg>=t(idx_sitphases(4)))=4;

%% Identify noise
%Option 1) amplitude threshold

% %Check amplitudes to define threshold
% Amp = cellfun(@(x) max(abs(x)),LFPseg);
% 
% for idx_subj = 1:size(Amp,3)
% amp = reshape(Amp(2,:,idx_subj),[],1);
% figure,histogram(amp,1000)
% xlim([0 5])
% % title(subjNames{idx_subj})
% title(mad(amp))
% pause(0)
% end

clear Idx_noise
disp('Finding noisy segments')
thr_noise = 0.5;
Idx_noise = cellfun(@(x) any(isnan(x)) | max(abs(x))>thr_noise,LFPseg);

%% Identify noise
%Option 2) Manual spectra sorting

%% Identify noise
%Option 3) semiautomatic scores and period inspection

if ~exist('Idx_noise') %create new Idx_noise if nonexistant
Idx_noise = false(size(LFPseg,1),size(LFPseg,2),size(LFPseg,3));
end

for idx_subj = 6%1:size(LFPseg,3)
    for coi = 1%:size(LFPseg,1)

        fprintf('\n')
        fprintf('Finding noise of subject %d/%d, region %s \n',idx_subj,size(LFPseg,3),regionlabels{coi})
        lfp = cell2mat(LFPseg(coi,:,idx_subj)')';

        if all(isnan(lfp))
            continue
        end

        %Sets predefined noise
        idx_noise = Idx_noise(coi,:,idx_subj)';
        lfp(:,idx_noise) = NaN;

        %Find noise
        idx_noise = findlfpnoiseepochs(lfp,Fs,[1 50],[ 7],true);

        Idx_noise(coi,:,idx_subj) = idx_noise';

    end
end


%% Power spectrum density
flim = [1 50];
nfft = 2^nextpow2(L); %2^nextpow2(L);
[f,Nf] = getfvector(nfft,Fs,flim);

PSD = nan(Nf,size(LFPseg,2),size(LFPseg,1),size(LFPseg,3));

for idx_subj = 1:size(LFPseg,3)

    for coi = 1:size(LFPseg,1)
        fprintf('Computing PSD of subject %d, region %d \n',idx_subj,coi)

        for j = 1:size(LFPseg,2)

            lfp = LFPseg{coi,j,idx_subj}';

            if ~Idx_noise(coi,j,idx_subj)
                % [psd,f] = periodogram(lfp,hamming(L),[],Fs);
                [psd,f] = pmtm(lfp,2,nfft,Fs);
                    psd = psd(f>=flim(1) & f<=flim(2));
                    f = f(f>=flim(1) & f<=flim(2));
            else
                continue
            end

            PSD(:,j,coi,idx_subj) = psd;

        end
    end
end
fprintf('\n')

%% Relative power
clear PSDrel
for idx_subj = 1:size(PSD,4)
    for coi = 1:size(PSD,3)

        psd = PSD(:,:,coi,idx_subj);
        
        %Relative to each segment
        psdrel = psd./repmat(nansum(psd,1),size(psd,1),1);

        %Relative to specific period
        % psdrel = reltf(psd,[],f,[1 50]); %spectrum sum of whole time
        % psdrel = reltf(psd,Idx_sitphases==1,f,[1 50]); %spectrum sum at hab
        % psdrel = reltf(psd,Idx_sitphases==2,f,[1 50]); %spectrum sum at cage

        PSDrel(:,:,coi,idx_subj) = psdrel;

    end
end

%% Baseline-normalized power
clear normPSD
for idx_subj = 1:size(PSD,4)
    for coi = 1:size(PSD,3)

        psd = PSD(:,:,coi,idx_subj);
        
        % normpsd = specgramnorm(psd,Idx_sitphases==1,5); %habituation
        normpsd = specgramnorm(psd,Idx_sitphases==2,5); %cage

        normPSD(:,:,coi,idx_subj) = normpsd;

    end
end


%% Functional connectivity
rng(1)
flim = [1 50];
[f,Nf] = getfvector(2^nextpow2(L),Fs,flim);

Nchch = nchoosek(size(LFPseg,1),2);
FC = nan(Nf,size(LFPseg,2),Nchch,size(LFPseg,3));
MPD = nan(Nf,size(LFPseg,2),Nchch,size(LFPseg,3));
FCsurr = nan(Nf,size(LFPseg,2),Nchch,size(LFPseg,3));

for idx_subj = 1:size(LFPseg,3)

    iter=0;
    c1c2=[];
    for c1 = 1:size(LFPseg,1)
    for c2 = c1+1:size(LFPseg,1)
        fprintf('Computing coherence of subject %d, regions %d x %d \n',idx_subj,c1,c2)

        iter=iter+1;
        c1c2(iter,:)=[c1 c2];

        % for j = 1:size(LFPseg,2)

            % lfp1 = LFPseg{c1,j,idx_subj}';
            % lfp2 = LFPseg{c2,j,idx_subj}';

            lfp1 = cell2mat(cellfun(@transpose,LFPseg(c1,:,idx_subj),'un',0));
            lfp2 = cell2mat(cellfun(@transpose,LFPseg(c2,:,idx_subj),'un',0));

            if ~Idx_noise(c1,j,idx_subj) && ~Idx_noise(c2,j,idx_subj)
                [fc,mpd,f] = fcspectrum(lfp1,lfp2,[3 5],[],[],Fs,'dwpli',flim);
                [fcsurr] = fcspectrum(lfp1,lfp2,[3 5],[],[],Fs,'dwpli',flim,[],true);
            else
                continue
            end

            FC(:,:,iter,idx_subj) = fc;
            try FCsurr(:,:,iter,idx_subj) = fcsurr; end
            MPD(:,:,iter,idx_subj) = mpd;

        % end
    end
    end
end
FC(FC<0)=0; %for dwPLI

%% Baseline-normalized coherence
clear normFC
for idx_subj = 1:size(FC,4)
    for coi = 1:size(FC,3)

        fc = FC(:,:,coi,idx_subj);
        
        normfc = specgramnorm(fc,Idx_sitphases==2,5); %baseline

        normFC(:,:,coi,idx_subj) = normfc;

    end
end

%% Spectrogram per subject
for idx_subj = 1:size(PSD,4)
    % continue
    for coi = 4%1:size(PSD,3)

    tf = 10*log10(PSD(:,:,coi,idx_subj));
    % tf = PSDrel(:,:,coi,idx_subj);
    % tf = normPSD(:,:,coi,idx_subj);
    % tf = FCsurr(:,:,coi,idx_subj);
    % tf = normCOH(:,:,coi,idx_subj);

    figure,
    % imagesc(tseg,f,tf);
    imagesc(tseg,f,smoothn(tf,[]));
    axis xy
    colormap jet

    ylim([1 50])
    xlim([-20 30])
    caxis([-60 -30])
    % caxis([0 .025])
    % caxis([-2 2])
    % caxis([0 1])
    % caxis([0 .5])

    ylabel('Frequência (Hz)')
    xlabel('Tempo (min)')
    title([subjNames{idx_subj} ' ' regionlabels{coi}])

    xline(t(idx_sitphases(2:end)),'--w')
    text(t(idx_sitphases(1:end))+1,45*ones(1,Nsitphases),sitphaselabels, ...
        'color','w')

    figurapadraoas(1)
    % colorbarpadrao('Potência (mV^{2}/Hz)',2);
    colorbarpadrao('Potência (dB)',2);

    end
end

%% Manually define included/excluded subjects
idx_incl = true(size(PSD,4),1); %all included
idx_excl = []; %exclude none

% idx_excl = [1 5 12 15 16 17];   idx_incl(idx_excl)=false;
% idx_excl = [5 6];               idx_incl(idx_excl)=false;
% idx_excl = [5];                 idx_incl(idx_excl)=false;
% idx_excl = [4 9 6];                 idx_incl(idx_excl)=false;
% idx_excl = [7];                 idx_incl(idx_excl)=false;

%Plot average spectrum
tf = squeeze(nanmean(PSD(:,:,4,:),2));
figure,
h=plot(f,tf,'linewidth',2); legend
set(h(idx_incl),'linestyle','-','color','k'); legend off
set(h(idx_excl),'linestyle','--','color','r')


%% Average spectrogram
clear sex

% for sex = SexNames
for group = GroupNames
    soi = groupNames==group & idx_incl; %non-sex loop
    % soi = groupNames==group & sexNames==sex & idx_incl; %for sex loop

    for coi = 1:size(PSD,3)
    
        % tf = nanmean((PSD(:,:,coi,soi)),4); YLabel = 'Power';
        % tf = nanmean(10*log10(PSD(:,:,coi,soi)),4); YLabel = 'Power (dB)';
        % tf = nanmean(PSDrel(:,:,coi,soi),4); YLabel = 'Relative power (u.a.)';
        tf = nanmean(normPSD(:,:,coi,soi),4); YLabel = 'Norm. power (Z-score)';
        % tf = nanmean(FC(:,:,coi,soi),4); YLabel = 'dwPLI';
        % tf = nanmean(normFC(:,:,coi,soi),4); YLabel = 'Norm. dwPLI (Z-score)';
    
        % tf = smoothn(tf,1);
    
        figure,
        imagesc(tseg,f,tf);
        imagesc(tseg,f,smoothn(tf,12));
        % contourfpadrao(tseg,f,smoothn(tf,10),60);
        axis xy
        colormap jet

        ylim([1 50])

        % xlim([-20 30])
        xlim([-10 10])
        
        % caxis([-60 -30])
        % caxis([0 .02])
        % caxis([-2 2])
        % caxis([0 1])
        % caxis([0 .5])
        caxis([-2 2])
    
        ylabel('Frequency (Hz)')
        xlabel('Time (min)')

        if exist('sex')
            % title([char(group) ' ' char(sex) ' ' regionlabels{coi}])
            % title([char(group) ' ' char(sex) ' ' regionlabels{c1c2(coi,1)} '-' regionlabels{c1c2(coi,2)}])
        else
            % title([char(group) ' ' regionlabels{coi}])
            % title([char(group) ' ' regionlabels{c1c2(coi,1)} '-' regionlabels{c1c2(coi,2)}])
        end

        xline(t(idx_sitphases(2:end)),'--k','linewidth',1)
        text(t(idx_sitphases(1:end))+1,45*ones(1,Nsitphases),sitphaselabels, ...
            'color','k')
    
        figurapadraoas(1)
        colorbarpadrao(YLabel,2);
        figresize(1,1.5)
    
    % copygraphics(gcf),pause
    % figname = lower(sprintf('avgspectrogram_%s_%s-%s_dwpli',group,regionlabels{c1c2(coi,1)},regionlabels{c1c2(coi,2)}));
    figname = lower(sprintf('avgspectrogram_%s_%s_normpower',group,regionlabels{coi}));
    fprintf('Saving %s \n',figname)
    savefigpadraoas(figname)
    end
% end %sex
end %group

%% Average PSD
toi = Idx_sitphases==3; %specific phase
toi = tseg>=0 & tseg<=5; %OR specific time

% figure
% for idx_sex = 1:length(SexNames)
    figure
for idx_group = 1:length(GroupNames)
    % soi = groupNames==GroupNames(idx_group) & sexNames==SexNames(idx_sex) & idx_incl; 
    soi = groupNames==GroupNames(idx_group) & idx_incl; idx_sex=1;

    for coi = 1:size(PSD,3)
    
        tf = nanmean(10*log10(PSD(:,toi,coi,soi)),2); YLabel = 'Power (dB)';
        % tf = nanmean(PSDrel(:,toi,coi,soi),2); YLabel = 'Relative power (u.a.)';
        % tf = nanmean(normPSD(:,toi,coi,soi),2); YLabel = 'Norm. power (Z-score)';
        % tf = nanmean(FC(:,toi,coi,soi),2); YLabel = 'dwPLI';
        % tf = nanmean(FCsurr(:,toi,coi,soi),2); YLabel = 'dwPLI';
        % tf = nanmean(normCOH(:,toi,coi,soi),2); YLabel = 'Norm. coherence (Z-score)';

        tf = permute(tf,[4 1 2 3]);
        
        hold on,plot(f,tf','linewidth',2,'color',groupcolor2(idx_group,:,idx_sex))
        % plotmeandev(f,tf,groupcolor2(idx_group,:,idx_sex),1);

        xlim([2 40])

        xlabel('Frequency (Hz)')
        ylabel(YLabel)
        % title(regionlabels{ch})
       
        figurapadraoas
    
    end
% end %sex
end %group

%% Average band power
clear bandlims bandlabels
bandlims(1,:) = [2 5]; bandlabels{1} = 'delta';
bandlims(2,:) = [5 10]; bandlabels{2} = 'theta';
bandlims(3,:) = [10 30]; bandlabels{3} = 'alfa/beta';
bandlims(4,:) = [30 50]; bandlabels{4} = 'gamma';

clear bandpsd bandpsdrel bandnormpsd
for boi = 1:size(bandlims,1)

    foi = f>=bandlims(boi,1) & f<=bandlims(boi,2);

    bandpsd(boi,:,:,:) = 10*log10(nanmean(PSD(foi,:,:,:),1));
    bandpsdrel(boi,:,:,:) = nanmean(PSDrel(foi,:,:,:),1);
    bandnormpsd(boi,:,:,:) = nanmean(normPSD(foi,:,:,:),1);

    % bandpsd(boi,:,:,:) = 10*log10(nanmax(PSD(foi,:,:,:),[],1));
    % bandpsdrel(boi,:,:,:) = nanmax(PSDrel(foi,:,:,:),[],1);
    % bandnormpsd(boi,:,:,:) = nanmax(normPSD(foi,:,:,:),[],1);

    try
    bandfc(boi,:,:,:) = nanmean(FC(foi,:,:,:),1);
    bandnormfc(boi,:,:,:) = nanmean(normFC(foi,:,:,:),1);
    bandmpd(boi,:,:,:) = nanmean(MPD(foi,:,:,:),1);
    end
end

%% Average band peak power
clear bandpklims bandpklabels
bandpklims(1,:) = [1 4]; bandlabels{1} = 'delta';
bandpklims(2,:) = [5 10]; bandlabels{2} = 'teta';

clear bandpsdpk bandpsdrelpk bandnormpsdpk
for boi = 1:size(bandpklims,1)

    foi = find(f>=bandpklims(boi,1) & f<=bandpklims(boi,2));

    for idx_subj = 1:size(PSD,4)
        for coi = 1:size(PSD,3)
            for idx_seg = 1:size(PSD,2) 
                
                psd = 10*log10(PSD(foi,idx_seg,coi,idx_subj));
                
                [psdpk,idx_pk] = findpeaks(psd,'NPeaks',1);

                if isempty(psdpk)
                    bandpsdpk(boi,idx_seg,coi,idx_subj) = NaN;
                    bandpsdrelpk(boi,idx_seg,coi,idx_subj) = NaN;
                    bandnormpsdpk(boi,idx_seg,coi,idx_subj) = NaN;
                else
                    bandpsdpk(boi,idx_seg,coi,idx_subj) = psdpk;
                    bandpsdrelpk(boi,idx_seg,coi,idx_subj) = PSDrel(foi(idx_pk),idx_seg,coi,idx_subj);
                    bandnormpsdpk(boi,idx_seg,coi,idx_subj) = normPSD(foi(idx_pk),idx_seg,coi,idx_subj);
                end


            end
        end
    end

end

%% Average band power over time
toi = tseg>-10 & tseg<10; %time of interest
boi = 4; %band of interest

%Number of blocks
Nbl = 12; %if L=5 ~1min

results = [];
    for coi = 1:1:size(PSD,3)
% for idx_sex = 1:length(SexNames)
    figure
for idx_group = 1:length(GroupNames)
    % soi = groupNames==GroupNames(idx_group) & sexNames==SexNames(idx_sex) & idx_incl; %per sex
    soi = groupNames==GroupNames(idx_group) & idx_incl; idx_sex=1; %no sex
    
    
        % avgband = bandpsd(boi,toi,coi,soi); YLabel = ['Potência de ' bandlabels{boi} ' (dB)'];
        % avgband = bandpsdrel(boi,toi,coi,soi); YLabel = ['Potência relativa de ' bandlabels{boi} ' (a.u.)'];
        avgband = bandnormpsd(boi,toi,coi,soi); YLabel = ['Potência norm. de ' bandlabels{boi} ' (Z-score)'];
        % avgband = bandfc(boi,toi,coi,soi); YLabel = [bandlabels{boi} ' dwPLI'];
        % avgband = bandnormfc(boi,toi,coi,soi); YLabel = ['Norm. ' bandlabels{boi} ' dwPLI (Z-score)'];
        % avgband = bandmpd(boi,toi,coi,soi); YLabel = ['Coerência norm. de ' bandlabels{boi} ' (Z-score)'];

        % avgband = bandpsdpk(boi,toi,coi,soi); 
        % avgband = bandpsdrelpk(boi,toi,coi,soi);
        % avgband = bandnormpsdpk(boi,toi,coi,soi); YLabel = ['Potência do pico norm. de ' bandlabels{boi} ' (Z-score)'];

        avgband = permute(avgband,[4 2 1 3]);

        %Average over blocks
        tblock = cellfun(@nanmean,segmentation(tseg(toi)',Nbl));
        avgband = cellfun(@nanmean,segmentation(avgband,Nbl));
        results{idx_sex,idx_group} = avgband;
        
        % hold on,plot(tblock,avgband','linewidth',3,'color',groupcolor2(idx_group,:,idx_sex))
        plotmeandev(tblock,avgband,groupcolor2(idx_group,:,idx_sex),1);
        xlim(get(gca,'XLim'))
    
        ylabel(YLabel)
        xlabel('Time (min)')
        % title(regionlabels{coi})
        % title(sprintf('%s-%s',regionlabels{c1c2(coi,1)},regionlabels{c1c2(coi,2)}))

        xline(t(idx_sitphases(2:end)),'--k')
        text(t(idx_sitphases(1:end))+1,45*ones(1,Nsitphases),sitphaselabels, ...
            'color','k')
    
        figurapadraoas
   
    end %channel
    
    % copygraphics(gcf),pause
    % figname = lower(sprintf('avgbandxtime_salvsvpa_%s-%s_%s_dwpli', ...
    %     regionlabels{c1c2(coi,1)},regionlabels{c1c2(coi,2)},bandlabels{boi}));
    % fprintf('Saving %s \n',figname)
    % savefigpadraoas(figname)

% end %sex
end %group

%Concatenate M  F
results = {cell2mat(results(:,1))  cell2mat(results(:,2)) };

% open results
[h,p,tstat] = ttest2(results{1},results{2});
p=p(tblock>0)

%Transpose
results = cellfun(@transpose,results,'un',0);
% close all
open results

%% Compare average band per SIT phase
coi = 1; %channel
boi = 2; %band of interest

results = [];
id_results = []; id=0;
for idx_group = 1:length(GroupNames)
for idx_sex = 1:length(SexNames)
    soi = groupNames==GroupNames(idx_group) & sexNames==SexNames(idx_sex) & idx_incl;
    id=id+1;
oi = subjNames(soi), continue
    for sitoi = 1:4
        toi = Idx_sitphases==sitoi; %each phase
        if sitoi == 3
        toi = tseg>=0 & tseg<=5; %specific time of interest
        end

        % avgband = nanmean(bandpsd(boi,toi,coi,soi),2); YLabel = ['Potência de ' bandlabels{boi} ' (dB)'];
        % avgband = nanmean(bandpsdrel(boi,toi,coi,soi),2); YLabel = ['Potência rel. de ' bandlabels{boi} ' (a.u.)'];
        avgband = nanmean(bandnormpsd(boi,toi,coi,soi),2); YLabel = ['Potência norm. de ' bandlabels{boi} ' (Z-score)'];
        
        avgband = permute(avgband,[4 2 1 3]);
        
        results{idx_sex,idx_group}(:,sitoi) = avgband; %{sex,group}(obs,sitphase)
        id_results{idx_sex,idx_group}(:,1) = idx_group*ones(size(avgband)); %group
        id_results{idx_sex,idx_group}(:,2) = idx_sex*ones(size(avgband)); %sex
        id_results{idx_sex,idx_group}(:,3) = id*ones(size(avgband)); %group/sex
        if idx_group==1 && idx_sex==1
            id_results{idx_sex,idx_group}(:,4) = 1*ones(size(avgband)); %SAL F
        elseif idx_group==2 && idx_sex==1
            id_results{idx_sex,idx_group}(:,4) = 2*ones(size(avgband)); %VPA F
        elseif idx_group==1 && idx_sex==2
            id_results{idx_sex,idx_group}(:,4) = 3*ones(size(avgband)); %SAL M
        elseif idx_group==2 && idx_sex==2
            id_results{idx_sex,idx_group}(:,4) = 4*ones(size(avgband)); %VPA M
        end

    end
end
end

open id_results
open results

%Concatenate M  F
results = {cell2mat(results(:,1))  cell2mat(results(:,2)) };
id_results = {cell2mat(id_results(:,1))  cell2mat(id_results(:,2)) };
% open results

fixscale = 0;

%Compare within group
disp('Comparison VEH social vs cage:')
[h,p,~,tstat] = ttest(results{1}(:,3),results{1}(:,2)); %VEH social vs cage
p,disp(tstat)

figure,
% b3=boxchart([results{1}(:,2) results{1}(:,3)], ...
%     'BoxFaceAlpha',.3,'BoxEdgeColor','k','BoxFaceColor',groupcolor(2,:));
b3=boxchart([id_results{1}(:,1) ; id_results{1}(:,1)+1],[results{1}(:,2) ; results{1}(:,3)], ...
    'BoxFaceAlpha',.3,'BoxEdgeColor','k','BoxFaceColor',groupcolor(2,:));
% hold on,b1=bar([1],[nanmean(results{1}(:,2),1)]+fixscale,'facecolor',groupcolor(1,:),'linewidth',1);
% hold on,b2=bar([2],[nanmean(results{1}(:,3),1)]+fixscale,'facecolor',groupcolor(2,:),'linewidth',1);
hold on,plot([results{1}(:,2) results{1}(:,3)]'+fixscale,'o-','linewidth',1 ...
    ,'color','k','markerfacecolor',groupcolor(2,:),'markeredgecolor','k')
ylabel(YLabel); 

xticks([1 2]);
xticklabels({'Jaula','Social'})

title('VEH')
figurapadraoas

disp('Comparison VPA social vs cage:')
[h,p,~,tstat] = ttest(results{2}(:,3),results{2}(:,2)); %VPA social vs cage
p,disp(tstat)

figure,
% b3=boxchart([results{2}(:,2) results{2}(:,3)], ...
%     'BoxFaceAlpha',.3,'BoxEdgeColor','k','BoxFaceColor',groupcolor(2,:));
b3=boxchart([id_results{2}(:,1)-1 ; id_results{2}(:,1)-1+1],[results{2}(:,2) ; results{2}(:,3)], ...
    'BoxFaceAlpha',.3,'BoxEdgeColor','k','BoxFaceColor',groupcolor(2,:));
% hold on,b1=bar([1],[nanmean(results{2}(:,2),1)]+fixscale,'facecolor',groupcolor(1,:),'linewidth',1);
% hold on,b2=bar([2],[nanmean(results{2}(:,3),1)]+fixscale,'facecolor',groupcolor(2,:),'linewidth',1);
hold on,plot([results{2}(:,2) results{2}(:,3)]'+fixscale,'o-','linewidth',1 ...
    ,'color','k','markerfacecolor',groupcolor(2,:),'markeredgecolor','k')
ylabel(YLabel); 
xticks([1 2]); 
xticklabels({'Jaula','Social'})
title('VPA')
figurapadraoas

%Compare between groups
disp('Comparison between groups:')
[h,p,~,tstat] = ttest2(results{2},results{1});
p,disp(tstat)

figure,
b3=boxchart([id_results{1}(:,1) ; id_results{2}(:,1)],[results{1}(:,3) ; results{2}(:,3)], ...
    'BoxFaceAlpha',.75,'BoxEdgeColor','k','BoxFaceColor',groupcolor(2,:));

% hold on,b1=bar([1],[nanmean(results{1}(:,3),1)]+fixscale,'facecolor',groupcolor(1,:),'linewidth',1);
% hold on,eb1=errorbar([1],[nanmean(results{1}(:,3),1)]+fixscale, ...
%     dev2D(results{1}(:,3),1)+0,'color','k','linewidth',1);
% hold on,b2=bar([2],[nanmean(results{2}(:,3),1)]+fixscale,'facecolor',groupcolor(2,:),'linewidth',1);
% hold on,eb2=errorbar([2],[nanmean(results{2}(:,3),1)]+fixscale, ...
%     dev2D(results{2}(:,3),1)+0,'color','k','linewidth',1);

ylabel(YLabel); 
xticks([1 2]); 
xticklabels({'VEH','VPA'})
% xticklabels({'VEH F','VEH M','VPA F','VPA M'})
% xticklabels({'VEH F','VPA F','VEH M','VPA M'})
title(regionlabels{coi})
figurapadraoas

%% NESTED FUNCTIONS
function [subjName,groupName,sexName] = getsubjid(varName)
% Get rat name, group name, day number
    idx_ = findstr(varName,'_');

    subjName = varName(1:idx_(3)-1); 

    if contains(varName,'_M_')
        sexName = 'M';
    elseif contains(varName,'_F_')
        sexName = 'F';
    end

    if contains(varName,'_SAL_')
        groupName = 'SAL';
    elseif contains(varName,'VPA')
        groupName = 'VPA';
    end

end

function [groupName] = getgroupid(varName)
% Get group name
    idx_ = findstr(varName,'_');

    if numel(idx_)==1
        groupName = varName(idx_(1)+1:end);
    elseif numel(idx_)>1
        groupName = varName(idx_(1)+1:idx_(2)-1);
    end
end

% %get Fs-sampled behavioral occurrences
% function [idx_behavior] = classRObehavior(ts_obj,L,fs)
% end

function [] = figurapadraoas(box)

movegui('onscreen')

set(gcf,'Color','w')

set(gca,'xcolor','k','ycolor','k'...
    ,'FontName','Helvetica'...
    ,'FontSize',12 ...
    ,'FontWeight','normal' ...
    ,'Linewidth',1 ...
    ,'Box','off' ...
    ,'Layer','top')

if nargin==1
    if box
        set(gca,'box','on')
    end
end

end