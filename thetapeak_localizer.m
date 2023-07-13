%% THETA PEAK LOCALIZER
%addpath 'C:\Users\Hilde Althof\OneDrive\Universiteiten\Radboud Universiteit\BabyBRAIN Lab (BRC)\fieldtrip-20230118'
addpath 'C:\Users\luuks\OneDrive\Documents\BCRC\Fieldtrip Toolbox\fieldtrip-20221004'

clc
close all
clearvars
ft_defaults

% Introduce subject number
subject = 39;

% Read data and cut trials
cd(sprintf('//cnas.ru.nl//wrkgrp//STD-Donders-DCC-Hunnius-ReadyToLearnEEG//Ready_to_learn_Children_Study//Data//Raw data//Neural//sub-39',subject));

%up to participant 8, this should be changed to 'S 80'.
%participant >= 36 should also be 'S 80' again
eventval = {'S 80'};

cfg = [];
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.eventtype  = 'Stimulus'; % for participant <= 8 and >= 36 , this should be changed to 'Stimulus'
cfg.trialdef.eventvalue = eventval;
cfg.trialdef.prestim    = -0.5;
cfg.trialdef.poststim   = 2.5;
cfg.dataset = sprintf('subject%02d.eeg',subject);
cfg = ft_definetrial(cfg);
data_trials=ft_preprocessing(cfg);

cfg=[];
%cfg.channel  = {'4','3','6','7','28','29','30','31','2','10','21'}; %for participants with workspace 'error'
channels_selection =ft_selectdata (cfg,data_trials);

% Quick artefact rejection
if ~exist(fullfile(sprintf('dataclean_s%02d.mat',subject)))
    cfg=[];
    cfg.method='summary';
    dataclean=ft_rejectvisual(cfg,channels_selection);
    dataclean.goodchannels=dataclean.label(find(~isnan(dataclean.trial{1,1}(:,1))));
    save(sprintf('dataclean_s%02d.mat',subject),'dataclean');
else
    load (fullfile(sprintf('dataclean_s%02d.mat',subject)));
end

% Re-reference
cfg             = [];
cfg.derivative = 'yes'; %'Turn it on if peak is not clear without it'
cfg.reref       = 'yes';
cfg.implicitref = 'FCz';
cfg.refchannel  = {'TP9', 'TP10'};
data_rereference    = ft_preprocessing(cfg, dataclean);

% Power anaylsis
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.foilim = [2 15];
cfg.taper = 'hanning';
freq_outcome = ft_freqanalysis(cfg, data_rereference);

cfg=[];
cfg.avgoverchan = 'yes';
cfg.nanmean     ='yes';
TFR_averag_channels = ft_selectdata(cfg, freq_outcome);

% Quick plot
plot(TFR_averag_channels.freq,TFR_averag_channels.powspctrm);

%%
% Possibilities and correspoding freq code
% 6Hz - 20
% 5.5Hz - 22
% 5Hz - 24
% 4.5Hz - 27
% 4Hz - 30
% 3.5Hz - 34
% 3Hz - 40