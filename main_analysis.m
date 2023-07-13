%% Main analysis for Theta Flicker: Children
% 12-07-2023

clear all
close all
clc

%% Restore FT defaults
restoredefaultpath
addpath 'W:\Ready_to_learn_Children_Study\Experiment_files\fieldtrip-20230118'
ft_defaults

cd 'W:\Ready_to_learn_Children_Study\Analysis'

%% General Cfg interpolating
if ~exist(fullfile(cd, 'my_neighbours.mat'))

    %prepare layout
    cfg         = [];
    cfg.layout  = 'ActiCap_64Ch_DCC_customized.mat';
    layout = ft_prepare_layout(cfg);

    %prepare neighborhood file
    cfg               = [];
    cfg.method        = 'triangulation';
    cfg.layout        = layout;
    neighbours = ft_prepare_neighbours(cfg);

    %edit neighbours
    cfg = [];
    cfg.neighbours=neighbours;
    cfg.enableedit='yes';
    cfg.layout  = 'ActiCap_64Ch_DCC_customized.mat';
    my_neighbours=ft_neighbourplot(cfg,neighbours);

    save(fullfile(cd, 'my_neighbours'), 'my_neighbours');
else
    load (fullfile(cd, 'my_neighbours.mat'));
end

%% Read subject data
subject_file      = [cd filesep 'BIDS-v2' filesep 'participants.tsv'];
subject_table     = readtable(subject_file, 'FileType', 'text', 'Delimiter', '\t');
sub               = subject_table.participant_id;

power_table = table();
average_chan = table();

%% Start looping over subjects
for ii =  1 : numel(sub)
    if ismember(sub{ii}, excluded_ids)
        continue; 
    end

    subject_id  = sub{ii};
    folder_name = sprintf('%s', subject_id);

    eegfile           = [cd filesep 'BIDS-v2' filesep sub{ii} filesep 'eeg' filesep sub{ii} '_task-audiovisual_eeg.eeg'];
    eegdata           = [cd filesep 'BIDS-v2' filesep sub{ii} filesep 'eeg' filesep sub{ii} '_task-audiovisual_eeg.vhdr'];
    dirout            = [cd filesep 'Offline analysis' filesep 'Output' filesep folder_name];

    if exist(eegdata, 'file')      
        if ~exist(dirout, 'dir')
            mkdir(dirout);
        end

        event = ft_read_event(eegdata);
        hdr = ft_read_header(eegdata);
    end
    
    %% Cut into trials
    if ~exist(fullfile(dirout, 'data_baseline.mat'))

        %Apply filters
        cfg = [];
        cfg.dataset = eegfile;
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 0.5; 
        cfg.lpfreq = 30;
        data_continuous = ft_preprocessing(cfg);

        %Define trials
        cfg = [];
        cfg.dataset = eegdata;
        cfg.trialfun  = 'ft_trialfun_bids_hilde';
        cfg.trialdef.prestim    = 0;
        cfg.trialdef.poststim   = 6;
        cfg.trialdef.type  = {'theta', 'random'};
        cfg_trial = ft_definetrial(cfg);

        %Make data continuous
        cfg = [];
        cfg.trl = cfg_trial.trl;
        data_baseline = ft_redefinetrial(cfg, data_continuous);

        cfg = [];
        cfg.demean = 'yes';
        cfg.baselinewindow=[0 6];
        data_baseline = ft_preprocessing(cfg, data_baseline);

        save(fullfile(dirout, 'data_baseline'), 'data_baseline');
    else
        load (fullfile(dirout, 'data_baseline.mat'));
    end


    %% Reject exclusion behaviour artifacts
        if ~contains({event.type}, 'exclusion')
        data_clean = data_baseline;
    else 
        %select exclusion artifacts in trials
        cfg = [];
        cfg.dataset = eegdata;
        cfg.trialfun  = 'ft_trialfun_bids_hilde';
        cfg.trialdef.prestim    = 0;
        cfg.trialdef.poststim   = 6;
        cfg.trialdef.type  = {'exclusion'};
        cfg_artifact = ft_definetrial(cfg);

        %remove exclusion artifacts from data
        cfg = [];
        cfg.method = 'partial';
        cfg.artfctdef.exclusion.artifact = cfg_artifact.trl(:,[1 2]);
        data_clean = ft_rejectartifact(cfg, data_baseline);
     end

    %Remove first 500 ms after start trial
    trl = data_clean.sampleinfo;
    trl(:,1) = trl(:,1) + data_clean.fsample * 0.5;
    trl(:,3) = 0;
    sel = trl(:,1) > trl(:,2);
    trl(sel,:) = [];

    cfg = [];
    cfg.trl = trl;
    data_clean_notransient = ft_redefinetrial(cfg, data_clean);

    %% Segment data without overlap
    cfg = [];
    cfg.length = 1;
    cfg.overlap = 0;
    data_clean_segmented = ft_redefinetrial(cfg, data_clean);

    %% Reject visual (first time)
    % trials
    if ~exist(fullfile(dirout, 'data_noisytr.mat'))
        cfg = [];
        cfg.layout = 'ActiCap_64Ch_DCC_customized.mat';
        cfg.channel = 'eeg';
        cfg.ylim = [-1 1]*1e2; % limits the y-axis
        cfg.method = 'trial';
        data_noisytr = ft_rejectvisual(cfg, data_clean_segmented);
        save(fullfile(dirout, 'data_noisytr'), 'data_noisytr');
    else
        load (fullfile(dirout, 'data_noisytr.mat'));
    end

    % channels
    if ~exist(fullfile(dirout, 'data_noisych.mat'))
        cfg = [];
        cfg.layout = 'ActiCap_64Ch_DCC_customized.mat';
        cfg.channel = 'eeg';
        cfg.ylim = [-1 1]*1e2; % limits the y-axis
        cfg.method = 'channel';
        cfg.keepchannel = 'nan';
        data_noisych = ft_rejectvisual(cfg, data_noisytr);
        save(fullfile(dirout, 'data_noisych'), 'data_noisych');
    else
        load (fullfile(dirout, 'data_noisych.mat'));
    end

    %% Storing list of good and bad channels
    subjectdata.badchannels  =  data_noisych.label(find(isnan(data_noisych.trial{1,1}(:,1))));
    subjectdata.goodchannels =  data_noisych.label(find(~isnan(data_noisych.trial{1,1}(:,1))));

    %% ICA
    if ~exist(fullfile(dirout, 'data_cleanICA.mat'))
        cfg = [];
        cfg.channel = subjectdata.goodchannels;
        cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
        comp = ft_componentanalysis(cfg, data_noisych);

        figure
        cfg = [];
        cfg.component = 1:numel(subjectdata.goodchannels);
        cfg.layout    = 'ActiCap_64Ch_DCC_customized.mat';
        cfg.comment   = 'no';
        ft_topoplotIC(cfg, comp)

        cfg = [];
        cfg.layout = 'ActiCap_64Ch_DCC_customized.mat';
        cfg.viewmode = 'component';
        ft_databrowser(cfg, comp)

        prompt = 'Enter components to be removed:';
        toberemov_comp = input(prompt,'s');
        toberemov_comp= str2num(toberemov_comp);

        %reject components
        cfg = [];
        cfg.channel = subjectdata.goodchannels;
        cfg.component = toberemov_comp;
        cfg.demean = 'no';
        data_cleanICA = ft_rejectcomponent(cfg, comp);

        save(fullfile(dirout,'data_cleanICA'), 'data_cleanICA');

    else
        load (fullfile(dirout, 'data_cleanICA'));
    end

    %% Reject visual (second time)
    % trials
    cfg = [];
    cfg.layout = 'ActiCap_64Ch_DCC_customized.mat';
    cfg.channel = 'eeg';
    cfg.ylim = [-1 1]*1e2; % limits the y-axis
    cfg.method = 'trial';
    data_noisytr_afterICA = ft_rejectvisual(cfg, data_cleanICA);
    save(fullfile(dirout, 'data_noisytr_afterICA'), 'data_noisytr_afterICA');

    % channels
    cfg = [];
    cfg.layout = 'ActiCap_64Ch_DCC_customized.mat';
    cfg.channel = 'eeg';
    cfg.ylim = [-1 1]*1e2; % limits the y-axis
    cfg.method = 'channel';
    cfg.keepchannel = 'nan';
    data_noisych_afterICA = ft_rejectvisual(cfg, data_noisytr_afterICA);
    save(fullfile(dirout, 'data_noisych_afterICA'), 'data_noisych_afterICA');

    %% Storing list of good and bad channels
    subjectdata.badchannels  =  data_noisych_afterICA.label(find(isnan(data_noisych_afterICA.trial{1,1}(:,1))));
    subjectdata.goodchannels =  data_noisych_afterICA.label(find(~isnan(data_noisych_afterICA.trial{1,1}(:,1))));

    %% Perform interpolation of bad channels
    % Interpolate
    if ~exist(fullfile(dirout, 'interpol_data.mat'))
       if isempty(subjectdata.badchannels)
           interpol_data = data_cleanICA;
       else
        cfg                = [];
        cfg.method         = 'weighted';
        cfg.badchannel     = subjectdata.badchannels;
        cfg.missingchannel =[];
        cfg.layout         = 'ActiCap_64Ch_DCC_customized.mat';
        cfg.neighbours     = my_neighbours.neighbours;
        interpol_data = ft_channelrepair(cfg,data_cleanICA);
       end

        save(fullfile(dirout, 'interpol_data'), 'interpol_data');
    else
        load (fullfile(dirout, 'interpol_data.mat'));
    end

    %% Rereferencing
    if ~exist(fullfile(dirout, 'data_rereference.mat'))
        cfg = [];
        cfg.reref = 'yes';
        cfg.refchannel = 'all';
        data_rereference = ft_preprocessing(cfg, interpol_data);

        save(fullfile(dirout,'data_rereference'), 'data_rereference');

    else
        load (fullfile(dirout, 'data_rereference'));
    end

    %% Segment data with overlap
    cfg = [];
    cfg.continuous = 'yes';
    data_rereference_continuous = ft_redefinetrial(cfg, data_rereference);

    cfg = [];
    cfg.length = 2;
    cfg.overlap = 0.75;
    data_rereference_segmented = ft_redefinetrial(cfg, data_rereference_continuous);

    %% Power per condition                                                                                                                                                                                                                                                                                   
    % Theta condition
    if ~exist(fullfile(dirout, 'freq_outcome_theta.mat'))
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'pow';
        cfg.foilim = [3 6]; 
        cfg.taper = 'hanning';
        cfg.keeptrials = 'no';
        cfg.channel= {'O1', 'Oz', 'O2', 'P7', 'P3', 'Pz', 'P4', 'P8'}; 
        cfg.trials = find(strcmp(data_baseline.trialinfo.type, 'theta'));
        freq_outcome_theta = ft_freqanalysis(cfg, data_rereference_segmented);

        save(fullfile(dirout,'freq_outcome_theta'), 'freq_outcome_theta');

    else
        load (fullfile(dirout, 'freq_outcome_theta'));
    end      
                         
    %Random condition
    if ~exist(fullfile(dirout, 'freq_outcome_random.mat'))
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output = 'pow';
        cfg.foilim = [3 6];
        cfg.taper = 'hanning';
        cfg.keeptrials = 'no';
        cfg.channel= {'O1', 'Oz', 'O2', 'P7', 'P3', 'Pz', 'P4', 'P8'}; 
        cfg.trials = find(strcmp(data_baseline.trialinfo.type, 'random'));
        freq_outcome_random = ft_freqanalysis(cfg, data_rereference_segmented);

        save(fullfile(dirout,'freq_outcome_random'), 'freq_outcome_random');

    else
        load (fullfile(dirout, 'freq_outcome_random'));
    end

    %% Average over channels of interest per condition
    %Theta
    if ~exist(fullfile(dirout, 'averag_channels_theta.mat'))
        cfg=[];
        cfg.avgoverchan = 'yes';
        averag_channels_theta = ft_selectdata(cfg, freq_outcome_theta);
        save(fullfile(dirout,'averag_channels_theta'), 'averag_channels_theta');

    else
        load (fullfile(dirout, 'averag_channels_theta'));
    end

    %Random
    if ~exist(fullfile(dirout, 'averag_channels_random.mat'))
        cfg=[];
        cfg.avgoverchan = 'yes';
        averag_channels_random = ft_selectdata(cfg, freq_outcome_random);
        save(fullfile(dirout,'averag_channels_random'), 'averag_channels_random');

    else
        load (fullfile(dirout, 'averag_channels_random'));
    end

    close all
    figure;
    hold on;
    plot(averag_channels_theta.freq, (averag_channels_theta.powspctrm), 'linewidth', 1.5)
    plot(averag_channels_random.freq, (averag_channels_random.powspctrm),'linewidth', 1.5)
    legend('Theta condition', 'Random condition')
    xlabel('Frequency (Hz)')
    ylabel('Power (\mu V^2)')

    %wait for button press
    disp('Have a look at figure')
    input('')

    %safe average of channels per conditions in table for all subjects
    average_chan.sub{ii} = sprintf('%02d', ii);
    average_chan.aver_theta{ii} = mean(averag_channels_theta.powspctrm);
    average_chan.aver_random{ii} = mean(averag_channels_random.powspctrm);

   save (fullfile('W:\Ready_to_learn_Children_Study\Analysis\Offline analysis\Output','average_chan'), 'average_chan', '-v7');
    
   %% Extract value for individual theta frequency per condition
    cfg = [];
    cfg.frequency = subject_table.thetafreq(ii); 
    IndivHz_theta = ft_selectdata(cfg, averag_channels_theta);

    cfg = [];
    cfg.frequency = subject_table.thetafreq(ii); 
    IndivHz_random = ft_selectdata(cfg, averag_channels_random);

end

%% Difference between conditions (theta-random)
cfg = [];
cfg.parameter    = 'powspctrm';
cfg.operation    = '(x1-x2)';
Cond_difference = ft_math(cfg,IndivHz_theta,IndivHz_random);

%PLEASE NOTE: you need to change the directory below to your own directory

if ~exist(fullfile('W:\Ready_to_learn_Children_Study\Analysis\Offline analysis\Output', 'power_diff_all.mat'))
    power_diff_all = zeros(30,2);
else load (fullfile('W:\Ready_to_learn_Children_Study\Analysis\Offline analysis\Output', 'power_diff_all.mat'))
    power_diff_all(ii,1) = ii;
    power_diff_all(ii,2) = Cond_difference.powspctrm;
end

save (fullfile('W:\Ready_to_learn_Children_Study\Analysis\Offline analysis\Output','power_diff_all'), 'power_diff_all');

%% Log transform power
cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'log10';
difference_logpow    = ft_math(cfg, Cond_difference);

%PLEASE NOTE: you need to change the directory below to your own directory

if ~exist(fullfile('W:\Ready_to_learn_Children_Study\Analysis\Offline analysis\Output', 'power_diff_all_log.mat'))
    power_diff_all_log = zeros(30,2);
else load (fullfile('W:\Ready_to_learn_Children_Study\Analysis\Offline analysis\Output', 'power_diff_all_log.mat'))
    power_diff_all_log(ii,1) = ii;
    power_diff_all_log(ii,2) = difference_logpow.powspctrm;
end

save (fullfile('W:\Ready_to_learn_Children_Study\Analysis\Offline analysis\Output','power_diff_all_log'), 'power_diff_all_log');
