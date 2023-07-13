% this script converts the fiorst version of the BIDS dataset such that all
% recordings are at 250Hz. Furthermore, it ensures that the channel names are correct
% (which is not the case in some of the original recordings and in v1 of the BIDS
% dataset)

%%

% these directories will be different on another computer
bids_v1 = 'W:\Ready_to_learn_Children_Study\Analysis\BIDS-v1';
bids_v2 = 'W:\Ready_to_learn_Children_Study\Analysis\BIDS-v2';

%%

participants_tsv = ft_read_tsv(fullfile(bids_v1, 'participants.tsv'));

clear data_raw data_resampled data_resampled_renamed

for i=17:size(participants_tsv,1)

  sub = participants_tsv.participant_id{i}; % this includes "sub-"
  task = 'audiovisual';

  %% READ THE ORIGINAL DATA

  cfg = [];
  cfg.dataset = fullfile(bids_v1, sub, 'eeg', sprintf('%s_task-%s_eeg.vhdr', sub, task));
  data_raw = ft_preprocessing(cfg);

  %% FIX THE SAMPLING RATE

  if data_raw.fsample==250
    % keep as it is
    data_resampled = data_raw;

  elseif data_raw.fsample==500
    % resample to 250 Hz
    cfg = [];
    cfg.resamplefs = 250;
    data_resampled = ft_resampledata(cfg, data_raw);

  else
    error('Unexpected sampling rate')
  end

  %% FIX THE CHANNEL NAMES

  incorrect = {'Fp1' 'Fp2' 'F3' 'F4' 'C3' 'C4' 'P3' 'P4' 'O1' 'O2' 'F7' 'F8' 'T7' 'T8' 'P7' 'P8' 'Fz' 'Cz' 'Pz' 'FC1' 'FC2' 'CP1' 'CP2' 'FC5' 'FC6' 'CP5' 'CP6' 'TP9' 'TP10' 'Eog' 'Ekg1' 'Ekg2'};
  correct = {'Fp1' 'Fz' 'F3' 'F7' 'FT9' 'FC5' 'FC1' 'C3' 'T7' 'TP9' 'CP5' 'CP1' 'Pz' 'P3' 'P7' 'O1' 'Oz' 'O2' 'P4' 'P8' 'TP10' 'CP6' 'CP2' 'Cz' 'C4' 'T8' 'FT10' 'FC6' 'FC2' 'F4' 'F8' 'Fp2'};

  assert(length(incorrect)==32);
  assert(length(correct)==32);

  montage = [];
  montage.tra = eye(32);
  montage.labelold = incorrect;
  montage.labelnew = correct;

  if strcmp(data_resampled.label{2}, correct{2})
    % this is the correct channel assignment
    data_resampled_renamed = data_resampled;

  elseif strcmp(data_resampled.label{2}, incorrect{2})
    % this is the incorrect channel assignment
    cfg = [];
    cfg.montage = montage;
    data_resampled_renamed = ft_preprocessing(cfg, data_resampled);

  else
    error('Unexpected channel names')
  end
 
  data_resampled_renamed.label = data_resampled_renamed.label(:);

  %% WRITE IT BACK TO DISK

  % there are a number of files with additional metadata
  channels_tsv = ft_read_tsv( fullfile(bids_v1, sub, 'eeg', sprintf('%s_task-%s_channels.tsv', sub, task)));
  events_tsv   = ft_read_tsv( fullfile(bids_v1, sub, 'eeg', sprintf('%s_task-%s_events.tsv', sub, task)));
  eeg_json     = ft_read_json(fullfile(bids_v1, sub, 'eeg', sprintf('%s_task-%s_eeg.json', sub, task)));

  % some of the metadata fields need to be updated
  eeg_json.SamplingFrequency = 250;

  channels_tsv.name = data_resampled_renamed.label;
  channels_tsv.type(:) = {'EEG'};
  channels_tsv.units(:) = {'uV'};
  channels_tsv.sampling_frequency(:) = 250;

  cfg = [];
  cfg.method = 'convert';
  cfg.bidsroot = bids_v2;
  cfg.sub = sub(5:end); % only the part following the "sub-"
  cfg.task = task;
  cfg.datatype = 'eeg';

  % this does not yet work, but I will improve data2bids so that this does work
  % cfg.channels = channels_tsv;

  % this is a work-around
  cfg.channels = table2struct(channels_tsv, 'ToScalar', true);

  cfg.events = events_tsv;
  cfg.eeg = eeg_json;

  data2bids(cfg, data_resampled_renamed);

end

%% COPY THE GENERAL METADATA FILES

copyfile(fullfile(bids_v1, 'dataset_description.json'), fullfile(bids_v2, 'dataset_description.json'))
copyfile(fullfile(bids_v1, 'participants.tsv'), fullfile(bids_v2, 'participants.tsv'))
% copyfile(fullfile(bids_v1, 'README'), fullfile(bids_v2, 'README'))
% copyfile(fullfile(bids_v1, 'CHANGES'), fullfile(bids_v2, 'CHANGES'))
