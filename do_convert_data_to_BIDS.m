%Convert_data_to_BIDS_ThetaFlickerC
%12-07-2023

clear all
close all
clc

%% Restore FT defaults
restoredefaultpath
addpath 'W:\Ready_to_learn_Children_Study\Experiment_files\fieldtrip-20230118'
ft_defaults

%% Set the path for these scripts and for the data

scripts     = 'W:\Ready_to_learn_Children_Study\Analysis\BIDS-v1';
sourcedata  = 'W:\Ready_to_learn_Children_Study\Data\Raw data'; 
bidsroot    = 'W:\Ready_to_learn_Children_Study\Analysis\BIDS-v1'; 
results     = 'W:\Ready_to_learn_Children_Study\Analysis\BIDS-v1'; 
derivatives = 'W:\Ready_to_learn_Children_Study\Analysis\BIDS-v1'; 

%% Conversion of the original "sourcedata" into BIDS
%% 1.1 Specification of folders

cd(sourcedata)

% Delete the current BIDS folder if it already exists to have create the latest version
if exist(bidsroot, 'dir')
  rmdir(bidsroot, 's');
end

%% 1.2 Subject information

% Read the excel file containing subject info
subject_file            = [sourcedata filesep 'S_EEG_ThetaFlicker_participant_info_all'];
[num, text]             = xlsread(subject_file);
sub                     = text(2:end, 1);
sex                     = text(2:end, 3);
age                     = num(:,1);
thetafreq               = num(:,3);              

%% 1.3 General information for the data2bids function

% Start looping over subjects 
excluded_ids = {'sub-04','sub-12', 'sub-27', 'sub-29'}; 

for ii =  1 :length(sub)
    if ismember(sub{ii}, excluded_ids)
        continue; 
    end
 
  cfg                                         = [];
  cfg.method                                  = 'copy';
  cfg.bidsroot                                = bidsroot;
  cfg.datatype                                = 'eeg';
  cfg.writejson                               = 'replace';

  %% 1.4 Create the metafile in the dataset: dataset_description.json
  
  cfg.dataset_description.Name                = 'Theta-induced effects on memory encoding in 4-year-olds';
  cfg.dataset_description.DatasetType         = 'raw';
  cfg.dataset_description.BIDSVersion         = '1.2.0'; 
  cfg.dataset_description.Authors             = {'Hilde Althof'; 'Marlene Meyer'; 'Robert Oostenveld'; 'Sabine Hunnius'};
  cfg.dataset_description.License             = 'ODC-ODbL-1.0';

  %% 1.5 Create a file with project level participant information: participants.tsv
  
  % All subjects
  cfg.sub                               = sub{ii};
  cfg.sub(1:4)                          = []; % removes 'sub-' from the subject identifier
  
  cfg.participants.age                  = age(ii);
  cfg.participants.sex                  = sex{ii};
  cfg.participants.thetafreq            = thetafreq(ii);

  %% 1.6 Select the indiviudal EEG dataset per participant
 
  cfg.dataset                           = ['Neural' filesep sub{ii} filesep 'subject' sub{ii}(end-1:end) '_encoding.vhdr'];
  if cfg.dataset
    hdr                                 = ft_read_header(cfg.dataset);
  end

  %% 1.7 Create metadata information on the EEG recordings
  
  % Describing the task
  cfg.TaskName                                       = 'audiovisual';
  cfg.TaskDescription                                = 'Subjects watched a flickering flower followed by a picture of a cartoon animal in an certain environment. The flickerus luminosity was modulated with an individual theta frequency in half of the trials and at random, i.e. at no fixed frequency, in the rest';
  cfg.Instructions                                   = 'Subjects were asked to watch the flickering flower and animal carefully, and to try to sit still';
  
  % Describing the recording setup
  cfg.InstitutionName                                = 'Baby and Child Research Center';
  cfg.InstitutionAddress                             = 'Thomas van Aquinostraat 4, 6525 GD Nijmegen, the Netherlands';
  cfg.InstitutionalDepartmentName                    = 'BabyBRAIN group';
  
  cfg.Manufacturer                                   = 'Brain Products GmbH';
  cfg.ManufacturersModelName                         = 'BrainAmp Standard';
  cfg.eeg.CapManufacturer                            = 'Brain Products GmbH';
  cfg.eeg.CapManufacturersModelName                  = 'actiCAP 32Ch';
  cfg.eeg.EEGPlacementScheme                         = '10-20';
  cfg.eeg.EEGReference                               = 'FCz'; 
  cfg.eeg.EEGGround                                  = 'AFz';  
 
  if str2num(sub{ii}(end-1:end))>=23 || str2num(sub{ii}(end-1:end))<=35
      cfg.eeg.PowerLineFrequency                         = 50;
      cfg.eeg.HardwareFilters.LowCutoff.Frequency        = 0.1; 
      cfg.eeg.HardwareFilters.HighCutoff.Frequency       = 1000;
      cfg.eeg.SoftwareFilters.LowCutoff.Frequency        = 0.3; 
      cfg.eeg.SoftwareFilters.HighCutoff.Frequency       = 70;
      cfg.eeg.EEGChannelCount                            = 32;
  else
      cfg.eeg.PowerLineFrequency                         = 50;
      cfg.eeg.HardwareFilters.LowCutoff.Frequency        = 0.1; 
      cfg.eeg.HardwareFilters.HighCutoff.Frequency       = 1000;
      cfg.eeg.SoftwareFilters.LowCutoff.Frequency        = 0.1; 
      cfg.eeg.SoftwareFilters.HighCutoff.Frequency       = 125;
      cfg.eeg.EEGChannelCount                            = 32;
  end

  % Describing the recording
  cfg.eeg.RecordingType                              = 'continuous';
  
  %% 1.8 Create metadata on markers/triggers used in the EEG recordings: events.tsv

  % To do this, create segment data into events using ft_define_trial
  % Event code selection
  if str2num(sub{ii}(end-1:end))<=8 || str2num(sub{ii}(end-1:end))>=36
      event = ft_read_event(cfg.dataset);
      sel   = strcmp({event.value}, 'S 50') | strcmp({event.value}, 'S100') | strcmp({event.value}, 'S200') | strcmp({event.value}, 'S  1') | strcmp({event.value}, 'S  2') |...
              strcmp({event.value}, 'S  3') | strcmp({event.value}, 'S  4') | strcmp({event.value}, 'S  5') | strcmp({event.value}, 'S  6') | strcmp({event.value}, 'S  7') |...
              strcmp({event.value}, 'S  8') | strcmp({event.value}, 'S  9') | strcmp({event.value}, 'S 10') | strcmp({event.value}, 'S 11') | strcmp({event.value}, 'S 12')|...
              strcmp({event.value}, 'S 13') | strcmp({event.value}, 'S 14') | strcmp({event.value}, 'S 15') | strcmp({event.value}, 'S 16');
      event = event(sel);
      hdr   = ft_read_header(cfg.dataset);

      % Define animal names and animal locations
      animal_names = {'bear', 'giraffe', 'monkey', 'parrot', 'penguin', 'rabbit', 'seal', 'squirrel', ...
          'tiger', 'snake', 'hippo', 'zebra', 'ostrich', 'rhino', 'elephant', 'kangaroo'};

      animal_locations = containers.Map({'bear', 'giraffe', 'monkey', 'parrot', 'penguin', 'rabbit', 'seal', 'squirrel', ...
          'tiger', 'snake', 'hippo', 'zebra', 'ostrich', 'rhino', 'elephant', 'kangaroo'}, ...
          {'library', 'railway', 'pool', 'forest', 'street', 'museum', 'castle', 'beach', ...
          'bakery', 'hairdresser', 'cinema', 'classroom', 'gasstation', 'office', 'shop', 'disco'});

      % Name events as 'type' and specify which animal name corresponds to which event value
      type  = cell(size(event));
      animal = cell(size(event));
      for i = 1:numel(event)
          switch event(i).value
              case 'S 50'
                  type{i} = 'black';
                  animal{i} = 'n/a';
              case 'S100'
                  type{i} = 'theta';
                  animal{i} = 'n/a';
              case 'S200'
                  type{i} = 'random';
                  animal{i} = 'n/a';
              otherwise
                  type{i} = 'stimulus';
                  animal_idx = str2double(event(i).value(3:end));
                  animal{i} = animal_names{animal_idx};
          end
      end
  
  else
      event   = ft_read_event(cfg.dataset);
      sel     = strcmp({event.value}, 'R 50') | strcmp({event.value}, 'R100') | strcmp({event.value}, 'R200') | strcmp({event.value}, 'R  1') | strcmp({event.value}, 'R  2') |...
                strcmp({event.value}, 'R  3') | strcmp({event.value}, 'R  4') | strcmp({event.value}, 'R  5') | strcmp({event.value}, 'R  6') | strcmp({event.value}, 'R  7') |...
                strcmp({event.value}, 'R  8') | strcmp({event.value}, 'R  9') | strcmp({event.value}, 'R 10') | strcmp({event.value}, 'R 11') | strcmp({event.value}, 'R 12')|...
                strcmp({event.value}, 'R 13') | strcmp({event.value}, 'R 14') | strcmp({event.value}, 'R 15') | strcmp({event.value}, 'R 16');
      event   = event(sel);
      hdr     = ft_read_header(cfg.dataset);
    
     % Define animal names and animal locations
      animal_names = {'bear', 'giraffe', 'monkey', 'parrot', 'penguin', 'rabbit', 'seal', 'squirrel', ...
          'tiger', 'snake', 'hippo', 'zebra', 'ostrich', 'rhino', 'elephant', 'kangaroo'};

      animal_locations = containers.Map({'bear', 'giraffe', 'monkey', 'parrot', 'penguin', 'rabbit', 'seal', 'squirrel', ...
          'tiger', 'snake', 'hippo', 'zebra', 'ostrich', 'rhino', 'elephant', 'kangaroo'}, ...
          {'library', 'railway', 'pool', 'forest', 'street', 'museum', 'castle', 'beach', ...
          'bakery', 'hairdresser', 'cinema', 'classroom', 'gasstation', 'office', 'shop', 'disco'});

      % Name events as 'type' and specify which animal name corresponds to which event value
      type  = cell(size(event));
      animal = cell(size(event));
      for i = 1:numel(event)
          switch event(i).value
              case 'R 50'
                  type{i} = 'black';
                  animal{i} = 'n/a';
              case 'R100'
                  type{i} = 'theta';
                  animal{i} = 'n/a';
              case 'R200'
                  type{i} = 'random';
                  animal{i} = 'n/a';
              otherwise
                  type{i} = 'stimulus';
                  animal_idx = str2double(event(i).value(3:end));
                  animal{i} = animal_names{animal_idx};
          end
      end
  end

  % Assign corresponding locations to each animal name
  location = cell(size(animal));
  for i = 1:length(animal)
      if strcmp(animal{i}, 'n/a')
          location{i} = 'n/a';
      else
          location{i} = animal_locations(animal{i});
      end
  end

  % Remembered during recall: Load excel file 
  behaviour_path   = ['W:\Ready_to_learn_Children_Study\Analysis\Offline analysis\Behaviour' filesep sub{ii} filesep 'sub-' sub{ii}(end-1:end) '_behavioural.xlsx'];
  behaviour_data   = readtable(behaviour_path);
  behaviour_data.remembered(isnan(behaviour_data.remembered)) = 2;

  % Initialize new column with 'n/a'
  remembered = nan(size(location));
  
  % Find the corresponding location in location_recall. If a match is found, use the value of 'remembered'
  for i = 1:length(location)
      if strcmp(location{i}, 'n/a')
          remembered(i) = NaN;
      else
          idx = find(strcmp(behaviour_data.location_recall, location{i}));
          remembered(i) = behaviour_data.remembered(idx);
      end
  end

  %Specify onset and duration in s
  onset    = ([event.sample]'-1)/hdr.Fs;  
  duration = zeros(size(event, 1), 1);
    for i = 1:numel(event)
        if strcmp(type{i}, 'black')
             duration(i) = 1;
         elseif strcmp(type{i}, 'theta') || strcmp(type{i}, 'random')
             duration(i) = 6;
         elseif strcmp(type{i}, 'stimulus')
             duration(i) = 5.5;
         else
             error(['Unknown event type: ' event.type{i}]);
         end
    end

  % Transpose cell arrays
  type = type(:);
  animal = animal(:);
  location = location(:);
  remembered = remembered(:);
  duration = duration(:); 

  %Include exclusion behaviour during videocoding 
  % Prepare event table and remove first black screen
  eventtable = struct2table(event);
  eventtable(:,1) = [];

  % remove the duration, timestamp, offset column
  eventtable(:,3:5) = [];
  eventtable.onset = (eventtable.sample - 1)/hdr.Fs;
   
  %Read videocoding table 
  videocoding = ['W:\Ready_to_learn_Children_Study\Analysis\Offline analysis\Behaviour' filesep sub{ii} filesep 'sub-' sub{ii}(end-1:end) '_videocoding.tsv'];
  videotable = readtable(videocoding, 'Delimiter', 'tab', 'FileType', 'text', 'TreatAsEmpty', 'n/a');
  videotable.Properties.VariableNames =  {'type'  'empty'  'onset'  'duration'  'value'};

  %Define EEG triggers and videocoding triggers and put them side-by-side
  if str2num(sub{ii}(end-1:end))<=8 || str2num(sub{ii}(end-1:end))>=36
      triggercode = {'S  1', 'S  2', 'S  3', 'S  4', 'S  5', 'S  6', 'S  7', 'S  8', 'S  9', 'S 10', 'S 11', 'S 12', 'S 13', 'S 14', 'S 15', 'S 16'};
      triggercode = triggercode(:);
  else
      triggercode = {'R  1', 'R  2', 'R  3', 'R  4', 'R  5', 'R  6', 'R  7', 'R  8', 'R  9', 'R 10', 'R 11', 'R 12', 'R 13', 'R 14', 'R 15', 'R 16'};
      triggercode = triggercode(:);
  end

  stimuluscode = [    1;    2;    3;    4;    5;    6;    7;    8;    9;    10;    11;    12;    13;    14;    15;    16];
  mapping = table(triggercode, stimuluscode);

  % find the video frame at which each stimulus starts
  sel = strcmp(videotable.type, 'trial_stimulus');
  subset_videotable = videotable(sel,:);
  [~, ~, ib1] = intersect(mapping.stimuluscode, subset_videotable.value, 'stable');

  %find the onset in seconds at which the trigger happens
  sel = ismember(eventtable.value, triggercode);
  value_event = {[]};
  onset_event = nan;  
  subset_eventtable = table(value_event, onset_event);
  j = 1;

  for i = 2:numel(sel)
      if str2num(sub{ii}(end-1:end))<=8 || str2num(sub{ii}(end-1:end))>=36
          if sel(i) == 1 && (strcmp(eventtable.value(i-1), 'S100') || strcmp(eventtable.value(i-1), 'S200'))
              subset_eventtable.value_event(j) = eventtable.value(i);
              subset_eventtable.onset_event(j) = eventtable.onset(i-1);
              j = j + 1;
          end
      else
          if sel(i) == 1 && (strcmp(eventtable.value(i-1), 'R100') || strcmp(eventtable.value(i-1), 'R200'))
              subset_eventtable.value_event(j) = eventtable.value(i);
              subset_eventtable.onset_event(j) = eventtable.onset(i-1);
              j = j + 1;
          end
      end
  end

  assert(j==17) 

  [~, ~, ib2] = intersect(mapping.triggercode, subset_eventtable.value_event, 'stable');

  % ib1 and ib2 should be the same length and order
  assert(isequal(ib1, ib2))

  %add seconds in videotable
  y = subset_eventtable.onset_event;
  x = subset_videotable.onset;
  p = polyfit(x,y, 1);

  videotable.seconds = polyval(p, videotable.onset);

  %%Add exclusion behaviour to cfg.events
  %add empty exclusion column 
  exclusion = nan(size(onset));
  overview_events = table(onset, duration, type, animal, location, remembered, exclusion);
  
  %select rows and columns from videotable
  sel = strcmp(videotable.type, 'exclusion_behaviour');
  videotable = videotable(sel,:);
  videotable(:,2) = [];
  videotable.duration = videotable.duration/1000;
  videotable.type(:) = {'exclusion'};
  videotable.onset = videotable.seconds;
  videotable.exclusion = videotable.value;
  videotable(:,[4 5]) = [];
  videotable.animal(:) = {'n/a'};
  videotable.location(:) = {'n/a'};
  videotable.remembered(:) = nan;
  videotable = videotable(:,[2 3 1 5 6 7 4]);

  cfg.events = cat(1, overview_events, videotable);

  [~,idx] = sort(cfg.events.onset);
  cfg.events = cfg.events(idx,:);

  %% 1.9 Create metadata on the channels used in the EEG recordings: channels.tsv
  
  % Double info with eeg.tsv --> here only fill it out if it is channel specific
  cfg.channels.name               = hdr.label;
  cfg.channels.type               = repmat({'EEG'}, 32, 1);  % Type of channel
  cfg.channels.units              = repmat({'uV'}, 32, 1); % Physical unit of the data values recorded by this channel in SI
  cfg.channels.sampling_frequency = repmat(250, 32, 1); % Sampling rate of the channel in Hz.
  cfg.channels.low_cutoff         = repmat(0.1, 32, 1); % Frequencies used for the hardware high-pass filter applied to the channel in Hz
  cfg.channels.high_cutoff        = repmat(1000, 32, 1); % Frequencies used for the hardware low-pass filter applied to the channel in Hz.
  cfg.channels.notch              = repmat(nan, 32, 1); % Frequencies used for the notch filter applied to the channel, in Hz. If no notch filter applied, use n/a.

  %% Call data2bids
  
  data2bids(cfg);

  %% 1.10 Create a description of the markers/triggers used in the EEG recordings: events.json
  
  events_json                                 = [];
  events_json.onset.description               = 'Onset of the event';
  events_json.onset.units                     = 'seconds';
  events_json.duration.description            = 'Duration of the event';
  events_json.duration.units                  = 'seconds';
  events_json.begsample.description           = 'Sample at which event begins (measured from start of recording)';
  events_json.begsample.units                 = 'sample number';
  events_json.endsample.description           = 'Sample at which event ends (measured from start of recording)';
  events_json.endsample.units                 = 'sample number';
  events_json.offset.description              = 'Offset from begsample till start of the trial';
  events_json.offset.units                    = 'sample number';
  events_json.marker.description              = 'Marker number corresponding to this event as indicated in the .vmrk file';
  events_json.stimulus.description            = 'Type of stimulus presented to the subject';
  events_json.stimulus.black                  = 'Presentation of black screen before flickering flower. Marker = R50 '; 
  events_json.stimulus.theta                  = 'Presentation of a flickering flower in theta frequency. Marker = S100 or R100';
  events_json.stimulus.random                 = 'Presentation of a flickering flower in random frequency. Marker = S200 or R200';
  events_json.stimulus.animal                 = 'Presentation of a cartoon animal image. Marker: S1-S16 or R1-16';
  events_json.stimulus.location               = 'Presentation of location linked to cartoon animal';
  
  foldername                                  = [bidsroot filesep 'sub-' cfg.sub filesep 'eeg'];
  filename                                    = [foldername filesep 'sub-' cfg.sub '_task-' cfg.TaskName '_events.json'];
  
  ft_write_json(filename, events_json);

end 

%% 1.11 Create description of metadata on participants: participants.json

participants_json.participant_id.description    = 'Subject identifier';
participants_json.age.description               = 'age of each subject';
participants_json.age.units                     = 'monts.days';
participants_json.sex.description               = 'gender of each subject';
participants_json.sex.levels                    = {'female', 'male'};

filename                                        = [bidsroot filesep 'participants.json'];

ft_write_json(filename, participants_json);

%% Add the matlab code used to generate BIDS to a subfolder

destination                                     = [bidsroot filesep 'code'];
this_script                                     = '\\cnas.ru.nl\wrkgrp\STD-Donders-DCC-Hunnius-ReadyToLearnEEG\Ready_to_learn_Children_Study\Analysis\Offline analysis\scripts\do_convert_data_to_BIDS.m';

mkdir(destination);
copyfile(this_script, destination);