function [bhv_times, SessionData, curr_trial_list] = ReadBpodShiftTime2WS(bhv_fullfilename, WS_data)
% function [bhv_times, SessionData, curr_trial_list] = ReadBpodShiftTime2WS(bhv_fullfilename, WS_data)
% Convert all Bpod behavior into Ephus timestamps
%
% INPUT:
%   bhv_fullfilename: filenames of behavior file
%      will be prompted to select if not entered; can also be data structure contains .RawEvents.Trial
%
% OUTPUT:
%   bhv_times: behavior file from bpod, with all times converted into ephus recording clock (sampling)
% 
% modified from AP's. HL 2016-2-15
%  calls function: AP_ReadBitCodeWS

% prompt file selection:
if ~exist('bhv_fullfilename','var')
    [bhv_filename bhv_path] = uigetfile('*.mat','Pick behavior file','Multiselect','off');
    bhv_fullfilename = [bhv_path bhv_filename];
end
% Turn warning off for loading, otherwise get message about % StateMachineAssembler class
if ischar(bhv_fullfilename) % if input file name
    warning off
    load(bhv_fullfilename,'-MAT');
    warning on
    bhv = SessionData.RawEvents.Trial;
    bhv_times = bhv;
else % input is data after data_raw reading: data_raw.data_bpod
    bhv = bhv_fullfilename.RawEvents.Trial;
    SessionData = bhv_fullfilename;
    bhv_times = bhv;
end
  
% Get WS sampling rate
WS_sample_rate = WS_data.sr;

% Get trials in raw samples since started
WS_data.ch_names = cellfun(@(x) strip(x), WS_data.ch_names, 'UniformOutput', false);
trial_channel = cellfun(@(x) any(strfind(x,'Bicode_DI')),WS_data.ch_names);
curr_trial_list = AP_ReadBitCodeWS(WS_data.ch_data(:,trial_channel),WS_data.sr);
% bitcode is the same btw dispatcher and Bpod
% time clock is ephus recording N/samping rate

% delete trials the Bitcode is not properly recorded in WS
if length(bhv_times) >  size(curr_trial_list,1)
    bhv_times = bhv_times(1:size(curr_trial_list,1));
end

% HL 2021-8-10: add 2nd layer of correction as now record continuesly. 
% if the curr trial list trial number is not consecutive, fix it
if length(unique(diff(curr_trial_list(:,2))))>1
    disp(curr_trial_list(:,2)')
    warning('Trial number not consecutive. Auto-correcting...')
    curr_trial_list_auto = curr_trial_list;
    curr_trial_list_auto(:,2) = 1:size(curr_trial_list,1);
    disp((curr_trial_list_auto(:,2) - curr_trial_list(:,2))')
    disp('Diff between auto and original trial number (above)')
    curr_trial_list = curr_trial_list_auto;
end
% loop through those trials, find the offsets
for curr_trial_indx = 1:size(curr_trial_list,1)
    
    curr_trial = curr_trial_list(curr_trial_indx,2);
    
    % skip if it's the last trial and not completed in behavior
    if curr_trial > length(bhv) || curr_trial < 1
        continue
    end
        
    % the start time is the rise of the first bitcode, in Bpod data
    % structure, event or state time is recorded in relate to the trial
    % start (bitcode sync pulse)
    % this is different from dispatcher whose zero is the session start
    % curr_bhv_start =  bhv{curr_trial}.states.bitcode(1); %
    curr_xsg_bhv_offset = 0 - curr_trial_list(curr_trial_indx,1);
    % apply the offset to all numbers within the trial
    % Find all fields in overall structure of trial
    curr_fieldnames = fieldnames(bhv_times{curr_trial});
    
    for curr_field = 1:length(curr_fieldnames)
        % a) get subfields
        curr_subfields = fieldnames(bhv_times{curr_trial}.(curr_fieldnames{curr_field}));
        % b) find which subfields are numeric
        curr_numeric_subfields = structfun(@isnumeric,bhv_times{curr_trial}.(curr_fieldnames{curr_field}));
        % c) subtract offset from numeric fields and convert to frames
        for subfield_offset_fix = find(curr_numeric_subfields)'
            % get current values of the subfield
            curr_bhv_times = bhv_times{curr_trial}.(curr_fieldnames{curr_field}).(curr_subfields{subfield_offset_fix});
            % compensate for offset
            curr_bhv_times = curr_bhv_times - curr_xsg_bhv_offset;
            % convert to frames (get the closest frame from frame times) 
            bhv_times{curr_trial}.(curr_fieldnames{curr_field}).(curr_subfields{subfield_offset_fix}) = ...
                curr_bhv_times;
        end
    end
end