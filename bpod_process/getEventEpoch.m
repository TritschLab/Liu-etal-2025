function [event_frames] = getEventEpoch (bhv_frames, eventname)
% function [event_frames] = getEventEpoch (bhv_frames, eventname)
% extracts the event_frames from bhv_frames according to the eventname
% specified. It also works for getting the time stamp from bhv_ts (replacing bhv_frames)
%
% INPUT
%   bhv_frames: from data_processed in routine procedure
%   eventname: -string subfield in the bhv_frames 
% OUTPUT
%   event_frames: nx2 matrix for epoches (states) in behavioral data
%
% HL 2016-2-28 modified from getEventframes.m
% 
%%
event_frames = []; 
disp(eventname);
% get subfield names
subfieldnames = fieldnames(bhv_frames{1}); %
% length(bhv_frames)
for i_field = 1:length(subfieldnames)
    
    for i_trial = 1:length(bhv_frames)
        if isfield(bhv_frames{i_trial}.(subfieldnames{i_field}), eventname) 
            if  ~isnan(bhv_frames{i_trial}.(subfieldnames{i_field}).(eventname)(:,1))% excluded NaN state
            event_frames = cat(1,event_frames,bhv_frames{i_trial}.(subfieldnames{i_field}).(eventname));
            end
        end
    end
    if isempty(event_frames)
        disp(['not in field: ', subfieldnames{i_field}]);
    else
        disp([' in ',subfieldnames{i_field}, ' break']);
        break;
    end
end

if isempty(event_frames)
   warning('No such event found, RETURNed empty matrix')
end