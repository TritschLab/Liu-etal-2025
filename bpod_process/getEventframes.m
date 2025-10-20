function [event_frames] = getEventframes (bhv_frames, eventname)
% function [event_frames] = getEventframes (bhv_frames, eventname)
% extract the event_frames from bhv_frames according to the eventname
% specified
% for states, get the first (ON)
% also works for getting the time stamp from bhv_ts
% INPUT
%   bhv_frames: from data_processed in routine procedure
%   eventname: -string subfield in the bhv_frames 
% OUTPUT
%   event_frames: n x 1 vector of frame numbers
% Haixin Liu 2016/1/13

%%
event_frames = []; disp(eventname);
% get subfield names
subfieldnames = fieldnames(bhv_frames{1}); %
% length(bhv_frames)
for i_field = 1:length(subfieldnames)
    
    for i_trial = 1:length(bhv_frames)
        if isfield(bhv_frames{i_trial}.(subfieldnames{i_field}), eventname) 
            % given Bpod only have .State and .Events
            if  ~isnan(bhv_frames{i_trial}.(subfieldnames{i_field}).(eventname)(:,1))% excluded NaN state
                switch subfieldnames{i_field}
                    % for Event, it is vector
                    % for State, it is Nx2 matrix
                    case 'States'
                        event_frames = cat(1,event_frames,bhv_frames{i_trial}.(subfieldnames{i_field}).(eventname)(:,1));
                    case 'Events'
                        event_frames = cat(1,event_frames,bhv_frames{i_trial}.(subfieldnames{i_field}).(eventname)(:));
                    otherwise
                        warning('Fields not in Bpod (States or Events):');
                        disp(subfieldnames{i_field});
                end
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
   warning('No such event found, RETURN empty matrix')
end