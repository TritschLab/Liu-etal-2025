function TrialNumber = AP_ReadBitCodeWS(BitCode_Ch, sr)
% TrialNumber = AP_ReadBitCodeWS(BitCode_Ch, sr)
%
% INPUT:
%   BitCode_Ch: Binary vector representing the bitcode
%   sr: Sampling rate in Hz
% OUTPUT: 
%   TrialNumber: trialnumber x 2, TrialNumber col 2 is trial number
%                                 TrialNumber col 1 is time (in seconds)
%                       precision determined by sample rate, which is 0.1 ms (default): bitcode_sync(curr_bitcode_sync)/xsg_sample_rate;
% HL 2020-3-5 adapted for WS recording
%
% This reads the bitcode from Dispatcher
% This is not to read the default trialnum_indicator with Dispatcher,
% because that doesn't work on our equipment.
% Instead, this reads a bitcode which has a sync signal followed by 12 bits
% for the trial number, which all have 5ms times with 5ms gaps in between.
% Bitcode is most significant bit first (2048 to 1).
% Total time: 5ms sync, 5ms*12 gaps + 5ms*12 bits = 125ms
% The temporal resolution of the linux state machine gives >0.1 ms loss
% per 5 ms period - keep in mind that this causes ~2.6ms to be
% lost over the course of 24 states.
%
% Note: the start of the trial is defined as the START of the bitcode
%
%%
num_bits = 12;

BinaryThreshold = BitCode_Ch; % as in Reach rig, Bitcode is DI channel,
ShiftBinaryThreshold = [NaN; BinaryThreshold(1:end-1)];
% Get raw times for rising edge of signals
rising_bitcode = find(BinaryThreshold==1 & ShiftBinaryThreshold==0);

% detect artifact: where the pulse is less than the designed pulse 5ms
bit_dur_design = 0.005*sr;
fall_bit = find(BinaryThreshold==0 & ShiftBinaryThreshold==1);
bit_dur = fall_bit - rising_bitcode;

% assume artifact is short and rare
[uni_bit_dur,ia,ic] = unique(bit_dur);
fprintf('Bit durations (sampling gap = %.2f ms): \n', 1/sr*1000)
for ii = 1:length(uni_bit_dur)
   fprintf('%i : %i, number - %i\n', ii, uni_bit_dur(ii), length(find(ic==ii))) ;
end
disp('auto artifact exclusion uses 1.1ms > or < than designed 5ms pulse')
temp_idx2delete = find(bit_dur < bit_dur_design-0.0011*sr | bit_dur > bit_dur_design+0.0011*sr); % use 1.1 ms around desinged duration as threshold 
rising_bitcode(temp_idx2delete) = [];

% Set up the possible bits, 12 values, most significant first
bit_values = [num_bits-1:-1:0];
bit_values = 2.^bit_values;

% Find the sync bitcodes: anything where the difference is larger than the
% length of the bitcode (16 ms - set as 20 ms to be safe)
bitcode_time_samples = 125*(sr/1000); % 125 ms
bitcode_sync = find(diff(rising_bitcode) > bitcode_time_samples);

% Assume that the first rising edge is a sync signal
if isempty(rising_bitcode)    
    TrialNumber = [];    
else    
    bitcode_sync = rising_bitcode([1;bitcode_sync + 1]); 
    % add the first one back and shift back to the rising pulse; get the bitcode index in time
    % Initialize TrialNumber output - col2 = trial number, col1 = time
    TrialNumber = zeros(length(bitcode_sync),2);
    % for each bitcode sync, check each bit and record as hi or low
    for curr_bitcode_sync = 1:length(bitcode_sync)
        curr_bitcode = zeros(1,num_bits);
        for curr_bit = 1:num_bits
            % boundaries for bits: between the half of each break
            % (bitcode_sync+5ms+2.5ms = 7.5ms)
            bit_boundary_min = bitcode_sync(curr_bitcode_sync) + 7.5*(sr/1000) + ...
                (curr_bit-1)*10*(sr/1000);
            bit_boundary_max = bitcode_sync(curr_bitcode_sync) + 7.5*(sr/1000) + ...
                (curr_bit)*10*(sr/1000);
            if any(rising_bitcode > bit_boundary_min & rising_bitcode < bit_boundary_max)
                curr_bitcode(curr_bit) = 1;
            end
        end
        curr_bitcode_trial = sum(curr_bitcode.*bit_values);
        % TrialNumber col 2 is trial number
        TrialNumber(curr_bitcode_sync,2) = curr_bitcode_trial;
        % TrialNumber col 1 is time (in seconds)
        TrialNumber(curr_bitcode_sync,1) = bitcode_sync(curr_bitcode_sync)/sr;
        
        % Catch the rare instance of the xsg file cutting out before the end of the bitcode 
        if bit_boundary_max > length(BinaryThreshold)
           TrialNumber(curr_bitcode_sync,:) = []; 
        end 
    end
    
    % Check here if anything fishy is going on, and warn user
    if ~all(diff(TrialNumber(:,2)))
        warning('TRIAL NUMBER WARNING: Nonconsecutive trials');
    end
    
    % Check here if increasing consecutive for continues recording, warn
    % user
    if any(find(diff(TrialNumber(:,2)) ~= 1))
        warning('TRIAL NUMBER WARNING: NOT consecutively increasing');
    end    
    
end

