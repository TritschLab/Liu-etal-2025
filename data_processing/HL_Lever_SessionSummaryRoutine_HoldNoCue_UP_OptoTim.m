%% this is a routime pipeline to process raw data and save into a specific format for further analysis, also have options to plot initial analysis result for quick examination

% function [Data] = HL_Lever_SessionSummaryRoutine_HoldNoCue_UP_OptoTim(mp4_name, flag_plot, params)
% INPUT:
%       mp4_name: 'HLXXX_date_sessionname', e.g 'HL112_200619_Lever_train'
%       flag_plot: some plot features swtich 
%       params: params structure for file names and channel for running routine without GUI input
% OUTPUT:
%       Data: structure containing session summary data
% 
% Note: for lever press with signal increase => flip the sign to make it negative
% and going down

function [Data] = HL_Lever_SessionSummaryRoutine_HoldNoCue_UP_OptoTim(mp4_name, flag_plot, params)
    
if nargin < 2
    flag_plot = false;
end

%% default folder names and parameters
default_folders;
% threshold to determine movement in lever trace analysis
movethresh = 0.0002; %
movethresh_20 = 0.00003;

%% initiate paramters names
if nargin < 1
    mp4_name = input('Type in the session name: HLXXX_date_sessionname e.g.: HL112_200619_Lever_train\n', 's');
end       

WS_fn = [];
Bpod_fn = [];

% tracking data skipped
% mat_fn = []; % processed mat file
% ses_fn = []; % tracking result mat file

% fetch Bpod file name and WS file name, load in data, and parse mp4 name
fprintf('Processing... %s\n',mp4_name);
tokens = regexp(mp4_name, ...
    '(?<animalname>\w\w\d\d\d)_(?<datename>\d\d\d\d\d\d)_(?<proc_name>\w+)*', 'names');
animalname = tokens.animalname;
ses_date = tokens.datename;
proc_name = tokens.proc_name;
fprintf('To save: \nanimalname - %s, ses_date - %s, proc_name - %s\n', ...
    animalname , ses_date , proc_name );

%%
if nargin<3 % use the input name
    check_temp = input('Please correct proc_name, if needed, otherwise press ENTER to continue');
    if ~isempty(check_temp)
        proc_name = check_temp;
    end
    %% get WS file
    % select the WS file if any, parse it or generate trial types
    % directly input ws_fn name
    [ws_fn,ws_f_path] = uigetfile(fullfile(WS_path,ses_date,animalname,'WS','*.h5'),...
        ['Select WaveServer data for ' animalname,' ',ses_date, ' ', proc_name, 'ESC if no WS data']);
    if any(ws_fn~=0) % find WS file
        WS_fn = [ws_f_path ws_fn];
    else % no WS data
        WS_fn= 0;
    end
    
    % get Bpod file
    [bpod_fn,bpod_f_path] = uigetfile(fullfile(WS_path,ses_date,animalname,'Bpod','*.mat'),...
        ['Select Bpod data for ' animalname,' ',ses_date, ' ', proc_name, 'ESC if no Bpod data']);
    if any(bpod_fn~=0) % find WS file
        Bpod_fn = [bpod_f_path bpod_fn];
    else % no WS data
        Bpod_fn = 0;
    end
else   % get file names from the input params structure
    WS_fn = params.WS_fn;
    Bpod_fn = params.Bpod_fn;
end

%% Load WS file and parse the data
% read in WS file
fprintf(2, 'WS file: %s \n', WS_fn);
if any(WS_fn~=0) % having WS file
    WS_data = HL_FP_loadWS_parseData(WS_fn);
    [WS_trial, ~, ~] = HL_FP_parseWSStiLib(WS_data.StiLib);
else
    warning('No WS data, just report Bpod result');
end

% load Bpod
% if having WS, register the time using bitcode; 
% if not, just load Bpod data
if any(WS_fn~=0)
    [bhv_times, SessionData] = ReadBpodShiftTime2WS(Bpod_fn,WS_data);
    while bhv_times{end}.States.BitcodeByte0(1) == 0 % for somereason not converted
        bhv_times(end) = [];
        fprintf(2, 'caution! last trial bitcode time start with 0')
    end
else
    bhv = load(Bpod_fn);
    SessionData = bhv.SessionData; clear bhv;
end

%% Bpod only situation: report Bpod result reward rate, then save and return
if WS_fn==0 % Not having WS file
    fprintf(2, 'No WS data, not implemented yet. RETURNED \n');
    return;

    n_trial = length(bhv_times);
    Outcome = zeros(n_trial,1);

    for i_trial = 1:length(bhv_times)
        % if rewarded
        Outcome(i_trial) = ~isnan(SessionData.RawEvents.Trial{i_trial}.States.Reward(1));
    end

    Reward_rate = sum(Outcome)/length(union(trial_ind_Cue_only, trial_ind_MaskCue));

    fprintf('Trial N=%d: Cue Only-%d, MaskLight and Cue -%d, MaskLight Only -%d. \n', ...
        n_trial, length(trial_ind_Cue_only), length(trial_ind_MaskCue), length(trial_ind_Mask_only));
    fprintf('Rewarded Rate (in cue trials)=%.3f %%\n', ...
        Reward_rate*100);

        % plot session summary
    figure; hold on;
    ylim([min(SessionData.TrialTypes)-0.5 max(SessionData.TrialTypes)+0.5]); 
    set(gca,'YTick', [1 2 3], 'YTickLabel', {'CueOnly', 'MaskLitCue', 'MaskLitOnly'});
    ylabel('Trial type')
    plot(1:n_trial, SessionData.TrialTypes, 'ob');
    % reward 
    plot(find(Outcome==1), SessionData.TrialTypes(Outcome==1), 'o', 'markerfacecolor','g', 'markeredgecolor', 'none');
    xlabel('Trial #');
    legend('All', 'Rewarded');
    title([animalname ' ' ses_date ' ' proc_name])

    % save and exit
    if exist(fullfile(save_path, animalname, ses_date),'dir') ~=7
    mkdir(fullfile(save_path, animalname, ses_date)) ;
    end

    save(fullfile(save_path, animalname, ses_date, [proc_name '.mat']),...
        'mp4_name','animalname', 'ses_date', 'proc_name', 'Bpod_fn',  'Reward_rate', 'Outcome', 'n_trial', 'trial_ind_*', ...
        'SessionData', 'WS_fn', '-v7.3')
    % 'bhv_times','lever*', 'movethresh', 'Lever*','idx_*',
    disp('Saved');
    Data.mp4_name           = mp4_name;
    Data.animalname         = animalname;
    Data.ses_date           = ses_date;
    Data.proc_name          = proc_name;
    Data.Bpod_fn            = Bpod_fn;
    Data.WS_fn              = WS_fn;
    Data.Outcome            = Outcome;
    Data.n_trial            = n_trial;
    Data.trial_ind_Cue_only = trial_ind_Cue_only;
    Data.trial_ind_MaskCue  = trial_ind_MaskCue;
    Data.trial_ind_Mask_only= trial_ind_Mask_only;
    Data.Reward_rate        = Reward_rate;
    Data.SessionData        = SessionData;

    disp('Saved');
    return;
end

%% Now having WS data, process lever signal and get behavior measures 
% Process lever signal
% get ch index for lever and others
WS_data.ch_names = cellfun(@(x) strip(x), WS_data.ch_names, 'UniformOutput', false);
idx_Lever_ch = find(cellfun(@(x) strcmp(x,'Lever'), WS_data.ch_names));
idx_cam_ch = find(cellfun(@(x) strcmp(x,'CamEXP_DI'), WS_data.ch_names));
if isempty(idx_cam_ch) % not reach rig, check LeftFP rig
    idx_cam_ch = find(cellfun(@(x) strcmp(x,'CamExp'), WS_data.ch_names));
    RigName = 'LeftFP';
else
    RigName = 'ReachRig';
end
idx_Lick_ch = find(cellfun(@(x) strcmp(x,'Lick'), WS_data.ch_names));


% flip the signal (-> negative) to have signal go down
[lever_active, lever_force_smooth, lever_velocity_envelope_smooth] ...
    = HL_parseLeverMovement(-WS_data.ch_data(:,idx_Lever_ch), WS_data.sr, movethresh);
% process lever trace for speed 20Hz 
[lever_active_20, lever_force_smooth_20,...
    lever_force_resample, lever_velocity_resample_smooth_20, movethresh_20] ...
    = HL_parseLeverMovement_Spd(-WS_data.ch_data(:,idx_Lever_ch), WS_data.sr, movethresh_20);

% 10Hz smoothed -> downsample time to 1k Hz
WS_data.ts_1k = downsample(WS_data.ts, WS_data.sr/1000);

% plot whole session: lever, cue, mvoement, reward
ITI = getEventEpoch (bhv_times, 'ITI'); % fixed amount <=2s
Reward = getEventEpoch (bhv_times, 'Reward');

if flag_plot

    figure;
    hist(diff(Reward(:,1)),[0:1:50])
    xlabel('Time (s)')
    ylabel('Reward count');
    title('Distribution of inter-reward distribution')
    set(gca,'XTick',[0 1  3  5 7 10 15 20 50])



    figure;
    hold on;
    plot(WS_data.ts_1k, lever_force_smooth, 'b-');
    plot(WS_data.ts_1k, lever_active*0.05 - 2, 'm-');

    for ii = 1:size(ITI,1)
        plot(ITI(ii,:), [1.3 1.3], '-r','linewidth',3);
    end

    for ii = 1:size(Reward,1)
        plot(Reward(ii,1), [1.33], '*g');    
    end

    title({'Session Overview'; 'g*-reward r-ITI m-lever active c-masklight'});
    xlabel('Time (s)');
    ylabel('Smoothed Lever (V, blue)');

end

%% get some behavior measures: correct rate respond rate and align movement to take a quick look 
n_trial = length(bhv_times);
Outcome = zeros(n_trial,1);
Resp = zeros(n_trial,1); % lever moved during opp win (reward window)


for i_trial = 1:length(bhv_times)
    % if rewarded
    Outcome(i_trial) = ~isnan(bhv_times{i_trial}.States.Reward(1));
    
    % rewwin are the same for all trials, ITI are different
    switch SessionData.TrialTypes(i_trial)
        case {1,3} % normal
            temp_win = bhv_times{i_trial}.States.OppWin;
            temp_idx = find(WS_data.ts_1k >=temp_win(1) & WS_data.ts_1k <=temp_win(2)); % time in lever active
            Resp(i_trial) = sum(lever_active(temp_idx))>0; % if lever is active during cue
                        
        case 2 % OptoTime Opto trigger given at Reward window onset 
            if isnan(bhv_times{i_trial}.States.CuePort1In(1)) % a miss trial
                temp_win = [bhv_times{i_trial}.States.OppWin(1) bhv_times{i_trial}.States.Miss(1)];
                temp_idx = find(WS_data.ts_1k >=temp_win(1) & WS_data.ts_1k <=temp_win(2)); % time in lever active
                Resp(i_trial) = sum(lever_active(temp_idx))>0; % if lever is active during cue            
            elseif isnan(bhv_times{i_trial}.States.Miss(1)) % a rewarded trial
                temp_win = [bhv_times{i_trial}.States.OppWin(1) bhv_times{i_trial}.States.CuePort1In(1)];
                temp_idx = find(WS_data.ts_1k >=temp_win(1) & WS_data.ts_1k <=temp_win(2)); % time in lever active
                Resp(i_trial) = sum(lever_active(temp_idx))>0; % if lever is active during cue
            else
                error('No such trial condition, please check code');
            end
        
        otherwise
            error('Not implemented yet or no this trial type'); 
    end
    
   clear temp_*   
   
end


Reward_rate = sum(Outcome)/n_trial; %length(cat(2,trial_ind_Cue_only, trial_ind_MaskCue, trial_ind_Cue_only_Opto, trial_ind_MaskCue_Opto));
Resp_rate = sum(Resp)/n_trial; %length(cat(2,trial_ind_Cue_only, trial_ind_MaskCue, trial_ind_Cue_only_Opto, trial_ind_MaskCue_Opto));

fprintf('Trial N=%d, rewarded %d: Respond Rate (moved during RewardWindow)=%.3f %%, Rewarded Rate=%.3f %%\n', ...
     n_trial, sum(Outcome), Resp_rate*100, Reward_rate*100);

%% rewarded time: from cue onset to reward delivery (success)
trial_ind_rew = find(Outcome == 1);
Time_rew = nan(n_trial,1);
Time_rew (trial_ind_rew)= arrayfun(@(x) bhv_times{x}.States.Reward(1) - bhv_times{x}.States.OppWin(1), trial_ind_rew);
% note: theres is a delay from press to reward ... just to be consistent
% with reward time
%% movement time: movement onset after open of opp window onset or before if moving at cue onset
% this movement time is an estimate, not use to define move before cue
% onset ...
Time_mov = nan(n_trial,1); % Time of movement before or after Reward window onset

% get lever on and off time and idx
lever_active_shift = [0;  lever_active(1:end-1)];
lever_ON_idx = find(lever_active == 1 &  lever_active_shift == 0);
lever_OFF_idx = find(lever_active == 0 &  lever_active_shift == 1);
% match the on and off
if length(lever_ON_idx) > length(lever_OFF_idx) % when terminate, moving => get rid of the last movement, not quite useful
    lever_ON_idx(end) = [];
elseif length(lever_ON_idx) < length(lever_OFF_idx)
    error ('OFF number should NOT be larger than ON number')
end

lever_ON_t = WS_data.ts_1k(lever_ON_idx);
lever_OFF_t = WS_data.ts_1k(lever_OFF_idx);

% show the inter movement duration distribution 
MV_Dur =  lever_OFF_t - lever_ON_t;
MV_interval = [lever_ON_t(2:end)-lever_OFF_t(1:end-1); NaN]; % interval is current to next movement, last one use NaN to make the legnth the same

if flag_plot
    figure;
    subplot(2,2,1); 
    cdfplot(MV_Dur);
    title('CDF: MV Dur (s)');
    subplot(2,2,3);
    hist(MV_Dur, [0:0.5:5]);
    title('Hist: MV DUr (s)');
    ylabel('MV Count');

    subplot(2,2,2); 
    cdfplot(MV_interval);
    title('CDF: inter-movement interval (s)');
    subplot(2,2,4);
    hist(MV_interval, [0:0.5:10]);
    title('Hist: inter-movement interval (s)');
    ylabel('MV interval Count');
    xlim([0 10]);
end

%
for i_trial = 1:length(bhv_times)
    temp_t = bhv_times{i_trial}.States.OppWin;
    if ~Resp(i_trial) % no movement in this trial, NaN
        Time_mov (i_trial) = NaN;
    else
        % find the first movement onset during this OppWin
        tmp_idx = find(lever_ON_t >= temp_t(1) & lever_ON_t<= temp_t(2));
        if isempty(tmp_idx) % no movement 
            Time_mov (i_trial) = NaN;
        else
            Time_mov (i_trial) = lever_ON_t(tmp_idx(1)) - temp_t(1);
        end
        clear temp_* tmp_idx
    end
end

% fprintf('Among responded trials, %.3f %% had movement before Reward window onset\n', length(find(Time_mov<0))/length(find(~isnan(Time_mov)))*100 )
Time_mov2rew = Time_rew - Time_mov;

if flag_plot
    figure;
    cdfplot(Time_mov)
    title('cdf plot of movement time');
    xlabel('Movement Onset to RewardWindow ON (nearest movement) [negatie means before window on]');

    % plot time from mvoement onset to reward for trials initiate after cue
    figure;
    cdfplot(Time_mov2rew)
    title('Reward to MoveOnset time After the open of RewardWindow')
end

%% Plot Movement around Opp window onset
if flag_plot
    pre = -10000; % ms
    post = 10000; % ms
    x_plot = (pre:post)/1000;

    LT = nan(n_trial, length(x_plot));
    LT_m = LT; % nan(length(use_lever_ON_t), length(x_plot));

    for i_t = 1:n_trial
        temp_ts_idx = round((bhv_times{i_t}.States.OppWin(1))*1000);
        
        temp_idx_plot = (temp_ts_idx+pre):(temp_ts_idx+post);
        
        if temp_idx_plot(end) > length(lever_force_smooth) ||  temp_idx_plot(1) <= 0% out of bound, throw the last trial or 1st
            
        else
        LT(i_t,:) = lever_force_smooth(temp_idx_plot);
        
        LT_m(i_t, :) = lever_force_smooth(temp_idx_plot).*lever_active(temp_idx_plot);
        end
    end
    LT_m(LT_m == 0) = NaN;
 
    figure;
    a = [];
    a(1)=subplot(2,2,1); hold on; title('Lever Trace Cue only trials aligned to CueONset: r-Opto')
    plot(x_plot, LT, 'b');
    a(2)=subplot(2,2,2);hold on; title('Lever Trace MovOnly Cue only trials: r-Opto')
    plot(x_plot, LT_m, 'b');
    a(3)=subplot(2,2,3);
    plot_jbfill_mean_se(LT, x_plot, 'b'); hold on;
    a(4)=subplot(2,2,4);
    plot_jbfill_mean_se(LT_m, x_plot, 'b'); hold on;

    linkaxes(a(3:4), 'y')
    linkaxes(a(1:2), 'y')
    linkaxes(a, 'x')
    xlim([pre post]/1000);

    % move or not plot
    figure;
    imagesc(x_plot, 1:size(LT,1), LT_m, [0.9 1.3]); colorbar;
    title('Lever Trace Cue only trials');
    colormap(hot);

    % velocity
    figure;
    a = [];
    a(1)=subplot(2,2,1); hold on; title('Lever Velocity aligned to RewardWindow Onset')
    plot(x_plot(1:end-1), diff(LT,1,2), 'b');
    a(2)=subplot(2,2,2);hold on; title('Lever Velocity MovOnly')
    plot(x_plot(1:end-1), diff(LT_m,1,2), 'b');
    a(3)=subplot(2,2,3);
    plot_jbfill_mean_se(diff(LT,1,2), x_plot(1:end-1), 'b'); hold on;
    a(4)=subplot(2,2,4);
    plot_jbfill_mean_se(diff(LT_m,1,2), x_plot(1:end-1), 'b'); hold on;
    linkaxes(a(3:4), 'y')
    linkaxes(a(1:2), 'y')
    linkaxes(a, 'x')
    xlim([pre post]/1000);


    figure;
    imagesc(x_plot(1:end-1), 1:size(LT,1), diff(LT_m,1,2), [-6 6]*10^-4); colorbar;
    title('Lever Velocity aligned to RewardWindow Onset');
    colormap(jet);
end

%% report some measure to get a sense of the behavior:

OppWin = getEventEpoch(bhv_times, 'OppWin');
% ITI length
figure;
a = [];
a(1)=subplot(4,1,1); 
plot(1:n_trial, OppWin(:,1)- ITI(:,1));
ylabel('Time (s)');
title(sprintf('%s %s %s\nITI length (s)', animalname, ses_date, proc_name));

a(2)=subplot(4,1,2); 
plot(1:n_trial, Time_mov);
ylabel('Time (s)');
title('Time from Reward window open to first movement onset (s)');

a(3)=subplot(4,1,3); 
plot(1:n_trial, Time_rew);
ylabel('Time (s)');
title('Time from Reward window open to reward (s)');

a(4)=subplot(4,1,4); 
plot(1:n_trial, Time_mov2rew);
ylabel('Time (s)');
title('Time from movement onset to reward (s)');
xlabel('Trial #');

linkaxes(a, 'x');

if flag_plot
    figure;
    subplot(2,2,1);
    cdfplot(OppWin(:,1)- ITI(:,1));
    title('CDF: ITI length (s)');
    subplot(2,2,2);
    cdfplot(Time_mov);
    title('CDF: T RewardWinON to 1st Mov');
    subplot(2,2,3);
    cdfplot(Time_rew);
    title('CDF: T RewardWinON to Rew (s)');
    subplot(2,2,4);
    cdfplot(Time_mov2rew);
    title('CDF: T MovON to Rew (s)');
end

figure; a= [];
a(1)=subplot(2,2,1);
hist(OppWin(:,1)- ITI(:,1), [0:1:10]);
title(sprintf('%s %s %s\nHist: ITI length (s)', animalname, ses_date, proc_name));
ylabel('Trial Count')
a(2)=subplot(2,2,2);
hist(Time_mov, [0:1:10]);
title('Hist: T RewardWinON to 1st Mov');
a(3)=subplot(2,2,3);
hist(Time_rew, [0:1:10]);
title('Hist: T RewardWinON to Rew (s)');
a(4)=subplot(2,2,4);
hist(Time_mov2rew, [0:1:10]);
title('Hist: T MovON to Rew (s)');
linkaxes(a, 'x');
xlim([-0.5 10.5])


% save and exit
if exist(fullfile(save_path, animalname, ses_date),'dir') ~=7
   mkdir(fullfile(save_path, animalname, ses_date)) ;
end

save(fullfile(save_path, animalname, ses_date, [proc_name '.mat']),...
    'mp4_name','animalname', 'ses_date', 'proc_name', 'Bpod_fn',  'WS_fn', ...
    'Outcome', 'Resp', 'n_trial', 'Time_*', 'lever_*', 'idx_*', '*_rate' ,...
    'SessionData', 'bhv_times', 'WS_data', 'WS_trial', 'movethresh*', 'ITI', 'Reward', 'MV_*', '-v7.3');
% 'lever*',  'Lever*',,

disp(['Saved: ' fullfile(save_path, animalname, ses_date, [proc_name '.mat'])]);
Data.mp4_name           = mp4_name;
Data.animalname         = animalname;
Data.ses_date           = ses_date;
Data.proc_name          = proc_name;
Data.Bpod_fn            = Bpod_fn;
Data.WS_fn              = WS_fn;
Data.Outcome            = Outcome;
Data.n_trial            = n_trial;

Data.Reward_rate        = Reward_rate;
Data.SessionData        = SessionData;

Data.Resp               = Resp;
Data.n_trial            = n_trial;
Data.trial_ind_rew      = trial_ind_rew;
Data.Time_mov           = Time_mov;
Data.Time_mov2rew       = Time_mov2rew;
Data.Time_rew           = Time_rew;
Data.lever_active       = lever_active;
Data.lever_active_shift = lever_active_shift;
Data.lever_force_smooth = lever_force_smooth;
Data.lever_OFF_idx      = lever_OFF_idx;
Data.lever_OFF_t        = lever_OFF_t;
Data.lever_ON_idx       = lever_ON_idx;
Data.lever_ON_t         = lever_ON_t;
Data.lever_velocity_envelope_smooth = lever_velocity_envelope_smooth;
Data.MV_Dur = MV_Dur;
Data.MV_interval = MV_interval;

Data.idx_cam_ch         = idx_cam_ch;
Data.idx_Lever_ch       = idx_Lever_ch;

Data.SessionData        = SessionData;
Data.bhv_times          = bhv_times;
Data.WS_data            = WS_data;
Data.WS_trial           = WS_trial;
Data.movethresh         = movethresh;

Data.lever_active_20    = lever_active_20;
Data.lever_force_smooth_20 = lever_force_smooth_20;
Data.lever_force_resample = lever_force_resample;
Data.lever_velocity_resample_smooth_20 = lever_velocity_resample_smooth_20;
Data.movethresh_20 = movethresh_20;
return;