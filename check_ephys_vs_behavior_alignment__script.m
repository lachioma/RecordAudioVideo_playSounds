
clear

folder_root = 'Y:\Users\ariadna\ephys\h-tester-align-data\h-tester_micLn9\2024-03-22';
% folder_root = 'Y:\Users\ariadna\ephys\h-tester-align-data\h-test_mic9_\2024-03-18';
% folder_root = 'Y:\Users\ariadna\ephys\h-tester-align-data\h-test_ttlcam_rest\2024-03-18';
% folder_root = 'Y:\Users\ariadna\ephys\h-tester-align-data\h-test_ttlcam\2024-03-18';

%% Load data.mat file

datamat_ls = dir([folder_root filesep '*_data.mat']);

D = load( fullfile(datamat_ls.folder, datamat_ls.name));

%%

% fmt = '.wav';
fmt = '.flac';
% fmt = '.mic'; % binary file

ls_mic = dir([folder_root filesep '*_audiorec' fmt]);
% ls_mic = dir([folder_root filesep '*_TTLcam' fmt]);
if length(ls_mic) == 1
    path_mic = ls_mic(1).folder;
    filename_mic = ls_mic(1).name;
else
    [filename_mic, path_mic] = uigetfile([folder_root filesep '*_audiorec' fmt]);
end
filename_datetime_tag = filename_mic(1:end-7-length(fmt));

if strcmp(fmt, '.mic')
    fid = fopen(fullfile(path_mic,filename_mic), 'r');
    % Read all data from the audio file in a matrix
    nAudioChannels_input_toSave = 4;
    % Read the whole data matrix:
    y = fread(fid,[nAudioChannels_input_toSave Inf],'single=>single')';
    % n_Samps = size(y,1);
    % Just get the nr. data points from the filesize:
    fseek(fid, 0, 'eof');
    filesize = ftell(fid);
    n_Samps = filesize/ 4 / nAudioChannels_input_toSave; % single precision takes 4 bytes per value
    % Close the file
    fclose(fid);
    Fs_mic = 192000;
else
    [y,Fs_mic] = audioread(fullfile(path_mic,filename_mic));
    n_Samps = length(y);
end

data_mic = y(:,1);


%% Get timestamps of microphone audio data (PTB clock)

if contains(filename_mic, 'audiorec')
    MicNrSamples  = D.audio_rec.MicNrSamples;
    MicTimeStamps = D.audio_rec.MicTimeStamps;
elseif contains(filename_mic, 'TTLcam')
    MicNrSamples  = D.audio_rec.ttl.MicNrSamples;
    MicTimeStamps = D.audio_rec.ttl.MicTimeStamps;
end
MicSamps_all = [1 : n_Samps]';
% % MicSamps = cumsum(MicNrSamples) - MicNrSamples(1) + 1;
% Now, MicTimeStamps gives you the timestamps (according to Behavior PC
% clock) corresponding to the audio file samples indicated in MicSamps.
% % T_all = interp1(MicSamps, MicTimeStamps, MicSamps_all, 'linear', 'extrap');

% % MicSamps = [1; cumsum(MicNrSamples(1:end-1)) + 1];
% % T_all = interp1(MicSamps, MicTimeStamps, MicSamps_all, 'linear', 'extrap');


tCaptureStart = MicTimeStamps(1);
% T_all = (MicSamps_all-1) / Fs_mic;
T_all = (MicSamps_all-1) / Fs_mic + tCaptureStart;

%% Load TTL ephys and get timestamp 

ttlephys_ls = dir([folder_root filesep '*_TTLephys.flac']);

ttl_ephys.fullpath_flac = fullfile(ttlephys_ls.folder, ttlephys_ls.name);

disp(['Extracting TTL ephys timestamps from audio channel '  '...']) 

[recaudiodata,Fs] = audioread(ttl_ephys.fullpath_flac);
% Fs = D.audio_rec.ttl_ephys.SampleRate;
ttl_eph = recaudiodata(:);
clear recaudiodata
                    
val_low  =  D.audio_rec.ttl_ephys.val_low;
val_high =  D.audio_rec.ttl_ephys.val_high;
dval_min =  D.audio_rec.ttl_ephys.dval_min;
val_min  =  val_high;
                    
                    
% Extract timestamps of White Matter Ephys TTLs (in microphone samples)

min_dist = 1; % Minimum distance of ttl onsets, in Hz
% Find peak (% locs contains the timestamps of the
% TTL onsets (in microphone samples))
[pks,locs_down] = findpeaks(ttl_eph, "MinPeakDistance",Fs/min_dist,...
    "MinPeakProminence",dval_min, "MinPeakHeight",val_min); 
nr_detectedPulses = length(locs_down);

if nr_detectedPulses < 1
    fprintf('Warning: No TTL Ephys pulses found! \n');
elseif nr_detectedPulses > 1
    fprintf('!!! Error: Nr. %d TTL Ephys pulses found.\n', nr_detectedPulses);
else
    fprintf('TTL Ephys pulse found, good. \n');
end

MicSamples     = locs_down;
MicTimeStamps  = D.audio_rec.ttl_ephys.MicTimeStamps;
MicNrSamples   = D.audio_rec.ttl_ephys.MicNrSamples;
ttlEphysTimeStamps = nan(nr_detectedPulses,1);
for f = 1 : nr_detectedPulses
    samp = locs_down(f); % frame onset in sound card samples
    ind_MicNrSamples = find( (samp - MicNrSamples)> 0, 1,"first");
    d_samps = samp - MicNrSamples(ind_MicNrSamples);
    ttlEphysTimeStamps(f) = MicTimeStamps(ind_MicNrSamples) + d_samps/D.audio_rec.ttl_ephys.SampleRate;
end
D.audio_rec.ttl_ephys.TimeStamps = ttlEphysTimeStamps;
D.audio_rec.ttl_ephys.MicSamples = MicSamples;
    
figure
k = 1000;
inds_toPlot = MicSamples+[-k:k];
plot(inds_toPlot, ttl_eph(inds_toPlot))

%% Load ephys data 

%%% The original ephys data .bin file was separately saved into .mat file
% ephysdata = load('Y:\Users\ariadna\ephys\h-tester-align-data\h-tester_micLn9\2024-03-22\HSW_2024_03_22__16_41_25__01min_06sec__hsamp_64ch_25000sps.mat');
% data_ephys = ephysdata.data(1,:)';
% Fs_ephys   = double(ephysdata.sr);
% clear ephysdata


%%% Use directly the ephys data .bin file %%%%%%%%

% ephysdata_binfile = "Y:\Users\ariadna\ephys\h-tester-align-data\h-tester_micLn9\2024-03-22\HSW_2024_03_22__16_41_25__01min_06sec__hsamp_64ch_25000sps.bin";
ephysdata_binfile_ls = dir([folder_root filesep '*sps.bin']);
ephysdata_binfile    = fullfile(ephysdata_binfile_ls.folder, ephysdata_binfile_ls.name);

sr_str = char(regexp(ephysdata_binfile,'_\d{4,5}sps', 'match'));
sr = str2double(sr_str(2:strfind(sr_str,'sps')-1));
Fs_ephys = sr;

ch_str = char(regexp(ephysdata_binfile,'_\d{1,4}ch_', 'match'));
n_channels = str2double(ch_str(2:strfind(ch_str,'ch')-1));

fid = fopen(ephysdata_binfile, 'r');
fseek(fid, 8, 'bof'); % Skip the first 8 bytes
yephys = fread(fid,[n_channels Inf],'int16=>double')';
fclose(fid);
fprintf('Ephys data .bin file loaded!\n')
% Convert values to microvolts:
yephys = int16(yephys*6.25e3/32768);

% figure;
% hold on;
% plot(data_ephys)
% plot(yephys(:,1))
% 
% figure;
% plot(data_ephys-yephys(:,1))

data_ephys = yephys(:,1);


sampvec_ephys  = [1 : length(data_ephys)]';
tvec_ephys     = [sampvec_ephys-1] / Fs_ephys;
tvec_ephys_ptb = tvec_ephys + ttlEphysTimeStamps;


%%

trange_ephys = [26.5 26.57];
trange_ephys_samp = trange_ephys * Fs_ephys;
trange_ephys_ptb  = trange_ephys + ttlEphysTimeStamps;


%%

figure
p = plot(sampvec_ephys, data_ephys);
p.DataTipTemplate.DataTipRows(1).Format = '%d'; % set precision of any datatip you add to 1 unit (by default it could be 5 units)
ylabel('Ephys trace')
xlabel('Samples (ephys samps)')
xlim(trange_ephys_samp)

figure
plot(tvec_ephys, data_ephys)
ylabel('Ephys trace')
xlabel('Time (s, ephys clock)')
xlim(trange_ephys)

figure
plot(tvec_ephys_ptb, data_ephys)
ylabel('Ephys trace')
xlabel('Time (s, ptb clock)')
xlim(trange_ephys_ptb)

fprintf('\nUse the plot with Ephys data as a function of Samples to choose one sample point as you like.\n') 
fprintf('The same point will be found in the behavior data stream (mic or ttlcam) and used for checking alignment.')
fprintf('Enter the sample number (x value) of this point in the variable ''ev_samp''\n\n');

%%

ev_samp   = 662791;
ev_samp   = 1563899;
% ev_samp   = 405482;
% ev_samp   = 139705;
% ev_samp   = 230350;
ev_tephys = tvec_ephys(ev_samp);
ev_tptb   = tvec_ephys_ptb(ev_samp);

trange = [-0.1 0.1];
trange_tephys = trange + ev_tephys;
trange_tptb   = trange + ev_tptb;

figure;
hold on;
plot(tvec_ephys_ptb, data_ephys)
plot(ev_tptb       , data_ephys(ev_samp), 'v')
ylabel('Ephys data')
xlabel('Time (s, ptb clock)')
xlim(trange_tptb)



%%

[d1, s1] = min(abs(T_all-trange_tptb(1)));
[d2, s2] = min(abs(T_all-trange_tptb(2)));
trange_mic_samp = [s1 s2];
trange_mic_ptb = T_all(trange_mic_samp);

figure
plot(T_all, data_mic);
ylabel('Behavior data (mic or ttlcam)')
xlabel('Time (s, ptb clock)')
xlim(trange_mic_ptb)

figure
p = plot(MicSamps_all, data_mic);
ylabel('Behavior data (mic or ttlcam)')
xlabel('Samples (samps)')
xlim(trange_mic_samp)
p.DataTipTemplate.DataTipRows(1).Format = '%d'; % set precision of any datatip you add to 1 unit (by default it could be 5 units)

fprintf('\nUse the plot with Behavior data as a function of Samples to find the same point you chose before.\n') 
fprintf('Enter the sample number (x value) of this point in the variable ''ev_mic_samp''\n\n');

%%

ev_mic_samp = 7998175;
ev_mic_samp = 14918705;
% ev_mic_samp = 4804690;
ev_mic_tptb = T_all(ev_mic_samp);

ev_dt = ev_tptb - ev_mic_tptb;
fprintf('dt_ev = %5.4f msec\n\n', ev_dt*1000);

%%
ttlpulse_dur = abs(29102 - 29296)/D.audio_rec.ttl_ephys.SampleRate;
fprintf('ttlpulse_dur = %5.4f msec\n\n', ttlpulse_dur*1000);

