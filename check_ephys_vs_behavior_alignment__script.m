
% Place behavior and ephys files in the same folder 'folder_root'

clear

% folder_root = 'Y:\Users\ariadna\ephys\h-tester-align-data\h-tester_micLn9\2024-03-22';
% folder_root = 'Y:\Users\ariadna\ephys\h-tester-align-data\h-test_mic9_\2024-03-18';
% folder_root = 'Y:\Users\ariadna\ephys\h-tester-align-data\h-test_ttlcam_rest\2024-03-18';
% folder_root = 'Y:\Users\ariadna\ephys\h-tester-align-data\h-test_ttlcam\2024-03-18';
% folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\trial03';
% folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-25';
% folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-26\trial001'; % Misalignment: ev_dt = -2.1617 msec
% folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-26\trial002'; % Misalignment: ev_dt = -2.1308 msec
% folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-26\trial003'; % Misalignment: ev_dt = -2.1931 msec
% folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-29\trial001'; % Misalignment: ev_dt = -2.1573 msec
% folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-29\trial002';
% folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-30\trial002';
% folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-30\trial003';
folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-30\trial004_buf512';  % Misalignment: ev_dt = -2.1725 msec
folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-30\trial005_buf8192'; % Misalignment: ev_dt = -1.0617 msec
folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-30\trial006_buf128';  % Misalignment: ev_dt = 0.2756 msec
% folder_root = 'Z:\Users\Alessandro La Chioma\Ariadna\ale_tmp\smart_mouse\2024-04-30\trial007_buf128';  % Misalignment: ev_dt = ~-1.5 msec



%% Load data.mat file

datamat_ls = dir([folder_root filesep '*_data.mat']);

D = load( fullfile(datamat_ls.folder, datamat_ls.name));

%%

% fmt = '.wav';
fmt = '.flac';
fmt = '.mic'; % binary file

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

MicSamps_all = [1 : n_Samps]';


% if contains(filename_mic, 'audiorec')
%     MicNrSamples  = D.audio_rec.MicNrSamples;
%     MicTimeStamps = D.audio_rec.MicTimeStamps;
% elseif contains(filename_mic, 'TTLcam')
%     MicNrSamples  = D.audio_rec.ttl.MicNrSamples;
%     MicTimeStamps = D.audio_rec.ttl.MicTimeStamps;
% end
% tCaptureStart = MicTimeStamps(1);
% T_all = (MicSamps_all-1) / Fs_mic + tCaptureStart;


% % MicSamps = cumsum(MicNrSamples) - MicNrSamples(1) + 1;
% Now, MicTimeStamps gives you the timestamps (according to Behavior PC
% clock) corresponding to the audio file samples indicated in MicSamps.
% % T_all = interp1(MicSamps, MicTimeStamps, MicSamps_all, 'linear', 'extrap');


% MicSamps = [1; cumsum(MicNrSamples(1:end-1)) + 1];
% T_all = interp1(MicSamps, MicTimeStamps, MicSamps_all, 'linear', 'extrap');


T_all = (MicSamps_all-1) / Fs_mic;

%% Load TTL ephys and get timestamp 


%     ttl_ephys_type = 'dedicated';
    ttl_ephys_type = 'mic';

switch ttl_ephys_type
    case 'dedicated'

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
        %%% If ephys TTL is fed to one of the four mic channels:
            MicTimeStamps = D.audio_rec.MicTimeStamps;
            MicNrSamples  = D.audio_rec.MicNrSamples;
        %%%
        ttlEphysTimeStamps = nan(nr_detectedPulses,1);
        for f = 1 : nr_detectedPulses
            samp = locs_down(f); % frame onset in sound card samples
            ind_MicNrSamples = find( (samp - MicNrSamples)> 0, 1,"first");
            d_samps = samp - MicNrSamples(ind_MicNrSamples);
            ttlEphysTimeStamps(f) = MicTimeStamps(ind_MicNrSamples) + d_samps/D.audio_rec.ttl_ephys.SampleRate;
        end
        D.audio_rec.ttl_ephys.TimeStamps = ttlEphysTimeStamps;
        D.audio_rec.ttl_ephys.MicSamples = MicSamples;
    
    case 'mic'

        ttl_eph = y(:,2);
        
        figure;plot(ttl_eph)
        Fs = Fs_mic;
        val_low  = 0;
        val_high = max(ttl_eph);
        dval_min = val_high/2;
        val_min  =  dval_min;
        
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
        % MicTimeStamps  = D.audio_rec.ttl_ephys.MicTimeStamps;
        % MicNrSamples   = D.audio_rec.ttl_ephys.MicNrSamples;
        %%% If ephys TTL is fed to one of the four mic channels:
        % MicTimeStamps = D.audio_rec.MicTimeStamps;
        % MicNrSamples  = D.audio_rec.MicNrSamples;
        %%%
        ttlEphysTimeStamps = nan(nr_detectedPulses,1);
        for f = 1 : nr_detectedPulses
            samp = locs_down(f); % frame onset in sound card samples
            ttlEphysTimeStamps(f) = T_all(samp);
        end

end


figure
k = Fs/4; % 1000;
inds_toPlot = MicSamples+[-k:k];
p = plot(inds_toPlot, ttl_eph(inds_toPlot));
p.DataTipTemplate.DataTipRows(1).Format = '%d'; % set precision of any datatip you add to 1 unit (by default it could be 5 units)

% Find TTL ephys offset, so to calculate pulse duration:
[~,locs_down2] = findpeaks(-ttl_eph(samp:samp+2*Fs), "MinPeakDistance",Fs/min_dist,...
    "MinPeakProminence",dval_min, "MinPeakHeight",val_min*0.7); 
samp_offset = locs_down2(1)+samp-1;
% % Pulse duration (manually enter the sample nr of the offset):
% samp_offset = 3843411;
hold on;
plot([samp samp_offset], [ttl_eph(samp) ttl_eph(samp_offset)], 'v');

ttlpulse_dur = abs(samp - samp_offset)/Fs;
fprintf('ttlpulse_dur = %5.4f msec\n\n', ttlpulse_dur*1000);

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
% t = fread(fid,1,'uint64=>uint64'); % read timestamp from start of file
yephys = fread(fid,[n_channels Inf],'int16=>single')';
fclose(fid);
fprintf('Ephys data .bin file loaded!\n')
% Convert values to microvolts:
yephys = yephys*6.25e3/32768;
% To save memory, you can also convert to int16 (2 bytes per digit, vs 4
% bytes for single and 8 for double):
% yephys = int16(yephys);

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

figure;
plot(tvec_ephys, data_ephys)
ylabel('Ephys trace')
xlabel('Time (s, ephys clock)')

%%

trange_ephys = [26.5 26.57];
trange_ephys = [0.05 0.1];
trange_ephys = [56.23 56.25];
trange_ephys = [5.9 6];

trange_ephys_samp = trange_ephys * Fs_ephys;
trange_ephys_ptb  = trange_ephys + ttlEphysTimeStamps;


%%



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

figure
p = plot(sampvec_ephys, data_ephys);
p.DataTipTemplate.DataTipRows(1).Format = '%d'; % set precision of any datatip you add to 1 unit (by default it could be 5 units)
ylabel('Ephys trace')
xlabel('Samples (ephys samps)')
xlim(trange_ephys_samp)
title(sprintf('Use this plot to choose one sample point as you like.'));

fprintf('\nUse the plot with Ephys data as a function of Samples to choose one sample point as you like.\n') 
fprintf('The same point will be found in the behavior data stream (mic or ttlcam) and used for checking alignment.\n')
fprintf('Enter the sample number (x value) of this point in the variable ''ev_samp''\n\n');

%%

% ev_samp   = 662791;
% ev_samp   = 1563899;
% ev_samp   = 405482;
% ev_samp   = 139705;
% ev_samp   = 230350;
ev_samp   = 148727;

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
title(sprintf('Use this plot to find the same point you chose before.'));

fprintf('\nUse the plot with Behavior data as a function of Samples to find the same point you chose before.\n') 
fprintf('Enter the sample number (x value) of this point in the variable ''ev_mic_samp''\n\n');

%% Get misalignment between ephys and behavior

ev_mic_samp = 7998175;
ev_mic_samp = 14918705;
ev_mic_samp = 3175352;

ev_mic_tptb = T_all(ev_mic_samp);

ev_dt = ev_tptb - ev_mic_tptb;
fprintf('Misalignment: ev_dt = %5.4f msec\n\n', ev_dt*1000);


%% Plot everything together

data_mic = y(:,1);
ttl_eph  = y(:,2);

MicSamps_all = [1 : n_Samps]';
tvec_mic = (MicSamps_all-1) / Fs_mic;

% Get range of data_mic to rescale ephys later on:
lims = prctile(data_mic, [0.1 99.9]);
lims = [-max(abs(lims)) max(abs(lims))];

figure;
hold on
p1 = plot(tvec_mic, data_mic);
% plot(tvec_mic, ttl_eph)
p1.DataTipTemplate.DataTipRows(1).Format = '%.5f'; % set precision of any datatip


data_ephys      = yephys(:,1);
% Rescale ephys data:
lims_ephys = prctile(data_ephys, [1 99]);
% lims_ephys = [min(data_ephys), max(data_ephys)];
data_ephys_norm = (data_ephys - lims_ephys(1)) / (lims_ephys(2)-lims_ephys(1)) * 2 - 1;
data_ephys_norm = data_ephys_norm * lims(2);

sampvec_ephys  = [1 : length(data_ephys)]';
tvec_ephys     = [sampvec_ephys-1] / Fs_ephys;
tvec_ephys2mic = tvec_ephys + ttlEphysTimeStamps;

plot([ttlEphysTimeStamps ttlEphysTimeStamps],[-1 1]*lims(2)*10,'k--')
p2 = plot(tvec_ephys2mic, data_ephys_norm);

p2.DataTipTemplate.DataTipRows(1).Format = '%.5f'; % set precision of any datatip

% ylim(lims)
ylim(lims*2)


