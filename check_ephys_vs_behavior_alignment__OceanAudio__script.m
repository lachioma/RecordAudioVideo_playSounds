
% Script for when you do NOT use RecordAudioVideo_playSounds, but you use
% OceanAudio to record mics and ephys TTL and copy of ephys test signal.
% 
% Alessandro La Chioma - 2024-10-18


% Place behavior and ephys files in the same folder 'folder_root'

clear

% OceanAudio
folder_root = 'Y:\Users\ariadna\ephys\h-tester\2024-10-18\test1'; % Misalignment: ev_dt = -2.1433 msec
folder_root = 'Y:\Users\ariadna\ephys\h-tester\2024-10-18\test2'; % Misalignment: ev_dt = -2.2090 msec
folder_root = 'Y:\Users\ariadna\ephys\h-tester\2024-10-18\test3'; % Misalignment: ev_dt = -2.1975 msec
% Play_BatterySounds_20241022.mlapp
folder_root = 'Y:\Users\ariadna\ephys\h-tester\2024-10-22\test1'; % Misalignment: ev_dt = -2.2010 msec
folder_root = 'Y:\Users\ariadna\ephys\h-tester\2024-10-22\test2'; % Misalignment: ev_dt = -2.2137 msec
folder_root = 'Y:\Users\ariadna\ephys\h-tester\2024-10-22\test3'; % Misalignment: ev_dt = -2.1237 msec



%%

fmt = '.wav';
fmt = '.flac';

ls_mic = dir([folder_root filesep '*' fmt]);

path_mic = ls_mic(1).folder;
filename_mic = ls_mic(1).name;

filename_datetime_tag = filename_mic(1:end-7-length(fmt));

[y,Fs_mic] = audioread(fullfile(path_mic,filename_mic));
n_Samps = length(y);

data_mic = y(:,1);


%% Get timestamps of microphone audio data (PTB clock)

MicSamps_all = [1 : n_Samps]';

tCaptureStart = 0;
T_all = (MicSamps_all-1) / Fs_mic + tCaptureStart;


%% Load TTL ephys and get timestamp 


    ttl_ephys_type = 'dedicated'; % _TTLephys.flac
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
            % MicTimeStamps = D.audio_rec.MicTimeStamps;
            % MicNrSamples  = D.audio_rec.MicNrSamples;
        %%%
        ttlEphysTimeStamps = nan(nr_detectedPulses,1);
        for f = 1 : nr_detectedPulses
            samp = locs_down(f); % frame onset in sound card samples
            % ind_MicNrSamples = find( (samp - MicNrSamples)> 0, 1,"first");
            % d_samps = samp - MicNrSamples(ind_MicNrSamples);
            % ttlEphysTimeStamps(f) = MicTimeStamps(ind_MicNrSamples) + d_samps/D.audio_rec.ttl_ephys.SampleRate;
            ttlEphysTimeStamps(f) = MicTimeStamps(1) + (samp-1)/D.audio_rec.ttl_ephys.SampleRate;
        end
        D.audio_rec.ttl_ephys.TimeStamps = ttlEphysTimeStamps;
        D.audio_rec.ttl_ephys.MicSamples = MicSamples;
    
    case 'mic'

        ttl_eph = y(:,6);
        
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
if length(ephysdata_binfile_ls) > 1
    warning('Multiple ephys bin files found!');
    keyboard;
    ephysdata_binfile_ls = ephysdata_binfile_ls(1);
elseif length(ephysdata_binfile_ls) < 1
    error('No ephys bin files found!');
end
ephysdata_binfile    = fullfile(ephysdata_binfile_ls.folder, ephysdata_binfile_ls.name);

sr_str = char(regexp(ephysdata_binfile,'_\d{4,5}sps', 'match'));
sr = str2double(sr_str(2:strfind(sr_str,'sps')-1));
Fs_ephys = sr;

ch_str = char(regexp(ephysdata_binfile,'_\d{1,4}ch_', 'match'));
n_channels = str2double(ch_str(2:strfind(ch_str,'ch')-1));


bytesToSkip = 8; % Skip the first 8 bytes
fid = fopen(ephysdata_binfile, 'r');

fseek(fid, 0, 'eof');  % Move to the end of the file
fileSize_bytes  = ftell(fid);  % Get the current position, which is the size in bytes
fileSize_Gbytes = fileSize_bytes/10^9;
n_bytes_perInt16 = 2;
n_data_points = (fileSize_bytes - bytesToSkip)/n_channels/n_bytes_perInt16;
assert(mod(n_data_points,1)==0, 'The number of data points has decimals, so something is wrong.')

fseek(fid, bytesToSkip, 'bof'); % Skip the first 8 bytes
% t = fread(fid,1,'uint64=>uint64'); % read timestamp from start of file

% n_channels_toRead = 5; % Load only the first 'n_channels_toRead' channels
n_channels_toRead = n_channels; % Load all channels

if fileSize_Gbytes > 10 && n_channels_toRead > 30
    warning('Are you sure you want to proceed? This might take long to load and even result in out of memory.')
    n_channels_toRead = 1;
    keyboard;
end

yephys = fread(fid,[n_channels_toRead n_data_points],'int16=>single', n_bytes_perInt16*(n_channels-n_channels_toRead))';
% yephys = fread(fid,[n_channels Inf],'int16=>single')';

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

% trange_ephys = [26.5 26.57];
trange_ephys = [0.05 0.1];
% trange_ephys = [56.23 56.25];
trange_ephys = [15.1 15.3];
trange_ephys = [20.1 20.2];

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

% ev_samp   = 190790;
ev_samp   = 284232;
ev_samp   = 378942;
ev_samp   = 640766;
ev_samp   = 1064094;
ev_samp   = 503144;

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

ev_mic_samp = 2641716;
ev_mic_samp = 3855044;
ev_mic_samp = 4702094;
ev_mic_samp = 5178640;
ev_mic_samp = 6532205;
ev_mic_samp = 10669890;
ev_mic_samp = 6866719;

ev_mic_tptb = T_all(ev_mic_samp);

ev_dt = ev_tptb - ev_mic_tptb;
fprintf('Misalignment: ev_dt = %5.4f msec\n\n', ev_dt*1000);


%% Plot everything together


figure;
hold on
% Mic data with test ephys data 
m1 = prctile(data_mic,99);
p = plot(T_all, data_mic/m1);
p.DataTipTemplate.DataTipRows(1).Format = '%.4f'; % set precision of any datatip you add to 1 unit (by default it could be 5 units)
% Test ephys data 
m2 = prctile(data_ephys,99);
plot(tvec_ephys_ptb, data_ephys/m2)
% Ephys TTL pulse
plot(T_all, ttl_eph)
title('Run the next section if you want to measure the dT between 2 points!')
disp('Run the next section if you want to measure the dT between 2 points!')

% % data_mic = y(:,1);
% % % ttl_eph  = y(:,2);
% % 
% % MicSamps_all = [1 : n_Samps]';
% % tvec_mic = (MicSamps_all-1) / Fs_mic;
% % 
% % % Get range of data_mic to rescale ephys later on:
% % lims = prctile(data_mic, [0.1 99.9]);
% % lims = [-max(abs(lims)) max(abs(lims))];
% % 
% % figure;
% % hold on
% % p1 = plot(tvec_mic, data_mic);
% % % plot(tvec_mic, ttl_eph)
% % p1.DataTipTemplate.DataTipRows(1).Format = '%.5f'; % set precision of any datatip
% % 
% % 
% % data_ephys      = yephys(:,1);
% % % Rescale ephys data:
% % lims_ephys = prctile(data_ephys, [1 99]);
% % % lims_ephys = [min(data_ephys), max(data_ephys)];
% % data_ephys_norm = (data_ephys - lims_ephys(1)) / (lims_ephys(2)-lims_ephys(1)) * 2 - 1;
% % data_ephys_norm = data_ephys_norm * lims(2);
% % 
% % sampvec_ephys  = [1 : length(data_ephys)]';
% % tvec_ephys     = [sampvec_ephys-1] / Fs_ephys;
% % tvec_ephys2mic = tvec_ephys + ttlEphysTimeStamps;
% % 
% % plot([ttlEphysTimeStamps ttlEphysTimeStamps],[-1 1]*lims(2)*10,'k--')
% % p2 = plot(tvec_ephys2mic, data_ephys_norm);
% % 
% % p2.DataTipTemplate.DataTipRows(1).Format = '%.5f'; % set precision of any datatip
% % 
% % % ylim(lims)
% % ylim(lims*2)

%% Measure Time between 2 selected points

% Get two points interactively
[x, y] = ginput(2);
% Calculate the Euclidean distance between the two points
% distance = sqrt((x(2) - x(1))^2 + (y(2) - y(1))^2);
distance = x(2) - x(1);
% Display the distance
disp(['Time between the 2 selected points: ', num2str(distance*1000), ' ms']);


