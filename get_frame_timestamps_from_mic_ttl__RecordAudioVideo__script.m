
%% Select and load audio file with FLIR camera TTL

clear

folder_root = 'Y:\Users\ariadna\behavior_PP\25681\25681_3_sounds\2023-11-24';

fmt = '.wav';
% fmt = '.flac';
% fmt = '.ttlcam'; % binary file

ls_mic = dir([folder_root filesep '*_TTLcam' fmt]);
if length(ls_mic) == 1
    path_mic = ls_mic(1).folder;
    filename_mic = ls_mic(1).name;
else
    [filename_mic, path_mic] = uigetfile([folder_root filesep '*_TTLcam' fmt]);
end
filename_datetime_tag = filename_mic(1:end-7-length(fmt));

if strcmp(fmt, '.ttlcam')
    fid = fopen(fullfile(path_mic,filename_mic), 'r');
    % Read all data from the audio file in a vector
    nAudioChannels_input_toSave = 1;
    y = fread(fid,[nAudioChannels_input_toSave Inf],'single=>single')';
    % Close the file
    fclose(fid);
    Fs = 192000/50;
else
    [y,Fs] = audioread(fullfile(path_mic,filename_mic));
end

fprintf('File loaded: %s\n', fullfile(path_mic,filename_mic))

%% Load .csv camera timestamps

% path_vidcsv     = path_mic;
% filename_vidcsv = [filename_datetime_tag '_cam1.csv'];

% Column 1: camera timestamps of each frame according to the FLIR camera internal clock (in nanoseconds)
% Column 2: camera timestamps (in Bonsai clock) of every time Bonsai received a frame from the camera (in seconds)
% Column 3: frame IDs as given by the FLIR camera processor (this might be useful to identify dropped frames)

% T = readmatrix(fullfile(path_vidcsv,filename_vidcsv), ...
%         'OutputType','double', 'Delimiter',',');

path_data     = path_mic;
filename_data = [filename_datetime_tag '_data.mat'];

D = load( fullfile(path_data, filename_data));

T(:,1) = D.cam.camflir.Timestamp_camera_clock_us;
T(:,2) = D.cam.camflir.Timestamp_bonsai_clock_s;
T(:,3) = D.cam.camflir.FrameID;


%% Extract timestamps of camera TTLs (in microphone samples)

% First assess the shape of the TTL, i.e. the values of the HIGH and LOW
% levels. These values change depending on the camera frame rate and also
% on the resistance applied by the potentiometer.

ttl = single(y(:,1));

% Take a piece from the entire ttl trace, that is 2 seconds long:
nr_samps = Fs*2;
% Take this piece from the end of the trace, because in the beginning the
% ttl changes shape. Exclude the very last 1 second because there should be
% no ttl at the very end:
final_samps_toExclude = Fs*10;

ttl_piece = ttl(end-final_samps_toExclude-nr_samps+1:end-final_samps_toExclude);

% val_low  = mean(ttl_piece(ttl_piece<0));
val_low  = prctile(ttl_piece(ttl_piece<0), 25);
val_high = mean(ttl_piece(ttl_piece>0));
% val_low  = -0.02;
% val_high =  0.01;
dval_min = (val_high - val_low)/2;
val_min  = val_high;


%%% Extract timestamps of camera TTLs (in microphone samples)

min_dist = 200; % Minimum distance of ttl onsets, in Hz (set reasonably > than actual camera frame rate)

[pks,locs] = findpeaks(ttl, "MinPeakDistance",Fs/min_dist,...
    "MinPeakProminence",dval_min ...
    , "MinPeakHeight",val_min...
    );

% locs contains the timestamps of the TTL onsets (in microphone samples)

figure;
hold on;
plot(ttl)
plot(locs,pks, 'v');
xlabel('Audio samples')
% xlim([0 Fs*2]) % first 2 seconds
% xlim([length(ttl)-[Fs*2, 0]]) % last 2 seconds

figure
histogram(diff(locs)/Fs*1000)
xlabel('Frame time (ms)')


%% Get duration, frame time, and frame rate according to all clocks


dur_snd = (locs(end) - locs(1)) / Fs;
dur_cam = (T(end,1) - T(1,1)) / 1e9;
dur_bon = (T(end,2) - T(1,2));

frametime_snd = mean(diff(locs)/Fs);
frametime_cam = mean(diff(T(1:end,1))) / 1e9;
frametime_bon = mean(diff(T(1:end,2)));

framerate_snd = 1/frametime_snd;
framerate_cam = 1/frametime_cam;
framerate_bon = 1/frametime_bon;

%% Get timestamps of FLIR camera TTL (PTB clock)

MicNrSamples  = D.audio_rec.ttl.MicNrSamples;
MicTimeStamps = D.audio_rec.ttl.MicTimeStamps;
MicSamps = cumsum(MicNrSamples) - MicNrSamples(1) + 1;
% Now, MicTimeStamps gives you the timestamps (according to Behavior PC
% clock) corresponding to the audio file samples indicated in MicSamps.
MicSamps_all = [1 : length(y)]';

T_all = interp1(MicSamps, MicTimeStamps, MicSamps_all, 'linear', 'extrap');

camflirTimeStamps_all = T_all(locs);


%% Get timestamps of microphone audio data (PTB clock)

% % MicNrSamples  = D.audio_rec.MicNrSamples;
% % MicTimeStamps = D.audio_rec.MicTimeStamps;
% % MicSamps = cumsum(MicNrSamples) - MicNrSamples(1) + 1;
% % % Now, MicTimeStamps gives you the timestamps (according to Behavior PC
% % % clock) corresponding to the audio file samples indicated in MicSamps.
% % MicSamps_all = [1 : length(y)]';
% % 
% % T_all = interp1(MicSamps, MicTimeStamps, MicSamps_all, 'linear', 'extrap');


%% Get number of video frames and check consistency with TTLs

path_vid     = path_mic;
filename_vid_root = [filename_datetime_tag '_CamFlir1_*.avi'];
ls_vid       = dir( fullfile(path_vid, filename_vid_root) );
filename_vid = ls_vid(1).name;

fprintf('Loading the video file to get the nr. of video frames (this will take a minute for long videos or when loading from the server)...\n');
v = VideoReader(fullfile(path_vid,filename_vid));
fprintf('Video file loaded, nr. of video frames: %d \n', v.NumFrames);


nr_dropped_frames = length(locs) - v.NumFrames;

if nr_dropped_frames ~= 0
    fprintf('\n ! ! ! ! ! \n');
    fprintf(' There are %d dropped frames ! \n', nr_dropped_frames);
    fprintf(' The number of video frames is %d, the number of TTL onsets extracted is %d \n', v.NumFrames, length(locs));
    fprintf(' The TTL onset extraction could have missed some TTLs, or some frames could have been dropped.\n')
    fprintf(' N.B. dropped frames do have a TTL (=camera sensor has been exposed) but they are not saved in the video file.\n')
    fprintf(' ! ! ! ! ! \n\n');
else
    fprintf('\n Great, no dropped frames! \n')
end


thr_dt = 0.001; % how many sec a frame has to be offset to detect a dropped frame
d_locs_sec = diff(locs)/Fs;
frametime_snd = mean(d_locs_sec);
problematic_ttls =  find(d_locs_sec > (frametime_snd+thr_dt) | d_locs_sec < (frametime_snd-thr_dt)) ;
nr_problematic_ttls = length( problematic_ttls );
if nr_problematic_ttls
    fprintf('\n ! ! ! ! ! \n');
    fprintf(' Based on inter-TTL interval, %d TTLs were too close or too far apart compared to expected! \n', nr_problematic_ttls);
    fprintf(' There might be something wrong with the TTL onset extraction (but not necessarily). \n')
    fprintf(' ! ! ! ! ! \n\n');
end



%% Check for dropped frames

if nr_dropped_frames ~= 0

    % nr_dropped_frames_id = (T(end,3)+1) - v.NumFrames;
    inds_dropped_frames_id = find(diff(T(:,3)) > 1);
    nr_dropped_frames_id = length( inds_dropped_frames_id );
    if nr_dropped_frames_id
        fprintf('\n ! ! ! ! ! \n');
        fprintf(' Based on the frame IDs, %d frames were dropped ! \n', nr_dropped_frames_id);
        fprintf(' ! ! ! ! ! \n\n');
    end
    dt_dropped_frames_id = (T(inds_dropped_frames_id+1,1)-T(inds_dropped_frames_id,1))/1e9;
    
    
    thr_dt = 0.005; % how many sec a frame has to be offset to detect a dropped frame
    dt_cam_sec = diff(T(:,1)/1e9);
    frametime_cam = mean(dt_cam_sec);
    inds_dropped_cam = find( (dt_cam_sec > (frametime_cam+thr_dt)) | (dt_cam_sec < (frametime_cam-thr_dt)));
    nr_dropped_frames_cam = length( inds_dropped_cam );
    if nr_dropped_frames_cam
        fprintf('\n ! ! ! ! ! \n');
        fprintf(' Based on inter-frame interval of camera timestamps, %d frames were dropped ! \n', nr_dropped_frames_cam);
        fprintf(' ! ! ! ! ! \n\n');
    end
    dt_dropped_frames_cam = (T(inds_dropped_cam+1,1)-T(inds_dropped_cam,1))/1e9;
    
    
    dt_bon_sec = diff(T(:,2));
    frametime_bon = mean(dt_bon_sec);
    inds_dropped_frames_bon =  find(dt_bon_sec > frametime_bon*1.9);
    nr_dropped_frames_bon   = length( inds_dropped_frames_bon );
    % if nr_dropped_frames_bon
    %     fprintf('\n ! ! ! ! ! \n');
    %     fprintf(' Based on inter-frame interval of bonsai timestamps, %d frames were dropped ! \n', nr_dropped_frames_bon);
    %     fprintf(' ! ! ! ! ! \n\n');
    % end



    if ~isempty(inds_dropped_frames_id) && isequal(inds_dropped_frames_id, inds_dropped_cam)
        fprintf('Dropped frames based on frame IDs and based on inter-frame interval of camera timestamps are matching. This is good... \n')
    else
        fprintf('\n ! ! ! ! ! \n');
        fprintf('Dropped frames based on frame IDs and based on inter-frame interval of camera timestamps are empty or NOT matching. \n')
        fprintf('You need to decide which dropped frames to take. By default we take Dropped frames based on frame IDs \n')
        fprintf(' ! ! ! ! ! \n\n');
    end


    if length(inds_dropped_frames_id) == nr_dropped_frames
        inds_dropped_frames = inds_dropped_frames_id;

    % If inds_dropped_frames_id has one frame less than the actual
    % nr_dropped_frames, take also the very last ttl as dropped frame:
    elseif length(inds_dropped_frames_id) == nr_dropped_frames-1
        inds_dropped_frames = [inds_dropped_frames_id; length(camflirTimeStamps_all)];

    else
        fprintf('\n ! ! ! ! ! \n');
        fprintf('Dropped frames based on frame IDs are too few compared to nr_dropped_frames \n')
        fprintf('Something is wrong (more than usual), check manually what is going on... \n')
        fprintf(' ! ! ! ! ! \n\n');

        % inds_dropped_frames = inds_dropped_frames_bon;
    end

else

    inds_dropped_frames = [];

end



%% Remove dropped frames from TTL timestamps

camflirTimeStamps = camflirTimeStamps_all;
% camflirTimeStamps = D.cam.camflir.TimeStamps; % as extracted at the end of trial acquisition

if nr_dropped_frames ~= 0

    camflirTimeStamps(inds_dropped_frames) = [];

end

assert( length(camflirTimeStamps) == v.NumFrames, 'The number of video frames is different from the number of TTL onsets extracted !')


%% Get frames corresponding to sound onsets

soundTimeStamps_eachType = D.sounds.soundTimeStamps;
n_sound_types = length(soundTimeStamps_eachType);

soundOnsets_frameNr_eachType = cell(n_sound_types, 1);

thr_dt = frametime_snd + 0.005;
for s = 1 : n_sound_types
    soundTimeStamps_thisType = soundTimeStamps_eachType{s};
    n_sound_events_thisType = length(soundTimeStamps_thisType);
    for ss = 1 : n_sound_events_thisType
        % N.B. soundTimeStamps_thisType already includes the
        % PredictedOutputLatency, so do not add it again!
        % sound_onset = soundTimeStamps_thisType(ss) + D.sounds.PredictedOutputLatency;
        sound_onset = soundTimeStamps_thisType(ss);
        % % Get closest frame to sound onset (it can be the frame right 
        % after or right before the sound onset):
        % [dt, ix_min] = min( abs(sound_onset - camflirTimeStamps) );
        % Get first frame after sound onset:
        ix_min = find( (camflirTimeStamps - sound_onset) > 0, 1,'first' );
        dt = camflirTimeStamps(ix_min) - sound_onset;
        if dt > thr_dt
            warning(sprintf('The difference between sound onset and frame time is beyond threshold, something is wrong! Sound type %d, sound event nr. %d', s,ss));
            soundOnsets_frameNr_eachType{s}(ss,1) = nan;
        else
            soundOnsets_frameNr_eachType{s}(ss,1) = ix_min;
        end
    end
end


%% Save camflirTimeStamps and soundOnsets_frameNr_eachType

D.sounds.soundOnsets_frameNr_eachType = soundOnsets_frameNr_eachType;
D.cam.camflir.TimeStamps_corr         = camflirTimeStamps;
sounds = D.sounds;
cam    = D.cam;
save(fullfile(path_data, filename_data), 'sounds', 'cam', '-append');

fprintf(' .mat file has been updated: %s \n\n', fullfile(path_data, filename_data))


%% get_frameID_from_frames_via_OCR__script

get_frameID_from_frames_via_OCR__script.m

