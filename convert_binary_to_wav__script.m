
clear

folder_root = 'Y:\Users\ariadna\behavior_PP\25679\25679_3_sounds\2023-11-24';

%% Load _data.mat file
ls_data = dir([folder_root filesep '*_data.mat']);

if length(ls_data) == 1
    path_data = ls_data(1).folder;
    filename_data = ls_data(1).name;
elseif length(ls_data) > 1
    [filename_data, path_data] = uigetfile([folder_root filesep '*_data.mat']);
else
    disp('No _audiorec.mic file found')
end

filename_datetime_tag = filename_data(1:end-9);

D = load( fullfile(path_data, filename_data));


%% Load and convert audiorec.mic into .flac (4 microphones)


if exist('path_data','var') && ~isempty(path_data)
    path_mic = path_data;
    filename_mic = [filename_datetime_tag '_audiorec.mic'];
else
    ls_mic = dir([folder_root filesep '*_audiorec.mic']);
    if length(ls_mic) == 1
        path_mic = ls_mic(1).folder;
        filename_mic = ls_mic(1).name;
    elseif length(ls_mic) > 1
        [filename_mic, path_mic] = uigetfile([folder_root filesep '*_audiorec.mic']);
    else
        disp('No _audiorec.mic file found')
    end
end

audio_rec = [];
audio_rec.fullpath_bin  = fullfile(path_mic, filename_mic);
audio_rec.fullpath_wav  = fullfile(path_mic, [D.audio_rec.filename '.wav']);
audio_rec.fullpath_flac = fullfile("C:\Users\dailyuser\Desktop\New folder", [D.audio_rec.filename '.flac']);

if ~isfile(audio_rec.fullpath_wav) && ~isfile(audio_rec.fullpath_flac)

    disp('Loading .mic file...')
    tic_wav_conversion = tic;
    fid = fopen(audio_rec.fullpath_bin, 'r');
    % Read all data from the audio file in a vector
    recaudiodata = fread(fid,[D.audio_rec.nAudioChannels_input_toSave Inf],'single=>single')';
    % Close the file
    fclose(fid);
    % Save .wav file:
    disp('...converting into .flac file...')
    audiowrite(audio_rec.fullpath_flac, recaudiodata, D.audio_rec.SampleRate, 'BitsPerSample',24);
    disp(['  .mic file saved as .wav! Time elapsed: ' num2str(toc(tic_wav_conversion)) ' s']) 

end

%% Load and convert TTLcam.ttlcam into .flac (FLIR camera TTL)


if exist('path_data','var') && ~isempty(path_data)
    path_cam = path_data;
    filename_cam = [filename_datetime_tag '_TTLcam.ttlcam'];
    if ~isfile(fullfile(path_cam, filename_cam))
        disp('No _TTLcam.ttlcam file found')
    end
else
    ls_mic = dir([folder_root filesep '*_TTLcam.ttlcam']);
    if length(ls_mic) == 1
        path_cam = ls_mic(1).folder;
        filename_cam = ls_mic(1).name;
    elseif length(ls_mic) > 1
        [filename_cam, path_cam] = uigetfile([folder_root filesep '*_TTLcam.ttlcam']);
    else
        disp('No _TTLcam.ttlcam file found')
    end
end

ttl = [];
ttl.fullpath_bin  = fullfile(path_cam, filename_cam);
ttl.fullpath_wav  = fullfile(path_cam, [D.audio_rec.filename '.wav']);
ttl.fullpath_flac = fullfile(path_cam, [D.audio_rec.ttl.filename '.flac']);

if ~isfile(ttl.fullpath_wav) && ~isfile(ttl.fullpath_flac)

    disp('Loading .ttlcam file...')
    tic_wav_conversion = tic;
    fid = fopen(ttl.fullpath_bin, 'r');
    % Read all data from the audio file in a vector
    recaudiodata = fread(fid,[D.audio_rec.ttl.nAudioChannels_input_toSave Inf],'single=>single')';
    % Close the file
    fclose(fid);
    % Save .wav file:
    disp('...converting into .flac file...')
    audiowrite(ttl.fullpath_flac, recaudiodata, D.audio_rec.ttl.SampleRate, 'BitsPerSample',24);
    disp(['  .mic file saved as .flac! Time elapsed: ' num2str(toc(tic_wav_conversion)) ' s']) 

end


