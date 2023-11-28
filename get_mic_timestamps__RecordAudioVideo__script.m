
%% Select and load audio file (4 microphones)

folder_root = 'Y:\Users\ariadna\behavior_PP\25679\25679_3_sounds\2023-11-24';

fmt = '.wav';
fmt = '.flac';
fmt = '.mic'; % binary file

ls_mic = dir([folder_root filesep '*_audiorec' fmt]);
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
    % y = fread(fid,[nAudioChannels_input_toSave Inf],'single=>single')';
    % n_Samps = size(y,1);
    % Just get the nr. data points from the filesize:
    fseek(fid, 0, 'eof');
    filesize = ftell(fid);
    n_Samps = filesize/ 4 / nAudioChannels_input_toSave; % single precision takes 4 bytes per value
    % Close the file
    fclose(fid);
    % Fs = 192000;
else
    [y,Fs] = audioread(fullfile(path_mic,filename_mic));
    n_Samps = length(y);
end



%% Get timestamps of microphone audio data (PTB clock)

MicNrSamples  = D.audio_rec.MicNrSamples;
MicTimeStamps = D.audio_rec.MicTimeStamps;
MicSamps = cumsum(MicNrSamples) - MicNrSamples(1) + 1;
% Now, MicTimeStamps gives you the timestamps (according to Behavior PC
% clock) corresponding to the audio file samples indicated in MicSamps.
MicSamps_all = [1 : n_Samps]';

T_all = interp1(MicSamps, MicTimeStamps, MicSamps_all, 'linear', 'extrap');
