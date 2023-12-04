# RecordAudioVideo_playSounds

![image](https://github.com/lachioma/RecordAudioVideo_playSounds/assets/29898879/f8d87686-0dfb-4f97-a40f-b3818db4a8ca)


## Usage
### To start recording
* Start the recording pressing the Record button.
* Then start the Bonsai program to record the FLIR camera video. Start this within 20 seconds from the Matlab RecordAudioVideo_playSounds if you want the video file to be automatically moved to the Matlab data folder from the Bonsai data folder.
### To stop
* First stop Bonsai by pressing F8.
* Then stop the Matlab RecordAudioVideo_playSounds.

## Other info

* You can play up to 5 distinct sounds on demand (press the button). Just place the sound files in the subfolder /Sounds, with filenames 'sound1*.wav' etc. (* means you can have any other text).

* When the sound is played, a TTL is also recorded by sound card. This can be used to extract sound timestamps in sound card clock.

* The TTL signal has 1 pulse for sound 1, 2 pulses for sound 2, etc. so that they can be distinguished even by just using that signal if needed.
Sounds are also timestamped in Matlab clock slightly less accurate than sound card clock by less than 1-3 msec in most cases (so pretty negligible).

* Option to try find the corresponding FLIR video file to move it and rename it. For this to work I set that the Bonsai recording must start <20 seconds after the trial start in Matlab.
  
* I've also added a feature to extract on the fly the FLIR frames timestamps in both sound card clock and Matlab clock.
  
* So the system allows time-stamping of FLIR cameras, microphones, and sounds in both clocks, and orange cameras in Matlab clock.

* I've added a "LogDiary" (v1.2): for each trial, all messages displayed in the command window are saved as a txt file.

* You can play a Looming Disc on demand (press the button) (v1.2). At button press, the disc appears at its onset size, it expands linearly until final size, it stays at final size for a post-stim interval, then it disappears. Three timestamps (PTB clock) are taken at onset, final size, disappearance.


## Data

The Matlab file `_data.mat` contains the structures `settings`, `audio_rec`, `cam`, `sounds`.

### settings

Some basic parameters of the recorded trial.

### audio_rec

Settings and timestamps of sound card in input.

* `MicTimeStamp`s: PTB timestamps taken when every chunk of audio samples is acquired, with nr samples of each chunk indicated in MicNrSamples. The corresponding audio data is the file `_audiorec.wav`. 

* `MicNrSamples`: nr. audio samples of each chunk that is acquired. The corresponding audio data is the file _audiorec.wav.


#### audio_rec.ttl

Settings and timestamps of the 5th input channel of the sound card, carrying the exposure of the FLIR camera.
The sampling of this input channel is downsampled by a factor of `ttl.downsampling` with respect to audio_rec to obtain a sampling rate of `ttl.SampleRate`.

* `MicTimeStamps`: PTB timestamps taken when every chunk of audio samples is acquired, with nr samples of each chunk indicated in MicNrSamples. The corresponding audio data is the file `_TTLcam.wav`. 
* `MicNrSamples`: nr. audio samples of each chunk that is acquired. The corresponding audio data is the file `_audiorec.wav`.


### audio_rec.ttl_sound

Timestamps and audio samples of sound output recorded back as input. N.B. This is not done wiring the sound output as audio input, but using the [slave mode 64 of PsychPortAudio](http://psychtoolbox.org/docs/PsychPortAudio-OpenSlave#:~:text=The%20slave%2Donly%20mode%20flag%2064).

The corresponding audio data is not saved entirely, but only the events (TTL value changes,`dValues_TTL`) and their timings (`MicTimeStamps_TTL`, `MicNrSamples_TTL`) are stored in the binary file `_TTLsound.ttlsnd`, which is then deleted once the .mat file has been saved.

* `MicTimeStamps`: PTB timestamps taken when every chunk of audio samples is acquired, with nr samples of each chunk indicated in MicNrSamples. The corresponding audio data is actually deleted, so this variable is in fact useless.
* `MicNrSamples`: nr. audio samples of each chunk that is acquired. The corresponding audio data is actually deleted, so this variable is in fact useless.
* `MicTimeStamps_TTL`: PTB timestamps of every TTL change (from 0 to 1 or from 1 to 0). These timestamps come from audio samples converted into PTB timestamps.
* `MicNrSamples_TTL`: audio samples of every TTL change (from 0 to 1 or from 1 to 0).
* `dValues_TTL`: value changes in the TTL signal (+1 when the TTL goes from 0 to 1, -1 when the TTL goes from 1 to 0).


### cam

Settings and timestamps of IP ("orange") cameras and FLIR camera.

* `camTimeStamps`: PTB timestamps of each frame of every IP camera. Keep in mind that these timestamps are not very accurate and should have a delay of ~0.1 seconds.

#### cam.camflir

* `Timestamp_camera_clock_us`: timestamps of each frame, according to camera internal clock (in microseconds). This variable is just taken from the .csv file produced by Bonsai.
* `Timestamp_bonsai_clock_s`: timestamps of when each frame is acquired by Bonsai, according to Bonsai clock (in seconds). This variable is just taken from the .csv file produced by Bonsai.
* `FrameID`: frame ID as given by the camera, with numbering starting from 0 (the first frame ID is subtracted). This variable is just taken from the .csv file produced by Bonsai.
* `TimeStamps`: PTB timestamps of each frame, as extracted post-hoc from the TTL signal in the 5th audio channel (these timestamps come from audio samples converted into PTB timestamps). 

The corresponding audio data is the file `_TTLcam.wav`. The sampling of this input channel is downsampled by a factor of `audio_rec.ttl.downsampling` with respect to audio_rec to obtain a sampling rate of `ttl.SampleRate`. 

This extraction is done after the trial recording is stopped and it might be wrong. If the number of timestamps is different from the number of frames in the video file, something likely went wrong.

* `MicSamples`: audio samples of each frame, as extracted post-hoc from the TTL signal in the 5th audio channel. 

The corresponding audio data is the file `_TTLcam.wav`. The sampling of this input channel is downsampled by a factor of `audio_rec.ttl.downsampling` with respect to audio_rec to obtain a sampling rate of `ttl.SampleRate`.

This extraction is done after the trial recording is stopped and it might be wrong. If the number of timestamps is different from the number of frames in the video file, something likely went wrong.


### sounds

Settings and timestamps of sounds that are played back
* `soundTimeStamps`: PTB timestamps of every sound played, taken by the PTB sound playback command PsychPortAudio(`Start`).

### About the binary files

Binary files are used to save the data on the fly during recording. After each trial, these files are converted or used, and if all went well they are then deleted.

* _audiorec.mic: 4 microphones audio file; _audiorec.mic converted into _audiorec.flac, then the binary file is deleted.
  
* _TTLcam.ttlcam: TTL from FLIR camera; _TTLcam.ttlcam converted into _TTLcam.flac, then the binary file is deleted. From this file, cam FLIR timestamps are extracted.

* _TTLsound.ttlsnd: TTL from sounds; this file is not converted, but values are extracted and then the ttlsnd is deleted.
