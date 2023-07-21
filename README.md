# RecordAudioVideo_playSounds

![image](https://github.com/lachioma/RecordAudioVideo_playSounds/assets/29898879/5d5760f3-dec1-443a-abf9-bc1da8cda41b)


* You can play up to 5 distinct sounds on demand (press the button). Just place the sound files in the subfolder /Sounds, with filenames 'sound1*.wav' etc. (* means you can have any other text).

* When the sound is played, a TTL is sent to the 6th channel of the sound card. This can be used to extract sound timestamps in sound card clock.

* The TTL signal has 1 pulse for sound 1, 2 pulses for sound 2, etc. so that they can be distinguished even by just using that signal if needed.
Sounds are also timestamped in Matlab clock slightly less accurate than sound card clock by less than 1-3 msec in most cases (so pretty negligible).

* Option to try find the corresponding FLIR video file to move it and rename it. For this to work I set it that no more than 15 seconds of difference there must be between trial start in Matlab and bonsai start.
  
* I've also added a feature to extract on the fly to e FLIR frames timestamps in both sound card clock and Matlab clock.
  
* So the system allows time-stamping of FLIR cameras, microphones, and sounds in both clocks, and orange cameras in Matlab clock.
