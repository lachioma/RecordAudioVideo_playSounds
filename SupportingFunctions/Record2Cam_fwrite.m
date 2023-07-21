function Record2Cam_fwrite(camIP1, saveName1, camIP2, saveName2, saveDir, camUsername, camPWD, TimeStamps_saveName, nFrames)
% Record2Cam(camIP1, saveName1, camIP2, saveName2, saveDir, camUsername, camPWD, D, nFrames)
%
% USAGE
% Record 2 cameras:
%   Record2Cam(camIP1, saveName1, camIP2, saveName2, saveDir, camUsername, camPWD, D)
% Record only 1 camera:
%   Record2Cam(camIP1, saveName1, [], [], saveDir,  camUsername, camPWD, D)
% Use nFrames to record for a fixed number of frames.
% 
%  IMPORTANT! The camObj must be created inside this function. So, never
% change this function such that you run ipcam outside this function and
% then give camObj as input to this function, it will not work! If not,
% when executing the function via parfeval you get a strange error.
% 

% Alessandro La Chioma ..... 2020/07

% camObj = ipcam(camIP, camUsername, camPWD);
% preview(camObj)
% closePreview(camObj)

if ~isempty(camIP1) && ~isempty(camIP2) % 2 cams
    camObj1 = ipcam(camIP1, camUsername, camPWD);
    camObj2 = ipcam(camIP2, camUsername, camPWD);
    vidWriterObj1 = VideoWriter([saveDir, filesep, saveName1]);
    vidWriterObj2 = VideoWriter([saveDir, filesep, saveName2]);
    open(vidWriterObj1);
    open(vidWriterObj2);
    TimeStamps_fid = fopen(TimeStamps_saveName,'w'); % this line is in
%     the main script
    if nargin < 9 || isempty(nFrames) 
        while 1
%             pause(.025)
            [im, ~]  = snapshot(camObj1);
            [im2, ~] = snapshot(camObj2);
            ts = GetSecs(); % GetSecs is part of Psychtoolbox
            writeVideo(vidWriterObj1, im);
            writeVideo(vidWriterObj2, im2);
%             send(D,ts);
            fwrite(TimeStamps_fid, ts, 'double');
        end
    else
        for f = 1 : nFrames
%             pause(.05)
            [im, ~] = snapshot(camObj1);
            [im2, ~] = snapshot(camObj2);
            ts = GetSecs();  % GetSecs is part of Psychtoolbox
            writeVideo(vidWriterObj1, im);
            writeVideo(vidWriterObj2, im2);
%             send(D,ts);
            fwrite(TimeStamps_fid, ts, 'double');
        end
    end
    
else  % 1 cam
    camObj1 = ipcam(camIP1, camUsername, camPWD);
    vidWriterObj1 = VideoWriter([saveDir, filesep, saveName1]);
    open(vidWriterObj1);
    TimeStamps_fid = fopen(TimeStamps_saveName,'w'); % this line is in
%     the main script
    if nargin < 9 || isempty(nFrames) 
        while 1
%             pause(.05)
            [im, ~]  = snapshot(camObj1);
            ts = GetSecs(); % GetSecs is part of Psychtoolbox
            writeVideo(vidWriterObj1, im);
%             send(D,ts);
            fwrite(TimeStamps_fid, ts, 'double');
        end
    else
        for f = 1 : nFrames
%             pause(.05)
            [im, ~] = snapshot(camObj1);
            ts = GetSecs(); % GetSecs is part of Psychtoolbox
            writeVideo(vidWriterObj1, im);
%             send(D,ts);
            fwrite(TimeStamps_fid, ts, 'double');
        end
    end
end