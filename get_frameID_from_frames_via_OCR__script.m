

path_vid     = path_mic;
filename_vid_root = [filename_datetime_tag '_CamFlir1_*.avi'];
ls_vid       = dir( fullfile(path_vid, filename_vid_root) );
filename_vid = ls_vid(1).name;

fprintf('Loading the video file to get the nr. of video frames (this will take a minute for long videos or when loading from the server)...\n');
v = VideoReader(fullfile(path_vid,filename_vid));
fprintf('Video file loaded, nr. of video frames: %d \n', v.NumFrames);


% frame = readFrame(v);
% Load last frame (because it contains the largest frame ID number, so you
% can make sure that the ROI will contain the whole frame ID):
frame = read(v,5249);

figure;image(frame);


%% Set ROI and test OCR

rois = [3 5 100 25 ]; % [x top-left corner, y top-left corner, width, height]

n = 1;
box = frame([rois(n,2):rois(n,2)+rois(n,4)],[rois(n,1):rois(n,1)+rois(n,3)],1);

figure; imshow(box);



txt = ocr(frame,rois, 'CharacterSet',"0123456789", 'TextLayout',"word")


frameIDocr = str2double(txt.Text)


%% Run OCR through all frames

tic
frameIDocr = nan(v.NumFrames,1);

v.CurrentTime = 0;
cnt = 0;

fprintf('Processing this video should take ~%3.1f min (speed is ~22 fps)\n', v.NumFrames/22/60)

disp_sameline([], 0)

while hasFrame(v)
    cnt = cnt + 1;
    disp_sameline( sprintf('Processing frame %1.0f/%1.0f\n', cnt,v.NumFrames) );

    frame = readFrame(v);
    % figure; imshow(frame);

    txt = ocr(frame,rois, 'CharacterSet',"0123456789", 'TextLayout',"word");
    
    frameIDocr(cnt) = str2double(txt.Text);
end

toc


%%

frameIDcam = T(:,3);

isequal(frameIDcam, frameIDocr)

figure;
plot(frameIDcam - frameIDocr)


%% Supporting function disp_sameline
function disp_sameline(msg, cnt)
% Display and update message printed in command window, keeping the same
% line.
% Do not print any other line between one call of disp_sameline and the
% next one.
% It is better to run disp_sameline([], 0) in the caller function before
% the actual first message printing.
% Otherwise the persistent revereseStr will remember its value from
% previous runs of the same caller function amd disp_sameline will wrongly
% delete the last line printed in the command window.
% If cnt == 0, initialize persistent variable reverseStr to ''
%
% Demo:
%     disp_sameline([], 0)
%     for i = 1 : iters
%         disp_sameline( sprintf('This is iteration %1.0f\n', i) );
%         % /!\ don't use disp('hello') or similars between one call of 
%         % disp_sameline and the next one!
%     end
%
% Alessandro La Chioma - 2017-08-01

persistent reverseStr

if nargin < 2 
    if isempty(reverseStr)
        reverseStr = '';
    end
else
    if cnt == 0
        reverseStr = '';
        return
    end
end
        
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));
end