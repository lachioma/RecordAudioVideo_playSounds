function LoomingDisc_demo()
%
tic_start = tic;

%% Set parameters

% Select screen for output window:
screenid = max(Screen('Screens'));
% screenid = 1;

screenSize_cm_horz = 50; % cm
screenSize_cm_vert = 35; % cm
mouse_dist_cm      = 37; % cm

discSize_init_deg    =  2.6; % deg
discSize_end_deg     = 47.2; % deg; if [], it will set as screen width
expansionSpeed_degps = 11.2; % deg/s, linear expansion

% Set disc color:
% discColor = [1 1 1]; % white
discColor = [0 0 0]; % black
% Set background color
BackgroundLuminance = 0.5;
% Once disc reached final size, wait post_stim_interval_sec before
% displaying only the background screen:
post_stim_interval_sec = 1; %sec

%%

% Setup defaults and unit color range:
PsychDefaultSetup(2);
% Disable synctests for this quick demo:
oldSyncLevel = Screen('Preference', 'SkipSyncTests', 2); %(use 1 to continue with a warning or 2 to remove the warning completely)

% white = WhiteIndex(screenid);
% black = BlackIndex(screenid);
% gray  = round((white+black)/2);

% Open a fullscreen, onscreen window with gray background. Enable 32bpc
% floating point framebuffer via imaging pipeline on it, if this is possible
% on your hardware while alpha-blending is enabled. Otherwise use a 16bpc
% precision framebuffer together with alpha-blending. We need alpha-blending
% here to implement the nice superposition of overlapping of discs. The demo will
% abort if your graphics hardware is not capable of any of this.
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

[win, screenRect] = PsychImaging('OpenWindow', screenid, BackgroundLuminance);

% Retrieve size of window in pixels, need it later to make sure that our
% moving discs don't move out of the visible screen area:
[width_px, height_px] = RectSize(screenRect);

% Query frame duration: We use it later on to time 'Flips' properly for an
% animation with constant framerate:
ifi = Screen('GetFlipInterval', win);

% Enable alpha-blending
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

toc(tic_start)

pixel_cm      = screenSize_cm_horz/screenRect(3); % calculates the size of a pixel in cm
pixel_cm_vert = screenSize_cm_vert/screenRect(4); % calculates the size of a pixel in cm
degperpix = (2*atan(pixel_cm./(2*mouse_dist_cm(1)))).*(180/pi);  %%true only for small angles, more generally would be: atan(pix./animal_dist_cm)  %*(180/pi) is necessary because atan in matlab gives result in radians
pixperdeg = 1./degperpix;

width_deg  = width_px * degperpix;
height_deg = height_px * degperpix;

% Radius and diameter of the disc (they must be integer)
radius_init_px   = round( discSize_init_deg/2 * pixperdeg );
discSize_init_px = radius_init_px * 2;

expansionIncrement_degpf = expansionSpeed_degps * ifi ;
expansionIncrement_pxpf  = expansionIncrement_degpf * pixperdeg;

if isempty(discSize_end_deg)
    discSize_end_deg = max([width_deg width_deg]);
end
discSize_end_px  = discSize_end_deg * pixperdeg;

% smoothing sigma in pixel
sigma = 3;
% use alpha channel for smoothing edge of disc?
useAlpha = true;
% smoothing method: cosine (0), smoothstep (1), inverse smoothstep (2)
smoothMethod = 1;

% Build a procedural disc
discTexture = CreateProceduralSmoothedDisc(win, discSize_init_px, discSize_init_px, [0 0 0 0], radius_init_px, sigma, ...
                                           useAlpha, smoothMethod);

% Preallocate array with destination rectangles:
texRect = Screen('Rect', discTexture);
dstRects(:,1) = CenterRectOnPointd(texRect, width_px/2, height_px/2);

rotAngles = 0;
myAlpha = 1; % 1 no transparency, 0 fully transparent
discColor = [discColor, myAlpha];

% Initially sync us to VBL at start of animation loop.
count = 0;
Screen('FillRect', win, BackgroundLuminance);
vbl = Screen('Flip', win);

discSize_next_px  = discSize_init_px;
discSize_next_deg = discSize_init_deg;

toc(tic_start)

% Animation loop: Run until any keypress:
while (discSize_next_deg < discSize_end_deg)
    % Step one: Batch-Draw all discs at the positions (dstRects) and
    % orientations (rotAngles) and colors (colours)
    % and with the stimulus parameters 'discParameters'
    Screen('DrawTextures', win, discTexture, [], dstRects, rotAngles, [], [], discColor);

    % Mark drawing ops as finished, so the GPU can do its drawing job while
    % we can compute updated parameters for next animation frame. This
    % command is not strictly needed, but it may give a slight additional
    % speedup, because the CPU can compute new stimulus parameters in
    % Matlab, while the GPU is drawing the stimuli for this frame.
    % Sometimes it helps more, sometimes less, or not at all, depending on
    % your system and code, but it only seldomly hurts.
    % performance...
    Screen('DrawingFinished', win);

    % Done. Flip one video refresh after the last 'Flip', ie. try to
    % update the display every video refresh cycle if you can.
    % This is the same as Screen('Flip', win);
    % but the provided explicit 'when' deadline allows PTB's internal
    % frame-skip detector to work more accurately and give a more
    % meaningful report of missed deadlines at the end of the script. Not
    % important for this demo, but here just in case you didn't know ;-)
    vbl = Screen('Flip', win, vbl + 0.5 * ifi);

    count = count + 1;

    discSize_prev_px  = discSize_next_px;
    discSize_prev_deg = discSize_next_deg;
    discSize_next_px  = discSize_prev_px  + expansionIncrement_pxpf;
    discSize_next_deg = discSize_prev_deg + expansionIncrement_degpf;

    scale = discSize_next_px / discSize_prev_px;
%     texRect = texRect * scale;
    texRect([3,4]) = texRect([3,4]) + expansionIncrement_pxpf;

    dstRects(:, 1) = CenterRectOnPointd(texRect, width_px/2, height_px/2)';

    if KbCheck
        break
    end
end

tic_post = tic;
while (toc(tic_post) < post_stim_interval_sec)
    if KbCheck
        break
    end
end

Screen('FillRect', win, BackgroundLuminance);
Screen('Flip', win);

KbWait; % just wait and exit upon key press
while ~KbCheck
    % just wait
end
% Close onscreen window, release all ressources:
sca;

% Restore old settings for sync-tests:
Screen('Preference', 'SkipSyncTests', oldSyncLevel);
