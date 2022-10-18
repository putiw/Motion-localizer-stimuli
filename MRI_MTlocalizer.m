clc
clear
tbUse vistadisp
Screen('Preference','TextRenderer', 0);
%--------------%
% MT localizer %
%--------------%
doEyelink = 0;
a = cd;
if a(1)=='/'
    a = PsychHID('Devices');
    for i = 1:length(a), d(i) = strcmp(a(i).usageName, 'Keyboard'); end
    keybs = find(d);
else
    keybs = [];
end

mykeys = {'1';'2';'3';'4';'6';'7';'8';'9'};
myKbCheck = [];
mypress = [];


% Syncobx?
syncbox = 0;
commandwindow
PsychDebugWindowConfiguration(1,1)
% Close (eventually) open connections and PTB screens
IOPort('CloseAll');
Screen('CloseAll');

% Trick suggested by the PTB authors to avoid synchronization/calibration
% problems
figure(1)
plot(sin(0:0.1:3.14));
% Close figure with sin plot (PTB authors trick for synchronization)
close Figure 1
timetopress = 2;
% Synchronization tests procedure - PTB
% Do you want to skipsync tests (1) or not (0) ?


% KbName will switch its internal naming
% scheme from the operating system specific scheme (which was used in
% the old Psychtoolboxes on MacOS-9 and on Windows) to the MacOS-X
% naming scheme, thereby allowing to use one common naming scheme for
% all operating systems
KbName('UnifyKeyNames');
% Code to identify "escape" key
escapekeycode = KbName('ESCAPE');

% open syncbox connection
if syncbox
    syncbox_handle = IOPort('OpenSerialPort', 'COM2', 'BaudRate=57600 DataBits=8 Parity=None StopBits=1 FlowControl=None');
    IOPort('Flush',syncbox_handle);
end


%----------------------------------------------------%
%%
fprintf('\n')
initials = input('Please enter subjct initials: ', 's');
sesNum = input('Please enter session number: ', 's');
sesNum = str2double(sesNum);

sesFileName = sprintf('%s%d', initials, sesNum);

while exist(sprintf('%s.edf',sesFileName), 'file')
    
    fprintf('\nFilename %s exists. Please re-enter subj ID and session number.\n', sesFileName)
    initials = input('Please enter subjct initials: ', 's');
    sesNum = input('Please enter session number: ', 's');
    sesNum = str2double(sesNum);
    sesFileName = sprintf('%s%d%s', initials, sesNum);
    
end


%%

AssertOpenGL;
%%
d = dir(fullfile('./', '*.png'));

load params

assert(numel(d)>0)
[x,y] = meshgrid(linspace(-1, 1, 30));
maskblob = 255*(x.^2+y.^2<1);

for ii = 1:length(d)
    im = imread(fullfile(d(ii).name));
    im = im(7:121, 2:116,:);
    im = imresize(im, [1 1]* 30);
    im(:,:,end+1) = maskblob;
    fixStimulus.images{ii} = im;
end


%%

try
    
    % ------------------------
    % set dot field parameters
    % ------------------------
    
    nframes     = 18000; % number of animation frames in loop
%     nframes     = 181; % number of animation frames in loop

    mon_width   = 36.2;   % horizontal dimension of viewable screen (cm)
    v_dist      = 83.5;   % viewing distance (cm)
    dot_speed   = 5; %7;    % dot speed (deg/sec)
    f_kill      = 0.00; % fraction of dots to kill each frame (limited lifetime)
    ndots       = 100; %2000; % number of dots
    max_d       = 5;%15;   % maximum radius of  annulus (degrees)
    min_d       = 0.3; %1;    % minumum (degrees)
    dot_w       = 0.1;  % width of dot (deg)
    fix_r       = 0.09; %0.15; % radius of fixation point (deg)
    waitframes  = 1; % Show new dot-images at each waitframes'th monitor refresh.  'waitframes' Number of video refresh intervals to show each image before updating the dot field. Defaults to 1 if omitted.
    
        
    fixation_train = {ones(1,460);ones(1,600);ones(1,560)}
    stim.fixation_sequence = [];
    
    while length(stim.fixation_sequence) < nframes
        
        
        p = randperm(3);
        multi = randperm(2);
        stim.fixation_sequence = [stim.fixation_sequence fixation_train{p(1)}.*multi(1)'];
    end
    
    stim.fixation_sequence = stim.fixation_sequence(1:nframes);
    stim.smiley_change_sec = find((abs(diff(stim.fixation_sequence)))==1)/60;
    
    %     smiley_change_frame = [fixation_change fixation_change+480];
    %     smiley_change_frame = smiley_change_frame/60;
    % ---------------
    % My Parameters
    % ---------------
    
    stim.offset_left=10;  % (degrees)
    stim.offset_right=10; % (degrees)
    stim.offset_center=0;  % (degrees)
    
    efr=60; % estimated-target frame rate (Hz)
    TR=1;   % MRI TR (seconds)
    
    stim.my_protocol=zeros(1,nframes);
    
    % rest
    repeat = 4;
    conditions = repmat([1 -1 2 -2 3 -3],[1 4]);
    
    framesperblock = 720;
    for c = 1 : length(conditions)
        
        stim.my_protocol((1:framesperblock)+(c-1)*framesperblock) = conditions(c);
        
    end
    
    
    
    % ---------------
    % open the screen
    % ---------------
    %%
    screens=Screen('Screens');
    screenNumber=max(screens);
    [w, rect] = Screen('OpenWindow', screenNumber, 0);
    clear params
    params.display.windowPtr = w;
    fixationStimulus = makeTextures(params.display, fixStimulus);
    
%%
if doEyelink
    fprintf('\n[%s]: Setting up Eyelink..\n',mfilename)
    
    Eyelink('SetAddress','192.168.1.5');
    el = EyelinkInitDefaults(params.display.windowPtr);
    EyelinkUpdateDefaults(el);
    %
    % %     Initialize the eyetracker
    Eyelink('Initialize', 'PsychEyelinkDispatchCallback');
    % %     Set up 5 point calibration
    s = Eyelink('command', 'calibration_type=HV5');
    %
    % %     Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    %
    % %     Throw an error if calibration failed
    if s~=0
        error('link_sample_data error, status: ', s)
    end
    
    el = prepEyelink(params.display.windowPtr);
    ELfileName = sprintf('%s.edf', sesFileName);
    edfFileStatus = Eyelink('OpenFile', ELfileName);
    
    if edfFileStatus ~= 0, fprintf('Cannot open .edf file. Exiting ...');
        try
            Eyelink('CloseFile');
            Eyelink('Shutdown');
        end
        return;
    else
        fprintf('\n[%s]: Succesfully openend Eyelink file..\n',mfilename)
    end
    
    cal = EyelinkDoTrackerSetup(el);
    
end



    % Enable alpha blending with proper blend-function. We need it
    % for drawing of smoothed points:
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [center(1), center(2)] = RectCenter(rect);
%     fps=Screen('FrameRate',w);      % frames per second
%     ifi=Screen('GetFlipInterval', w);
    fps = 60;
%     if fps==0
%         fps=1/ifi;
%     end;
    ifi = 1/fps;
    white = WhiteIndex(w);
    HideCursor;	% Hide the mouse cursor
    Priority(MaxPriority(w));
    
    % Do initial flip...
    vbl=Screen('Flip', w);
    
    % ---------------------------------------
    % initialize dot positions and velocities
    % ---------------------------------------
    
%     ppd = pi * (rect(3)-rect(1)) / atan(mon_width/v_dist/2) / 360;    % pixels per degree
    
    ss.ScreenWidth      = 60;
    ss.ScreenHeight     = 36.2;
    ss.ViewDistance     = 85;
    ss.ScreenSizePixels = [0 0 1920 1080];
     
    ss.VisAngle = (2*atand(ss.ScreenHeight/(2*ss.ViewDistance)));
    ppd=round(ss.ScreenSizePixels(4)/ss.VisAngle);
    pfs = dot_speed * ppd / fps;                            % dot speed (pixels/frame)
    s = dot_w * ppd; % dot size (pixels)
    smin=s;
    smax=s+s;
    fix_cord = [center-fix_r*ppd center+fix_r*ppd];
    
    rmax = max_d * ppd;	% maximum radius of annulus (pixels from center)
    rmin = min_d * ppd; % minimum
    
    ms=(smax-smin)/(rmax-rmin);
    
    % IN #####################
    r_IN = rmax * sqrt(rand(ndots,1));	% r
    r_IN(r_IN<rmin) = rmin;
    t_IN = 2*pi*rand(ndots,1);                     % theta polar coordinate
    cs_IN = [cos(t_IN), sin(t_IN)];
    xy_IN = [r_IN r_IN] .* cs_IN;   % dot positions in Cartesian coordinates (pixels from center)
    
    %mdir = 2 * floor(rand(ndots,1)+0.5) - 1;    % motion direction (in or out) for each dot
    mdirIN=ones(ndots,1)-2;
    drIN = pfs * mdirIN;                            % change in radius per frame (pixels)
    dxdyIN = [drIN drIN] .* cs_IN;                       % change in x and y per frame (pixels)
    
    
    
    % OUT #####################
    r_EXT = rmax * sqrt(rand(ndots,1));	% r
    r_EXT(r_EXT<rmin) = rmin;
    t_EXT = 2*pi*rand(ndots,1);                     % theta polar coordinate
    cs_EXT = [cos(t_EXT), sin(t_EXT)];
    xy_EXT = [r_EXT r_EXT] .* cs_EXT;   % dot positions in Cartesian coordinates (pixels from center)
    
    %mdir = 2 * floor(rand(ndots,1)+0.5) - 1;    % motion direction (in or out) for each dot
    mdirEXT=ones(ndots,1);
    drEXT = pfs * mdirEXT;                            % change in radius per frame (pixels)
    dxdyEXT = [drEXT drEXT] .* cs_EXT;                       % change in x and y per frame (pixels)
    
    
    % Create a vector with different colors for each single dot, if
    % requested:
    colvect=white;
    
    %%
    Screen('FillRect',w,[0 0 0]);

    Screen('DrawTexture', w, fixationStimulus.textures(1))
    vbl=Screen('Flip', w);
    
    
    
        if doEyelink
        Eyelink('StartRecording');
        WaitSecs(2)
        end
    
    
    wait4T(keybs);    
    t = 0;
    t0 = GetSecs;
   
    %%
    
    % --------------
    % animation loop
    % --------------
    
    
    

    while t <= nframes/fps;
        
        
        t = GetSecs - t0;
        i = ceil(t/(1/fps));
%         
%         if mod(i,720) == 0;
%             myt(ct) = t;
%             t
%             ct = ct + 1;
%         end
                
        if (i>1) && i<=nframes
            %             Screen('FillOval', w, uint8(white), fix_cord);	% draw fixation dot (flip erases it)
        
            
            
            Screen('FillRect',w,[0 0 0]);
            Screen('DrawTexture', w, fixationStimulus.textures(1))
            
            if stim.fixation_sequence(i) == 1
                
                %                 Screen('FillOval', w, uint8([255 0 0]), fix_cord);
                Screen('DrawTexture', w, fixationStimulus.textures(1))
                
            else
                Screen('DrawTexture', w, fixationStimulus.textures(2))
                
                %                 Screen('DrawTexture', w, params.display.fixationStimulus.textures(2), ...
                %             params.display.fixationStimulus.srcRect, params.display.fixationStimulus.destRect);
                %                 Screen('FillOval', w, [0 255 0], fix_cord);
                
                % draw fixation dot (flip erases it)
                
            end
            
            %%
            mydiff = t - stim.smiley_change_sec;
            mydiff(mydiff>timetopress) = NaN;
            mydiff(mydiff<0) = NaN;
            lowerthanone = find(mydiff>0);
            
            
             if ~isempty(lowerthanone)
                [pressed,~,keyCode] = KbCheck;
                
                if pressed
                    
                    keypressed = KbName(keyCode);
                    Index = find(contains(mykeys,keypressed(1)));
                    
                    if ~isempty(Index)
                        mypress = [mypress t];
                        myKbCheck = [myKbCheck KbCheck];
                    end
                    
                else
                    
                    mypress = [mypress t];
                    myKbCheck = [myKbCheck -1];
                    
                end
            else
                
                [pressed,~,keyCode] = KbCheck;
                
                if pressed
                    
                    keypressed = KbName(keyCode);
                    Index = find(contains(mykeys,keypressed));
                    
                    if ~isempty(Index)
                        
                        mypress = [mypress t];
                        myKbCheck = [myKbCheck KbCheck];
                        
                    end
                else
                    
                    mypress = [mypress t];
                    myKbCheck = [myKbCheck KbCheck];
                end
             end
            
            %%
            
            if stim.my_protocol(i)==0
                % do nothing
                Screen('DrawTexture', w, fixationStimulus.textures(1))

            elseif stim.my_protocol(i)==1
                my_xymatrix(1,:)=xymatrix(1,:)+stim.offset_center*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN + dxdyIN;						% move dots
                r_IN = r_IN + drIN;							% update polar coordinates too
                xy_EXT = xy_EXT + dxdyEXT;				% move dots
                r_EXT = r_EXT + drEXT;					% update polar coordinates too
                
                
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip
                
            elseif stim.my_protocol(i)==-1
                my_xymatrix(1,:)=xymatrix(1,:)+stim.offset_center*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN;  						% move dots
                r_IN = r_IN;							% update polar coordinates too
                xy_EXT = xy_EXT;      				% move dots
                r_EXT = r_EXT;   		    			% update polar coordinates too
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            elseif stim.my_protocol(i)==2
                my_xymatrix(1,:)=xymatrix(1,:)-stim.offset_left*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN + dxdyIN;						% move dots
                r_IN = r_IN + drIN;							% update polar coordinates too
                xy_EXT = xy_EXT + dxdyEXT;				% move dots
                r_EXT = r_EXT + drEXT;					% update polar coordinates too
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            elseif stim.my_protocol(i)==-2
                my_xymatrix(1,:)=xymatrix(1,:)-stim.offset_left*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN;						% move dots
                r_IN = r_IN;						% update polar coordinates too
                xy_EXT = xy_EXT;	     			% move dots
                r_EXT = r_EXT;					% update polar coordinates too
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            elseif stim.my_protocol(i)==3
                my_xymatrix(1,:)=xymatrix(1,:)+stim.offset_right*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN + dxdyIN;						% move dots
                r_IN = r_IN + drIN;							% update polar coordinates too
                xy_EXT = xy_EXT + dxdyEXT;				% move dots
                r_EXT = r_EXT + drEXT;					% update polar coordinates too
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            elseif stim.my_protocol(i)==-3
                my_xymatrix(1,:)=xymatrix(1,:)+stim.offset_right*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN;						% move dots
                r_IN = r_IN;							% update polar coordinates too
                xy_EXT = xy_EXT;				% move dots
                r_EXT = r_EXT;					% update polar coordinates too
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            end
            
        end
        
        
        
        [keyIsDown, secs, keyCode] = KbCheck;
        
        % The user asked to exit the program
        if keyIsDown==1 && keyCode(escapekeycode)
            
            % Close PTB screen and connections
            Screen('CloseAll');
            IOPort('CloseAll');
            ShowCursor;
            Priority(0);
            
            % Launch window with warning of early end of program
            warndlg('The run was terminated with ''Esc'' before the end!','Warning','modal')
            
            return % abort program
        end
        
        
        % check to see which dots have gone beyond the borders of the annuli
        r_out = find(r_EXT > rmax | r_EXT < rmin | rand(ndots,1) < f_kill);	% dots to reposition
        nout = length(r_out);
        if nout
            
            % choose new coordinates
            r_EXT(r_out) = rmin; %* sqrt(rand(nout,1));
            r_EXT(r_EXT<rmin) = rmin;
            t_EXT(r_out) = 2*pi*(rand(nout,1));
            
            % now convert the polar coordinates to Cartesian
            cs_EXT(r_out,:) = [cos(t_EXT(r_out)), sin(t_EXT(r_out))];
            xy_EXT(r_out,:) = [r_EXT(r_out) r_EXT(r_out)] .* cs_EXT(r_out,:);
            
            % compute the new cartesian velocities
            dxdyEXT(r_out,:) = [drEXT(r_out) drEXT(r_out)] .* cs_EXT(r_out,:);
            
        end;
        xymatrix_EXT = transpose(xy_EXT);
        
        
        
        
        
        % check to see which dots have gone beyond the borders of the annuli
        r_out = find(r_IN > rmax | r_IN < rmin | rand(ndots,1) < f_kill);	% dots to reposition
        nout = length(r_out);
        if nout
            
            % choose new coordinates
            r_IN(r_out) = rmax; %* sqrt(rand(nout,1));
            r_IN(r_IN<rmin) = rmax;
            t_IN(r_out) = 2*pi*(rand(nout,1));
            
            % now convert the polar coordinates to Cartesian
            cs_IN(r_out,:) = [cos(t_IN(r_out)), sin(t_IN(r_out))];
            xy_IN(r_out,:) = [r_IN(r_out) r_IN(r_out)] .* cs_IN(r_out,:);
            
            % compute the new cartesian velocities
            dxdyIN(r_out,:) = [drIN(r_out) drIN(r_out)] .* cs_IN(r_out,:);
        end
        
        xymatrix_IN = transpose(xy_IN);
        
        
        
        my_s(1:ndots)=ms*r_IN'+smin;
        my_s(ndots+1:2*ndots)=ms*r_EXT'+smin;
        
        % RICARDO - record video
        %imageArrayN=Screen('GetImage', w, [], [], [], []);
        %imwrite(imageArrayN,[num2str(1000+i),'.png'],'png');
        % RICARDO - record video
        
        % merge IN and OUT dots
        xymatrix=[xymatrix_IN xymatrix_EXT];
        vbl=Screen('Flip', w, vbl + (waitframes)*ifi);

        toc
        
    end
    
    if  doEyelink
        Eyelink('StopRecording');
        Eyelink('ReceiveFile', ELfileName, fileparts(vistadispRootPath) ,1);
        
        Eyelink('CloseFile');
        
        Eyelink('Shutdown');
    end
    
    
    K = strfind(myKbCheck,[-1 1]); %find good clicks
    N = strfind(myKbCheck,[0 1]); % find bad clicks
    mygoodclicks = mypress(K+1); % find if good click happened couple times
    myavg = find(diff(mygoodclicks)<1); % correct for that
    K(myavg+1) = []; %remove repetitions
    mygoodclicks_corr = mypress(K+1); % assing to new var
    
    stim.clickaccuracy = length(K)/(length(stim.smiley_change_sec))*100; %calculate accuracy
    sca
    save(sprintf('./stim_params/run%i_%s',r,subj),'stim');
    fprintf('\n\n\n Accuracy = %.2f%% \n',stim.clickaccuracy)
    
    
    
    
    Priority(0);
    ShowCursor
    Screen('CloseAll');
    
catch ME
    
    Priority(0);
    ShowCursor
    Screen('CloseAll');
    
end



function wait4T(keybs)

ch = '';
while ~strcmp(ch,'5');
    [ keyIsDown, timeSecs, keyCode ] = KbCheck(keybs);
    keyPressed= KbName(keyCode);
    if ~isempty(keyPressed)
        ch = keyPressed(1);
    end
end

end