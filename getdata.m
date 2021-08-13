function [mocapdata,treadmilldata] = getdata(name, detail)
    % the name must be, for example: 'Par7_PRE\mocap0001.txt'
    
    % settings
    threshold = 50;     % threshold (N) for heelstrike detection
    cutoff_freq = 6;    % GRF filter cutoff frequency for heelstrike detection
   
    % find the path to the shared data folders
    path = getpath();  

    % import the mocap data
    mocapdata = importdata([path name]);
    mocapdata = cleanup(mocapdata);   % interpolate missing markers, and insert NaN for large gaps
    mocapdata.name = name;  % add the short filename to the data struct
    % make it suitable for Latex interpreter
    latexname = strrep(name,'\','/');
    mocapdata.latexname = strrep(latexname,'_','\_');

    % make new time stamps (we know that Cortex sends data at exactly 100 Hz
    mocapdata.data(:,1) = mocapdata.data(1,1) + 0.01 * (0:(size(mocapdata.data,1)-1));
   
    % construct treadmill filename from the trial number and import that file too
    if nargout > 1   % do this only if the treadmill data is needed
        num = extract(name, digitsPattern);
        i = strfind(name,'\'); % find the backslashes in the mocap file name
        name = [name(1:i(end)) 'Treadmill' num{end} '.txt']; % make the treadmill file name
        treadmilldata = importdata([path name]);
        treadmilldata.name = name;
    end
    
    % finally, determine the heelstrikes (we need them everywhere)
    
    % extract the vertical GRF and do some smoothing
    left_GRF  = getcolumn(mocapdata, 'FP1.ForY');
	right_GRF = getcolumn(mocapdata, 'FP2.ForY');    
    time = mocapdata.data(:,1) - mocapdata.data(1,1);  % time stamps, relative to start of file
	[b, a] = butter(2, (cutoff_freq/(100/2)));
	right_GRF = filtfilt(b, a, right_GRF);
	left_GRF = filtfilt(b, a, left_GRF);    
    
    % detect the heelstrikes
    sign_L = left_GRF > threshold;  % sign is 1 when above threshold, zero otherwise
    sign_R = right_GRF > threshold;
    d_signL = diff(sign_L);         % determine the changes in sign
    d_signR = diff(sign_R);
    Lhs = find(d_signL == 1);      % heelstrike is when sign changes from 0 to 1
    Rhs = find(d_signR == 1);
    
    % remove any heelstrikes before 2 seconds (walking is not quite
    % periodic yet)
    Lhs = Lhs(time(Lhs) > 2);
    Rhs = Rhs(time(Rhs) > 2);
    
    % store them in the mocapdata
    mocapdata.Lhs = Lhs;
    mocapdata.Rhs = Rhs;
  
    if (detail)
        % make a plot for verification
        close all
		figure(1)
        screen = get(0,'screensize');
        set(gcf,'Position',[1 floor(0.5*screen(4)) 1280 floor(0.5*screen(4))]); % use top of the screen
		
        % break it up into 3 subplots
        duration = time(end)/3;  % duration of each plot
        for i = 1:3
            subplot(3,1,i);
            plot(time, left_GRF,'r', time, right_GRF, 'b');
    		ylabel('GRF (N)')
            hold on
            plotheelstrikes(time,Lhs,'r');
            plotheelstrikes(time,Rhs,'b');
            set(gca,'XLim',duration*[(i-1) i]);
        end
        legend('Left', 'Right');
        xlabel('time(s)');

        disp('Check Figure 1 for correct heelstrike detection in this trial');
        disp('(problems in the last 30 seconds can be ignored)');
        disp('Hit ENTER to continue, CTRL-C to quit');
        pause
    end
end
%==================================================
function plotheelstrikes(time,hs,color)
    ylim = get(gca,'YLim');
    for i=1:numel(hs)
        t = time(hs(i));
        plot([t t], ylim, color, 'LineWidth',2);
    end
end

    
    