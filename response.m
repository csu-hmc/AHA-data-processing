function result = response(filename)
% quantification of the perturbation response
% 
% Input:    filename........name of the mocap data file, must start with "Mocap"
% Output:   result..........
    global testing      % a logical flag that is true when testing

    % if no file is specified, we use one particular file for testing
    if nargin < 1
        testing = 1;
        shortname = 'Par2_PRE\Mocap0001.txt';
        computer = getenv('COMPUTERNAME');
        if strcmp(computer,'LRI-102855')
            filename = ['C:\Users\Ton\Cleveland State University\Hala E Osman - Hala data\' shortname];    
        elseif strcmp(computer,'DESKTOP-0HN0T6U')
            filename = ['C:\Users\hallo\OneDrive - Cleveland State University\Hala data\' shortname];    
        else
            fprintf('Your computername is: %s\n', computer);
            fprintf('Please edit response.m by changing yourcomputername and yourdatafolder in lines 15-16\n');
            error('exiting now');
        end
        close all
        if exist('tmp.ps')
            delete('tmp.ps');
        end
    else
        testing = 0;
    end
    
    % import the data
    mocapdata = importdata(filename);
    mocapdata = cleanup(mocapdata);   % interpolate missing markers, and insert NaN for large gaps
    treadmillfile = strrep(filename,'Mocap','Treadmill');
    treadmilldata = importdata(treadmillfile);
    
    % resample the treadmill data at the same time stamps as the mocap data
    treadmilltimes = treadmilldata.data(:,1);  % column 1 of treadmill file
    mocaptimes     = mocapdata.data(:,1);      % column 1 of mocap file
    % Note: currently we don't actually use the resampled treadmill data
    tmdata = interp1(...
        treadmilltimes, treadmilldata.data(:,2:end), mocaptimes);
    
    % find the frame and time when the perturbation happened
	speedcol  = findcolumn(treadmilldata, 'Actual belt speed');
    speed = treadmilldata.data(:,speedcol); 
    acc = diff(speed)./diff(treadmilltimes);
    perturbtime = treadmilltimes(find(acc > 5.0, 1));  % find first time when acceleration exceeded 5 m/s2
    fprintf('Perturbation happened at %8.3f seconds into the trial\n', perturbtime-mocaptimes(1));  
    [~,perturbframe] = min(abs(mocaptimes-perturbtime));
    
    % find the heelstrikes using vertical GRF
    [Lhs, Rhs] = heelstrikes(mocapdata);
    PLhs = find(Lhs < perturbframe, 1, 'last' );
    PRhs = find(Rhs < perturbframe, 1, 'last' );
    
    % find out which foot was in stance phase during perturbation
    % and calculate gait deviation   
    % the gaitdeviation function also estimates missing mocapdata
    if Lhs(PLhs) > Rhs(PRhs)    % which side was the last heelstrike before perturbation?
        fprintf('Left foot was in stance phase during perturbation\n');
        [t1,t2,mocapdata] = gaitdeviation(mocapdata, Lhs, perturbtime);
    else
        fprintf('Right foot was in stance phase during perturbation\n');
        [t1,t2,mocapdata] = gaitdeviation(mocapdata, Rhs, perturbtime);
    end
    
    if (testing)
        % change backslashes and underscores, they mess up the plot labels
        % because they are Latex comments
        shortname = strrep(shortname,'\','/');
        shortname = strrep(shortname,'_','\_');
     
        % calculate some descriptive variables
        Rkneeangle = 180/pi * kneeangles(mocapdata, 'R');
        Lkneeangle = 180/pi * kneeangles(mocapdata, 'L');
        Rhipangle = 180/pi * hipangles(mocapdata, 'R');
        Lhipangle = 180/pi * hipangles(mocapdata, 'L');

        % anterior position of heel relative to sacrum
        Rheelpos = footplacement(mocapdata, 'R');
        Lheelpos = footplacement(mocapdata, 'L'); 

        % plot average and SD of normal gait, superimpose the gait just
        % before and after the perturbation
        avcycles(mocaptimes, Rkneeangle, Rhs, PRhs, perturbtime, [shortname 'RIGHT KNEE ANGLE']);
        avcycles(mocaptimes, Lkneeangle, Rhs, PRhs, perturbtime, [shortname 'LEFT KNEE ANGLE']);
        avcycles(mocaptimes, Rhipangle, Rhs, PRhs, perturbtime, [shortname 'RIGHT HIP ANGLE']);
        avcycles(mocaptimes, Lhipangle, Rhs, PRhs, perturbtime, [shortname 'LEFT HIP ANGLE']);
        avcycles(mocaptimes, -Rheelpos(:,3), Rhs, PRhs, perturbtime, [shortname 'R anterior foot pos']);
        avcycles(mocaptimes, -Lheelpos(:,3), Rhs, PRhs, perturbtime, [shortname 'L anterior foot pos']);
        avcycles(mocaptimes, Rheelpos(:,1), Rhs, PRhs, perturbtime, [shortname 'R lateral foot pos']);
        avcycles(mocaptimes, -Lheelpos(:,1), Rhs, PRhs, perturbtime, [shortname 'L lateral foot pos']);

        % also take a look at sacrum XYZ trajectories
        % (we can also use Center of Mass, as defined in the MOS analysis)
        Sacrum = marker(mocapdata, 'SACR');
        avcycles(mocaptimes, Sacrum(:,1), Rhs, PRhs, perturbtime, [shortname 'SACRUM X(right)']);
        avcycles(mocaptimes, Sacrum(:,2), Rhs, PRhs, perturbtime, [shortname 'SACRUM Y(up)']);
        avcycles(mocaptimes, Sacrum(:,3), Rhs, PRhs, perturbtime, [shortname 'SACRUM Z(posterior)']);
    end
    
    % store result of the calculations in the result structure
    result.t1 = t1;
    result.t2 = t2;
    result.mocapdata = mocapdata;
    
end
%========================================================
function [avtime, avvar, sdvar] = avcycles(times, var, hs, phs, ptime, titletext)
    % average variable "var" across all the cycles (heelstrike to heelstrike)
    
    % split and time-normalize all gait cycles
    ncycles = numel(hs)-1;
    nsam = 500;  % number of samples per cycle
    allvar = zeros(ncycles, nsam);
    cycletimes = zeros(ncycles,1);
    for i = 1:ncycles
        cycletimes(i) = times(hs(i+1)) - times(hs(i));
        % we will resample nsam+1 samples, the first and last are at
        % heelstrike
        newtimes = times(hs(i)) + (0:nsam-1)*cycletimes(i)/nsam;
        allvar(i,:) = interp1(times(hs(i):hs(i+1)),var(hs(i):hs(i+1)),newtimes);
    end
    
    % create an average of the gait cycles before the perturbation
    cycletime = mean(cycletimes(1:phs-1));
    avvar = mean(allvar(1:phs-1,:),'omitnan')';   % NaN values are not included in average
    sdvar = std(allvar(1:phs-1,:),'omitnan')';
    avtime = (0:nsam-1)' / nsam;
    
    % create shaded area for several normal gait cycles
    ncycles = 10;
    tcycles = [];
    for i = (-ncycles):(ncycles-1)
        tcycles = [tcycles ; avtime+i];
    end
    sdx = [tcycles; flipud(tcycles)];
    sdy = [repmat(avvar-sdvar,2*ncycles,1); flipud(repmat(avvar+sdvar,2*ncycles,1))];
    color = [0.6 0.6 0.6];
    figure();
    set(gcf,'Position',[8 610 1232 488]);
    subplot(2,1,1);
    fill(sdx,sdy,color,'LineStyle','none');
    % ylabel('knee flexion (deg)');
    title(titletext);
    % set(gca,'YLim',[0 80]);
    
    % superimpose individual (time-normalized) gait cycles before and after the perturbation
    pvar = reshape(allvar(phs-ncycles:phs+ncycles-1,:)', 2*ncycles*nsam, 1);
    hold on
    plot(tcycles, pvar,'r');
    
    % draw a vertical dashed line when the perturbation happened
    ptime = (ptime - times(hs(phs))) / mean(cycletimes);
    ylim = get(gca,'YLim');
    plot([ptime ptime], ylim,'r--');
    text(ptime+0.3, (ylim(1)+9*ylim(2))/10, 'perturbation time','Color','r')
    
    % in a subplot, show deviation normalized to SD
    deviation = abs(repmat(avvar,2*ncycles,1) - pvar);
    deviation = deviation ./ repmat(sdvar,2*ncycles,1);
    subplot(2,1,2)
    plot(tcycles,deviation);
    set(gca,'YLim',[0 10]);
    xlabel('time (gait cycles, zero is at heelstrike before perturb)');
    ylabel('perturbation effect (SD)');
       
end
%=======================================================
function [angles] = kneeangles(data, side)
    Hip   = marker(data, [side 'GTRO']);
    Knee  = marker(data, [side 'LEK']);
    Ankle = marker(data, [side 'LM']);
    
    % construct unit vectors from hip to knee and from knee to ankle
    thigh = Knee - Hip;
    shank = Ankle - Knee;
    crossnorm = vecnorm(cross(thigh,shank),2,2);  % norm of the cross product
    thighnorm = vecnorm(thigh,2,2);               % norm of the thigh vector
    shanknorm = vecnorm(shank,2,2);               % norm of the shank vector
    
    % angle is inverse sine of (cross product divided by both vector norms)
    angles = asin(crossnorm ./ thighnorm ./ shanknorm);
    
end
%==========================================================
function [pos] = footplacement(data, side)
    % calculate position of heel relative to sacrum
    Heel   = marker(data, [side 'HEE']);
    Sacrum  = marker(data, 'SACR');
    pos = Heel - Sacrum;   
end
%======================================================
function xyz = marker(data, name)
% extract a 3D marker trajectory from data
	colX  = findcolumn(data, [name '.PosX']);
    % grab three successive columns from data.data
    xyz = data.data(:,colX+(0:2));
end
%=======================================================
function [angles] = hipangles(data, side)
    Hip   = marker(data, [side 'GTRO']);
    Knee  = marker(data, [side 'LEK']);
    Shoulder = marker(data, [side 'SHO']);
    
    % construct unit vectors from shoulder to hip and hip to knee
    trunk = Hip - Shoulder;
    thigh = Knee - Hip;
    crossnorm = vecnorm(cross(trunk,thigh),2,2);  % norm of the cross product
    thighnorm = vecnorm(thigh,2,2);               % norm of the thigh vector
    trunknorm = vecnorm(trunk,2,2);               % norm of the shank vector
    
    % angle is inverse sine of (cross product divided by both vector norms)
    angles = asin(crossnorm ./ thighnorm ./ trunknorm);
    
end
%======================================================
function [L_hs, R_hs] = heelstrikes(data)
% find the left and right heelstrikes
    
    % extract the vertical GRF from the data
	col_left_GRF  = findcolumn(data, 'FP1.ForY');
	col_right_GRF = findcolumn(data, 'FP2.ForY');
	time       = data.data(:, 1) - data.data(1,1);  % time from start of trial
	right_GRF  = data.data(:, col_right_GRF);
	left_GRF   = data.data(:, col_left_GRF);

    % filter the GRF to make sure that no rapid change at around 100 N.
	cutoff_freq = 6;
	[b, a] = butter(2, (cutoff_freq/(100/2)));
	right_GRF = filtfilt(b, a, right_GRF);
	left_GRF = filtfilt(b, a, left_GRF);
    
    % find the times when GRF goes above the threshold
    threshold = 50;   % newtons
    sign_L = left_GRF > threshold;  % sign is 1 when above threshold, zero otherwise
    sign_R = right_GRF > threshold;
    d_signL = diff(sign_L);         % determine the changes in sign
    d_signR = diff(sign_R);
    L_hs = find(d_signL == 1);      % heelstrike is when sign changes from 0 to 1
    R_hs = find(d_signR == 1);
    
    % remove any heelstrikes before 2 seconds (walking is not quite
    % periodic yet)
    L_hs = L_hs(time(L_hs) > 2.0);
    R_hs = R_hs(time(R_hs) > 2.0);

    % plot the gait cycle durations for checking that data is good
    % (we can perhaps automate this by detecting outliers)
    Ltimes = time(L_hs);  % subtract time of frame 1
    Rtimes = time(R_hs);
    Lcycletimes = diff(Ltimes);
    Rcycletimes = diff(Rtimes);
    plot(Ltimes(1:end-1),Lcycletimes,Rtimes(1:end-1),Rcycletimes);
    xlabel('time (s)');
    ylabel('gait cycle duration (s)');
    
    

end
%=========================================================================================
function [t1,t2,newdata] = gaitdeviation(data, hs, perturbtime)   
% calculate gait deviation timing (t1,t2) using Mahalanobis Distance
% and also estimate missing data
    global testing
    
    newdata = data;
    nsamples = size(data.data,1);
    t = data.data(:,1) - data.data(1,1);
    tperturb = perturbtime - data.data(1,1); % perturbation time relative to start of trial

    % find the columns with marker data and remove those with too many NaNs
    columns = find(contains(data.colheaders,'.Pos'))';
    remove = find(sum(isnan(data.data(:,columns))) > 0.5*nsamples);
    for i=1:numel(remove)
        fprintf('gait deviation analysis will ignore %s\n', data.colheaders{columns(remove(i))});
    end
    columns = setdiff(columns,columns(remove));  % these columns are kept
    ncolumns = numel(columns);

    % estimate the covariance matrix using frames from gait cycles before the
    % perturbation, and those after the perturbation
    disp('Estimating mean and covariance...');
    i1 = find(t(hs) < tperturb,    1, 'last');  % last heelstrike before perturbation
    i2 = find(t(hs) > tperturb+10, 1, 'first'); % first heelstrike ten seconds after perturbation
    normalframes = [hs(1):hs(i1) hs(i2):hs(end)];
    normalcycles = [hs(1:i1-1) hs(2:i1) ; hs(i2:end-1) hs(i2+1:end)];
    [mu,C] = ecmnmle(data.data(normalframes,columns));
    Cinv = inv(C);
   
    % impute missing data using idea from Rasmussen 2020
    disp('Estimating missing data...');
    tic
    for i = 1:nsamples
        x = data.data(i,columns)';
        m = find(isnan(x));  % elements of x that are missing data
        if ~isempty(m)
            if (toc>5)
                fprintf('Estimating %d coordinates in frame %d (%4.1f %% done)\n',numel(m),i,100*i/nsamples);
                tic
            end
            o = find(~isnan(x)); % elements of x that are observed
            Cmm = Cinv(m,m);
            Cmo = Cinv(m,o);
            x(m) = mu(m) - inv(Cmm)*Cmo*(x(o)-mu(o));
        end
        newdata.data(i,columns) = x';      % store in newdata
    end
    % we could now estimate mean and covariance again, but result is
    % extremely close to what we already have (I checked)

    % Calculate the Mahalanobis distance (squared) of each frame to the mean
    % We could also do this from the original incomplete data (see
    % notebook 3/23/2021) but result is exactly the same.
    % We could also do this after subtracting the average gait cycle, but
    % again, result is (almost) exactly the same.
    % This quantity is also known at the "Hotelling T-squared" statistic
    drel = newdata.data(:,columns)-repmat(mu',nsamples,1);
    T2 = dot(drel * Cinv, drel, 2);
            
    % low-pass filter to smooth the T2 curve
    Fs = 1000;
    Fc = 3;
    tnew = (t(hs(1)):(1/Fs):t(hs(end)))';  % resample to ensure constant sampling rate
    T2new = interp1(t,T2,tnew);
    [b,a] = butter(2,Fc/(Fs/2));
    T2f = filtfilt(b,a,T2new);
    
    % calculate p values from the Hotelling T-squared statistic
    % See https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
    np = ncolumns;
    n = nsamples;
    p = fcdf(T2f*(n-np)/np/(n-1),np,n-np,'upper');
    
    % gait deviation measure is the duration of the period when p<threshold
    threshold = 0.01;
    d = [0 ; diff(p<threshold)];  % diff detects when p went across the threshold
    rt1 = tnew(find(d==1  & tnew>tperturb,1,'first'));  % first time that p went below threshold
    rt2 = tnew(find(d==-1 & tnew>tperturb,1,'first'));  % first time that p went above threshold again
    t1 = rt1 - tperturb;
    t2 = rt2 - tperturb;
    fprintf('Start    of perturbation recovery: %.3f s. after perturbation\n',t1);
    fprintf('End      of perturbation recovery: %.3f s. after perturbation\n',t2);
    fprintf('Duration of perturbation recovery: %.3f s.\n',t2-t1);
    
    % make a figure to illustrate what was calculated
    % (only when testing)
    if (testing)
        xlim = [tperturb-2 tperturb+6]; % the time range we want to plot
        figure;
        subplot(2,1,1);
        plot(tnew,[T2new T2f]);
        title('Hotelling T^2 statistic');
        xlabel('time (s)');
        legend('T^2 unfiltered',num2str(Fc,'%.1f Hz double 2nd order'));
        set(gca,'XLim',xlim);

        subplot(2,1,2)
        ylim = [1e-6 1];
        semilogy(tnew,p);
        title('probability of normal gait');
        xlabel('time (s)');
        hold on
        semilogy([tperturb tperturb],get(gca,'YLim'),'r'); % show perturbation time
        text(tperturb+0.1,0.5,'perturbation','Color','r');
        semilogy([rt1 rt1],get(gca,'YLim'),'k--'); % show t1
        semilogy([rt2 rt2],get(gca,'YLim'),'k--'); % show t2
        semilogy(xlim,[threshold threshold],'k--');
        set(gca,'YLim',ylim);
        set(gca,'XLim',xlim);
        p = get(gca,'Position');  % we can use this to calculate coordinates relative to Figure window
        xrel1 = p(1) + p(3)*(tperturb+[0 t1]-xlim(1))/diff(xlim);
        xrel2 = p(1) + p(3)*(tperturb+[0 t2]-xlim(1))/diff(xlim);
        yrel1 = p(2) + p(4)*( (log(threshold)-log(ylim(1)) )/diff(log(ylim)) - 0.1);
        yrel2 = p(2) + p(4)*( (log(threshold)-log(ylim(1)) )/diff(log(ylim)) - 0.3);
        annotation('textarrow',xrel1,[yrel1 yrel1],'String', num2str(t1,'t1 = %.3f s.'));
        annotation('textarrow',xrel2,[yrel2 yrel2],'String', num2str(t2,'t2 = %.3f s.'));
        
        % append to the tmp.ps file
        print('-dpsc','-append','-bestfit','tmp.ps');
    end
end
