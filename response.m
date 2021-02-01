function result = response(filename)
% quantification of the perturbation response
% 
% Input:    filename........name of the mocap data file, must start with "Mocap"
% Output:   result..........

    % for testing
    if nargin < 1
        testing = 1;
        filename = 'C:\Users\Ton\Cleveland State University\Hala E Osman - Hala data\Par3_PRE\Mocap0001.txt';
        close all
        delete('tmp.ps');
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
    tmdata = interp1(...
        treadmilltimes, treadmilldata.data(:,2:end), mocaptimes);
    
    % find the frame when the perturbation happened, and which foot is in
    % stance phase at that time
	speedcol  = findcolumn(treadmilldata, 'Actual belt speed');
    speed = tmdata(:,speedcol-1);  % -1 because time stamps were removed
    acc = diff(speed)./diff(mocaptimes);
    perturbframe = find(acc > 5.0, 1);  % find first frame where acceleration exceeded 5 m/s2
    perturbtime = mocaptimes(perturbframe);
    fprintf('Perturbation happened at %8.3f seconds into the trial\n', perturbtime-mocaptimes(1));  
    
    % find the heelstrikes before the perturbation, using vertical GRF
    [Lhs, Rhs] = heelstrikes(mocapdata);
    PLhs = max(find(Lhs < perturbframe));
    PRhs = max(find(Rhs < perturbframe));
    if Lhs(PLhs) > Rhs(PRhs)    % which side was the last heelstrike before perturbation?
        fprintf('Left foot was in stance phase during perturbation\n');
    else
        fprintf('Right foot was in stance phase during perturbation\n');
    end
    
    % calculate knee angle (or whatever else we want to use to quantify the
    % perturbation response)
    Rkneeangle = 180/pi * kneeangles(mocapdata, 'R');
    Lkneeangle = 180/pi * kneeangles(mocapdata, 'L');
    
    % plot all cycles before perturbation
    % and superimpose the cycle in which perturbation happened
    avcycles(mocaptimes, Rkneeangle, Rhs, PRhs, perturbtime, 'Par3 PRE Mocap0001 RIGHT KNEE ANGLE');
	avcycles(mocaptimes, Lkneeangle, Rhs, PRhs, perturbtime, 'Par3 PRE Mocap0001 LEFT KNEE ANGLE');
    
    % also take a look at sacrum XYZ trajectories
    % (we can also use Center of Mass, as defined in the MOS analysis)
    Sacrum = marker(mocapdata, 'SACR');
    avcycles(mocaptimes, Sacrum(:,1), Rhs, PRhs, perturbtime, 'Par3 PRE Mocap0001 SACRUM X');
    avcycles(mocaptimes, Sacrum(:,2), Rhs, PRhs, perturbtime, 'Par3 PRE Mocap0001 SACRUM Y');
    avcycles(mocaptimes, Sacrum(:,3), Rhs, PRhs, perturbtime, 'Par3 PRE Mocap0001 SACRUM Z');

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
    
    % append to the tmp.ps file
    print('-dpsc','-append','-bestfit','tmp.ps');
       
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
%======================================================
function xyz = marker(data, name)
% extract a 3D marker trajectory from data
	colX  = findcolumn(data, [name '.PosX']);
    % grab three successive columns from data.data
    xyz = data.data(:,colX+(0:2));
end
%======================================================
function [L_hs, R_hs] = heelstrikes(data)
% find the left and right heelstrikes
    
    % extract the vertical GRF from the data
	col_left_GRF  = findcolumn(data, 'FP1.ForY');
	col_right_GRF = findcolumn(data, 'FP2.ForY');
	time       = data.data(:, 1);
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
    
    % plot the gait cycle durations for checking that data is good
    % (we can perhaps automate this by detecting outliers)
    Ltimes = time(L_hs) - time(1);  % subtract time of frame 1
    Rtimes = time(R_hs) - time(1);
    Lcycletimes = diff(Ltimes);
    Rcycletimes = diff(Rtimes);
    plot(Ltimes(1:end-1),Lcycletimes,Rtimes(1:end-1),Rcycletimes);
    xlabel('time (s)');
    ylabel('gait cycle duration (s)');

end