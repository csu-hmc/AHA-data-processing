function result = mos(data, detail)
    % analysis of Margin of Stability
    % data should be a data structure from a mocap file, from getdata.m
    % use detail=1 to show and pause the results
    %
    % result is a table with one row and the following columns:
    % MOS_AP_left, MOS_AP_right, MOS_ML_right, MOS_ML_left
    % These are measured at each heelstrike, according to McAndrew Young et al 2012.

    % if no input is specified, we use one particular file for testing
    if nargin < 1
        detail = 1;
        data = getdata('Par3_PRE\Mocap0007.txt', 1);  % we only need the mocap data
    end
    
    fprintf('MoS analysis for %s\n', data.name);
    
    % find the heelstrikes and only keep those that occurred between 10 and 30 seconds
	% (NOTE: we could add additional frames after the perturbation, if we know when the perturbation happened)
	Lhs = data.Lhs;
    Rhs = data.Rhs;
    st = 10; % starting time
	ed = 30; % ending time
    time = data.data(:,1)-data.data(1,1);  % time relative to start of file
    Lhs = Lhs( (time(Lhs) >= st) & (time(Lhs) <= ed) );
    Rhs = Rhs( (time(Rhs) >= st) & (time(Rhs) <= ed) );
    
    % extract the required marker data   
    SACRx = getcolumn(data, 'SACR.PosX');
    SACRy = getcolumn(data, 'SACR.PosY');
    SACRz = getcolumn(data, 'SACR.PosZ');
	LHEEy = getcolumn(data, 'LHEE.PosY');
	LHEEz = getcolumn(data, 'LHEE.PosZ');
    RHEEy = getcolumn(data, 'RHEE.PosY');
    RHEEz = getcolumn(data, 'RHEE.PosZ');
	LTOEz = getcolumn(data, 'LTOE.PosZ');
    RTOEz = getcolumn(data, 'RTOE.PosZ');
	LMT5x = getcolumn(data, 'LMT5.PosX');
    RMT5x = getcolumn(data, 'RMT5.PosX');
    
    % calculate the mean distance from heel marker to Sacrum (in the sagittal YZ plane)
    Rdistance = sqrt( (SACRz-RHEEz).^2 + (SACRy-RHEEy).^2 );
    Ldistance = sqrt( (SACRz-LHEEz).^2 + (SACRy-LHEEy).^2 );
    L = mean([Rdistance ; Ldistance]); % everage leg length, used for XCoM calculation
    
    % calculate the Sacrum velocity in X and Z direction, and the extrapolated center of mass
    SACRvx = velocity(time, SACRx);
    SACRvz = velocity(time, SACRz);
    g = 9.81;
    w0 = sqrt(g/L);
    XCoMx = SACRx + SACRvx / w0;
    XCoMz = SACRz + SACRvz / w0;
    
    % find the MOS at the left and right heelstrikes
    % sign convention: MOS is positive when XCoM is inside BOS
    MOS_AP_left  = ( XCoMz(Lhs) - LTOEz(Lhs) );  % because Z is posterior
    MOS_AP_right = ( XCoMz(Rhs) - RTOEz(Rhs) );  % because Z is posterior
    MOS_ML_left  = ( XCoMx(Lhs) - LMT5x(Lhs) );  % because X is to the right
    MOS_ML_right = -( XCoMx(Rhs) - RMT5x(Rhs) );  % because X is to the right
        
    % calculate mean and SD for all variables, and store in result
    MOS_AP_left_mean  = mean(MOS_AP_left);
    MOS_AP_left_SD    = std(MOS_AP_left);
    MOS_AP_right_mean = mean(MOS_AP_right);
    MOS_AP_right_SD   = std(MOS_AP_right);
    MOS_ML_left_mean  = mean(MOS_ML_left);
    MOS_ML_left_SD    = std(MOS_ML_left);
    MOS_ML_right_mean = mean(MOS_ML_right);
    MOS_ML_right_SD   = std(MOS_ML_right);
    result = table(MOS_AP_left_mean, MOS_AP_left_SD, MOS_AP_right_mean, MOS_AP_right_SD, ...
                   MOS_ML_left_mean, MOS_ML_left_SD, MOS_ML_right_mean, MOS_ML_right_SD);
    
    if (detail)
		% plot the mos variables for all heelstrikes
        close all
		figure(1)
		subplot(2,1,1)
		plot(MOS_AP_left);hold on; plot(MOS_AP_right);
		ylabel('AP MOS (m)')
		subplot(2,1,2)
		plot(MOS_ML_left);hold on; plot(MOS_ML_right);
		xlabel('step number')
		ylabel('ML MOS (m)')
        legend('Right','Left');

        fprintf('Length L for XCoM calculation: %8.3f m.\n', L);
        fprintf('Left  AP MoS:   %8.3f %s %8.3f s.\n',MOS_AP_left_mean, char(177),MOS_AP_left_SD);
        fprintf('Right AP MoS:   %8.3f %s %8.3f s.\n',MOS_AP_right_mean,char(177),MOS_AP_right_SD);
        fprintf('Left  ML MoS:   %8.3f %s %8.3f s.\n',MOS_ML_left_mean, char(177),MOS_ML_left_SD);
        fprintf('Right ML MoS:   %8.3f %s %8.3f s.\n',MOS_ML_right_mean,char(177),MOS_ML_right_SD);
        disp('Check Figure 1, and the results printed above, for problems');
    	disp('Hit ENTER to Continue')
		pause
    end
    disp('MoS analysis (mos.m) completed.');

end
%==================================================
function fraction = missing(x)
    missingframes = find(isnan(x));
    fraction = numel(missingframes) / numel(x);
end
%==================================================
function v = velocity(t,x)
    % find the velocity v(t) for the position data x(t)
    cutoff_freq = 6.0;  % Hz, low pass filter to reduce noise in velocity
    Fs = 1/mean(diff(t)); % sampling frequency
    if std(diff(t)) > 0.001
        error('velocity calculation: sampling rate is not constant');
    end
    [b, a] = butter(2, (cutoff_freq/(Fs/2)));
	x = filtfilt(b, a, x);
    v = diff(x)./diff(t);
    v = [v(1) ; v];  % duplicate sample 1 so we have same number of samples in v and x
end