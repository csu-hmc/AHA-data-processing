function [SL_L, SL_R, ST_L, ST_R] = step_length(data, detail)
	
	% select frames during 20 seconds of normal walking, before the perturbation
	% (NOTE: we could add additional frames, if we know when the perturbation happened)
	st = 1000; % starting data frame
	ed = 3000; % ending data frame
	
% PART 1: step times, this only requires force plate data
    
    % set the threshold (Newtons) for the detection of heelstrike
    threshold = 50;         
	
	% extract the vertical GRF from the data
	col_left_GRF  = findcolumn(data, 'FP1.ForY');
	col_right_GRF = findcolumn(data, 'FP2.ForY');
	time       = data.data(st:ed, 1);
	right_GRF  = data.data(st:ed, col_right_GRF);
	left_GRF   = data.data(st:ed, col_left_GRF);
    
    if (detail)
		% plot GRF for checking
		figure(1)
		subplot(2,1,1)
		plot(time, [left_GRF right_GRF]);
		legend('Left', 'Right')
		ylabel('Ground Reaction Force (N)')
    end
    
   	% filter the GRF to make sure that no rapid change at around 100 N.
	% Heel Z is also filtered to make step length a little more
	% reliable
	cutoff_freq = 6;
	[b, a] = butter(2, (cutoff_freq/(100/2)));
	right_GRF = filtfilt(b, a, right_GRF);
	left_GRF = filtfilt(b, a, left_GRF);
    
    sign_L = left_GRF > threshold;  % sign is 1 when above threshold, zero otherwise
    sign_R = right_GRF > threshold;

    d_signL = diff(sign_L);         % determine the changes in sign
    d_signR = diff(sign_R);

    L_hs = find(d_signL == 1);      % heelstrike is when sign changes from 0 to 1
    R_hs = find(d_signR == 1);

    % force the left and right legs have the same number of steps, cut off by
    % the end steps
    num_step_min = min(length(L_hs), length(R_hs));
    L_hs = L_hs(1:num_step_min);
    R_hs = R_hs(1:num_step_min);    % extract the heel Z trajectory
    
    % calculate step times
    % definitions: Wall et al., Clinical Biomechanics 1987; 2: 119-125
    % NOTE: I changed the variable names, ST_R is the step time when the
    % right leg is stepping (and left is in stance phase)
    if L_hs(1) < R_hs(1)
        ST_R =  time(R_hs) - time(L_hs);
        ST_L =  time(L_hs(2:end)) - time(R_hs(1:end-1));
    else
        ST_R = time(R_hs(2:end)) - time(L_hs(1:end-1));
        ST_L = time(L_hs) - time(R_hs);
    end
    
% PART 2: determine step lengths, this requires marker data
    
	col_LHEEz     = findcolumn(data, 'LHEE.PosZ');
    if isempty(col_LHEEz)
        % markers were probably not labeled so we can't do step length
        disp('step_time: markers not labeled so only doing step time, not step length');
    end
    col_RHEEz     = findcolumn(data, 'RHEE.PosZ');
	right_HEEL  = data.data(st:ed, col_RHEEz);
	left_HEEL  = data.data(st:ed, col_LHEEz);

	% Heel Z is filtered to make step length a little more reliable
    % this only works if there are no remaining gaps in the data
    if isempty(find(isnan(right_HEEL)))
        right_HEEL = filtfilt(b, a, right_HEEL);
    end
    if isempty(find(isnan(left_HEEL)))
        left_HEEL  = filtfilt(b, a, left_HEEL); 
    end
    
	if (detail)
		% plot HEEL Z for checking
		figure(1)
		subplot(2,1,2)
		plot(time, [left_HEEL right_HEEL]);
		xlabel('Time (s)')
		ylabel('Heel Marker Z (m)')
	end

    % calculate step lengths, but only if less than 50% of each heel is missing
    % remember that negative Z is forward (towards the screen)
    % NOTE: I changed the variable names based on the definitions by Wall.
    % left step length is a step by the left leg
    if (missing(right_HEEL) < 0.5) && (missing(left_HEEL) < 0.5)
        SL_L = right_HEEL(L_hs) - left_HEEL(L_hs);  % at left heelstrike, how far left heel is in front of right heel
        SL_R =  left_HEEL(R_hs) - right_HEEL(R_hs);
    else
        SL_L = NaN(size(L_hs));
        SL_R = NaN(size(R_hs));
    end
  
    % remove the outliers
    SL_L = rmoutliers(SL_L);
    SL_R = rmoutliers(SL_R);
    ST_L = rmoutliers(ST_L);
    ST_R = rmoutliers(ST_R);
    
    if (detail)
    	disp('Click Continue to Continue')
		keyboard
    end
end
%==================================================
function fraction = missing(x)
    missingframes = find(isnan(x));
    fraction = numel(missingframes) / numel(x);
end