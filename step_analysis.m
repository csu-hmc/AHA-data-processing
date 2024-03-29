function result = step_analysis(data, detail)
    % analysis of step time and step length
    % data should be a data structure from a mocap file, from getdata.m
    % use detail=1 to show and pause the results
    % result is a table with one row and 8 columns: mean and SD of step time and step length (left and right)
        
    % if no input is specified, we use one particular file for testing
    if nargin < 1
        detail = 1;
        data = getdata('Par12_POST\Mocap0008.txt', 1);  % we only need the mocap file
    end

    fprintf('Step analysis for %s\n', data.name);
    
    % find the heelstrikes and only keep those that occurred between 10 and 30 seconds
	% (NOTE: we could add additional frames after the perturbation, if we know when the perturbation happened)
	L_hs = data.Lhs;
    R_hs = data.Rhs;
    st = 10; % starting time
	ed = 30; % ending time
    time = data.data(:,1)-data.data(1,1);  % time relative to start of file
    L_hs = L_hs( (time(L_hs) >= st) & (time(L_hs) <= ed) );
    R_hs = R_hs( (time(R_hs) >= st) & (time(R_hs) <= ed) );
    
    % find the left step times, which are RHS to LHS, as long as no RHS comes before the LHS
    ST_L = [];
    for i = 1:numel(R_hs)-1
        nextLHS = min(L_hs(L_hs > R_hs(i)));
        if nextLHS < R_hs(i+1)
            ST_L = [ST_L time(nextLHS)-time(R_hs(i))];
        end
    end
    
    % find the right step times, which are LHS to RHS, as long as no LHS comes before the RHS
    ST_R = [];
    for i = 1:numel(L_hs)-1
        nextRHS = min(R_hs(R_hs > L_hs(i)));
        if nextRHS < L_hs(i+1)
            ST_R = [ST_R time(nextRHS)-time(L_hs(i))];
        end
    end
        
    % step lengths, this requires marker data   
	left_HEEL  = getcolumn(data, 'LHEE.PosZ');
    right_HEEL = getcolumn(data, 'RHEE.PosZ');

	% Heel Z is filtered to make step length a little more reliable
    % this only works if there are no remaining gaps in the data
    cutoff_freq = 6;
	[b, a] = butter(2, (cutoff_freq/(100/2)));
    if isempty(find(isnan(right_HEEL)))
        right_HEEL = filtfilt(b, a, right_HEEL);
    else
        disp('right_HEEL not filtered because it has gaps.');
    end
    if isempty(find(isnan(left_HEEL)))
        left_HEEL  = filtfilt(b, a, left_HEEL); 
    else
        disp('left_HEEL not filtered because it has gaps.');
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
    
    if (0)     % use if(1) to see details of the step length analysis
        close all
        H = [right_HEEL left_HEEL];
        Hmax = max(max(H));
        Hmin = min(min(H));
        plot(time,H)
        hold on
        hs = [L_hs ; R_hs];
        for i=1:numel(hs)
            plot(time([hs(i) hs(i)]),[Hmin Hmax],'k')
        end
        set(gca,'XLim',[time(min(hs))-1 time(max(hs))+1]);
        legend('LHEE.PosZ','RHEE.PosZ')
        keyboard
    end
  
    % remove the outliers
    SL_L = rmoutliers(SL_L);
    SL_R = rmoutliers(SL_R);
    ST_L = rmoutliers(ST_L);
    ST_R = rmoutliers(ST_R);
    
    % calculate mean and SD for all variables, and store in result
    ST_left_mean  = mean(ST_L);
    ST_left_SD    = std(ST_L);
    ST_right_mean = mean(ST_R);
    ST_right_SD   = std(ST_R);
    SL_left_mean  = mean(SL_L);
    SL_left_SD    = std(SL_L);
    SL_right_mean = mean(SL_R);
    SL_right_SD   = std(SL_R);
    result = table(ST_left_mean,ST_left_SD,ST_right_mean,ST_right_SD, ...
                   SL_left_mean,SL_left_SD,SL_right_mean,SL_right_SD);
    
    if (detail)
		% plot the step times and step lengths
        close all
		figure(1)
		subplot(2,1,1)
		plot(ST_R);hold on; plot(ST_L);
		ylabel('step time (s)')
		subplot(2,1,2)
		plot(SL_R);hold on; plot(SL_L);
		xlabel('step number')
		ylabel('step length (m)')
        legend('Right','Left');

        fprintf('Left  step time:   %8.3f %s %8.3f s.\n',ST_left_mean,char(177),ST_left_SD);
        fprintf('Right step time:   %8.3f %s %8.3f s.\n',ST_right_mean,char(177),ST_right_SD);
        fprintf('Left  step length: %8.3f %s %8.3f s.\n',SL_left_mean,char(177),SL_left_SD);
        fprintf('Right step length: %8.3f %s %8.3f s.\n',SL_right_mean,char(177),SL_right_SD);
        disp('Check Figure 1, and the results printed above, for problems');
        disp('(affected side should have lower step length, can even be negative if foot drags behind)');
    	disp('Hit ENTER to Continue')
		pause
    end
    disp('Step analysis (step_analysis.m) completed.');
end
%==================================================
function fraction = missing(x)
    missingframes = find(isnan(x));
    fraction = numel(missingframes) / numel(x);
end