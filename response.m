function [result, mocapdata] = response(mocapdata, treadmilldata, options)
% quantification of the perturbation response
% 
% Input:    mocapdata.......mocap data, from getdata() function
%           treadmilldata...treadmill data, from getdata() function
%           options.........struct which can have the following fields:
%               filter....cutoff frequency (Hz) for smoothing of perturbation response (default: 1.0);    % use a 1 hz low pass filter for the Hotelling T-squared
%               pvalue....threshold for detecing abnormal gait with the T-squared (default 0.001)
%               markerset17...(logical) use the set of 17 markers, rather than all markers (default 0)
%               testing.......(logical) produce extra output for testing (default 1)
%
% Output:   result..........table with one row and the following columns:
%               t1............time when the perturbation response starts (relative to perturbation)
%               t2............time when the perturbation response ends
%               T2max.........maximum of the Hotelling T-squared after perturbation
%           mocapdata.....mocap data after PCA reconstruction (no missing markers)

    % use defaults for the options that were not specified
    if nargin < 3 , options = []; end
    if ~isfield(options,'filter') , options.filter = 1.0; end
    if ~isfield(options,'pvalue') , options.pvalue = 0.001; end
    if ~isfield(options,'testing') , options.testing = 0; end
    if ~isfield(options,'markerset17') , options.markerset17 = 1; end
    
    % if no file is specified, we use one particular file for testing
    if nargin < 1
        options.testing = 1;
        [mocapdata,treadmilldata] = getdata('Par3_PRE\Mocap0007.txt', 1);
    end
    
    fprintf('Perturbation response analysis for %s\n', mocapdata.name);
            
    % extract column 1 (time stamps) from both files
    treadmilltimes = treadmilldata.data(:,1);  % column 1 of treadmill file
    mocaptimes     = mocapdata.data(:,1);      % column 1 of mocap file
    
    % find the frame and time when the perturbation happened
	speed = getcolumn(treadmilldata, 'Actual belt speed');
    acc = diff(speed)./diff(treadmilltimes);
    perturbtime = treadmilltimes(find(acc > 5.0, 1));  % find first time when acceleration exceeded 5 m/s2
    fprintf('Perturbation happened at %8.3f seconds into the trial\n', perturbtime-mocaptimes(1));  
    [~,perturbframe] = min(abs(mocaptimes-perturbtime));  % find the closest mocap frame
    
    % find the heelstrikes of both feet before the perturbation
    Lhs = mocapdata.Lhs;
    Rhs = mocapdata.Rhs;
    PLhs = find(Lhs < perturbframe, 1, 'last' );
    PRhs = find(Rhs < perturbframe, 1, 'last' );
    
    % find out which foot was in stance phase during perturbation
    % and calculate gait deviation   
    % the gaitdeviation function also estimates missing mocapdata
    if Lhs(PLhs) > Rhs(PRhs)    % which side was the last heelstrike before perturbation?
        fprintf('Left foot was in stance phase during perturbation\n');
        [result,mocapdata] = gaitdeviation(mocapdata, Lhs, perturbtime, options);
    else
        fprintf('Right foot was in stance phase during perturbation\n');
        [result,mocapdata] = gaitdeviation(mocapdata, Rhs, perturbtime, options);
    end
    
    if options.testing
        % allow user to inspect the result, before proceeding
        disp('Please check Figure 1 for correct detection of the perturbation response');
        disp('Hit ENTER to continue');
        pause
    end
    
    % We currently don't do anything with the rest of the analysis, but the plots may
    % be useful to illustrate that the marker-based T2 is much better for seeing perturbation
    % responses than the single variables.
    % To turn on the code, use "if (1)" in the following line.
    if (0) 
        name = mocapdata.latexname;
        
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
        avcycles(mocaptimes, Rkneeangle, Rhs, PRhs, perturbtime, [name 'RIGHT KNEE ANGLE']);
        avcycles(mocaptimes, Lkneeangle, Rhs, PRhs, perturbtime, [name 'LEFT KNEE ANGLE']);
        avcycles(mocaptimes, Rhipangle, Rhs, PRhs, perturbtime, [name 'RIGHT HIP ANGLE']);
        avcycles(mocaptimes, Lhipangle, Rhs, PRhs, perturbtime, [name 'LEFT HIP ANGLE']);
        avcycles(mocaptimes, -Rheelpos(:,3), Rhs, PRhs, perturbtime, [name 'R anterior foot pos']);
        avcycles(mocaptimes, -Lheelpos(:,3), Rhs, PRhs, perturbtime, [name 'L anterior foot pos']);
        avcycles(mocaptimes, Rheelpos(:,1), Rhs, PRhs, perturbtime, [name 'R lateral foot pos']);
        avcycles(mocaptimes, -Lheelpos(:,1), Rhs, PRhs, perturbtime, [name 'L lateral foot pos']);

        % also take a look at sacrum XYZ trajectories
        % (we can also use Center of Mass, as defined in the MOS analysis)
        Sacrum = getcolumn(mocapdata, 'SACR');
        avcycles(mocaptimes, Sacrum(:,1), Rhs, PRhs, perturbtime, [name 'SACRUM X(right)']);
        avcycles(mocaptimes, Sacrum(:,2), Rhs, PRhs, perturbtime, [name 'SACRUM Y(up)']);
        avcycles(mocaptimes, Sacrum(:,3), Rhs, PRhs, perturbtime, [name 'SACRUM Z(posterior)']);
    end
    disp('Analysis of perturbation response (response.m) completed.')  
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
    Hip   = getcolumn(data, [side 'GTRO']);
    Knee  = getcolumn(data, [side 'LEK']);
    Ankle = getcolumn(data, [side 'LM']);
    
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
    Heel   = getcolumn(data, [side 'HEE']);
    Sacrum  = getcolumn(data, 'SACR');
    pos = Heel - Sacrum;   
end
%=======================================================
function [angles] = hipangles(data, side)
    Hip   = getcolumn(data, [side 'GTRO']);
    Knee  = getcolumn(data, [side 'LEK']);
    Shoulder = getcolumn(data, [side 'SHO']);
    
    % construct unit vectors from shoulder to hip and hip to knee
    trunk = Hip - Shoulder;
    thigh = Knee - Hip;
    crossnorm = vecnorm(cross(trunk,thigh),2,2);  % norm of the cross product
    thighnorm = vecnorm(thigh,2,2);               % norm of the thigh vector
    trunknorm = vecnorm(trunk,2,2);               % norm of the shank vector
    
    % angle is inverse sine of (cross product divided by both vector norms)
    angles = asin(crossnorm ./ thighnorm ./ trunknorm);
    
end
%=========================================================================================
function [result, newdata] = gaitdeviation(data, hs, perturbtime, options)   
% calculate gait deviation timing (t1,t2) using Mahalanobis Distance
% and also estimate missing data
    
    newdata = data;
    nsamples = size(data.data,1);
    t = data.data(:,1) - data.data(1,1);
    tperturb = perturbtime - data.data(1,1); % perturbation time relative to start of trial

    % find the columns with marker data
    if options.markerset17
        % use data from 17 defined markers
        markers = {'RASIS','LASIS','RPSIS','LPSIS', 'SACR', ...
                   'RGTRO','RLEK','RLM','RHEE','RMT5','RTOE', ...
                   'LGTRO','LLEK','LLM','LHEE','LMT5','LTOE'};
    else
        % use all marker data, which is recognized by '.Pos' in the label
        markers = {'.Pos'};
    end
    columns = find(contains(data.colheaders,markers))';
    
    % remove those with too many NaNs
    remove = [];
    for i=1:numel(columns)
        perc_missing = 100*sum(isnan(data.data(:,columns(i)))) / nsamples;
        if perc_missing > 50
            remove = [remove i];  % add i to the remove list            
            % generate a message once for each removed marker (but not S1,
            % it is not an actual marker)
            column_name = data.colheaders{columns(i)};
            if contains(column_name,'.PosX') && ~contains(column_name,'S1.Pos')
                marker_name = strrep(column_name,'.PosX','');
                fprintf('gait deviation analysis will ignore %s, it has %.1f%% missing.\n', marker_name, perc_missing);
            end
        end
    end
    columns = setdiff(columns,columns(remove));  % these columns are kept
    ncolumns = numel(columns);

    % estimate the covariance matrix using frames from gait cycles before the
    % perturbation, and those after the perturbation
    i1 = find(t(hs) < tperturb,    1, 'last');  % last heelstrike before perturbation
    i2 = find(t(hs) > tperturb+10, 1, 'first'); % first heelstrike ten seconds after perturbation
    normalframes = [hs(1):hs(i1) hs(i2):hs(end)];
    [mu,C] = ecmnmle_hash(data.data(normalframes,columns));
    Cinv = inv(C);
   
    % impute missing data using idea from Rasmussen 2020
    % for each frame, we reconstruct the missing markers by finding
    % coordinates that put this frame closest (using Mahalanobis Distance, MD)
    % to the mean (mu).
    % this is an optimization problem with an analytical solution
    % partition the inverse covariance matrix into parts corresponding to
    % observed coordinates (o) and missing coordinates (m)
    % Cinv = [Coo Com]  this is a symmetrical matrix, so Com = Cmo'
    %        [Cmo Cmm]
    % we minimize MD=(x-mu)'*Cinv*(x-mu), subject to x(o) = xo (observed data)
    % the solution is: x(m) = mu(m) - inv(Cmm)*Cmo*(xo-mu(o))
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

    % Calculate the Mahalanobis distance of each frame to the mean
    % We could also do this from the original incomplete data (see
    % notebook 3/23/2021) but result is exactly the same.
    % We could also do this after subtracting the average gait cycle, but
    % again, result is (almost) exactly the same.
    drel = newdata.data(:,columns)-repmat(mu',nsamples,1);
    T2 = dot(drel * Cinv, drel, 2);
            
    % low-pass filter to smooth the T2 curve
    Fs = 1000;
    Fc = options.filter;
    tnew = (t(hs(1)):(1/Fs):t(hs(end)))';  % resample to ensure constant sampling rate
    T2new = interp1(t,T2,tnew);
    [b,a] = butter(2,Fc/(Fs/2));
    T2f = filtfilt(b,a,T2new);
    
    % calculate p values from the T2, which has
    % a Chi-squared distribution with ncolumns degrees of freedom
    % See https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
    % we could calculate p directly from marker data: p = mvncdf(data, mu, C)
    % but Chi-squared is easier to explain, and we already have the T2f anyway
    p = chi2cdf(T2f, ncolumns, 'upper');
    
    % gait deviation measure is the duration of the period when p<threshold
    threshold = options.pvalue;
    d = [0 ; diff(p<threshold)];  % diff detects when p went across the threshold
    it1 = find(d==1  & tnew>(tperturb-1),1,'first');  % first time after the perturbation that p went below threshold
    it2 = find(d==-1 & tnew>tperturb,1,'first');  % first time after the perturbation that p went above threshold again
    rt1 = tnew(it1);
    rt2 = tnew(it2);
    t1 = rt1 - tperturb;
    t2 = rt2 - tperturb;
    T2max = max(T2f(it1:it2));
    pmin  = min(p(it1:it2));
    
    % create the result table
    if isempty(t1), t1 = 9999, end
    if isempty(t2), t2 = 9999, end
    if isempty(T2max), T2max = 9999, end
    result = table(t1,t2,T2max);
    
    % write results on screen, make a figure to illustrate what was calculated
    % (only when testing)
    if (options.testing)
        fprintf('Start (t1) of perturbation response:    %.3f s. after perturbation\n',t1);
        fprintf('End (t2) of perturbation response:      %.3f s. after perturbation\n',t2);
        fprintf('Magnitude (T2max) of response:          %.3f\n',T2max);
        fprintf('Minimum p value (pmin) during response: %.3e\n',pmin);
        xlim = [tperturb-2 tperturb+10]; % the time range we want to plot
        close all  % remove figures from screen

        figure(1);
        subplot(2,1,1);
        plot(tnew,[T2new T2f]);
        title([data.latexname ': Hotelling T^2']);
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
        pos = get(gca,'Position');  % we can use this to calculate coordinates relative to Figure window
        xrel1 = pos(1) + pos(3)*(tperturb+[0 t1]-xlim(1))/diff(xlim);
        xrel2 = pos(1) + pos(3)*(tperturb+[0 t2]-xlim(1))/diff(xlim);
        yrel1 = pos(2) + pos(4)*( (log(threshold)-log(ylim(1)) )/diff(log(ylim)) - 0.1);
        yrel2 = pos(2) + pos(4)*( (log(threshold)-log(ylim(1)) )/diff(log(ylim)) - 0.3);
        if ~isnan(t1)
            annotation('textarrow',xrel1,[yrel1 yrel1],'String', num2str(t1,'t1 = %.3f s.'));
        end
        annotation('textarrow',xrel2,[yrel2 yrel2],'String', num2str(t2,'t2 = %.3f s.'));
        
        % append to the tmp.ps file
        % print('-dpsc','-append','-bestfit','tmp.ps');
    end
end
%========================================================
function [mu,C] = ecmnmle_hash(data)
% estimate the mean and covariance of a data matrix
% use a hash method to see if we already have done this before
    % Convert data into a byte array called B...
    B = typecast(data(:),'uint8');
    if (isempty(B))
        error('ecmnmle_hash: B is empty, you may not have any marker data');
    end
    % Or, as suggested, using the undocumented function getByteStreamFromArray:
    % B = getByteStreamFromArray(data);

    % Create an instance of a Java MessageDigest with the desired algorithm:
    md = java.security.MessageDigest.getInstance('SHA-1');
    md.update(B);

    % Properly format the computed hash as an hexadecimal string:
    hash = reshape(dec2hex(typecast(md.digest(),'uint8'))',1,[]);
    
    % Construct the filename and see if the file exists
    filename = ['cache/' hash '.mat'];
    if exist(filename)
        % we have seen this data before, so simply reload the existing mean
        % and covariance
        disp('Loading mean and covariance data from cache...');
        load(filename);
    else
        % we have not seen this data before, so we need to estimate
        % the mean and covariance
        disp('Estimating mean and covariance of marker data...');
        disp('  (can take up to 6 min for full marker set, 30 seconds for 17 markers)');
        [mu,C] = ecmnmle(data);
        save(filename,'mu','C');  % save the result for later
    end
end
