function processdata(folder, detail)
    % process data from a series of trials in one folder
    % Inputs:
    %   folder............Folder where trial files are stored
    %   detail............(optional) use 1 when you want to inspect all results
    %
    % Results are stored in files in the same folder:
    %   0<foldername>.xlsx..........Discrete variables from each trial
    %   0step_analysis.png..........Step variables plotted over the session
    
    % use the folder name to create the name of the Excel file
    folderpath = [getpath() folder '\'];
    excelfile = [folderpath '0' folder '.xlsx'];
    tmppsfile = [folderpath 'tmp.ps'];
    
    % delete any files that may have been created previously
    if exist(excelfile), delete(excelfile); end
    if exist(tmppsfile), delete(tmppsfile); end
    
    % turn detail off if not specified
    if nargin < 2
        detail = 0;
    end

    % load the list of files that must be processed
    filelist = importdata([folderpath '0FileList.txt']);
    for trial_num = 1:numel(filelist)
        name = strtrim(filelist{trial_num});
        name = [folder '\' name];
        fprintf('   ...%s\n', name);
        TrialNum(trial_num,1) = trial_num;
        FileName{trial_num,1} = name;
        [mocapdata,treadmilldata] = getdata(name);
        
        % do the analysis of the perturbation response and store result
        options.testing = detail;
        options.markerset17 = 1;  % use the smaller marker set
        result = response(mocapdata, treadmilldata, options);
        T2max(trial_num,1) = result.T2max;
        t2(trial_num,1)    = result.t2;
        mocapdata = result.mocapdata;  % use the PCA-filled mocap data for the other measures
        
        % do the step analysis and store result
        result = step_analysis(mocapdata, detail);
        ST_right_mean(trial_num,1) = result.ST_right_mean;
        ST_right_SD(trial_num,1)   = result.ST_right_SD;
        ST_left_mean(trial_num,1)  = result.ST_left_mean;
        ST_left_SD(trial_num,1)    = result.ST_left_SD;
        SL_right_mean(trial_num,1) = result.SL_right_mean;
        SL_right_SD(trial_num,1)   = result.SL_right_SD;
        SL_left_mean(trial_num,1)  = result.SL_left_mean;
        SL_left_SD(trial_num,1)    = result.SL_left_SD;
        
        % do the Margin of Stability analysis during normal walking (20
        % seconds)
        % ...
        
    end

    % make a table with all results and write it on the Excel file
    T = table(TrialNum, FileName, ...
        T2max, t2, ...         
        ST_left_mean,ST_left_SD,ST_right_mean,ST_right_SD, ...
        SL_left_mean,SL_left_SD,SL_right_mean,SL_right_SD);     
    writetable(T, excelfile);

    %% plot results of the step analysis, showing trend over the trials
    if detail
        fig2 = figure(2);
        subplot(2,2,1)
        errorbar(T.SL_left_mean, T.SL_left_SD)
        ylabel('Left Step Length (m)')
        xlabel('Trial')
        ylim([0 0.5])
        xlim([0, 16])

        subplot(2,2,2)
        errorbar(T.SL_right_mean, T.SL_right_SD)
        ylabel('Right Step Length (m)')
        xlabel('Trial')
        ylim([0 0.5])
        xlim([0, 16])

        subplot(2, 2, 3)
        errorbar(T.ST_left_mean, T.ST_left_SD)
        ylabel('Left Step Time (s)')
        xlabel('Trial')
        ylim([0.4 1.4])
        xlim([0, 16])

        subplot(2, 2, 4)
        errorbar(T.ST_right_mean, T.ST_right_SD)
        ylabel('Right Step Time (s)')
        xlabel('Trial')
        ylim([0.4 1.4])
        xlim([0, 16])

        saveas(fig2, strcat(folderpath, '0step_analysis.png'))

        disp('Check Figure 2 for problems.');
        disp('Hit ENTER to continue');
    end
end
