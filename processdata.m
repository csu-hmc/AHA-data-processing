function processdata(folder, detail)
    % process data from a series of trials in one folder
    % Inputs:
    %   folder............Folder where trial files are stored
    %   detail............(optional) use 1 when you want to inspect all results
    %
    % Results are stored in files in the same folder:
    %   results.xlsx..........Discrete variables from each trial
    %   plots.pdf.............Plots to help with quality control
    
    % use the folder name to create the name of the Excel file
    folderpath = [getpath() folder '\'];
    excelfile = [folderpath folder '.xlsx'];
    tmppsfile = [folderpath 'tmp.ps'];
    
    % delete any files that may have been created previously
    if exist(excelfile), delete(excelfile); end
    if exist(tmppsfile), delete(tmppsfile); end
    
    % turn detail off if not specified
    if nargin < 2
        detail = 0;
    end

    % create an empty table to store the results from all trials
    % (this table will go into the Excel file)
    T = table;

    % load the list of files that must be processed
    filelist = importdata([folderpath '0FileList.txt']);
    for trial_num = 1:numel(filelist)
        name = strtrim(filelist{trial_num});
        name = [folder '\' name];
        fprintf('   ...%s\n', name);
        T.TrialNum(trial_num) = trial_num;
        T.FileName{trial_num} = name;
        [mocapdata,treadmilldata] = getdata(name);
        
        % do the analysis of the perturbation response and store result in table T
        options.testing = detail;
        options.markerset17 = 1;  % use the smaller marker set
        result = response(mocapdata, treadmilldata, options);
        T.T2max(trial_num) = result.T2max;
        T.t2(trial_num)    = result.t2;
        mocapdata = result.mocapdata;  % use the PCA-filled mocap data for the other measures
        
        % do the step analysis and store result in table T
        result = step_analysis(mocapdata, detail);
        T.ST_right_mean(trial_num) = result.ST_right_mean;
        T.ST_right_SD(trial_num)   = result.ST_right_SD;
        T.ST_left_mean(trial_num)  = result.ST_left_mean;
        T.ST_left_SD(trial_num)    = result.ST_left_SD;
        T.SL_right_mean(trial_num) = result.ST_right_mean;
        T.SL_right_SD(trial_num)   = result.ST_right_SD;
        T.SL_left_mean(trial_num)  = result.ST_left_mean;
        T.SL_left_SD(trial_num)    = result.SL_left_SD;

        % do the Margin of Stability analysis during normal walking (20
        % seconds)
        % ...
        
    end

    % write the table T on the Excel file
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

        saveas(fig2, strcat(folderpath, 'step_analysis.png'))

        disp('Check Figure 2 for problems.');
        disp('Hit ENTER to continue');
    end
end
