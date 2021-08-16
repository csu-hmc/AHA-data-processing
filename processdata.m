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
    if exist(excelfile) 
        error('Could not delete %s. It may be open in another application. Close it and retry.', excelfile);
    end
    if exist(tmppsfile) 
        error('Could not delete %s. It may be open in another application. Close it and retry.', tmppsfile);
    end
    
    % turn detail off if not specified
    if nargin < 2
        detail = 0;
    end

    % load the list of files that must be processed
    filelist = importdata([folderpath '0FileList.txt']);
    
    % process all files in the list, and create a table with all results
    T = table;  % start with an empty table
    for TrialNum = 1:numel(filelist)
        FileName = strtrim(filelist{TrialNum});  % make sure file name has no leading or trailing spaces
        FileName = [folder '\' FileName];
        fprintf('   ...%s\n', FileName);
        [mocapdata,treadmilldata] = getdata(FileName, detail);
        
        % do the analysis of the perturbation response and store result
        options.testing = detail;
        options.markerset17 = 1;  % use the smaller marker set
        [result1, mocapdata] = response(mocapdata, treadmilldata, options);  % mocapdata now has PCA-filed data!
        
        % do the step analysis and store result
        result2 = step_analysis(mocapdata, detail);
        
        % do the Margin of Stability analysis during normal walking (20
        % seconds)
        result3 = mos(mocapdata, detail);
        
        % add results from this trial to the table T
        T = [T ; table(TrialNum) table({FileName}) result1 result2 result3];
        
    end

    % convert the table to an Excel file   
    writetable(T, excelfile);

    %% plot results of the step analysis and MOS, showing trend over the trials
    fig2 = figure(2);
    latexfolder = strrep(folder,'\','/');
    latexfolder = strrep(latexfolder,'_','\_');  % make foldername suitable for plot title
    trials = T.TrialNum;

    subplot(2,2,1)
    errorbar([trials trials], ...
        [T.SL_left_mean T.SL_right_mean], [T.SL_left_SD T.SL_right_SD])
    set(gca,'XLim',[min(trials)-0.5 max(trials)+0.5]);
    ylabel('Step Length (m)')
    xlabel('Trial')
    title(latexfolder)

    subplot(2,2,2)
    errorbar([trials trials], ...
        [T.ST_left_mean T.ST_right_mean], [T.ST_left_SD T.ST_right_SD])
    set(gca,'XLim',[min(trials)-0.5 max(trials)+0.5]);
    ylabel('Step Time (s)')
    xlabel('Trial')
    title(latexfolder)

    subplot(2,2,3)
    errorbar([trials trials], ...
        [T.MOS_AP_left_mean T.MOS_AP_right_mean], [T.MOS_AP_left_SD T.MOS_AP_right_SD])
    set(gca,'XLim',[min(trials)-0.5 max(trials)+0.5]);
    ylabel('MOS AP (m)')
    xlabel('Trial')
    title(latexfolder)

    subplot(2,2,4)
    errorbar([trials trials], ...
        [T.MOS_ML_left_mean T.MOS_ML_right_mean], [T.MOS_ML_left_SD T.MOS_ML_right_SD])
    set(gca,'XLim',[min(trials)-0.5 max(trials)+0.5]);
    ylabel('MOS ML (m)')
    xlabel('Trial')
    title(latexfolder)

    legend('left','right')
    saveas(fig2, strcat(folderpath, '0step_analysis.png'))
    
    if (detail)
        disp('Check Figure 2 for problems.');
        disp('Hit ENTER to continue');
    end
end
