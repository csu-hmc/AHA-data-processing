function result = c3dtotxt(c3d_filename, txt_filename)
% Marker data is extracted from a C3D file, synchronized with a previously
% recorded mocap file from D-Flow.  A new mocap file is generated
% containing the original force plate data recorded in D-Flow, and the marker
% data from the C3D file.
% The name of the new mocap file is <original mocap filename>_edited.txt
%
% Output is a struct with the following fields:
%   result.nmissing_before  number of missing marker coordinates in original mocap fole
%   result.nmissing_after:   number of missing marker coordinates in new mocap fole
%   result.info:
%       0: success
%       1: number of force samples on c3d file was not 10x number of marker samples 
%       2: could not write the _edited.txt file
%
% The C3D file reading is done using the ezc3d toolbox (if installed) or with the
% Opensim C3D tools: https://simtk-confluence.stanford.edu:8443/display/OpenSim/C3D+(.c3d)+Files

    testing = (nargin == 0);    % if no file names are specified, we're testing
    if (testing)
        % here we define where the test file is
        % edit this code to match your computer
        computer = getenv('COMPUTERNAME');
        if strcmp(computer, 'LRI-102855')   % Ton's computer
            datapath = 'C:\Users\Ton\Cleveland State University\Hala E Osman - Hala data\';
        elseif strcmp(computer, 'DESKTOP-0HN0T6U')   % Hala's computer
            datapath = 'C:\Users\hallo\OneDrive - Cleveland State University\Hala data\';
        else
            fprintf('Your computer name is: %s\n', computer);
            fprintf('Please configure c3dtotxt.m for your computer.\n');
        end
        trial = 'Par12_POST\Mocap0001';
        txt_filename = [datapath trial '.txt'];
        c3d_filename = [datapath trial '.c3d'];
    end

    % open the log file with 'a' (append)
    logfile = 'c3d.log';
    log = fopen(logfile,'a');
    if (log < 0)
        error('cannot write %s', logfile);
    end
    fprintf(log,'------------------------------------------------------\n');
    fprintf(log,'Processing %s...\n', c3d_filename);
    
    % read the C3D file with the ezc3d toolbox, if it is installed
    tic;
    if which('ezc3dRead')
        fprintf('Reading %s...\n   (will take 1-2 seconds for a 90-second trial)...\n', c3d_filename);
        c3d = ezc3dRead(c3d_filename);
        c3dMarkerNames = c3d.parameters.POINT.LABELS.DATA;
        c3dMarkerData = 0.001 * c3d.data.points;    % convert from mm to m
        c3dFy = -c3d.data.analogs(:,[3 9]);  % vertical force in force plate 1 and 2
    else
        fprintf('ezc3d toolbox is not installed, continuing with OpenSim...\n');
        import org.opensim.modeling.*
        fprintf('Reading %s...\n   (will take 5-10 minutes for a 90-second trial)...\n', c3d_filename);
        c3d = osimC3D(c3d_filename,0);
        [markerStruct, forceStruct] = c3d.getAsStructs();
        c3dMarkerNames = fieldnames(markerStruct);
        c3dMarkerNames = c3dMarkerNames(1:end-1);  % remove the last one, which is time
        for i = 1:numel(c3dMarkerNames)
            c3dMarkerdata(:,i,:) = 0.001 * markerStruct.(c3dMarkerNames{i})';
        end
        % extract the two vertical forces
        c3dFy = [forceStruct.f1(:,2) forceStruct.f2(:,2)];
    end
    fprintf('...it took %.1f seconds\n',toc);

    c3dNframes = size(c3dMarkerData,3);
    c3dNmarkers = numel(c3dMarkerNames);  % the last "marker" in the c3d is actually time
    fprintf('C3D file has %d frames and %d markers\n', c3dNframes, c3dNmarkers);

    % downsample c3d Fy1 data from 1000 Hz to 100 Hz by averaging 10 sequential force
    % samples into 1
    if size(c3dFy,1) ~= 10*c3dNframes
        error('c3dtotxt.m: number of force samples in c3d is not 10x number of frames');
    end
    c3dFyr = zeros(c3dNframes,2); % to store resampled vertical force from C3D data
    for i = 1:c3dNframes
        j = 10*(i-1) + (1:10);  % the 10 samples that must be averaged
        c3dFyr(i,:)   = mean(c3dFy(j,:));  % take the average of the 10 samples
    end

    % import the original TXT file and generate more reliable time stamps;
    % we know that Cortex samples exactly every 10 ms
    data = importdata(txt_filename);
    txtNframes = size(data.data,1);
    fprintf('TXT file has %d frames\n', txtNframes);
    if txtNframes > 20000
        % this happened in Par12_POST: 300 Hz time stamps, more than 20000 frames
        % use only the first 9000 frames
        txtNframes = 9000;
        data.data = data.data(1:txtNframes,:);
        fprintf('TXT file truncated to %d frames\n', txtNframes);
    end
    data.data(:,1) = data.data(1,1) + (0:(txtNframes-1))*0.010;
%     else
%         newtimestamps = (data.data(1,1) : 0.01 : data.data(end,1))';
%         newframenumbers = (1:numel(newtimestamps))';
%         newdata = interp1(data.data(:,1), data.data(:,3:end), newtimestamps);
%         data.data = [ newtimestamps newframenumbers newdata ];
%         txtNframes = size(data.data,1);
%         fprintf('TXT file resampled to %d frames\n', txtNframes);
%     end

    % remove any spaces from the column headers and collect txtMarkerNames
    txtMarkerNames = {};
    for i = 1:numel(data.colheaders)
        data.colheaders{i} = strrep(data.colheaders{i},' ','');
        if contains(data.colheaders{i},'.PosX')
            txtMarkerNames = [ txtMarkerNames strrep(data.colheaders{i},'.PosX','') ];
        end
    end
    
    % If there is not a full markerset in the TXT file,
    % we add the standard 47-marker set and missing data for all markers and 
    % all frames. This allows processing to continue, and C3D marker data will be used
    % We also do this when the TXT file has unlabeled markers, 47 names all containing the letter M
    if numel(txtMarkerNames) < 47 || sum(contains(txtMarkerNames,'M')) == 47
        txtMarkerNames = {'LHEAD' 'THEAD' 'RHEAD' 'FHEAD' 'C7'   'T10'  'SACR' 'NAVE'  'XYPH'  'STRN'  ...
                          'BBAC'  'LSHO'  'LDELT' 'LLEE'  'LMEE' 'LFRM' 'LMW'  'LLW'   'LFIN'  'RSHO'  ...
                          'RDELT' 'RLEE'  'RMEE'  'RFRM'  'RMW'  'RLW'  'RFIN' 'LASIS' 'RASIS' 'LPSIS' ...
                          'RPSIS' 'LGTRO' 'FLTHI' 'LLEK'  'LATI' 'LLM'  'LHEE' 'LTOE'  'LMT5'  'RGTRO' ...
                          'FRTHI' 'RLEK'  'RATI'  'RLM'   'RHEE' 'RTOE' 'RMT5'};
        for i = 1:numel(txtMarkerNames)
            markername = txtMarkerNames{i};
            columnnames = { [markername '.PosX'] [markername '.PosY'] [markername '.PosZ'] };
            data.colheaders = [data.colheaders columnnames];
        end
        data.data = [ data.data zeros(txtNframes, 3*numel(txtMarkerNames)) ];
    end
             
    % C3D sometimes only has the first 4 characters of the marker name
    % so we try to match them.  
    % also remove 'FullBodyef:' which sometimes occurs in a C3D marker name
    % also remove 'FULLBODY47:' which sometimes occurs in a C3D marker name
    for i = 1:numel(c3dMarkerNames)
        c3dMarkerNames{i} = strrep(c3dMarkerNames{i},'FullBodyef:','');
        c3dMarkerNames{i} = strrep(c3dMarkerNames{i},'FULLBODY47:','');
        matches = find(strncmp(txtMarkerNames,c3dMarkerNames{i},4));
        if numel(matches) == 0
            error('Marker %s was not found in the original TXT file',c3dMarkerNames{i});
        elseif numel(matches) > 1
            error('C3D marker %s has multiple matches in the TXT file')
        else
            c3dMarkerNames(i) = txtMarkerNames(matches);
        end
    end
        
    % get the Fy1 signal from the original D-Flow TXT file
    iFy1 = strcmp(data.colheaders,'FP1.ForY');  % channel number for FP1.ForY
    Fy1 = data.data(:,iFy1);                          % FP1.ForY signal from the TXT file

    % find the time shift using peak of the cross correlation on Fy of force plate 1
    [c,lags] = xcorr(Fy1,c3dFyr(:,1),1000);  % max lag is 1000 frames (10 seconds)
    [~,imax] = max(c);
    lag = lags(imax);
    if (abs(lag) > 500)
        fprintf('WARNING: lag is very large: %.2f seconds\n', 0.01*lag);
        fprintf(log,'WARNING: lag is very large: %.2f seconds\n', 0.01*lag);
    end

    % check the agreement between Fy1 from TXT and C3D in the first 10 seconds
    if (lag<0)
        RMSdiff = rms(c3dFyr(-lag+(1:1000),1)-Fy1(1:1000));
    else
        RMSdiff = rms(c3dFyr(1:1000,1)-Fy1(lag+(1:1000)));
    end    
    maxdiff = 20.0;   % we allow a difference of 10 N
    if (RMSdiff > maxdiff)
        if (lag<0)
            plot([Fy1(1:1000) c3dFyr(-lag+(1:1000),1)]);
        else
            plot([Fy1(lag+(1:1000)) c3dFyr(1:1000,1)]);
        end 
        legend('Fy1 TXT','Fy1 C3D')
        fprintf('ERROR: Vertical force of forceplate 1 does not agree within %.1f N\n',maxdiff);
        fprintf(log,'ERROR: Based on force plate data, it seems that TXT and C3D are not the same trial');
        disp('type dbcont if you want to continue anyway');
        keyboard
    end
    
    % set up the copying of marker data from c3d to txt
    if (lag >= 0)
        fprintf('C3D file started %.3f s later than the D-Flow TXT file\n', 0.01*lag);
        fprintf('   the new TXT file will have the first %d frames of the original TXT file\n', lag);
        fprintf(log,'C3D file started %.3f s later than the D-Flow TXT file\n', 0.01*lag);
        fprintf(log,'   the new TXT file will have the first %d frames of the original TXT file\n', lag);
        txt_start = lag;   % start here when copying C3D data into the new txt file
        c3d_start = 0;       % start here when copying the c3d data
    else
        fprintf('C3D file started %.3f ms earlier than the TXT file\n', -0.01*lag);
        fprintf(log,'C3D file started %.3f ms earlier than the TXT file\n', -0.01*lag);
        txt_start = 0;     
        c3d_start = -lag; 
    end

    mis = txtNframes - c3dNframes - lag;  % number of C3D frames missing at the end
    if (mis >= 0)
        fprintf('C3D file stopped %.3f s earlier than the D-Flow TXT file\n', 0.01*mis);
        fprintf('   the new TXT file will have the last %d frames of the original TXT file\n', mis);
        fprintf(log,'C3D file stopped %.3f s earlier than the D-Flow TXT file\n', 0.01*mis);
        fprintf(log,'   the new TXT file will have the last %d frames of the original TXT file\n', mis);
        ncopy = c3dNframes - c3d_start;  % how many frames to copy from C3D to TXT
    else
        fprintf('D-Flow TXT file stopped %.3f s earlier than the C3D file\n', -0.01*mis);
        fprintf('   WARNING: last %d frames of new TXT file will have incomplete force data (vertical only)\n', -mis);
        fprintf(log,'D-Flow TXT file stopped %.3f s earlier than the C3D file\n', -0.01*mis);
        fprintf(log,'   WARNING: last %d frames of new TXT file will have incomplete force data (vertical only)\n', -mis);
        ncopy = c3dNframes - c3d_start;
    end

    % find the relevant columns in the original TXT data
    colnames = {};
    for i=1:c3dNmarkers
        colnames = [ colnames [c3dMarkerNames{i} '.PosX'] [c3dMarkerNames{i} '.PosY'] [c3dMarkerNames{i} '.PosZ'] ];
    end
    fpmnames = {'FP1.ForX','FP1.ForY','FP1.ForZ','FP1.CopX','FP1.CopY','FP1.CopZ','FP1.MomX','FP1.MomY','FP1.MomZ', ...
                'FP2.ForX','FP2.ForY','FP2.ForZ','FP2.CopX','FP2.CopY','FP2.CopZ','FP2.MomX','FP2.MomY','FP2.MomZ'};
    colnames = [colnames fpmnames];  % add the force plate data column names
    txtcol = [];
    for i = 1:numel(colnames)
        if isempty(find(strcmp(data.colheaders,colnames{i})))
            fprintf('Column %s was not found in the original TXT file; C3D data will be used',colnames{i});
        end
        txtcol = [txtcol find(strcmp(data.colheaders,colnames{i}))];
    end

    % first copy the relevant columns from original TXT data
    newdata = data.data(:,txtcol);
    
    % if the C3D file continued past the TXT file, use more data from C3D
    if (mis < 0)
        newdata = [newdata ; zeros(-mis,numel(txtcol))];
        % only copy vertical force (others are complicated)
        newdata(txtNframes+(1:(-mis)),3*c3dNmarkers+[2 11]) = c3dFyr(txtNframes-txt_start+(1:(-mis)),:);
        txtNframes = size(newdata,1);
    end

    % insert ncopy frames from C3D markerStruc, markers only; also convert mm to m
    for i = 1:c3dNmarkers
        columns = 3*(i-1) + (1:3);
        c3dData = squeeze(c3dMarkerData(:,i,c3d_start+(1:ncopy)))';  % extract ncopy x 3 marker coordinates
        missingc3d = find(isnan(c3dData(:,1)));
        c3dData(missingc3d,:) = 0.0;  % use 0.0 as missing marker indicator (D-Flow convention) 
        missingtxt = find(newdata(:,columns(1))==0);
        if numel(missingc3d) > numel(missingtxt)
            fprintf('WARNING: C3D has more missing data than TXT for %s; C3D data not used.\n', c3dMarkerNames{i});
            fprintf(log,'WARNING: C3D has more missing data than TXT for %s; C3D data not used.\n', c3dMarkerNames{i});
        else
            newdata(txt_start+(1:ncopy),columns) = c3dData;
        end
    end

    % convert NaN to zero (indicates missing marker data)
    markercol = (1:3*c3dNmarkers);
    for i = markercol
        k = isnan(newdata(:,i));
        newdata(k,i) = 0.0;             % put zeros where the NaN were
    end

    % find out how many times a marker was missing in original TXT and new TXT
    result.missing_before = 0;
    result.missing_after = 0;
    for i = 1:c3dNmarkers
        col = find(strcmp(data.colheaders, [c3dMarkerNames{i} '.PosX']));
        if isempty(col)     % this marker was not found in the original TXT file
            nmissing = nFramesTxt;   % all frames are missing for this marker
        else
            nmissing = numel(find(~data.data(:,col)));  % find how many zeros in column
        end
        nmissingnew = numel(find(~newdata(:,3*(i-1)+1)));  % find how many zeros in x of this marker
        fprintf(log,'%6s: %5d missing in TXT, %5d missing in edited TXT (out of %5d)\n', ...
            c3dMarkerNames{i},nmissing,nmissingnew,txtNframes);
        % fprintf('%s in TXT has %d missing data\n',data.colheaders{col},nmissing);
        result.missing_before = result.missing_before + nmissing;
        result.missing_after = result.missing_after + nmissingnew;
    end

    % Write the marker and force data to the new .txt file
    % marker data will come from the C3D file
    % time stamp, frame number, and force data will be copied from the original TXT file
    filename = strrep(txt_filename, '.txt', '_edited.txt');  % change .c3d to .txt in the file name
    fprintf('Writing %s...\n', filename);
    fprintf(log,'Writing %s...\n', filename);
    fid = fopen(filename,'w');
    if (fid < 0)
        error('cannot write %s', filename);
        result.info = 2;
        return
    end

    % write the header line
    fprintf(fid,'TimeStamp\tFrameNumber');
    for i = 1:numel(colnames)
        fprintf(fid,'\t%s',colnames{i});
    end
    fprintf(fid,'\n');  % end of line

    % write the data
    for iframe = 1:txtNframes
        fprintf(fid,'%f\t%d',data.data(1,1)+0.01*(iframe-1), data.data(1,2)+(iframe-1));    % write time stamp and frame number
        fprintf(fid,'\t%f', newdata(iframe,:));  % write new data
        fprintf(fid,'\n');  % end of line
    end

    fclose(fid);
    fprintf(log,'Done.\n');
    fclose(log);

    % info = 0 to indicate success
    result.info = 0;
    
    % if we are not testing, we are done and we exit this function
    if ~testing
        return
    end
    
    % if testing, compare the original file to the new file
    data1 = importdata(txt_filename);  % the original data
    newtxt_filename = strrep(txt_filename, '.txt', '_edited.txt');  % this is the name of the new .txt file
    data2 = importdata(newtxt_filename);  % the data converted from c3d

    % remove spaces from column headers in the original file
    for i = 1:numel(data1.colheaders)
        data1.colheaders{i} = strrep(data1.colheaders{i},' ','');
    end

    % compare each channel in the new file to the corresponding channel in the original file
    close all
    t1 = data1.data(:,1);  % timestamp is in column 1
    t2 = data2.data(:,1);
    t1 = t1(1) + 0.01 * (0:numel(t1)-1);  % make new time stamps for data1 (original TXT)
    for col2 = 1:numel(data2.colheaders)
        varname = data2.colheaders{col2};  % name of variable i in data 2
        col1 = find(strcmp(data1.colheaders, varname));  % column number in data1
        if isempty(col1)
            fprintf('%s is not available in original TXT file\n', varname);
        else
            d1 = data1.data(:,col1);
            d2 = data2.data(:,col2);

            % if variable name includes 'Pos', replace zeros by NaN so they are not plotted
            if findstr(varname,'Pos')
                d1(~d1) = nan;  % replace zeros (missing markers) by NaN, so plot does not show those
                d2(~d2) = nan;
            end

            plot(t2,d2,t1,d1);
            xlabel('timestamp');
            legend('TXT converted from C3D','original TXT file');
            title(varname);
            disp('Hit ENTER to continue or CTRL-C to stop');
            pause
        end
    end

end