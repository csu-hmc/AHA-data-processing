function result = c3dtotxt(c3d_filename, txt_filename)
% Marker data is extracted from a C3D file, synchronized with a previously
% recorded mocap file from D-Flow.  A new mocap file is generated
% containing the original force plate data, and the marker data from the
% C3D file.
% The name of the new mocap file is <original mocap filename>_edited.txt
% result.nmissing_before  number of missing marker coordinates in original mocap fole
% result.nmissing_after:   number of missing marker coordinates in new mocap fole
% result.info:
%   0: success
%   1: number of force samples on c3d file was not 10x number of marker samples 
%   2: could not write the _edited.txt file
%
% The C3D file reading is done using the Opensim C3D tools: https://simtk-confluence.stanford.edu:8443/display/OpenSim/C3D+(.c3d)+Files

    % open the log file with 'a' (append)
    logfile = 'c3d.log';
    log = fopen(logfile,'a');
    if (log < 0)
        error('cannot write %s', logfile);
    end
    fprintf(log,'------------------------------------------------------\n');
    fprintf(log,'Processing %s...\n', c3d_filename);
    
    % read the C3D file with the ezc3d toolbox, if it is installed
    if which('ezc3dRead')
        c3d = ezc3dRead(c3d_filename);
        c3dMarkerNames = c3d.parameters.POINT.LABELS.DATA;
        c3dMarkerData = 0.001 * c3d.data.points;    % convert from mm to m
        c3dFy1 = -c3d.data.analogs(:,3);  % vertical force in force plate 1
    else

        % During testing, we sometimes read the C3D data (as structures) from tmp.mat because
        % extracting data from the C3D file is very slow
        if exist('tmp.mat')
            load('tmp.mat');
        else
            % load OpenSim libs
            import org.opensim.modeling.*

            % Construct an opensimC3D object with input c3d path
            % Constructor takes full path to c3d file and an integer for forceplate representation (1 = COP).
            fprintf('Opening %s...\n', c3d_filename);
            c3d = osimC3D(c3d_filename,0);

            % Get the c3d data as Matlab Structures
            fprintf('Extracting data... (will take approx. 5 minutes for a 90-second trial)\n');
            tic
            [markerStruct, forceStruct] = c3d.getAsStructs();
            fprintf('...it took %.1f minutes\n',toc/60);
            % save('tmp.mat','markerStruct','forceStruct');  % use this only for testing
            c3dMarkerNames = fieldnames(markerStruct);
            c3dMarkerNames = c3dMarkerNames(1:end-1);  % remove the last one, which is time
            c3dNframes = size(markerStruct.time, 1);
            for i = 1:numel(c3dMarkerNames)
                c3dMarkerdata(:,i,:) = 0.001 * markerStruct.(markernames{i})';
            end
            c3dFy1 = forceStruct.f1(:,2);
        end
    end

    c3dNframes = size(c3dMarkerData,3);
    c3dNmarkers = numel(c3dMarkerNames);  % the last "marker" in the c3d is actually time
    fprintf('C3D file has %d frames and %d markers\n', c3dNframes, c3dNmarkers);

    % downsample c3d Fy1 data from 1000 Hz to 100 Hz by averaging 10 sequential force
    % samples into 1
    if numel(c3dFy1) ~= 10*c3dNframes
        error('c3dtotxt.m: number of force samples in c3d is not 10x number of frames');
    end
    Fy1c3dr = zeros(c3dNframes,1); % to store vertical force from forceplate 1 in C3D data
    for i = 1:c3dNframes
        j = 10*(i-1) + (1:10);  % the 10 samples that must be averaged
        Fy1c3dr(i)   = mean(c3dFy1(j));  % take the average of the 10 samples
    end

    % import the original TXT file and generate more reliable time stamps
    data = importdata(txt_filename);
    txtNframes = size(data.data,1);
    fprintf('TXT file has %d frames\n', txtNframes);
    data.data(:,1) = data.data(1,1) + (0:(txtNframes-1))*0.010;

    % remove any spaces from the column headers and collect txtMarkerNames
    txtMarkerNames = {};
    for i = 1:numel(data.colheaders)
        data.colheaders{i} = strrep(data.colheaders{i},' ','');
        if contains(data.colheaders{i},'.PosX')
            txtMarkerNames = [ txtMarkerNames strrep(data.colheaders{i},'.PosX','') ];
        end
    end

    % ezc3d only has the first 4 characters of the marker name
    % copy the marker names from the TXT file, and flag error if
    % the matching is ambigious
    % also remove 'FullBodyef:' when it occurs in the C3D marker name
    for i = 1:numel(c3dMarkerNames)
        c3dMarkerNames{i} = strrep(c3dMarkerNames{i},'FullBodyef:','');
        matches = find(strncmp(txtMarkerNames,c3dMarkerNames{i},4));
        if numel(matches) == 0
            error('Marker %s was not found in the original TXT file',c3dMarkerNames{i});
        elseif numel(matches) > 1
            error('C3D marker %s has multiple matches in the TXT file')'
        else
            c3dMarkerNames(i) = txtMarkerNames(matches);
        end
    end
        
    % get the Fy1 signal from the original D-Flow TXT file
    iFy1 = strcmp(data.colheaders,'FP1.ForY');  % channel number for FP1.ForY
    Fy1 = data.data(:,iFy1);                          % FP1.ForY signal from the TXT file

    % find the time shift using peak of the cross correlation on Fy of force
    % plate 1
    [c,lags] = xcorr(Fy1,Fy1c3dr,1000);  % max lag is 1000 frames (10 seconds)
    [~,imax] = max(c);
    lag = lags(imax);
    if (lag > 500)
        fprintf('WARNING: lag is larger than 5 seconds\n');
        fprintf(log,'WARNING: lag is larger than 5 seconds\n');
    end

    % check the agreement between Fy1 from TXT and C3D in the first 10 seconds
    if (lag<0)
        RMSdiff = rms(Fy1c3dr(-lag+(1:1000))-Fy1(1:1000));
    else
        RMSdiff = rms(Fy1c3dr(1:1000)-Fy1(lag+(1:1000)));
    end    
    maxdiff = 10.0;   % we allow a difference of 10 N
    if (RMSdiff > maxdiff)
        if (lag<0)
            plot([Fy1(1:1000) Fy1c3dr(-lag+(1:1000))]);
        else
            plot([Fy1(lag+(1:1000)) Fy1c3dr(1:1000)]);
        end 
        legend('Fy1 TXT','Fy1 C3D')
        fprintf('ERROR: Vertical force of forceplate 1 does not agree within %.1f N\n',maxdiff);
        fprintf(log,'ERROR: Based on force plate data, it seems that TXT and C3D are not the same trial');
    end
    
    % set up the copying of marker data from c3d to txt
    if (lag >= 0)
        fprintf('C3D file started %d ms later than the D-Flow TXT file\n', 10*lag);
        fprintf('   the new TXT file will have the first %d frames of the original TXT file\n', lag);
        fprintf(log,'C3D file started %d ms later than the D-Flow TXT file\n', 10*lag);
        fprintf(log,'   the new TXT file will have the first %d frames of the original TXT file\n', lag);
        txt_start = lag;   % start here when copying C3D data into the new txt file
        c3d_start = 0;       % start here when copying the c3d data
    elseif (lag < 0)
        fprintf('C3D file started %d ms earlier than the TXT file\n', -10*lag);
        fprintf(log,'C3D file started %d ms earlier than the TXT file\n', -10*lag);
        txt_start = 0;     
        c3d_start = -lag; 
    end

    mis = txtNframes - c3dNframes - lag;  % number of C3D frames missing at the end
    if (mis >= 0)
        fprintf('C3D file stopped %d ms earlier than the D-Flow TXT file\n', 10*mis);
        fprintf('   the new TXT file will have the last %d frames of the original TXT file\n', mis);
        fprintf(log,'C3D file stopped %d ms earlier than the D-Flow TXT file\n', 10*mis);
        fprintf(log,'   the new TXT file will have the last %d frames of the original TXT file\n', mis);
        ncopy = c3dNframes - c3d_start;  % how many frames to copy from C3D to TXT
    elseif(mis < 0)
        fprintf('C3D file stopped %d ms later than the D-Flow TXT file\n', 10*mis)
        fprintf(log,'C3D file stopped %d ms later than the D-Flow TXT file\n', 10*mis)
        ncopy = txtNframes - txt_start;
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
            error('Column %s was not found in the original TXT file',colnames{i});
        end
        txtcol = [txtcol find(strcmp(data.colheaders,colnames{i}))];
    end

    % first copy the relevant columns from original TXT data
    newdata = data.data(:,txtcol);

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

    % Write the marker and force data to .txt file
    % marker data will come from the C3D file
    % time stamp, frame number, and force data will be copied from the original TXT file
    filename = strrep(txt_filename, '.txt', '_edited.txt');  % change .c3d to .txt in the file name
    fprintf('Writing %s...\n', filename);
    fprintf(log,'Writing %s...\n', filename);
    fid = fopen(filename,'w');
    if (fid < 0)
        warning('cannot write %s', filename);
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
        fprintf(fid,'%f\t%d',data.data(iframe,1:2));    % write time stamp fand frame number
        fprintf(fid,'\t%f', newdata(iframe,:));  % write new data
        fprintf(fid,'\n');  % end of line
    end

    fclose(fid);
    fprintf(log,'Done.\n');
    fclose(log);

    % info = 0 to indicate success
    result.info = 0;

end