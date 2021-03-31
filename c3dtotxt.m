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
    fprintf('Extracting data... (will take approx. %.1f minutes for a 90-second trial)\n',10);
    tic
    [markerStruct, forceStruct] = c3d.getAsStructs();
    fprintf('...it took %.1f minutes\n',toc/60);
    % save('tmp.mat','markerStruct','forceStruct');  % use this only for testing
end

markernames = fieldnames(markerStruct);
nFrames = size(markerStruct.time, 1);
nMarkers = numel(markernames) - 1;  % the last "marker" in the c3d is actually time

% downsample Fy1 data from 1000 Hz to 100 Hz by averaging 10 sequential force
% samples into 1
nforcesamples = size(forceStruct.f1,1);
if nforcesamples ~= 10*nFrames
    error('c3dtotxt.m: number of force samples in c3d is not 10x number of frames');
end
fpmdata = zeros(nFrames,18); % to store force,CoP,moment from the two force plates
fpmnames = {'FP1.ForX','FP1.ForY','FP1.ForZ','FP1.CopX','FP1.CopY','FP1.CopZ','FP1.MomX','FP1.MomY','FP1.MomZ', ...
            'FP2.ForX','FP2.ForY','FP2.ForZ','FP2.CopX','FP2.CopY','FP2.CopZ','FP2.MomX','FP2.MomY','FP2.MomZ'};
for i = 1:nFrames
    j = 10*(i-1) + (1:10);  % the 10 samples that must be averaged
    % also convert mm and Nmm to m and Nm
    fpm = [ forceStruct.f1(j,:) forceStruct.p1(j,:)/1000 forceStruct.m1(j,:)/1000 ...
            forceStruct.f2(j,:) forceStruct.p2(j,:)/1000 forceStruct.m2(j,:)/1000 ];
    fpmdata(i,:)   = mean(fpm);  % take the average of the 10 samples
end

% import the original TXT file and generate more reliable time stamps
data = importdata(txt_filename);
nFramesTxt = size(data.data,1);
data.data(:,1) = data.data(1,1) + (0:(nFramesTxt-1))*0.010;

% remove spaces from the column headers
for i = 1:numel(data.colheaders)
    data.colheaders{i} = strrep(data.colheaders{i},' ','');
end

% get the Fy1 signal from the original D-Flow TXT file
iFy1 = strcmp(data.colheaders,'FP1.ForY');  % channel number for FP1.ForY
Fy1 = data.data(:,iFy1);                          % FP1.ForY signal from the TXT file

% find the time shift using peak of the cross correlation on Fy of force
% plate 1
[c,lags] = xcorr(Fy1,fpmdata(:,2),1000);  % max lag is 1000 frames (10 seconds)
[~,imax] = max(c);
lag = lags(imax);

% set up the copying of marker data from c3d to txt
if (lag > 0)
    fprintf(log,'C3D file started %d ms later than the D-Flow TXT file\n', 10*lag);
    fprintf(log,'   the new TXT file will have the first %d frames of the original TXT file\n', lag);
    txt_start = lag;   % start here when copying C3D data into the new txt file
    c3d_start = 0;       % start here when copying the c3d data
elseif (lag < 0)
    fprintf(log,'C3D file started %d ms earlier than the TXT file\n', 10*lag);
    txt_start = 0;     
    c3d_start = lag; 
end

mis = nFramesTxt - nFrames - lag;  % number of C3D frames missing at the end
if (mis > 0)
    fprintf(log,'C3D file stopped %d ms earlier than the D-Flow TXT file\n', 10*mis);
    fprintf(log,'   the new TXT file will have the last %d frames of the original TXT file\n', mis);
    ncopy = nFrames - c3d_start;  % how many frames to copy from C3D to TXT
elseif(mis < 0)
    fprintf(log,'C3D file stopped %d ms later than the D-Flow TXT file\n', 10*mis)
    ncopy = nFramesTxt - txt_start;
end

% find the relevant columns in the original TXT data
colnames = {};
for i=1:nMarkers
    colnames = [ colnames [markernames{i} '.PosX'] [markernames{i} '.PosY'] [markernames{i} '.PosZ'] ];
end
colnames = [colnames fpmnames];  % add the force plate data column names
txtcol = [];
for i = 1:numel(colnames)
    if isempty(find(strcmp(data.colheaders,colnames{i})))
        keyboard
    end
    txtcol = [txtcol find(strcmp(data.colheaders,colnames{i}))];
end

% first copy the relevant columns from original TXT data
newdata = data.data(:,txtcol);

% insert ncopy frames from C3D markerStruc, markers only; also convert mm to m
for i = 1:nMarkers
    columns = 3*(i-1) + (1:3);
    newdata(txt_start+(1:ncopy),columns) = 0.001 * markerStruct.(markernames{i})(c3d_start+(1:ncopy),:);
end

% convert NaN to zero (indicates missing marker data)
markercol = (1:3*nMarkers);
for i = markercol
    k = isnan(newdata(:,i));
    newdata(k,i) = 0.0;             % put zeros where the NaN were
end

% find out how many times a marker was missing in original TXT and new TXT
result.missing_before = 0;
result.missing_after = 0;
for i = 1:nMarkers
    col = find(strcmp(data.colheaders, [markernames{i} '.PosX']));
    if isempty(col)     % this marker was not found in the original TXT file
        nmissing = nFramesTxt;   % all frames are missing for this marker
    else
        nmissing = numel(find(~data.data(:,col)));  % find how many zeros in column
    end
    nmissingnew = numel(find(~newdata(:,3*(i-1)+1)));  % find how many zeros in x of this marker
    fprintf(log,'%6s: %5d missing in TXT, %5d missing in edited TXT (out of %5d)\n',markernames{i},nmissing,nmissingnew,nFramesTxt);
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
for iframe = 1:nFramesTxt
    fprintf(fid,'%f\t%d',data.data(iframe,1:2));    % write time stamp fand frame number
    for j = 1:size(newdata,2)
        fprintf(fid,'\t%f', newdata(iframe,j));  % write new data
    end
    fprintf(fid,'\n');  % end of line
end

fclose(fid);
fprintf(log,'Done.\n');
fclose(log);

% info = 0 to indicate success
result.info = 0;

end