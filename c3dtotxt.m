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
%   2: synchronization detected negative time lag.  ask Ton to modify the code
%   3: could not write the _edited.txt file
%
% The C3D file reading is done using the Opensim C3D tools: https://simtk-confluence.stanford.edu:8443/display/OpenSim/C3D+(.c3d)+Files

% During testing, we sometimes read the C3D data from tmp.mat because
% extracting data from the C3D file is very slow
if exist('tmp.mat')
    load('tmp.mat');
else
    % load OpenSim libs
    import org.opensim.modeling.*

    % Construct an opensimC3D object with input c3d path
    % Constructor takes full path to c3d file and an integer for forceplate representation (1 = COP).
    fprintf('Opening %s...\n', c3d_filename);
    fprintf('This will take several minutes.\n');
    c3d = osimC3D(c3d_filename,0);

    % Get the c3d data as Matlab Structures
    fprintf('Extracting data...\n');
    [markerStruct forceStruct] = c3d.getAsStructs();
    % save('tmp.mat','markerStruct','forceStruct');  % use this only for testing
end

markernames = fieldnames(markerStruct);
nFrames = size(markerStruct.time, 1);
nMarkers = numel(markernames) - 1;  % the last "marker" in the c3d is actually time

% downsample Fy1 data from 1000 Hz to 100 Hz by averaging 10 sequential force
% samples into 1
nforcesamples = size(forceStruct.f1,1);
if nforcesamples ~= 10*nFrames
    warning('c3dtotxt.m: number of force samples in c3d is not 10x number of frames');
    result.info = 1;
    return
end
Fy1c = zeros(nFrames,1); % to Fy, forceplate 1, from forceStruct
for i = 1:nFrames
    j = 10*(i-1) + (1:10);  % the 10 samples that must be averaged
    Fy1c(i)   = mean(forceStruct.f1(j,2));
end

% import the original TXT file and determine how much marker data was
% missing there
data = importdata(txt_filename);
nFramesTxt = size(data.data,1);
result.missing_before = 0;
for i = 1:nMarkers
    col = find(strcmp(data.colheaders, [markernames{i} '.PosX']));
    if isempty(col)     % this marker was not found in the original TXT file
        nmissing = nFramesTxt;   % all frames are missing for this marker
    else
        nmissing = numel(find(~data.data(:,col)));  % find how many zeros in column
    end
    % fprintf('%s in TXT has %d missing data\n',data.colheaders{col},nmissing);
    result.missing_before = result.missing_before + nmissing;
end

% get the same Fy1 signal from the original D-Flow TXT file
iFy1 = find(strcmp(data.colheaders,'FP1.ForY'));  % channel number for FP1.ForY
Fy1 = data.data(:,iFy1);                          % FP1.ForY signal from the TXT file

% find the time shift using peak of the cross correlation function
[c,lags] = xcorr(Fy1,Fy1c,1000);
[~,imax] = max(c);
lag = lags(imax);

% set up the copying of marker data from c3d to txt
if lag >= 0
    fprintf('C3D file was recorded %d ms later than the D-Flow TXT file\n', 10*lag);
    fprintf('The new TXT file will not include the first %d frames of the original TXT file\n', lag);
    txt_start = lag+1;   % start here when creating the new txt file
    c3d_start = 1;       % start here when copying the c3d data
else
    fprintf('C3D file was recorded %d ms earlier than the TXT file\n', 10*lag);
    txt_start = 1;       % start here when creating the new txt file
    c3d_start = lag+1;   % start here when copyng the c3d data
    return
end

missing_at_end = nFramesTxt - nFrames - lag;  % number of C3D frames missing at the end
if (missing_at_end >= 0)
    fprintf('C3D file stopped %d ms earlier than the D-Flow TXT file\n', 10*missing_at_end)
    fprintf('The new TXT file will not include the last %d frames of the original TXT file\n', missing_at_end);
    ncopy = nFrames + 1 - c3d_start;  % how many frames to copy from C3D to TXT
else
    ncopy = nFramesTxt + 1 - txt_start;
end

markerdata = zeros(ncopy,nMarkers,3);  % where the new marker data is stored
for i = 1:nMarkers
    % copy from markerStruct and also convert mm to m
    markerdata(:,i,:) = 0.001 * markerStruct.(markernames{i})(c3d_start-1+(1:ncopy),:);
end

% convert NaN to zero (indicates missing marker data)
result.missing_after = 0;
for i = 1:nMarkers
    d = squeeze(markerdata(:,i,:));  % make a 2d matrix with data from marker i
    [row,col] = find(isnan(d));
    % fprintf('%s in C3D has %d missing data\n',markernames{i},numel(row)/3);
    result.missing_after = result.missing_after + numel(row)/3;   % keep track of how much data is missing
    d(row,col) = 0.0;               % put zeros where the NaN were
    markerdata(:,i,:) = d;          % put it back into the 3d matrix
end

%% Write the marker and force data to .txt file
% marker data will come from the C3D file
% time stamp, frame number, and force data will be copied from the original TXT file
forcevar = {'FP1.CopX','FP1.CopY','FP1.CopZ','FP1.ForX','FP1.ForY','FP1.ForZ','FP1.MomX','FP1.MomY','FP1.MomZ', ...
            'FP2.CopX','FP2.CopY','FP2.CopZ','FP2.ForX','FP2.ForY','FP2.ForZ','FP2.MomX','FP2.MomY','FP2.MomZ'};
filename = strrep(txt_filename, '.txt', '_edited.txt');  % change .c3d to .txt in the file name
fprintf('Writing %s...\n', filename);
fid = fopen(filename,'w');
if (fid < 0)
    warning('cannot write %s', filename);
    result.info = 3;
    return
end

% write the header line
fprintf(fid,'TimeStamp\tFrameNumber');
for i = 1:nMarkers
    name = markernames{i};
    fprintf(fid,'\t%s.PosX\t%s.PosY\t%s.PosZ',name,name,name);
end
for i = 1:numel(forcevar)
    fprintf(fid,'\t%s', forcevar{i});
end
fprintf(fid,'\n');  % end of line

% write the data
for i = 1:ncopy
    iframe = txt_start + (i-1);
    fprintf(fid,'%f',data.data(iframe,1));          % write time stamp from original TXT file
    fprintf(fid,'\t%d',data.data(iframe,2));        % write frame number from original TXT file
    for j = 1:nMarkers
        fprintf(fid,'\t%f', markerdata(i,j,:));  % write x,y,z of this marker
    end
    for j = 1:numel(forcevar)
        col = strcmp(data.colheaders, forcevar{j});   % column number for this force variable
        fprintf(fid,'\t%f', data.data(iframe,col));   % write the value of this force variable in frame i
    end
    fprintf(fid,'\n');  % end of line
end

fclose(fid);

% info = 0 to indicate success
result.info = 0;

end