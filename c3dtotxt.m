function c3dtotxt(c3d_filename, txt_filename)
% converts a C3D file to txt
% c3d_filename: the file with C3D data
% txt_filename: the corresponding txt data, recorded with D-Flow, for synchronization
% uses the Opensim C3D tools: https://simtk-confluence.stanford.edu:8443/display/OpenSim/C3D+(.c3d)+Files

%% Load OpenSim libs
import org.opensim.modeling.*
 
%% Construct an opensimC3D object with input c3d path
% Constructor takes full path to c3d file and an integer for forceplate
% representation (1 = COP).
if exist('tmp.mat')
    load('tmp.mat')
else
    fprintf('Opening %s...\n', c3d_filename);
    fprintf('This will take several minutes.\n');
    c3d = osimC3D(filename,0);

    %% Get some stats...
    % Get the number of marker trajectories
    % nMarkers = c3d.getNumTrajectories();
    % Get the marker data rate
    % rMarkers = c3d.getRate_marker();
    % Get the number of forces
    % nForces = c3d.getNumForces();
    % Get the force data rate
    % rForces = c3d.getRate_force();

    % Get Start and end time
    % t0 = c3d.getStartTime();
    % tn = c3d.getEndTime();

    %% Rotate the data ?
    % c3d.rotateData('x',-90)

    %% Get the c3d data as Matlab Structures
    fprintf('Extracting data...\n');
    [markerStruct forceStruct] = c3d.getAsStructs();
    save('tmp.mat','markerStruct','forceStruct');
end

markernames = fieldnames(markerStruct);
nFrames = size(markerStruct.time, 1);

% downsample Fy1 data from 1000 Hz to 100 Hz by averaging 10 sequential force
% samples into 1
nforcesamples = size(forceStruct.f1,1);
if nforcesamples ~= 10*nFrames
    error('c3dtotxt.m: number of force samples in c3d is not 10x number of frames');
end
Fy1c = zeros(nFrames,1); % to Fy, forceplate 1, from forceStruct
for i = 1:nFrames
    j = 10*(i-1) + (1:10);  % the 10 samples that must be averaged
    Fy1c(i)   = mean(forceStruct.f1(j,2));
end

% get the same Fy1 signal from the original D-Flow TXT file
data = importdata(txt_filename);
nFramesTxt = size(data.data,1);
iFy1 = find(strcmp(data.colheaders,'FP1.ForY'));  % channel number for FP1.ForY
Fy1 = data.data(:,iFy1);                          % FP1.ForY signal from the TXT file

% find the time shift using peak of the cross correlation function
[c,lags] = xcorr(Fy1,Fy1c,1000);
[~,imax] = max(c);
lag = lags(imax);
if lag >= 0
    fprintf('C3D file was recorded %d ms later than the D-Flow TXT file\n', 10*lag);
    fprintf('The new TXT file will start with %d frames of missing markers\n', lag);
else
    fprintf('C3D file was recorded %d ms earlier than the TXT file\n', 10*lag);
    error('negative lag was not implemented yet');
end
missing = nFramesTxt - nFrames - lag;  % C3D data missing at the end
if (missing >= 0)
    fprintf('C3D file stopped %d ms earlier than the D-Flow TXT file\n', 10*missing)
    fprintf('The new TXT file will end with %d frames of missing markers\n', missing);
else
    missing = 0;
end

nMarkers = numel(markernames) - 1;  % the last "marker" is actually time
markerdata = zeros(nFramesTxt,nMarkers,3);  % where the new marker data is stored
for i = 1:nMarkers
    markerdata(1:lag, i, :) = 0.0;  % no marker data for the first lag frames
    ncopy = nFramesTxt-lag-missing;
    % copy from markerStruct and also convert mm to m
    markerdata(lag+1:end-missing,i,:) = 0.001 * markerStruct.(markernames{i})(1:ncopy,:);
end

% convert NaN to zero (indicates missing marker data)
for i = 1:nMarkers
    d = squeeze(markerdata(:,i,:));  % make a matrix with data from marker i
    [row,col] = find(isnan(d));
    d(row,col) = 0.0;               % put zeros where the NaN were
    markerdata(:,i,:) = d;          % put it back into the 3d matrix
end

% to do: handle negative lag
% to do: COP should be calculated from Moment (but we don't use it)

%% Write the marker and force data to .txt file
% marker data will come from the C3D file
% time stamp, frame number, and force data will be copied from the original TXT file
forcevar = {'FP1.CopX','FP1.CopY','FP1.CopZ','FP1.ForX','FP1.ForY','FP1.ForZ','FP1.MomX','FP1.MomY','FP1.MomZ', ...
            'FP2.CopX','FP2.CopY','FP2.CopZ','FP2.ForX','FP2.ForY','FP2.ForZ','FP2.MomX','FP2.MomY','FP2.MomZ'};
filename = strrep(c3d_filename, '.c3d', '.txt');  % change .c3d to .txt in the file name
if strcmp(filename, txt_filename)
    error('c3dtotxt.m: attempting to overwrite the original .txt file');
end
fprintf('Writing %s...\n', filename);
fid = fopen(filename,'w');
if (fid < 0)
    error('cannot write .txt file');
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
for i = 1:nFramesTxt
    fprintf(fid,'%f',data.data(i,1));          % write time stamp from original TXT file
    fprintf(fid,'\t%d',data.data(i,2));        % write frame number from original TXT file
    for j = 1:nMarkers
        fprintf(fid,'\t%f', markerdata(i,j,:));  % write x,y,z of this marker
    end
    for j = 1:numel(forcevar)
        col = find(strcmp(data.colheaders, forcevar{j}));  % column number for this force variable
        fprintf(fid,'\t%f', data.data(i,col));   % write the value of this force variable in frame i
    end
    fprintf(fid,'\n');  % end of line
end

fclose(fid);

end