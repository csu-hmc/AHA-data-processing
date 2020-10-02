function c3dtotxt
% converts a C3D file to txt
% uses the Opensim C3D tools: https://simtk-confluence.stanford.edu:8443/display/OpenSim/C3D+(.c3d)+Files

%% Load OpenSim libs
import org.opensim.modeling.*
 
%% Get the path to a C3D file
folder = 'C:\Users\Ton\Documents\CSU\Lab\People\Osman\software\';
filename = 'Hala Sub4-POST- trials1.cap_edited.c3d';
 
%% Construct an opensimC3D object with input c3d path
% Constructor takes full path to c3d file and an integer for forceplate
% representation (1 = COP).
if exist('tmp.mat')
    load('tmp.mat')
else
    c3d = osimC3D([folder filename],0);

    %% Get some stats...
    % Get the number of marker trajectories
    nMarkers = c3d.getNumTrajectories();
    % Get the marker data rate
    rMarkers = c3d.getRate_marker();
    % Get the number of forces
    nForces = c3d.getNumForces();
    % Get the force data rate
    rForces = c3d.getRate_force();

    % Get Start and end time
    t0 = c3d.getStartTime();
    tn = c3d.getEndTime();

    %% Rotate the data
    % c3d.rotateData('x',-90)

    %% Get the c3d data as Matlab Structures
    [markerStruct forceStruct] = c3d.getAsStructs();
end

markernames = fieldnames(markerStruct);
nframes = size(markerStruct.time, 1);

% downsample force data from 1000 Hz to 100 Hz by moving average

% here we should load the treadmill file and create time stamps
% after looking at Mx for the synchronization

%% Write the marker and force data to .txt file
filename = strrep(filename, '.c3d', '.txt');
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
fprintf(fid,'\tFP1.CopX\tFP1.CopY\tFP1.CopZ\tFP1.ForX\tFP1.ForY\tFP1.ForZ\tFP1.MomX\tFP1.MomY\tFP1.MomZ');
fprintf(fid,'\tFP2.CopX\tFP2.CopY\tFP2.CopZ\tFP2.ForX\tFP2.ForY\tFP2.ForZ\tFP2.MomX\tFP2.MomY\tFP2.MomZ');
fprintf(fid'\n');  % end of line

% write the data
for i = 1:nframes

end

fclose(fid);

end