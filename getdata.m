function [mocapdata,treadmilldata] = getdata(name)
    % the name must be, for example: 'Par7_PRE\mocap0001.txt'
   
    % find the path to the shared data folders
    path = getpath();  

    % import the mocap data
    mocapdata = importdata([path name]);
    mocapdata = cleanup(mocapdata);   % interpolate missing markers, and insert NaN for large gaps
    mocapdata.name = name;  % add the short filename to the data struct
    % make it suitable for Latex interpreter
    latexname = strrep(name,'\','/');
    mocapdata.latexname = strrep(latexname,'_','\_');

    % make new time stamps (we know that Cortex sends data at exactly 100 Hz
    mocapdata.data(:,1) = mocapdata.data(1,1) + 0.01 * (0:(size(mocapdata.data,1)-1));
   
    % construct treadmill filename from the trial number
    num = extract(name, digitsPattern);
    i = strfind(name,'\'); % find the backslashes in the mocap file name
    name = [name(1:i(end)) 'Treadmill' num{end} '.txt']; % make the treadmill file name
    treadmilldata = importdata([path name]);
    treadmilldata.name = name;
    
end
    
    