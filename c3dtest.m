% Testing the C3D to txt file conversion on Participant 4, POST test, trial 1
% This is only used during development. For normal use, c3dbatch.m will be
% more effective.

% edit the following two lines so they point to your test trial:
txt_filename = ['C:\Users\Ton\Cleveland State University\Hala E Osman - Hala data\Par4_POST\Mocap0001.txt'];
c3d_filename = ['C:\Users\Ton\Cleveland State University\Hala E Osman - Hala data\Par4_POST\Mocap0001.c3d'];

% run the conversion tool, to create ...\Mocap0001_edited.txt
result = c3dtotxt(c3d_filename, txt_filename);
result

% compare the original file to the new file
% they should be same, except for gaps filled in the c3d data
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
for col2 = 1:numel(data2.colheaders)
    varname = data2.colheaders{col2};  % name of variable i in data 2
    col1 = find(strcmp(data1.colheaders, varname));  % column number in data1
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
    disp('Hit ENTER to continue');
    pause
end
    
