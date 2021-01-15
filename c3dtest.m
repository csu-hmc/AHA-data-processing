% Testing the C3D to txt file conversion
folder = 'C:\Users\Ton\Documents\CSU\Lab\People\Osman\software\';

% convert c3d to txt
txt_filename = [folder 'Mocap0001.txt']; % this is the original txt file, used for synchronization
c3d_filename = [folder 'Hala Sub4-POST- trials1.cap_edited.c3d'];
c3dtotxt(c3d_filename, txt_filename);

% compare the original file to the new file
% they should be same, except for gaps filled in the c3d data
data1 = importdata(txt_filename);  % the original data
txt_filename = strrep(c3d_filename, '.c3d', '.txt');
data2 = importdata(txt_filename);  % the data converted from c3d

% compare each channel in the new file to the corresponding channel in the
% old file
close all
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

    plot([d1 d2]);
    legend('original TXT file','TXT converted from C3D');
    title(varname);
    disp('Hit ENTER to continue');
    pause
end
    
