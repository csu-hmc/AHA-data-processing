function c3dbatch
% Create edited mocap files from all mocap files in a series of folders
% Marker data is taken from c3d files in the same folder
% The original mocap files are named Mocapxxxx.txt, where xxxx is a 4-digit
% trial number.
% The c3d files are named Mocapxxxx.c3d
% The edited mocap files are named Mocapxxxx_edited.txt

    % here we define where the mocap folders are
    % edit this code to match your computer
    computer = getenv('COMPUTERNAME');
    if strcmp(computer, 'LRI-102855')   % Ton's computer
        datapath = 'C:\Users\Ton\Cleveland State University\Hala E Osman - Hala data\';
    elseif strcmp(computer, 'DESKTOP-0HN0T6U')   % Hala's computer
        datapath = 'C:\Users\hallo\OneDrive - Cleveland State University\Hala data\';
    else
        fprintf('Your computer name is: %s\n', computer);
        fprintf('Please configure c3dbatch.m for your computer.\n');
    end
    
    % here we define which folders must be processed
    % (eventually, this list should contain all folders, and we run this
    % code once to make sure that all C3D files were processed the same
    % way)
    folders = {'Par7_POST'};  % a cell array containing one or more folder names
    
    % create a log file, or append to an existing log file if it exists
    logfile = 'c3d.log';
    fid = fopen(logfile,'w');
    if (fid < 0)
        error('cannot write %s', logfile);
    end
    fprintf(fid,'======================================\n');
    fprintf(fid,'c3dbatch processing started on %s\n', datetime);
    fclose(fid);
    
    % go through the folders in the folder list
    forall = 0;
    for i = 1:numel(folders)
        folder = [datapath folders{i} '\']; 
        % get a list of all the c3d files in this folder
        files = dir([folder 'Mocap*.c3d']);
        if numel(files)==0
            error('there are no mocap files in %s', folder)
        end
        
        % go through all the trials for which c3d files were found
        for j = 1:numel(files)
            fid = fopen(logfile,'a');
            if (fid < 0)
                error('cannot write %s', logfile);
            end
            
            % generate the full file names
            c3d_filename = [folder files(j).name];
            txt_filename = strrep(c3d_filename, '.c3d', '.txt');
            edited_filename = strrep(c3d_filename, '.c3d', '_edited.txt');

            % if edited file already exists, ask user what to do (or use
            % their previous decision)
            if exist(edited_filename)
                if ~forall
                    [skip,forall] = askuser([edited_filename]);
                end
            else
                skip = 0;
            end
            if ~skip
                fprintf('Processing %s\n',c3d_filename);
                result = c3dtotxt(c3d_filename, txt_filename);
                if result.info == 0
                    fprintf(fid,'Done. Missing markers was %d, is now %d\n', result.missing_before, result.missing_after);
                else
                    fprintf(fid,'Error in c3dtotxt.m. Please report the problem.\n', c3d_filename);
                    result
                end
            else
                fprintf(fid,'Skipping %s\n', c3d_filename);
            end
            fclose(fid);
        end
    end
end
%=================================================
function [skip, forall] = askuser(filename)
    skip = 1;
    forall = 0;
    d = dialog('Position',[200 200 500 150],'Name','File exists');
    uicontrol(d,'Style','text','Position',[20 100 460 40],'String',filename, ...
        'FontSize',12,'FontWeight','bold', 'HorizontalAlignment','left');
    uicontrol(d,'Style','text','Position',[20 50 300 40],'String','already exists.', ...
        'FontSize',12,'FontWeight','bold', 'HorizontalAlignment','left');
       
    uicontrol(d,'Position',[15 20 80 25],'String','skip file','Callback',@button1, ...
        'FontSize',12,'FontWeight','bold');
    uicontrol(d,'Position',[110 20 80 25],'String','redo file','Callback',@button2, ...
        'FontSize',12,'FontWeight','bold');
    uicontrol(d,'Position',[205 20 170 25],'Style','checkbox','String','do this for all files', ...
        'HorizontalAlignment','left','Callback',@box, ...
        'Value',0,'FontSize',12,'FontWeight','bold');
    uiwait(d);
    
    function button1(src,event)
        skip = 1;
        delete(gcf);
    end
    function button2(src,event)
        skip = 0;
        delete(gcf);
    end
    function box(src,event)
        forall = src.Value;
    end
end

