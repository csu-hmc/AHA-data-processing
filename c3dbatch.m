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
    elseif strcmp(computer, '...hala computer name here...')   % Hala's computer
        datapath = '...path to data folders on hala computer';
    else
        fprintf('Your computer name is: %s\n', computer);
        fprintf('Please configure c3dbatch.m for your computer.\n');
    end
    
    % here we define which folders must be processed
    folders = {'Par4_POST'};  % a cell array containing one or more folder names
    
    % create a log file, or append to an existing log file if it exists
    logfile = 'c3dbatch.log';
    fid = fopen(logfile,'a');
    if (fid < 0)
        error('cannot write %s', logfile);
    end
    fprintf(fid,'c3dbatch processing started on %s\n', datetime);
    fclose(fid);
    
    % go through the folders in the folder list
    for i = 1:numel(folders)
        folder = [datapath folders{i} '\']; 
        % get a list of all the c3d files in this folder
        files = dir([folder 'Mocap*.c3d']);
        if numel(files)==0
            error('there are no mocap files in %s', folder)
        end
        
        % go through all the trials for which c3d files were found
        for j = 1:numel(files)
            if ~strcmp(files(j).name, 'Mocap0001.c3d')   % for testing only
                break
            end
            % write something on the log file
            fid = fopen(logfile,'a');
            if (fid < 0)
                error('cannot write %s', logfile);
            end
            fprintf(fid,'   %s: ',[folders{i} '\' strrep(files(j).name,'.c3d','')]);

            % if edited file does not exist yet, run the conversion tool, to create ..._edited.txt
            c3d_filename = [folder files(j).name];
            txt_filename = strrep(c3d_filename, '.c3d', '.txt');
            edited_filename = strrep(c3d_filename, '.c3d', '_edited.txt');
            if exist(edited_filename)
                fprintf(fid,'_edited.txt already exists.  Skipping.\n');
            else    
                result = c3dtotxt(c3d_filename, txt_filename);
                if result.info == 0
                    fprintf(fid,'Done. Missing markers was %d, is now %d\n', result.missing_before, result.missing_after);
                else
                    fprintf(fid,'Error in c3dtotxt.m. Please report the problem.\n', c3d_filename);
                    result
                end
            end
            fclose(fid);
        end
    end
end
    
