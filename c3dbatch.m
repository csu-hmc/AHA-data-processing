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
    end
    
    % here we define which folders must be processed
    folders = {'Par4_POST','Par4_PRE'};  % a cell array containing folder names
    
    % create a log file
    logfile = 'c3dbatch.log';
    fid = fopen(logfile,'a');
    if (fid < 0)
        error('cannot write %s', logfile);
    end
    fprintf(fid,'c3dbatch processing started %s\n', datetime);
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
            c3d_filename = [folder files(j).name];
            txt_filename = strrep(c3d_filename, '.c3d', '.txt');
            edited_filename = strrep(c3d_filename, '.c3d', '_edited.txt');
            % if edited file does not exist yet, run the conversion tool, to create ..._edited.txt
            if exist(edited_filename)
                fprintf(fid,'%s already exist.  Skipping.\n', edited_filename);
            else    
                result = c3dtotxt(c3d_filename, txt_filename);
                if result.info == 0
                    fid = fopen(logfile,'a');
                    if (fid < 0)
                        error('cannot write %s', logfile);
                    end
                    fprintf(fid,'%s created.  Missing markers was %d, is now %d\n', ...
                        edited_filename, result.missing_before, result.missing_after);
                    fclose(fid);
                else
                    fprintf('Error in c3dtotxt.m while processing %s'.\n', c3d_filename);
                    result
                end
            end
        end
    end
end
    
