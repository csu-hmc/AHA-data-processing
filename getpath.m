function path = getpath()
% find the path to the shared data folder

    % location depends on the computer
    computer = getenv('COMPUTERNAME');
    if strcmp(computer,'LRI-102855')
        path = 'C:\Users\Ton\Cleveland State University\Hala E Osman - Hala data\';    
    elseif strcmp(computer,'DESKTOP-0HN0T6U')
        path = 'C:\Users\hallo\OneDrive - Cleveland State University\Hala data\';    
    else
        fprintf('Your computername (%s) was not recognized.\n', computer);
        fprintf('Customize getdata.m by adding two lines similar to lines 8-9\n');
        error('exiting now');
    end

end