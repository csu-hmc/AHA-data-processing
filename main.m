function main
% main program to process the data
close all

% set up what needs to be done
participants = [2 3 4];       % list of participant numbers, also used to construct folder names
% participants = [2];       % example, if you wanted to process just one participant
conditions = {'_PRE' '_POST'};  % conditions to be analyzed, also used to construct folder names
details = 0;                    % use 1 for testing, 0 for fully automated processing

% determine where the shared data folder is
computer = getenv('computername');
if strcmp(getenv('computername'), 'LRI-102855')
    % on Ton's computer the shared data folder is here:
    dataroot = 'C:\Users\Ton\Cleveland State University\Hala E Osman - Hala data';
else
    % on other computers, we start this program from a subfolder of the data folder
    dataroot = what('..').path;
end

% loop through the participant folders and use processdata.m to do the work
for i = 1:numel(participants)
	for j = 1:numel(conditions)
        datafolder = strcat(dataroot, '/Par', num2str(participants(i)), conditions{j});
        resultsfolder = strcat(dataroot, '/results/Par', num2str(participants(i)), conditions{j} );
        processdata(datafolder, resultsfolder, details);
    end
end

end