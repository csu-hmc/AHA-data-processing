function main
% main program to process the data
    close all

    % set up what needs to be done
    % participants = [2 3 4];       % list of participant numbers, also used to construct folder names
    participants = [7];       % example, if you wanted to process just one participant
    conditions = {'_PRE' '_POST'};  % conditions to be analyzed, also used to construct folder names
    details = 0;                    % use 1 for testing, 0 for fully automated processing

    % loop through the participant folders and use processdata.m to do the work
    for i = 1:numel(participants)
        for j = 1:numel(conditions)
            folder = strcat('Par', num2str(participants(i)), conditions{j});
            fprintf('Processing trials from %s\n', folder);
            processdata(folder, details);
        end
    end
end