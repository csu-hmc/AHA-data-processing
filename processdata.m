function processdata(datafolder, resultsfolder, detail)
    % process data from (up to) 15 trials of D-Flow data
    % Inputs:
    %   datafolder............Folder where D-Flow data is stored
    %   resultsfolder.........PNG and XLSX will be stored here
    %   detail................(optional) use 1 when you want to debug
    %
    % Outputs:
    %   XLSX and PNG files, will be stored in resultsfolder
    
    % first make sure that the data folder exists
    if ~exist(datafolder)
        error('processdata: data folder %s does not exist', datafolder);
    end
    
    % turn detail off if not specified
    if nargin < 3
        detail = 0;
    end
    
    % create the result folder if it does not exist
    if ~exist(resultsfolder, 'dir')
        mkdir(resultsfolder);
    end

    SL_L = {};
    SL_R = {};
    ST_L = {};
    ST_R = {};

    ntrials = numel(dir([datafolder '\Mocap*.txt']));  % determine number of trial files
    for trial_num = 1:ntrials 

        % load all the data
        filename = strcat(datafolder, '\Mocap000', num2str(trial_num), '.txt');
        fprintf('Processing %s\n', filename);
        data = importdata(filename);
        data = cleanup(data);   % interpolate missing markers, and insert NaN for large gaps
     
        % calculate step length and step time during normal walking
        if trial_num==7
            keyboard
        end
        [SL_L{trial_num}, SL_R{trial_num}, ST_L{trial_num}, ST_R{trial_num}] = ...
            step_length(data, detail);
        
        % do the Margin of Stability analysis during normal walking (20
        % seconds)
        % ...
        
        % determine the recovery time (how many cycles) from the perturbation
        % (we will need to load the treadmill file to determine when the
        % perturbation happened)
        % do this from the knee angle
        % see how they recover?
    end

%% save results
save_file_name = strcat(resultsfolder, '\normal_walking_SL_ST_20seconds.xlsx');
mean_SL_ST = zeros(ntrials, 5);
std_SL_ST = zeros(ntrials, 5);

for sheet = 1:ntrials

   StepLength_Left = SL_L{sheet};
   StepLength_Right = SL_R{sheet};
   
   StepTime_Left = ST_L{sheet};
   StepTime_Right = ST_R{sheet};
   
   min_raws = min([length(StepLength_Left), length(StepLength_Right), ...
        length(StepTime_Left), length(StepTime_Right)]);
    
   mean_SL_ST(sheet, 1) = mean(StepLength_Left(1:min_raws));
   mean_SL_ST(sheet, 2) = mean(StepLength_Right(1:min_raws));
   mean_SL_ST(sheet, 3) = mean(StepTime_Left(1:min_raws));
   mean_SL_ST(sheet, 4) = mean(StepTime_Right(1:min_raws));
   mean_SL_ST(sheet, 5) = mean(StepTime_Left(1:min_raws) + ...
                            StepTime_Right(1:min_raws));
                        
   std_SL_ST(sheet, 1)= std(StepLength_Left(1:min_raws));
   std_SL_ST(sheet, 2)= std(StepLength_Right(1:min_raws));
   std_SL_ST(sheet, 3)= std(StepTime_Left(1:min_raws));
   std_SL_ST(sheet, 4)= std(StepTime_Right(1:min_raws));
   std_SL_ST(sheet, 5)= std(StepTime_Left(1:min_raws) + ...
                            StepTime_Right(1:min_raws));
   
   
   T = table(StepLength_Left(1:min_raws), StepLength_Right(1:min_raws), ...
        StepTime_Left(1:min_raws), StepTime_Right(1:min_raws));
    
   T.Properties.VariableNames = {'StepLength_Left', 'StepLength_Right', ...
        'StepTime_Left', 'StepTime_Right'};
    
   writetable(T, save_file_name, 'Sheet', sheet, 'Range', 'A1');
   
end

T_mean_summary =  table(mean_SL_ST(:, 1), mean_SL_ST(:, 2), mean_SL_ST(:, 3),...
                        mean_SL_ST(:, 4), mean_SL_ST(:, 5));
T_mean_summary.Properties.VariableNames = {'Mean_StepLength_Left',...
    'Mean_StepLength_Right', 'Mean_StepTime_Left', 'Mean_StepTime_Right', ....
    'Mean_StrideTime'};

T_std_summary =  table(std_SL_ST(:, 1), std_SL_ST(:, 2), std_SL_ST(:, 3),...
                        std_SL_ST(:, 4), std_SL_ST(:, 5));
T_std_summary.Properties.VariableNames = {'STD_StepLength_Left',...
    'STD_StepLength_Right', 'STD_StepTime_Left', 'STD_StepTime_Right', ....
    'STD_StrideTime'};

writetable(T_mean_summary, save_file_name, 'Sheet', 16, 'Range', 'A1')
writetable(T_std_summary, save_file_name, 'Sheet', 16, 'Range', 'G1');




%% average the step lenght and step time ?



for trials = 1:15
    
    
    
    
end


%% plot results

fig2 = figure(2);
subplot(2,2,1)
errorbar(mean_SL_ST(:, 1), std_SL_ST(:, 1))
ylabel('Left Step Length (m)')
xlabel('Trial')
ylim([0 0.5])
xlim([0, 16])

subplot(2,2,2)
errorbar(mean_SL_ST(:, 2), std_SL_ST(:, 2))
ylabel('Right Step Length (m)')
xlabel('Trial')
ylim([0 0.5])
xlim([0, 16])

subplot(2, 2, 3)
errorbar(mean_SL_ST(:, 3), std_SL_ST(:, 3))
ylabel('Left Step Time (s)')
xlabel('Trial')
ylim([0.4 1.4])
xlim([0, 16])

subplot(2, 2, 4)
errorbar(mean_SL_ST(:, 4), std_SL_ST(:, 4))
ylabel('Right Step Time (s)')
xlabel('Trial')
ylim([0.4 1.4])
xlim([0, 16])

saveas(fig2, strcat(resultsfolder, '\result_fig.png'))


end
