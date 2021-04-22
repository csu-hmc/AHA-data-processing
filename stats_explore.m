function stats_explore
% exploring the statistical analysis
    close all

    % load the data
    datafile = '2021-04-16 data 8var 12part.xlsx';
    d = readcell(datafile);
    columns = [2 4:13];    % remove columns 1 (session) and 3 (initials)
    headers = d(1,columns);
    % extract 24 rows of data (12 participants pre and post)
    data = d(2:25,columns);
    data = cell2mat(data);  % convert to matrix
    
    % calculate the first two principal components
    d = data(:,end-7:end);  % these are the 8 outcome variables
    [coeff,score,latent] = pca(d);
    figure(1);
    subplot(2,2,1);
    plot(latent);
    xlabel('principal component');
    ylabel('eigenvalue');
    
    % calculate the mahalanobis distance
    C = cov(d);  % covariance matrix
    mu = mean(d);
    d = d - repmat(mu,size(d,1),1);  % subtract the mean from all rows
    MD = sqrt(diag(d * inv(C) * d'));  % mahalanobis distance
    subplot(2,2,2);
    plot(MD,'o');
    xlabel('data row');
    ylabel('Mahalanobis distance');
    
    % add PC1 and PC2 scores to the data matrix
    headers = [headers 'PC1' 'PC2'];
    data = [data score(:,1:2)];
    
    % time and group factors
    Time = [0 1]';          % pre and post time values
    Group = {'PT','Slip','Gaming'}';

    % sort the data, so PRE of all participants is followed by POST of the
    % same participants
    col_time         = find(strcmp(headers,'PRE/POST'));
    col_participant = find(strcmp(headers,'Participant'));
    col_intervention = find(strcmp(headers,'Intervention'));
    data = sortrows(data,[col_time col_participant]);
    nparticipants = size(data,1)/2;
    data_pre = data(1:nparticipants, :);
    data_post = data(nparticipants + (1:nparticipants),:);
    group = data_pre(:,col_intervention);
    if ~isequal(group,data_post(:,col_intervention))
        error('pre and post data are not paired correctly');
    end
    if ~isequal(data_pre(:,col_participant), (1:nparticipants)')
        error('PRE data does not have correct participant numbers');
    end
    if ~isequal(data_post(:,col_participant), (1:nparticipants)')
        error('POST data does not have correct participant numbers');
    end
    
    % make pre-post comparison plots for all variables
    figure(3)
    groupmarkers = 'ox+';
    groupcolors = 'rgb';
    for i = 4:numel(headers)
        subplot(3,4,i-3);
        for p = 1:nparticipants
           g = group(p);
           h = plot([1 2],[data_pre(p,i) data_post(p,i)],'Marker',groupmarkers(g),'Color',groupcolors(g),'LineWidth',1.5);
           hold on
        end
        set(gca,'XTick',1:2,'XTickLabel',{'PRE','POST'})
        set(gca,'XLim',[0.5 2.5])  
        if i==4
            legend(Group);
        end
        title(headers(i))
    end
    
    % list of dependent variables we want to explore
    depvars = {'MRP'};
    
    % do RANOVA for each dependent variable in the list
    % following the "Longitudinal Data" example in the help for "ranova"
    for i = 1:numel(depvars)
        depvar = depvars{i};
        % extract the relevant columns from data and make a table for rmfit
        outcome_pre  = data_pre(:,strcmp(headers,depvar));
        outcome_post = data_post(:,strcmp(headers,depvar));

        t = table(group,outcome_pre,outcome_post,...
        'VariableNames',{'Group','PRE','POST'});

        % fit the repeated measures model
        rm = fitrm(t,'PRE,POST ~ Group','WithinDesign',Time,'WithinModel','orthogonalcontrasts');
        ranovatbl = ranova(rm)
    end
    
end
