function stats_explore
% exploring the statistical analysis

    close all

    % load the data
    datafile = '2021-04-23 data 8var 12part.xlsx';
    d = readcell(datafile);
    columns = [2 4:13];    % remove columns 1 (session) and 3 (initials)
    headers = d(1,columns);
    % extract 24 rows of data (12 participants pre and post)
    data = d(2:25,columns);
    data(strcmp(data,'NaN')) = {NaN};   % put the value NaN where there is a string 'NaN'
    data = cell2mat(data);  % convert to matrix
    d = data(:,end-7:end);  % select the columns with the 8 outcome variables
    varnames = headers(end-7:end);
    
    % do the Matlab PCA
    nvar = size(d,2);
    mu = mean(d);  % mean values of all variables
    sd = std(d);
    n = size(d,1);  % number of observations
    drel = (d - repmat(mu,n,1)) ./ repmat(sd,n,1);  % data minus mean, normalized to SD
    [coeff, score, latent]  = pca(drel);
    figure(1);
    subplot(2,2,1);
    plot(latent);
    xlabel('principal component');
    ylabel('eigenvalue');
    
    % print the PCA loadings
    fprintf('Principal Component loadings:\n');
    fprintf('          ');
    for i = 1:nvar
        fprintf('     PC%d',i);
    end
    fprintf('\n');
    for i = 1:nvar
        fprintf('%10s',varnames{i});
        for j = 1:nvar
            fprintf('%8.3f',coeff(i,j));
        end
        fprintf('\n');
    end
 
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
    
    % group labels
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
    % for example: depvars = {'MRP','PC1'};
    depvars = headers(4:end);   % use all variables
    
    % do 1-way ANOVA for each dependent variable in the list
    for i = 1:numel(depvars)
        depvar = depvars{i};
        fprintf('===========================================================================\n');
        fprintf('Statistical analysis for variable: %s\n', depvar);
        
        % extract the relevant columns from data and make a (subj x group) matrix for ANOVA
        outcome_pre  = data_pre(:,strcmp(headers,depvar));
        outcome_post = data_post(:,strcmp(headers,depvar));
        effect = outcome_post - outcome_pre;
         for k = 1:3
            m(:,k) = effect(find(group==k));
         end
        fprintf('Post-pre differences (4 subjects x 3 groups):\n')
        m
        
        % Do the ANOVA
        fprintf('One-way ANOVA results:\n');
        [p,tbl,stats] = anova1(m,Group,'off');
        fprintf('---------------------------------------------------\n');
        fprintf('%7s %9s %5s %9s %7s %9s\n',tbl{1,:});
        fprintf('---------------------------------------------------\n');
        fprintf('%7s %9.4f %5d %9.5f %7.2f %9.4f\n',tbl{2,:});
        fprintf('%7s %9.4f %5d %9.5f\n',tbl{3,1:4});
        fprintf('%7s %9.4f %5d\n',tbl{4,1:3});
        fprintf('---------------------------------------------------\n');

        % do the pairwise comparisons, if one-way ANOVA had p<0.05
        if (p<0.05)
            [results,means] = multcompare(stats,'Display','off');
            for j = 1:size(results,1)
                group1 = Group{results(j,1)};
                group2 = Group{results(j,2)};
                fprintf('%6s vs. %6s: mean difference %8.3f (conf interval: %8.3f %8.3f) p=%8.3f\n', ...
                    group1, group2, results(j,3:6))
            end
        end
    end
end