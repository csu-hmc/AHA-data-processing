function stats_explore
% exploring the statistical analysis
    global Xrel group    % these variables are needed for the mymanova() function

    close all

    % load the data
    datafile = '2021-04-23 data 8var 12part.xlsx';
    d = readcell(datafile);
    columns = [2 4:13];    % remove columns 1 (session) and 3 (initials), we don't need them
    headers = d(1,columns);
    % extract 24 rows of data (12 participants pre and post)
    data = d(2:25,columns);
    data(strcmp(data,'NaN')) = {NaN};   % put the value NaN where there is a string 'NaN'
    data = cell2mat(data);  % convert to matrix
    d = data(:,end-7:end);  % select the columns with the 8 outcome variables
    varnames = headers(end-7:end);
    
    % do the Matlab PCA
    fprintf('=========================================================\n');
    fprintf('PRINCIPAL COMPONENT ANALYSIS\n');
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

    % sort the data, so PRE of all participants is first, followed by POST of the
    % participants in the same order
    col_time         = find(strcmp(headers,'PRE/POST'));
    col_participant = find(strcmp(headers,'Participant'));
    col_intervention = find(strcmp(headers,'Intervention'));
    data = sortrows(data,[col_time col_participant]);
    nparticipants = size(data,1)/2;
    data_pre = data(1:nparticipants, :);
    data_post = data(nparticipants + (1:nparticipants),:);
    data_diff = data_post-data_pre;
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
    
    % do one-way MANOVA using the differences in the 8 original variables
    fprintf('=========================================================\n');
    fprintf('ONE-WAY MANOVA\n');
    columns = 4:11;  % all 8 variables
    columns = [6 10]; % only NWS and RPS (variable 3 and 8)
    X = data_diff(:,columns);
    varnames = headers(columns);
    grouplabels = Group(group);
    % gplotmatrix(X,[],grouplabels,[],[],[],[],'variable',varnames)
    [d,p,stats] = manova1(X,grouplabels);
    fprintf('MANOVA result for c1: p = %8.4f\n',p(1));
    fprintf('MANOVA result for c2: p = %8.4f\n',p(2));
    
    % plot the first two canonical variables, by group
    c1 = stats.canon(:,1);
    c2 = stats.canon(:,2);
    figure
    gscatter(c2,c1,grouplabels,[],'oxs')
    title('Canonical variables c1,c2 after MANOVA');
    
    % show the transformation from original variables to c1
    showtransform(stats.eigenvec(:,1), varnames);
    
    % show the pairwise comparisons by doing anova on c1
    disp('Do not trust the p values for the pairwise comparisons')
    for k = 1:3
         m(:,k) = c1(find(group==k));
    end
    [p,tbl,stats] = anova1(m,Group,'off');
    results = multcompare(stats,'Display','off','CType','bonferroni');
    for j = 1:size(results,1)
        group1 = Group{results(j,1)};
        group2 = Group{results(j,2)};
        fprintf('%6s vs. %6s: mean difference %8.3f (conf interval: %8.3f %8.3f) p=%8.3f\n', ...
            group1, group2, results(j,3:6))
    end
    
    % do the MANOVA with a brute force method, where we constrain the signs of the transformation coefficients
    % and find the linear transformation that best separates the groups
    fprintf('=========================================================\n');
    fprintf('ONE-WAY MANOVA WITH SIGN CONSTRAINTS\n');
    Xrel = (X - mean(X));  % so we can interpret the coefficients in the same way as with manova1
    positive = [0 1 1 1 1 1 0 1];
    lb = zeros(nvar,1);
    ub = zeros(nvar,1);
    for i=1:nvar
        fprintf('coefficient of %9s is constrained to be: ', varnames{i});
        if positive(i)
            fprintf('positive\n');
            lb(i) = 0;
            ub(i) = 1;
        else
            fprintf('negative\n');
            lb(i) = -1;
            up(i) = 0;
        end
    end
    Aeq = ones(1,nvar);
    Beq = 1;
    opts = optimoptions('fmincon','Display','off');

    % look for a transformation c1 = a(1)*v(1) + a(2)*v(2) + ...
    % such that ANOVA on c1 has the lowest p value for the group effect
    % try this optimization with many initial guesses and keep the best result
    nattempts = 100;
    pbest = 1e10;
    for ii = 1:nattempts
        a0 = lb + rand(size(lb)).*(ub-lb);
        % lb=-ones(size(lb));ub=ones(size(ub));Aeq=[];Beq=[];   % uncomment this line to remove all constraints
        Aeq=[];Beq=[];
        [a,p] = fmincon(@mymanova, a0, [],[],Aeq,Beq,lb,ub,[],opts);
        fprintf('Attempt %3d -- one-way ANOVA result: p=%8.4f\n', ii, p);
        if p<pbest
            pbest = p;
            abest = a;
        end
    end
    % run mymanova one more time to get more information about the best result
    [p,tbl,results,c1] = mymanova(abest);

    fprintf('Sign-constrained MANOVA result for c1: p = %8.4f\n',p);
    showtransform(abest, varnames);
    figure
    gscatter(group,c1,grouplabels,[],'oxs')
    title('Canonical variable c1 after sign-constrained MANOVA');
    
    % report the pairwise comparisons on c1
    for j = 1:size(results,1)
        group1 = Group{results(j,1)};
        group2 = Group{results(j,2)};
        fprintf('%6s vs. %6s: mean difference %8.3f (conf interval: %8.3f %8.3f) p=%8.3f\n', ...
            group1, group2, results(j,3:6))
    end


    % also do the PCA on the post-pre differences in the original variables and add those scores to the data matrix
    [coeff, score, latent]  = pca(Xrel);
    figure(1);
    subplot(2,2,1);
    plot(latent);
    xlabel('principal component');
    ylabel('eigenvalue');
    
    % add PC1 and PC2 scores to the data matrix
    headers = [headers 'PCdiff1' 'PCdiff2'];
    data_diff = [data_diff score(:,1:2)];
    
    % list of dependent variables we want to explore
    % for example: depvars = {'MRP','PC1'};
    depvars = headers(4:end);   % use all variables we have, original and PC
    
    % do 1-way ANOVA for each dependent variable in the list
    fprintf('=========================================================\n');
    fprintf('ONE-WAY ANOVA ON ORIGINAL VARIABLES AND PC1 AND PC2\n');
    for i = 1:numel(depvars)
        depvar = depvars{i};
        fprintf('===========================================================================\n');
        fprintf('Statistical analysis for variable: %s\n', depvar);
        
        % extract the relevant columns from data and make a (subj x group) matrix for ANOVA
        column = find(strcmp(headers,depvar));  % column with variable i in it
        for k = 1:3
           m(:,k) = data_diff(find(group==k), column);
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
            [results,means] = multcompare(stats,'Display','off','CType','bonferroni');
            for j = 1:size(results,1)
                group1 = Group{results(j,1)};
                group2 = Group{results(j,2)};
                fprintf('%6s vs. %6s: mean difference %8.3f (conf interval: %8.3f %8.3f) p=%8.3f\n', ...
                    group1, group2, results(j,3:6))
            end
        end
    end
end
%==============================================================
function [p,tbl,results,c1] = mymanova(a)
    global Xrel group Group
    
    c1 = Xrel * a;  % transform the data using coefficients a
    % put the data in a 4x3 matrix for anova1
    for k = 1:3
         m(:,k) = c1(find(group==k));
    end
    [p,tbl,stats] = anova1(m,Group,'off');
    results = multcompare(stats,'Display','off','CType','bonferroni');
    
    % to minimize p1*p2*p3, uncomment the following line
    p = prod(results(:,6));
end
%===================================================================
function showtransform(a, varnames)
    % show how c1 depends on the original variables v
    % c1 = a(1)*v(1) + a(2)*v(2) + ...
    nvar = numel(varnames);
    fprintf('c1 =');
    for i = 1:nvar
        if i==1
            fprintf(' %0.3f*%s', a(i), varnames{i});
        else
            fprintf(' %0.3f*%s', abs(a(i)), varnames{i});
        end
        if i<nvar
            if a(i+1)<0
                fprintf(' -');
            else
                fprintf(' +');
            end
        else
            fprintf('\n');
        end
    end
    fprintf('NOTE: before transform, mean was subtracted from each original variable.\n');
end