function manovatest
% measure the type 1 error rate for Matlab MANOVA
% methods inspired by "MANOVA: Type I Error Rate Analysis" by CDW Ling
    close all
    ng = 3;                 % number of groups
    nvar = 8;               
    
    % construct the covariance matrix
    corr = 0.5;             % off diagonal elements in correlation matrix
    S = ones(nvar,1);       % standard deviations
    R = corr + zeros(nvar,nvar);
    for i=1:nvar
        R(i,i) = 1.0;
    end
    cov = S'*S*R;

    [n,nvar] = meshgrid(4:10, 2:8);
    nsig = zeros(size(n));
    for i = 1:size(n,1)
        for j = 1:size(n,2)
            n1 = n(i,j);
            nvar1 = nvar(i,j);
            nsig(i,j) = dotest(n1,ng,nvar1,cov(1:nvar1,1:nvar1));
        end
    end
    surface(n,nvar,nsig);
    view(3);
    xlabel('number of participants per group');
    ylabel('number of variables')
    title('empirical MANOVA type 1 error rate');
end
%==============================================================
function [type1rate] = dotest(n, ng, nvar, cov)
    % inputs:
    %       n       number of participants in each group
    %       ng      number of groups
    %       nvar    number of variables measured
    %       cov     covariance matrix for simulated data
    %
    % outputs:
    %       type1rate    type1 error measured with MANOVA Wilks test at p<0.05
    nsubj = n*ng;           % total number of subjects
    if nargin<4
        cov = eye(nvar);    % default covariance is identity matrix
    end
    mu = zeros(1,nvar);  % mean will be zero for all subjects and groups

    % generate group labels for each subject
    subj = 0;
    for i = 1:ng
        for j = 1:n
            subj = subj+1;
            group(subj) = i;   % group number for subj
        end
    end

    nrep = 5000;  % number of times the experiment is done
    nsig = 0;     % here we count how many times MANOVA said significant
    for k = 1:nrep
        % generate random data, assuming the null hypothesis 
        % (all subjects have the same probability distribution)
        d = mvnrnd(mu,cov,nsubj);

        % do the MANOVA
        [~,p] = manova1(d,group);
        if p(1)<0.008
            nsig = nsig+1;
        end
    end
    type1rate = nsig/nrep;
    fprintf('n=%d  ng=%d  nvar=%d  empirical type 1 error rate: %8.4f\n', n, ng, nvar, type1rate);
    pause(0.1);
end