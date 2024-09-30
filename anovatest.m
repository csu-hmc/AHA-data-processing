% measure the type 1 error rate for Matlab one-way ANOVA
close all
n = 4;                  % number of participants in each group
ng = 3;                 % number of groups
nsubj = n*ng;           % total number of subjects

% generate group labels for each subject
subj = 0;
for i = 1:ng
    for j = 1:n
        subj = subj+1;
        group(subj) = i;   % group number for subj
    end
end

nrep = 10000;
nsig = 0;
for k = 1:nrep
    % generate random data, assuming the null hypothesis (group means are all zero)
    d = randn(nsubj);  

    % do the ANOVA
    [p,tbl,stats] = anova1(d,group,'off');
    fprintf('ANOVA result: p = %8.4f\n',p);
    if p(1)<0.05
        nsig = nsig+1;
    end
end
fprintf('%d out of %d were significant\n', nsig, nrep);