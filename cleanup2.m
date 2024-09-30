function [d] = cleanup2(d)
% use PCA to reconstruct missing marker data
%
% Input:
%       d: data structure from importdata(mocapfile)
% Output:
%       d: data with gaps interpolated

    % find the column numbers of the marker data
    columns = find(contains(d.colheaders,'.Pos'))';

    % replace all zeros (missing data) by NaN 
    for i=1:numel(columns)
        d.data(d.data(:,columns(i))==0, columns(i)) = NaN;
    end

    % ignore those with more than 50% NaNs
    remove = [];
    nsamples = size(d.data,1);
    for i=1:numel(columns)
        perc_missing = 100*sum(isnan(d.data(:,columns(i)))) / nsamples;
        if perc_missing > 50
            remove = [remove i];  % add i to the remove list            
            % generate a message once for each removed marker (but not S1,
            % it is not an actual marker)
            column_name = d.colheaders{columns(i)};
            if contains(column_name,'.PosX') && ~contains(column_name,'S1.Pos')
                marker_name = strrep(column_name,'.PosX','');
                fprintf('gap filling will ignore %s, it has %.1f%% missing.\n', marker_name, perc_missing);
            end
        end
    end
    columns = setdiff(columns,columns(remove));  % these columns are kept
    
    % estimate mean and covariance of marker data
    [mu,C] = ecmnmle_hash(d.data(:,columns));
    Cinv = inv(C);
    
    % use Mahalanobis distance to remove outliers in any marker data
    disp('Removing outlier data...');
    newd = d;
    MD = zeros(nsamples,1);
    for i = 1:nsamples
        x = d.data(i,columns)';
        MD(i) = x' * Cinv * x;
    end
    plot(MD)
    
    
    % impute missing data using idea from Rasmussen 2020
    % for each frame, we reconstruct the missing markers by finding
    % coordinates that put this frame closest (using Mahalanobis Distance, MD)
    % to the mean (mu).
    % this is an optimization problem with an analytical solution
    % partition the inverse covariance matrix into parts corresponding to
    % observed coordinates (o) and missing coordinates (m)
    % Cinv = [Coo Com]  this is a symmetrical matrix, so Com = Cmo'
    %        [Cmo Cmm]
    % we minimize MD=(x-mu)'*Cinv*(x-mu), subject to x(o) = xo (observed data)
    % the solution is: x(m) = mu(m) - inv(Cmm)*Cmo*(xo-mu(o))
    disp('Estimating missing data...');
    for i = 1:nsamples
        x = newd.data(i,columns)';
        m = find(isnan(x));  % elements of x that are missing data
        if ~isempty(m)
            o = find(~isnan(x)); % elements of x that are observed
            Cmm = Cinv(m,m);
            Cmo = Cinv(m,o);
            x(m) = mu(m) - inv(Cmm)*Cmo*(x(o)-mu(o));
        end
        newd.data(i,columns) = x';      % store in newd
    end

    % plot newdata compared to original data
    if (1)
        clf
        for j = 1:numel(columns)
            i = columns(j);
            plot(d.data(:,1),[newd.data(:,i) d.data(:,i)]);
            legend('new data','original data');
            title(d.colheaders{i})
            disp('Hit ENTER to continue...');
            pause
        end
    end
end
%==================================================
function [mu,C] = ecmnmle_hash(data)
% estimate the mean and covariance of a data matrix
% use a hash method to see if we already have done this before
    % Convert data into a byte array called B...
    B = typecast(data(:),'uint8');
    if (isempty(B))
        error('ecmnmle_hash: B is empty, you may not have any marker data');
    end
    % Or, as suggested, using the undocumented function getByteStreamFromArray:
    % B = getByteStreamFromArray(data);

    % Create an instance of a Java MessageDigest with the desired algorithm:
    md = java.security.MessageDigest.getInstance('SHA-1');
    md.update(B);

    % Properly format the computed hash as an hexadecimal string:
    hash = reshape(dec2hex(typecast(md.digest(),'uint8'))',1,[]);
    
    % Construct the filename and see if the file exists
    filename = ['cache/' hash '.mat'];
    if exist(filename)
        % we have seen this data before, so simply reload the existing mean
        % and covariance
        disp('Loading mean and covariance data from cache...');
        load(filename);
    else
        % we have not seen this data before, so we need to estimate
        % the mean and covariance
        disp('Estimating mean and covariance of marker data...');
        disp('  (can take up to 10 min)');
        [mu,C] = ecmnmle(data);
        save(filename,'mu','C');  % save the result for later
    end
end
