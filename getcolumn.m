function d = getcolumn(data,labels)
    % extract one or more columns of data
    col = contains(data.colheaders, labels);
    if isempty(col)
        warning('data file did not have any column with label: %s\n', label)
        return
    end
    d = data.data(:,col);
end