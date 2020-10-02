function [colnum] = findcolumn(data, label)
    colnum = find(strcmp(data.colheaders, label));
    if isempty(colnum)
        warning('data file did not have a column with label: %s\n', label)
        return
    end
end
