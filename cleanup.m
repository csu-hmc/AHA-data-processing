function [d] = cleanup(d)
% interpolate missing markers when gaps are small
% gaps that are too large are replaced by NaN (not a number)
%
% Input:
%       d: data structure from importdata(mocapfile)
% Output:
%       d: data with gaps interpolated

    % we only interpolate when less than maxgap successive frames are
    % missing
    maxgap = 10;

    times = d.data(:,1);   % time values are in column 1
    [nframes,ncolumns] = size(d.data);
    for i = 1:ncolumns
        % see if it a marker channel
        colname = d.colheaders{i};
        if findstr(colname, '.PosX')
            markercolumns = i:i+2;  % PosX, PosY, PosZ
            % find the frames with good data (value not zero)
            goodframes = find(d.data(:,i));
            
            % if we don't have 2 frames, make this marker NaN for the whole trial
            if numel(goodframes) < 2
                d.data(:,markercolumns) = NaN(nframes,3);
                fprintf('%12s: COMPLETELY missing\n', colname(1:end-1));
            else
                % find large gaps, they should not be interpolated but replaced by NaN
                interpframes = setdiff(1:nframes, goodframes);   % normally, we interpolate all frames that are missing
                largegaps = find(diff(goodframes)>maxgap);
                for j = 1:numel(largegaps)
                    gapframes = goodframes(largegaps(j))+1:goodframes(largegaps(j)+1)-1;  % the frame numbers inside the gap
                    d.data(gapframes,markercolumns) = NaN;
                    interpframes = setdiff(interpframes, gapframes);  % remove these frames from the interpolation
                end
                
                % print some information
                fractionmissing = (nframes-numel(goodframes))/nframes;
                if (fractionmissing > 0)
                    fractioninterpolated = numel(interpframes)/nframes;
                    fprintf('%12s: %5.2f%% missing -- %5.2f%% interpolated\n', colname(1:end-1), 100*fractionmissing, 100*fractioninterpolated);
                end
                
                % interpolate across the gaps that are not too large
                d.data(interpframes,markercolumns) = interp1(times(goodframes), d.data(goodframes,i:i+2), times(interpframes), 'linear', 'extrap');
            end
        end
    end
end