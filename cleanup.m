function [d] = cleanup(d)
% interpolate missing markers when gaps are small
% gaps that are too large are replaced by NaN (not a number)
%
% Input:
%       d: data structure from importdata(mocapfile)
% Output:
%       d: data with gaps interpolated

    % replace zeros in mocapdata (missing markers) by NaN
    d.data( d.data == 0 ) = NaN;

    % first do gap filling on pelvis markers, using rigid body kinematics
%     col = contains(d.colheaders, {'RASIS','LASIS','RPSIS','LPSIS','SACR'} );
%     d.data(:,col) = rigidbodyfill(d.data(:,col));

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
            goodframes = find(~isnan(d.data(:,i)));
            
            % if we don't have 2 frames, make this marker NaN for the whole trial
            if numel(goodframes) < 2
                d.data(:,markercolumns) = NaN(nframes,3);
                % many TXT files contain a "S1" marker which is not real
                % and does not have data, no need to report it
                if ~contains(colname,'S1.Pos')
                    fprintf('%12s: COMPLETELY missing\n', colname(1:end-1));
                end
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
                d.data(interpframes,markercolumns) = interp1(times(goodframes), d.data(goodframes,i:i+2), times(interpframes), 'linear');
            end
        end
    end
end
%==================================================================
function d = rigidbodyfill(d)
% d is a N x 3M matrix with N frames of XYZ data for M markers
% try to fill the gaps using the assumption that all M markers are on the
% same rigid body
    N = size(d,1);
    M = size(d,2)/3;
    d1 = d;

    % find the frames where all markers are visible
    f = find(~isnan(sum(d,2)));
    
    % use the mean of these frames to define (arbitrary) local coordinates of all markers
    xyzloc = d(f(1),:);
    xyzloc = reshape(xyzloc,3,M)';
    
    % calculate local coordinates for the markers from all complete frames
    newxyzloc = [];
    for i = 1:numel(f)
        frame = f(i);
        xyz = reshape(d(frame,:),3,M)';
        [R,t,rms] = soder(xyz, xyzloc);
        if rms < 0.003
            newxyzloc = [newxyzloc ; reshape(R*xyz' + t,1,3*M)];
        end
    end
    xyzloc = reshape(mean(newxyzloc),3,M)';  % now we have more accurate local coordinates
    
    % go through all frames where something was missing
    f = find(isnan(sum(d,2)));
    for frame = 1:numel(f)
        i = f(frame);
        xyz = reshape(d(i,:), 3, M)';
        vis = find(~isnan(xyz(:,1)));  % which markers are visible
        SEMmin = 1e10;
        S = subsets(vis, 3);  % find subsets of vis with at least 3 markers
        % find the subset that gives the best fit
        for k = 1:numel(S)
            set = S{k};
            [R,t,rms] = soder(xyzloc(set,:), xyz(set,:));
            SEM = 3*rms / sqrt(3*numel(set)-6);  % unbiased SEM
            if SEM < SEMmin
                SEMmin = SEM;
                setmin = set;
                Rmin = R;
                tmin = t;
            end
        end
        % reconstruct the markers that are not in the subset
        if i == 4860
            keyboard
        end
        m = setdiff(1:M, setmin);
        for j = 1:numel(m)
            k = m(j);
            col = 3*(k-1) + (1:3);  % columns where XYZ of this marker are
            d(i,col) = ( Rmin*xyzloc(k,:)' + tmin )';
        end

    end
end
%======================================================================
function[R,d,rms] = soder(P1,P2)
% Estimates rigid body pose R,d from marker coordinates
%
% Input:
%	P1		(N x 3 matrix) coordinates of N markers in coordinate system 1
%	P2		(N x 3 matrix) coordinates of the same markers, now measured in coordinate system 2
%
% Output:
%	R		(3x3 matrix) 
%   d       (3x1 matrix) such that P2 = R*P1 + d)
%	rms		(scalar, optional) rms fit error of the rigid body model: P2 = T*P1 + error
%
% Method:
%	Soderkvist I. and Wedin P.A. (1993) Determining the movements of the skeleton using
%	well-configured markers.  Journal of Biomechanics, 26:1473-1477.    

	% error checking
	[nmarkers,ndim1]=size(P1);
	[nmarkers2,ndim2]=size(P2);
	if (ndim1 ~= 3) || (ndim2 ~= 3)
		error('hmc3_soder: matrix with marker coordinates must have 3 columns (x,y,z)');
	end
	if (nmarkers2 ~= nmarkers)
		error('hmc3+soder: matrices with marker coordinates must have same number of rows');
	end
		
	% construct matrices A and B, subtract the mean so there is only rotation
	m1=mean(P1);
	m2=mean(P2);
	for i=1:nmarkers
	  A(i,:)=P1(i,:)-m1;
	  B(i,:)=P2(i,:)-m2;
	end
	A = A';
	B = B';

	% The singular value decomposition to calculate rotation matrix R
	C=B*A';
	[P,T,Q]=svd(C);
	R=P*diag([1 1 det(P*Q')])*Q';

	% Calculate the translation vector from the centroid of all markers
	d=m2'-R*m1';

	% calculate RMS value of residuals
	sumsq = 0;
	for i=1:nmarkers
	  P2model = R*P1(i,:)' + d;
	  sumsq = sumsq + norm(P2model-P2(i,:)')^2;
	end
	rms = sqrt(sumsq/3/nmarkers);
end
%=======================================================================
function S = subsets(a,n)
    % find all subsets of a with at least n elements
    m = numel(a);
    S = { };
    for setsize = n:m
        s = nchoosek(a,setsize);
        for k = 1:size(s,1)
           S = [S s(k,:)];
        end
    end
end
