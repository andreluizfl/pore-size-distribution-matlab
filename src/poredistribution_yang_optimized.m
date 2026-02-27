function [C0, C1, Re] = poredistribution_yang_optimized(C, useParallel)
% poredistribution_yang_optimized
% Compute critical radius map (C0), propagated radii (C1),
% and the histogram/distribution of radii (Re) from a binary 3D volume.
%
%   [C0, C1, Re] = poredistribution_yang_optimized(C)
%   [C0, C1, Re] = poredistribution_yang_optimized(C, useParallel)
%
% INPUTS:
%   C           - logical 3D array representing the pore mask (true => pore).
%   useParallel - (logical, optional) request parallel processing for the
%                 radius propagation step. If omitted or empty, true is
%                 assumed and the function attempts to use a parpool.
%
% OUTPUTS:
%   C0 - double array of the same size as C containing the "critical radius"
%        measure for pore voxels. Values are of the form (k - 0.5) where k is
%        a positive integer; non-pore voxels are zero.
%   C1 - uint16 array of the same size as C containing propagated radius
%        labels. Each value represents the assigned radius bin (1..N), 0 if none.
%   Re - uint32 column vector containing the histogram counts for each
%        integer radius value. Index k contains the count for radius k.
%
% METHOD SUMMARY:
%   - C0 is derived from the Euclidean distance transform inside pore
%     regions: for pore voxels, C0 = ceil(distance - tol) - 0.5.
%   - For each integer radius s (descending), identify "centers" where
%     round(C0-0.5) == s and dilate them by s (using bwdist). The C1 label
%     for those voxels is set to s+1 (uint16). The parallel version
%     computes candidate volumes per radius and reduces by element-wise max.
%   - Re is the histogram of C1 values using bins centered at integers.
%
% NOTES:
%   - This function expects C to be logical. Non-logical inputs will be
%     coerced but a warning will be issued to inform the user.
%   - The implementation contains a fallback for older MATLAB releases
%     (uses histc when histcounts is not available).
%
% EXAMPLE:
%   [C0, C1, Re] = poredistribution_yang_optimized(C, true);
%
% Author: adapted and modularized
% Date:   (no date hard-coded here)
%

% -------------------- Input validation & defaults -----------------------
if nargin < 2 || isempty(useParallel)
    useParallel = true;
end

if ~islogical(C)
    warning('Input C is not logical. It will be converted to logical for processing.');
    C = logical(C);
end

% -------------------- Version detection --------------------
isLegacy = verLessThan('matlab', '9.1'); % compatibility for histcounts/histc

% -------------------- Parallel capability --------------------
if useParallel
    try
        hasParallel = license('test', 'Distrib_Computing_Toolbox');
        if hasParallel
            pool = gcp('nocreate');  % Verifica se já existe uma pool
            if isempty(pool)
				parpool('local');
                pool = gcp('nocreate');  % Atualiza a referência da pool
            end
            useParallel = ~isempty(pool);
        else
            useParallel = false;
        end
    catch
        useParallel = false;
    end
end

% -------------------- Compute C0 via distance transform -----------------
% Pad the volume to ensure border voxels have correct distances
pad = [1 1 1];
Cpad = padarray(C, pad, 0, 'both');
% Euclidean distance from background (~Cpad) within the padded volume
Dpad = bwdist(~Cpad, 'euclidean');
% Extract the central region corresponding to original volume
D_center = Dpad( (1+pad(1)):(end-pad(1)), (1+pad(2)):(end-pad(2)), (1+pad(3)):(end-pad(3)) );

tol = 1e-12;
C0 = zeros(size(C));        % default 0 for background voxels
C0(C) = ceil(D_center(C) - tol) - 0.5;


% -------------------- Compute C1 (propagated radii) ---------------------
% We assign integer radius labels: r = round(C0 - 0.5)
C1 = zeros(size(C), 'uint16');
r = double(C0) - 0.5;
r(~C) = -1;                % background marker
r = round(r);              % integer radii

r_values = unique(r(:));
r_values(r_values < 0) = [];    % drop background and negative entries
r_values = sort(r_values, 'descend');  % process larger radii first

if useParallel
    % Build candidate volumes for each radius in parallel, then reduce via max.
    K = numel(r_values);
    Cand = cell(K, 1);
    parfor kk = 1:K
        s = double(r_values(kk));
        if s < 0
            Cand{kk} = zeros(size(C), 'uint16');
            continue;
        end
        centers = (r == s);
        if ~any(centers(:))
            Cand{kk} = zeros(size(C), 'uint16');
            continue;
        end
        % distance (euclidean) from centers
        distToCenters = bwdist(centers);
        maskDil = (distToCenters <= s);
        tmp = zeros(size(C), 'uint16');
        % use label = s + 1 so that zero remains "unassigned"
        tmp(maskDil) = uint16(s + 1);
        Cand{kk} = tmp;
    end
    % Combine candidate volumes by taking the maximum per voxel
    for kk = 1:K
        C1 = max(C1, Cand{kk});
    end
else
    % Serial propagation: assign directly
    for s = r_values(:)'
        centers = (r == s);
        if ~any(centers(:)), continue; end
        distToCenters = bwdist(centers);
        maskDil = (distToCenters <= s);
        % assign radius label (s + 1) wherever mask applies (overwrites previous)
        % using max semantics could be chosen; here we set directly (larger s first)
        toAssign = maskDil & (C1 == 0);
        if any(toAssign(:))
            C1(toAssign) = uint16(s + 1);
        end
    end
end

% -------------------- Compute Re histogram (distribution) ----------------
% Maximum radius (in half-step units was stored as k - 0.5). Convert to integer bins.
dpm = max(C0(:)) - 0.5;
maxIndex = round(dpm + 1);
% Ensure histogram has at least 100 bins
Re = zeros( [max(maxIndex,100),1],'uint32');
vals = double(C1(:));
edges = 0.5:1:(maxIndex + 0.5);
% Central points of each bin (Centers are integers 1 to maxIndex)
% binCenters = edges(1:end-1) + diff(edges)/2
if isLegacy
    counts = histc(vals, edges);
    counts = counts(1:end-1);
else
    counts = histcounts(vals, edges)';
end
% Safely assign counts into Re (avoid dimension mismatch)
nValid = min(length(counts), maxIndex);
Re(1:nValid) = uint32(counts(1:nValid));

end
