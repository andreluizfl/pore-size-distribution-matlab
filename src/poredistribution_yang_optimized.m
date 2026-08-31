function [C0, C1, Re] = poredistribution_yang_merged_v2(C, useParallel)
% POREDISTRIBUTION_YANG_OPTIMIZED
% Compute critical radius map (C0), propagated radii (C1), and distribution (Re).
%
% CORE OPTIMIZATIONS IN THIS VERSION:
% 1. Serial by Default: Avoids worker overhead; maximizes legacy compatibility.
% 2. Bounding Box Confinement: Both serial and parallel modes calculate
%    distances strictly within localized sub-volumes.
% 3. Ultra-Low Memory Parallelism: Workers do NOT allocate full-size temporary
%    volumes. They return only localized bounding box results, which are
%    reduced into the main C1 matrix post-loop, avoiding RAM saturation.
% 4. Data Downcasting: Aggressive use of 'single' and 'int16' cuts memory by up to 75%.

% --- Input Validation & Defaults ---
if nargin < 2 || isempty(useParallel)
    useParallel = false; % Default to serial for absolute stability and low memory footprint
end

if ~islogical(C)
    warning('Input C is not logical. It will be converted to logical.');
    C = logical(C);
end

isLegacy = verLessThan('matlab','9.1');
sz = size(C);

% --- Parallel Verification ---
if useParallel
    try
        hasParallel = license('test','Distrib_Computing_Toolbox');
    catch
        hasParallel = false;
    end
    if ~hasParallel, useParallel = false; end
end

% --- 1. Compute C0 via Distance Transform ---
pad = [1 1 1];
Cpad = padarray(C, pad, 0, 'both');
Dpad = bwdist(~Cpad, 'euclidean');

% Extract the original center volume and clear heavy paddings immediately
D_center = Dpad((1+pad(1)):(end-pad(1)), (1+pad(2)):(end-pad(2)), (1+pad(3)):(end-pad(3)));
clear Dpad Cpad; 

tol = single(1e-12);
C0 = zeros(sz, 'single'); % Downcast to 4-byte float (single) to halve RAM
C0(C) = ceil(D_center(C) - tol) - 0.5;
clear D_center; % Free double-precision distance matrix

% --- 2. Radius Pre-processing ---
C1 = zeros(sz, 'uint16');
r = round(C0 - 0.5);
r(~C) = -1; % Mark background voxels
r = int16(r); % Downcast to 2-byte integer to free up memory

r_values = unique(r(:));
r_values(r_values < 0) = [];
r_values = sort(r_values, 'descend'); % Process largest pores first to respect geometry
K = numel(r_values);

% --- 3. Radius Propagation (Bounding Box Strategy) ---
if ~useParallel
    % SERIAL MODE: Strict overwriting control, zero worker overhead.
    for kk = 1:K
        s = double(r_values(kk));
        centers = (r == s);
        if ~any(centers(:)), continue; end

        % Locate centers and define Bounding Box limits to confine bwdist
        idx = find(centers);
        [y, x, z] = ind2sub(sz, idx);
        minY = max(1, min(y)-s); maxY = min(sz(1), max(y)+s);
        minX = max(1, min(x)-s); maxX = min(sz(2), max(x)+s);
        minZ = max(1, min(z)-s); maxZ = min(sz(3), max(z)+s);

        % Extract local sub-volumes
        subCenters = centers(minY:maxY, minX:maxX, minZ:maxZ);
        subC1      = C1(minY:maxY, minX:maxX, minZ:maxZ);

        % Calculate distance strictly within the localized box
        subDist = bwdist(subCenters);

        % Assign radius (s+1) only to unassigned voxels (subC1 == 0)
        maskToAssign = (subDist <= s) & (subC1 == 0);
        subC1(maskToAssign) = uint16(s + 1);

        % Write updated block back to the main volume
        C1(minY:maxY, minX:maxX, minZ:maxZ) = subC1;
    end
else
    % PARALLEL MODE: Ultra-low memory implementation. Threads return only small crops.
    Cand = cell(K, 1);
    parfor kk = 1:K
        s = double(r_values(kk));
        centers = (r == s);
        if ~any(centers(:)), continue; end

        idx = find(centers);
        [y, x, z] = ind2sub(sz, idx);
        minY = max(1, min(y)-s); maxY = min(sz(1), max(y)+s);
        minX = max(1, min(x)-s); maxX = min(sz(2), max(x)+s);
        minZ = max(1, min(z)-s); maxZ = min(sz(3), max(z)+s);

        subCenters = centers(minY:maxY, minX:maxX, minZ:maxZ);
        subDist = bwdist(subCenters);

        % Allocate subResult ONLY for the local bounding box size, avoiding OOM errors
        maskDil = (subDist <= s);
        subResult = zeros(size(subCenters), 'uint16');
        subResult(maskDil) = uint16(s + 1);

        % Store only coordinates and localized data
        Cand{kk} = {minY, maxY, minX, maxX, minZ, maxZ, subResult};
    end

    % Localized Serial Reduction: Merge threaded results back into the main C1 volume
    for kk = 1:K
        if isempty(Cand{kk}), continue; end
        item = Cand{kk};
        minY = item{1}; maxY = item{2};
        minX = item{3}; maxX = item{4};
        minZ = item{5}; maxZ = item{6};
        subResult = item{7};

        % Merge overlapping larger radii using max() locally
        subC1 = C1(minY:maxY, minX:maxX, minZ:maxZ);
        subC1 = max(subC1, subResult);
        C1(minY:maxY, minX:maxX, minZ:maxZ) = subC1;

        Cand{kk} = []; % Free cell immediately to keep RAM profile flat
    end
    clear Cand;
end

% --- 4. Histogram Generation (Re) ---
dpm = double(max(C0(:))) - 0.5;
maxIndex = round(dpm + 1);
Re = zeros([max(maxIndex, 100), 1], 'uint32');

vals = double(C1(:));
edges = 0.5:1:(maxIndex + 0.5);

if isLegacy
    counts = histc(vals, edges);
    counts = counts(1:end-1);
else
    counts = histcounts(vals, edges)';
end

nValid = min(length(counts), maxIndex);
Re(1:nValid) = uint32(counts(1:nValid));
end