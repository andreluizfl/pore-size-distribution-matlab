function [C0, C1, Re] = poredistribution_yang_optimized(C)
% poredistribution_yang_optimized
% Compute critical radius map (C0), propagated radii (C1),
% and the histogram/distribution of radii (Re) from a binary 3D volume.
%
%   [C0, C1, Re] = poredistribution_yang_legacy_optimized(C)
%
% INPUTS:
%   C  - logical 3D array representing the pore mask (true => pore).
%
% OUTPUTS:
%   C0 - single array of the same size as C containing the "critical radius"
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
%     Calculated using 'single' precision to reduce memory footprint by 50%.
%   - For each integer radius s (descending), identify "centers" where
%     round(C0-0.5) == s. To ensure high performance and maximum legacy 
%     compatibility, propagation is done in a strict serial loop.
%     The distance transform (bwdist) is confined spatially using a 
%     local Bounding Box around the centers. The C1 label for those 
%     voxels is set to s+1 (uint16).
%   - Re is the histogram of C1 values using bins centered at integers.
%
% NOTES:
%   - Parallel execution (parpool/parfor) is deliberately excluded from this 
%     version because 'bwdist' is not thread-safe in older MATLAB versions 
%     and causes severe crashes when executed inside workers.
%   - This function expects C to be logical. Non-logical inputs will be
%     coerced but a warning will be issued to inform the user.
%   - The implementation contains a fallback for older MATLAB releases
%     (uses histc when histcounts is not available).
%
% EXAMPLE:
%   [C0, C1, Re] = poredistribution_yang_legacy_optimized(C);
%
% Author: André Luiz Ferraz Lourenço
% Date:   01/09/2026
%

% -------------------- Input validation & defaults -----------------------
if ~islogical(C)
    warning('Input C is not logical. It will be converted to logical for processing.');
    C = logical(C);
end

% Check compatibility for histcounts/histc fallback
isLegacy = verLessThan('matlab','9.1');
sz = size(C);

% -------------------- 1. Compute C0 via distance transform --------------
% Pad the volume to ensure border voxels have correct distances
pad = [1 1 1];
Cpad = padarray(C, pad, 0, 'both');

% Euclidean distance from background (~Cpad) within the padded volume
Dpad = bwdist(~Cpad, 'euclidean');

% Extract the original center volume and clear heavy paddings immediately
D_center = Dpad((1+pad(1)):(end-pad(1)), (1+pad(2)):(end-pad(2)), (1+pad(3)):(end-pad(3)));
clear Dpad Cpad; 

% OPTIMIZATION: Downcast to 4-byte float (single) to halve RAM usage
tol = single(1e-12);
C0 = zeros(sz, 'single'); 
C0(C) = ceil(D_center(C) - tol) - 0.5;
clear D_center; % Free double-precision distance matrix early

% -------------------- 2. Radius Pre-processing --------------------------
% We assign integer radius labels: r = round(C0 - 0.5)
C1 = zeros(sz, 'uint16');
r = round(C0 - 0.5);
r(~C) = -1; % Mark background voxels

% OPTIMIZATION: Downcast to 2-byte integer (int16) to free up memory
r = int16(r); 

r_values = unique(r(:));
r_values(r_values < 0) = [];
r_values = sort(r_values, 'descend'); % Process larger radii first

% -------------------- 3. Radius Propagation (Bounding Box) --------------
% Serial propagation: absolute stability with spatial confinement
for s = double(r_values(:)')
    centers = (r == s);
    if ~any(centers(:)), continue; end

    % Locate centers
    idx = find(centers);
    [y, x, z] = ind2sub(sz, idx);
    
    % Confinement limits (Bounding Box) to avoid calculating bwdist 
    % on the entire 3D volume. min/max prevents indexing out of bounds.
    minY = max(1, min(y)-s); maxY = min(sz(1), max(y)+s);
    minX = max(1, min(x)-s); maxX = min(sz(2), max(x)+s);
    minZ = max(1, min(z)-s); maxZ = min(sz(3), max(z)+s);

    % Extract local sub-volumes
    subCenters = centers(minY:maxY, minX:maxX, minZ:maxZ);
    subC1      = C1(minY:maxY, minX:maxX, minZ:maxZ);

    % bwdist is executed strictly within the tiny localized box (ultra-fast)
    subDist = bwdist(subCenters);

    % Assign radius (s+1) only to unassigned voxels (subC1 == 0)
    maskToAssign = (subDist <= s) & (subC1 == 0);
    subC1(maskToAssign) = uint16(s + 1);

    % Write the updated block back to the main C1 volume
    C1(minY:maxY, minX:maxX, minZ:maxZ) = subC1;
end

% -------------------- 4. Compute Re histogram (distribution) ------------
% Maximum radius (in half-step units was stored as k - 0.5)
dpm = double(max(C0(:))) - 0.5;
maxIndex = round(dpm + 1);

% Ensure histogram has at least 100 bins
Re = zeros([max(maxIndex, 100), 1], 'uint32');

vals = double(C1(:));
edges = 0.5:1:(maxIndex + 0.5);

if isLegacy
    counts = histc(vals, edges);
    counts = counts(1:end-1); % Remove the last catch-all bin edge
else
    counts = histcounts(vals, edges)';
end

% Safely assign counts into Re (avoid dimension mismatch)
nValid = min(length(counts), maxIndex);
Re(1:nValid) = uint32(counts(1:nValid));

end