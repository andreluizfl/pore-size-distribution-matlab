function [C0, C1, Re] = poredistribution_yang_optimized(C)
% poredistribution_yang_optimized
% Compute critical radius map (C0), propagated radii (C1),
% and the histogram/distribution of radii (Re) from a binary 3D volume.
% This version features memory optimizations and a fast bounding-box 
% approach for radius propagation.
%
%   [C0, C1, Re] = poredistribution_yang_optimized(C)
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
%     regions: for pore voxels, C0 = ceil(distance - tol) - 0.5. Uses single
%     precision to halve memory usage compared to double.
%   - For each integer radius s (descending), identify "centers" where
%     round(C0-0.5) == s. Instead of full-volume dilation, it extracts a 
%     local bounding box sub-volume around the centers, dilates by s 
%     (using bwdist locally), and updates C1. This drastically improves speed.
%   - Re is the histogram of C1 values using bins centered at integers.
%
% NOTES:
%   - Memory optimized: uses single/int16 casting and clears large 
%     temporary variables immediately.
%   - This function expects C to be logical. Non-logical inputs will be
%     coerced but a warning will be issued.
%   - Includes a fallback for older MATLAB releases (uses histc instead of histcounts).
%   - Note: Legacy 'useParallel' checks remain in code for backward 
%     compatibility but it is no longer used in the main logic or signature.

% -------------------- Input validation & defaults -----------------------
if ~islogical(C)
    warning('Input C is not logical. It will be converted to logical for processing.');
    C = logical(C);
end

% -------------------- Version detection --------------------
% Checks if the MATLAB version is older than 9.1 (R2016b) to handle syntax compatibility
isLegacy = verLessThan('matlab', '9.1'); 

% -------------------- Compute C0 via distance transform -----------------
% Pad the array with 0s to prevent edge effects during the distance transform
pad = [1 1 1];
Cpad = padarray(C, pad, 0, 'both');

% Compute the Euclidean distance transform of the inverted padded matrix
% bwdist calculates the distance from every 0-voxel to the nearest non-zero voxel
Dpad = bwdist(~Cpad, 'euclidean');

% Remove the padding to get the distance values for the original center volume
D_center = Dpad( (1+pad(1)):(end-pad(1)), (1+pad(2)):(end-pad(2)), (1+pad(3)):(end-pad(3)) );

% OPTIMIZATION 2: Clear heavy variables from RAM immediately to free up memory
clear Dpad Cpad; 

tol = single(1e-12); % Tolerance for floating-point comparisons
% OPTIMIZATION 1: Downcast to float32 (single) to halve RAM usage compared to double
C0 = zeros(size(C), 'single'); 
C0(C) = ceil(D_center(C) - tol) - 0.5;

% OPTIMIZATION 2: Clear memory immediately
clear D_center; 

% -------------------- Compute C1 (propagated radii) ---------------------
C1 = zeros(size(C), 'uint16');
r = round(C0 - 0.5);       % 'r' inherits the 'single' data type from C0
r(~C) = -1;                

% OPTIMIZATION 1: Cast to a 2-byte integer (int16) to use 1/4 of the original RAM
r = int16(r);              

% Find all unique radius values present in the volume
r_values = unique(r(:));
r_values(r_values < 0) = [];    
r_values = sort(r_values, 'descend');  % Process from largest to smallest pores

% -------------------- Radius propagation (Bounding Box) -----------------
% Store the original volume size for index calculations
sz = size(C); 
for s = double(r_values(:)')
    
    centers = (r == s); % Logical mask of centers with radius 's'
    if ~any(centers(:))
        continue;
    end
    
    % Find the exact linear indices of the centers
    idx = find(centers);
    % Convert linear indices to X, Y, Z subscripts (coordinates)
    [y, x, z] = ind2sub(sz, idx);
    
    % Define the bounding box limits by expanding by the distance 's'
    % Using max() and min() ensures we do not attempt to crop outside the matrix boundaries
    minY = max(1, min(y) - s);
    maxY = min(sz(1), max(y) + s);
    minX = max(1, min(x) - s);
    maxX = min(sz(2), max(x) + s);
    minZ = max(1, min(z) - s);
    maxZ = min(sz(3), max(z) + s);
    
    % Extract only the sub-volume containing the centers (and the margin 's')
    subCenters = centers(minY:maxY, minX:maxX, minZ:maxZ);
    subC1      = C1(minY:maxY, minX:maxX, minZ:maxZ);
    
    % Calculate bwdist ONLY within the sub-volume (significantly faster than processing the whole volume)
    subDist = bwdist(subCenters);
    
    % Update only the voxels in the sub-volume that are within radius 's' and haven't been assigned yet (subC1 == 0)
    % This ensures smaller pores do not overwrite larger ones
    maskToAssign = (subDist <= s) & (subC1 == 0);
    subC1(maskToAssign) = uint16(s + 1);
    
    % Return the updated sub-volume to the original C1 matrix
    C1(minY:maxY, minX:maxX, minZ:maxZ) = subC1;
    
end

% -------------------- Compute Re histogram (distribution) ----------------
dpm = double(max(C0(:))) - 0.5;
maxIndex = round(dpm + 1);
Re = zeros( [max(maxIndex,100),1],'uint32');

vals = double(C1(:)); % Flatten the C1 matrix for histogram counting
edges = 0.5:1:(maxIndex + 0.5); % Define bin edges for the histogram

if isLegacy
    % Use histc for older MATLAB versions (deprecated in newer releases)
    counts = histc(vals, edges);
    counts = counts(1:end-1); % Remove the last edge catch-all bin
else
    % Use the modern histcounts function for newer MATLAB versions
    counts = histcounts(vals, edges)';
end

nValid = min(length(counts), maxIndex);
Re(1:nValid) = uint32(counts(1:nValid));
end