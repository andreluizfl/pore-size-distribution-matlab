function C = load_volume(imgDir, extType, volumetricSize, binThreshold, useParallel)
% LOAD_VOLUME  Read an image stack and produce a binary 3D volume.
%
%   C = LOAD_VOLUME(imgDir)
%   C = LOAD_VOLUME(imgDir, extType, volumetricSize, binThreshold, useParallel)
%
% INPUTS:
%   imgDir         - (char) path to the folder containing the image stack.
%   extType        - (char, optional) image extension to search for,
%                     e.g., '.bmp' (default), '.tif', '.png'.
%   volumetricSize - (optional) defines the sub-volume to read. Accepts:
%                     [] (default)         -> full image extents and all slices
%                     string or char 'fit' -> automatically fits a centered cubic volume 
%                     scalar N             -> NxNxN volume starting at (1,1,1)
%                     3x2 matrix           -> explicit ranges
%   binThreshold   - (numeric, optional) binarization threshold.
%                     -1 (default) => compute global Otsu threshold via sampling.
%                     Otherwise, use the provided numeric threshold (in [0,1]).
%   useParallel    - (logical, optional) request parallel processing.
%
% OUTPUT:
%   C - logical 3D array (rows x cols x slices) containing the binary pore mask.
%
% APPLIED OPTIMIZATIONS:
%   - Removed single-precision floating-point matrix allocation, saving memory.
%   - Stratified sampling for global Otsu threshold calculation.
%   - Read and binarize steps merged into a single I/O loop.
%   - Intelligent bypass detection for images that are natively logical (binary).
%   - Thread-safe I/O check: Preemptively prioritizes 'local' pool for TIFF/GIF formats.
%   - Enforced type consistency before cell2mat concatenation.

% -------------------- Defaults and input normalization --------------------
if nargin < 2 || isempty(extType)
    extType = '.*';
end
if nargin < 3
    volumetricSize = [];
end
if nargin < 4 || isempty(binThreshold)
    binThreshold = -1;
end
if nargin < 5 || isempty(useParallel)
    useParallel = true;
end

% Ensure imgDir ends with filesep for safe concatenation
if isempty(imgDir)
    error('imgDir must be a non-empty directory path.');
end
if imgDir(end) ~= filesep
    imgDir = [imgDir filesep];
end

% -------------------- Find and sort files --------------------------------
files = dir([imgDir '*' extType]);
files = files(~[files.isdir]);  % Remove directories

if isempty(files)
    error('No image files found in the specified directory: %s', imgDir);
end

% Numeric-aware sorting: extract trailing number groups and sort by the last group
fileNumbers = zeros(1, numel(files));
for ff = 1:numel(files)
    nums = regexp(files(ff).name, '\d+', 'match');
    if isempty(nums)
        fileNumbers(ff) = 0;
    else
        fileNumbers(ff) = str2double(nums{end});
    end
end
[~, sortIdx] = sort(fileNumbers);
files = files(sortIdx);
numImages = numel(files);

% -------------------- Validate consistent image geometry -----------------
firstInfo = imfinfo([imgDir files(1).name]);
origRows = firstInfo.Height;
origCols = firstInfo.Width;

for ff = 2:numel(files)
    infoFf = imfinfo([imgDir files(ff).name]);
    if infoFf.Height ~= origRows || infoFf.Width ~= origCols
        error('All images must have identical dimensions. File "%s" differs.', files(ff).name);
    end
end

% -------------------- Determine volumetric ranges ------------------------
if isempty(volumetricSize)
    rangeX = [1 origCols];
    rangeY = [1 origRows];
    rangeZ = [1 numImages];
elseif isstring(volumetricSize) || ischar(volumetricSize)
    if  string(volumetricSize) == "fit"
        N = min([origRows, origCols, numImages]);
        diffX = round((origCols - N) / 2);
        diffY = round((origRows - N) / 2);
        diffZ = round((numImages - N) / 2);
        rangeX = [diffX + 1, diffX + N];
        rangeY = [diffY + 1, diffY + N];
        rangeZ = [diffZ + 1, diffZ + N];
    else
        error("Invalid volumetricSize name. Must be 'fit'.");
    end
elseif isscalar(volumetricSize)
    N = round(volumetricSize);
    if N <= 0, error('Scalar volumetricSize must be positive.'); end
    rangeX = [1 N]; rangeY = [1 N]; rangeZ = [1 N];
elseif ismatrix(volumetricSize) && all(size(volumetricSize) == [3 2])
    rangeX = volumetricSize(1, :);
    rangeY = volumetricSize(2, :);
    rangeZ = volumetricSize(3, :);
else
    error('Invalid volumetricSize. Provide [] | scalar | 3x2 matrix.');
end

% Clamp user ranges to available image dimensions / slice count
rangeX(2) = min(rangeX(2), origCols);
rangeY(2) = min(rangeY(2), origRows);
rangeZ(2) = min(rangeZ(2), numImages);

rowsRange = rangeY(1):rangeY(2);
colsRange = rangeX(1):rangeX(2);
imgsRange = rangeZ(1):rangeZ(2);

numRows = numel(rowsRange);
numCols = numel(colsRange);
numImgs = numel(imgsRange);

fprintf('Using volume range: X=[%d %d], Y=[%d %d], Z=[%d %d]\n', ...
    rangeX(1), rangeX(2), rangeY(1), rangeY(2), rangeZ(1), rangeZ(2));

% -------------------- Parallel setup (best-effort) -----------------------
if useParallel
    ncores = feature('numcores');
    eff_ncores = floor(ncores*0.9); %avoid overhead and RAM excess
    try
        hasParallel = license('test', 'Distrib_Computing_Toolbox');
        if hasParallel
            % Determine if the file extension is safe for thread-based reading.
            % TIFF and GIF files use underlying C/C++ libraries that are not thread-safe.
            isThreadSafeIO = ~contains(extType, {'tif', 'tiff', 'gif'}, 'IgnoreCase', true);
            
            pool = gcp('nocreate');
            if isempty(pool)
                if isThreadSafeIO
                    % For safe formats (BMP, PNG, JPG), try threads first to save RAM
                    try
                        parpool('threads',eff_ncores);
                    catch
                        try
                            parpool('local',eff_ncores);
                        catch
                        end
                    end
                else
                    % For unsafe formats (TIFF), skip threads and go straight to process-based pool
                    try
                        parpool('local',eff_ncores);
                    catch
                    end
                end
                pool = gcp('nocreate'); % Update pool reference
            end
            
            % Safety handler: what if the user already had a thread pool open 
            % in the session before calling this function?
            if ~isempty(pool) && isa(pool, 'parallel.ThreadPool') && ~isThreadSafeIO
                warning(['A thread-based pool is already active, but the extension ' ...
                         '"%s" is not thread-safe. Parallelism will be disabled for this ' ...
                         'execution to prevent imread failures.'], extType);
                useParallel = false;
            else
                useParallel = ~isempty(pool);
            end
        else
            useParallel = false;
        end
    catch
        useParallel = false;
    end
end

% -------------------- Compute global Otsu threshold (Optimized via Sampling) --------
% Check if the first image is already natively logical (e.g., pure 1-bit BMP)
firstImageTest = imread([imgDir files(imgsRange(1)).name]);
isAlreadyLogical = islogical(firstImageTest);

if isAlreadyLogical
    % If it is already binary, we do not need Otsu thresholding
    globalLevel = 0.5; % Safety dummy value, will not be processed in imbinarize
elseif binThreshold == -1
    % Read samples distributed across the volume instead of all slices to save I/O and RAM
    numSamples = min(numImgs, 25); 
    sampleIndices = round(linspace(1, numImgs, numSamples));
    samplePixels = cell(numSamples, 1);
    
    % Define a reference data class based on the first sample to prevent cell2mat errors
    I_ref = imread([imgDir files(imgsRange(sampleIndices(1))).name]);
    if ndims(I_ref) == 3
        I_ref = rgb2gray(I_ref);
    end
    refClass = class(I_ref);
    
    for s = 1:numSamples
        idxFile = imgsRange(sampleIndices(s));
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3
            I = rgb2gray(I);
        end
        I = I(rowsRange, colsRange);
        
        % Force data type consistency before appending to cell array
        if ~isa(I, refClass)
            I = cast(I, refClass);
        end
        
        samplePixels{s} = I(:); 
    end
    
    sampledVol = cell2mat(samplePixels);
    globalLevel = graythresh(sampledVol);
else
    globalLevel = binThreshold;
end

% -------------------- Single-Pass Read and Binarize ----------------------
% Pre-allocate purely logical matrix.
C = false(numRows, numCols, numImgs);

if useParallel
    parfor ii = 1:numImgs
        idxFile = imgsRange(ii);
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3
            I = rgb2gray(I);
        end
        I = I(rowsRange, colsRange);
        
        % Direct assignment if already binary; otherwise, binarize
        if islogical(I)
            C(:, :, ii) = I;
        else
            C(:, :, ii) = imbinarize(I, globalLevel);
        end
    end
else
    for ii = 1:numImgs
        idxFile = imgsRange(ii);
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3
            I = rgb2gray(I);
        end
        I = I(rowsRange, colsRange);
        
        if islogical(I)
            C(:, :, ii) = I;
        else
            C(:, :, ii) = imbinarize(I, globalLevel);
        end
    end
end

end