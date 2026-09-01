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

% Check compatibility for legacy binarization
isLegacy = verLessThan('matlab', '9.0'); 

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

% Numeric-aware sorting
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

% Clamp user ranges
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
% (Mantido igual ao original, já está bem otimizado)
if useParallel
    ncores = feature('numcores');
    eff_ncores = floor(ncores*0.9);
    try
        hasParallel = license('test', 'Distrib_Computing_Toolbox');
        if hasParallel
            isThreadSafeIO = ~contains(extType, {'tif', 'tiff', 'gif'}, 'IgnoreCase', true);
            pool = gcp('nocreate');
            if isempty(pool)
                if isThreadSafeIO
                    try parpool('threads',eff_ncores); catch, try parpool('local',eff_ncores); catch; end; end
                else
                    try parpool('local',eff_ncores); catch; end
                end
                pool = gcp('nocreate');
            end
            if ~isempty(pool) && isa(pool, 'parallel.ThreadPool') && ~isThreadSafeIO
                warning(['A thread-based pool is already active, but the extension ' ...
                         '"%s" is not thread-safe. Parallelism will be disabled.'], extType);
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

% -------------------- Determine Binarization Strategy --------------------
firstImageTest = imread([imgDir files(imgsRange(1)).name]);
if ndims(firstImageTest) == 3
    firstImageTest = rgb2gray(firstImageTest);
end

isNativelyLogical = islogical(firstImageTest);
isPseudoBinary = false;
minVal = 0;

% Checagem de Imagens Pseudo-Binárias (ex: uint8, valores [0, 255])
if ~isNativelyLogical
    valores_unicos = unique(firstImageTest(:));
    if numel(valores_unicos) <= 2
        isPseudoBinary = true;
        minVal = min(valores_unicos);
        fprintf('Pseudo-binary images detected. Bypassing Otsu calculation.\n');
    end
end

% Cálculo Global Otsu (se for realmente escala de cinza)
if isNativelyLogical || isPseudoBinary
    globalLevel = 0.5; % Dummy value, não será usado.
elseif binThreshold == -1
    fprintf('Grayscale images detected. Computing global Otsu threshold via sampling...\n');
    numSamples = min(numImgs, 25); 
    sampleIndices = round(linspace(1, numImgs, numSamples));
    samplePixels = cell(numSamples, 1);
    
    refClass = class(firstImageTest);
    
    for s = 1:numSamples
        idxFile = imgsRange(sampleIndices(s));
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3, I = rgb2gray(I); end
        I = I(rowsRange, colsRange);
        
        if ~isa(I, refClass)
            I = cast(I, refClass);
        end
        samplePixels{s} = I(:); 
    end
    
    sampledVol = cell2mat(samplePixels);
    globalLevel = graythresh(sampledVol);
    fprintf('Otsu global threshold set to: %.4f\n', globalLevel);
else
    globalLevel = binThreshold;
end

% -------------------- Single-Pass Read and Binarize ----------------------
C = false(numRows, numCols, numImgs);

if useParallel
    parfor ii = 1:numImgs
        idxFile = imgsRange(ii);
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3, I = rgb2gray(I); end
        I = I(rowsRange, colsRange);
        
        % Nova Lógica de Atribuição Inteligente e Legada
        if isNativelyLogical
            C(:, :, ii) = I;
        elseif isPseudoBinary
            C(:, :, ii) = (I ~= minVal);
        else
            if isLegacy
                C(:, :, ii) = im2bw(I, globalLevel);
            else
                C(:, :, ii) = imbinarize(I, globalLevel);
            end
        end
    end
else
    for ii = 1:numImgs
        idxFile = imgsRange(ii);
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3, I = rgb2gray(I); end
        I = I(rowsRange, colsRange);
        
        % Nova Lógica de Atribuição Inteligente e Legada
        if isNativelyLogical
            C(:, :, ii) = I;
        elseif isPseudoBinary
            C(:, :, ii) = (I ~= minVal);
        else
            if isLegacy
                C(:, :, ii) = im2bw(I, globalLevel);
            else
                C(:, :, ii) = imbinarize(I, globalLevel);
            end
        end
    end
end

end