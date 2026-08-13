function C = load_volume2(imgDir, extType, volumetricSize, binThreshold, useParallel)
% LOAD_VOLUME  Read an image stack and produce a binary 3D volume.
%
%   C = LOAD_VOLUME(imgDir)
%   C = LOAD_VOLUME(imgDir, extType, volumetricSize, binThreshold, useParallel)
%
% INPUTS:
%   imgDir         - (char) path to the folder containing the image stack.
%   extType        - (char, optional) image extension to search for,
%                     e.g. '.bmp' (default), '.tif', '.png'.
%   volumetricSize - (optional) defines the sub-volume to read. Accepts:
%                     [] (default)         -> full image extents and all slices
%                     string or char 'fit' -> automatically fits a centered cubic volume 
%                     scalar N             -> NxNxN volume starting at (1,1,1)
%                     3x2 matrix           -> explicit ranges
%   binThreshold   - (numeric, optional) binarization threshold.
%                     -1 (default) => compute global Otsu threshold via sampling.
%                     Otherwise use the provided numeric threshold (in [0,1]).
%   useParallel    - (logical, optional) request parallel processing.
%
% OUTPUT:
%   C - logical 3D array (rows x cols x slices) containing the binary pore mask.
%
% OTIMIZAÇÕES APLICADAS:
%   - Remoção de alocação matriz em ponto flutuante (single), economizando memória.
%   - Amostragem estratificada para cálculo do limiar de Otsu global.
%   - Leitura e binarização fundidas em um único laço de I/O.
%   - Detecção inteligente (Bypass) para imagens que já são nativamente lógicas (binárias).

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
files = files(~[files.isdir]);  % remove directories

if isempty(files)
    error('No image files found in the specified directory: %s', imgDir);
end

% Numeric-aware sorting: extract trailing number groups and sort by last group
file_numbers = zeros(1, numel(files));
for ff = 1:numel(files)
    nums = regexp(files(ff).name, '\d+', 'match');
    if isempty(nums)
        file_numbers(ff) = 0;
    else
        file_numbers(ff) = str2double(nums{end});
    end
end
[~, sort_idx] = sort(file_numbers);
files = files(sort_idx);
num_images = numel(files);

% -------------------- Validate consistent image geometry -----------------
first_info = imfinfo([imgDir files(1).name]);
orig_rows = first_info.Height;
orig_cols = first_info.Width;

for ff = 2:numel(files)
    info_ff = imfinfo([imgDir files(ff).name]);
    if info_ff.Height ~= orig_rows || info_ff.Width ~= orig_cols
        error('All images must have identical dimensions. File "%s" differs.', files(ff).name);
    end
end

% -------------------- Determine volumetric ranges ------------------------
if isempty(volumetricSize)
    rangeX = [1 orig_cols];
    rangeY = [1 orig_rows];
    rangeZ = [1 num_images];
elseif isstring(volumetricSize) || ischar(volumetricSize)
    if  string(volumetricSize) == "fit"
        N = min([orig_rows,orig_cols,num_images]);
        diff_x = round((orig_cols-N)/2);
        diff_y = round((orig_rows-N)/2);
        diff_z = round((num_images-N)/2);
        rangeX = [diff_x+1 diff_x+N];
        rangeY = [diff_y+1 diff_y+N];
        rangeZ = [diff_z+1 diff_z+N];
    else
        error("Invalid volumetricSize. Provide a invalid name. Must be 'fit'");
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
rangeX(2) = min(rangeX(2), orig_cols);
rangeY(2) = min(rangeY(2), orig_rows);
rangeZ(2) = min(rangeZ(2), num_images);

rows_range = rangeY(1):rangeY(2);
cols_range = rangeX(1):rangeX(2);
imgs_range = rangeZ(1):rangeZ(2);

rows = numel(rows_range);
cols = numel(cols_range);
num_imgs = numel(imgs_range);

fprintf('Using volume range: X=[%d %d], Y=[%d %d], Z=[%d %d]\n', ...
    rangeX(1), rangeX(2), rangeY(1), rangeY(2), rangeZ(1), rangeZ(2));

% -------------------- Parallel setup (best-effort) -----------------------
if useParallel
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
                        parpool('threads');
                    catch
                        try
                            parpool('local');
                        catch
                        end
                    end
                else
                    % For unsafe formats (TIFF), skip threads and go straight to process-based pool
                    try
                        parpool('local');
                    catch
                    end
                end
                pool = gcp('nocreate'); % Update pool reference
            end
            
            % Safety handler: what if the user already had a thread pool open 
            % in the session before calling this function?
            if ~isempty(pool) && isa(pool, 'parallel.ThreadPool') && ~isThreadSafeIO
                warning(['Uma pool baseada em threads ja esta ativa, mas a extensao ' ...
                         '"%s" nao e thread-safe. O paralelismo sera desativado nesta ' ...
                         'execucao para evitar falhas no imread.'], extType);
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
I_first_test = imread([imgDir files(imgs_range(1)).name]);
isAlreadyLogical = islogical(I_first_test);

if isAlreadyLogical
    % If it is already binary, we do not need Otsu thresholding
    globalLevel = 0.5; % Valor de segurança, não será processado no imbinarize
elseif binThreshold == -1
    % Read samples distributed across the volume instead of all slices to save I/O and RAM
    num_samples = min(num_imgs, 25); 
    sample_indices = round(linspace(1, num_imgs, num_samples));
    sample_pixels = cell(num_samples, 1);
    
    for s = 1:num_samples
        idxFile = imgs_range(sample_indices(s));
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3
            I = rgb2gray(I);
        end
        I = I(rows_range, cols_range);
        sample_pixels{s} = I(:); % Keeps original type (uint8 or uint16)
    end
    
    sampled_vol = cell2mat(sample_pixels);
    globalLevel = graythresh(sampled_vol);
else
    globalLevel = binThreshold;
end

% -------------------- Single-Pass Read and Binarize ----------------------
% Pre-allocate purely logical matrix.
C = false(rows, cols, num_imgs);

if useParallel
    parfor ii = 1:num_imgs
        idxFile = imgs_range(ii);
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3
            I = rgb2gray(I);
        end
        I = I(rows_range, cols_range);
        
        % Direct assignment if already binary; otherwise, binarize
        if islogical(I)
            C(:, :, ii) = I;
        else
            C(:, :, ii) = imbinarize(I, globalLevel);
        end
    end
else
    for ii = 1:num_imgs
        idxFile = imgs_range(ii);
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3
            I = rgb2gray(I);
        end
        I = I(rows_range, cols_range);
        
        if islogical(I)
            C(:, :, ii) = I;
        else
            C(:, :, ii) = imbinarize(I, globalLevel);
        end
    end
end

end