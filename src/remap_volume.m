function C_sub = remap_volume(C_full, subvolSpec, alignment)
% ============================================================
% REMAP_VOLUME
% ============================================================
% Extracts a subvolume from a fully-loaded 3D binary or grayscale
% volume without reloading data from disk. This allows flexible
% slicing based on proportion, absolute size, or explicit ranges,
% with control over spatial alignment in X, Y, and Z directions.
%
% ============================================================
% INPUTS:
%   C_full       : 3D matrix (logical or numeric), the original volume
%   subvolSpec   : Specifies the size of the subvolume. Accepted forms:
%                  - Scalar <1: proportional cube (same ratio X,Y,Z)
%                  - 3-element vector <1: proportional per axis [X Y Z]
%                  - Scalar >=1: absolute cube size (same for all axes)
%                  - 3-element vector >=1: absolute size per axis [X Y Z]
%                  - 3x2 matrix: explicit coordinate ranges
%                    [x_min x_max; y_min y_max; z_min z_max]
%   alignment    : String specifying subvolume alignment along each axis:
%                  - X-axis: "left", "right", "center" (default)
%                  - Y-axis: "top", "bottom", "center" (default)
%                  - Z-axis: "front", "back", "center" (default)
%                  - Combinations allowed: e.g., "trf" (top-right-front)
%
% ============================================================
% OUTPUT:
%   C_sub        : Extracted subvolume, same class as C_full
%
% ============================================================
% OTIMIZAÇÕES APLICADAS:
%   - Remoção completa do overhead de I/O da função 'verLessThan'.
%   - Uso universal do 'strfind' nativo em C para compatibilidade imediata.
%   - Tratamento seguro de matrizes com a 3ª dimensão colapsada.
% ============================================================

%% -------------------- Input validation --------------------
% Default alignment is "center" if not specified
if nargin < 3 || isempty(alignment)
    alignment = 'center';
end

% Safe dimension extraction (avoids errors if Nz == 1)
Nx_full = size(C_full, 2);  % X-axis (columns)
Ny_full = size(C_full, 1);  % Y-axis (rows)
Nz_full = size(C_full, 3);  % Z-axis (slices)

if Nx_full == 0 || Ny_full == 0
    error('C_full cannot be empty.');
end

%% -------------------- Determine extraction mode --------------------
% Handle different forms of subvolSpec
if isscalar(subvolSpec)
    if subvolSpec > 0 && subvolSpec < 1
        % Proportional cube (same proportion in all axes)
        Nx = round(Nx_full * subvolSpec);
        Ny = round(Ny_full * subvolSpec);
        Nz = round(Nz_full * subvolSpec);
    else
        % Absolute cube
        Nx = round(subvolSpec);
        Ny = Nx;
        Nz = Nx;
    end

elseif isvector(subvolSpec) && numel(subvolSpec)==3
    if all(subvolSpec>0 & subvolSpec<1)
        % Proportional per axis
        Nx = round(Nx_full * subvolSpec(1));
        Ny = round(Ny_full * subvolSpec(2));
        Nz = round(Nz_full * subvolSpec(3));
    else
        % Absolute size per axis
        Nx = round(subvolSpec(1));
        Ny = round(subvolSpec(2));
        Nz = round(subvolSpec(3));
    end

elseif ismatrix(subvolSpec) && all(size(subvolSpec)==[3 2])
    % Explicit coordinate ranges
    rangeX = round(subvolSpec(1,:));
    rangeY = round(subvolSpec(2,:));
    rangeZ = round(subvolSpec(3,:));
    % Direct extraction and return
    C_sub = C_full(rangeY(1):rangeY(2), rangeX(1):rangeX(2), rangeZ(1):rangeZ(2));
    return;

else
    error('Invalid subvolSpec. Must be scalar, 3-element vector, or 3x2 matrix.');
end

%% -------------------- Clamp subvolume sizes --------------------
% Ensure subvolume does not exceed original volume
Nx = min(Nx, Nx_full);
Ny = min(Ny, Ny_full);
Nz = min(Nz, Nz_full);

%% -------------------- Parse alignment string --------------------
% Lowercase for consistency
alignment = lower(alignment);

% Unified fast string search.
% strfind is implemented in C at the core MATLAB level and is fully
% backward compatible, eliminating the need for slow verLessThan checks.
hasRight  = ~isempty(strfind(alignment,'right'));
hasLeft   = ~isempty(strfind(alignment,'left'));
hasTop    = ~isempty(strfind(alignment,'top'));
hasBottom = ~isempty(strfind(alignment,'bottom'));
hasCenter = ~isempty(strfind(alignment,'center')) || ~isempty(strfind(alignment,'c'));
hasFront  = ~isempty(strfind(alignment,'front'));
hasBack   = ~isempty(strfind(alignment,'back'));

%% -------------------- Compute start indices per axis --------------------
% X-axis (columns)
if hasRight
    x_start = Nx_full - Nx + 1; % Right-aligned
elseif hasCenter
    x_start = floor((Nx_full - Nx)/2) + 1; % Centered
else
    x_start = 1; % Left-aligned
end

% Y-axis (rows)
if hasBottom
    y_start = Ny_full - Ny + 1; % Bottom-aligned
elseif hasCenter
    y_start = floor((Ny_full - Ny)/2) + 1; % Centered
else
    y_start = 1; % Top-aligned
end

% Z-axis (slices)
if hasBack
    z_start = Nz_full - Nz + 1; % Back-aligned
elseif hasCenter
    z_start = floor((Nz_full - Nz)/2) + 1; % Centered
elseif hasFront
    z_start = 1; % Front-aligned
else
    z_start = 1; % Default to front if unspecified
end

%% -------------------- Clamp start indices --------------------
% Ensure indices are within valid bounds
x_start = max(1,min(x_start,Nx_full-Nx+1));
y_start = max(1,min(y_start,Ny_full-Ny+1));
z_start = max(1,min(z_start,Nz_full-Nz+1));

% Compute end indices
x_end = x_start + Nx - 1;
y_end = y_start + Ny - 1;
z_end = z_start + Nz - 1;

%% -------------------- Extract subvolume --------------------
C_sub = C_full(y_start:y_end, x_start:x_end, z_start:z_end);

%% -------------------- Debug / Info --------------------
fprintf('Remapped subvolume -> X=[%d %d], Y=[%d %d], Z=[%d %d], alignment="%s"\n', ...
    x_start, x_end, y_start, y_end, z_start, z_end, alignment);

end