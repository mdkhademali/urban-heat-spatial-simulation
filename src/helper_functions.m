function h = helper_functions()
%HELPER_FUNCTIONS  Collection of small reusable utilities (returned as handles).
%
%   h = helper_functions();
%   z = h.fractal_noise(500, 2.5, 7);
%
%   MATLAB only exposes the first function in a file to other files, so this
%   file returns a struct of function handles to its local functions. All
%   functions use base MATLAB only (no toolboxes).
%
%   Available handles
%     white_randn        deterministic, platform-independent N(0,1) noise
%     white_rand         deterministic, platform-independent U(0,1) noise
%     fractal_noise      spatially correlated (1/f^beta) noise field, z-scored
%     gaussian_smooth    separable Gaussian low-pass filter (replicate padding)
%     euclidean_distance exact Euclidean distance transform (metres)
%     percentile_value   percentile without the Statistics Toolbox
%     normalize01        rescale an array to [0, 1]
%     ensure_folder      create a folder if it does not exist
%     t_pvalue           two-sided p-value of a t statistic (base MATLAB)
%     make_colormap      interpolate a colormap from control colours

    h.white_randn        = @white_randn;
    h.white_rand         = @white_rand;
    h.fractal_noise      = @fractal_noise;
    h.gaussian_smooth    = @gaussian_smooth;
    h.euclidean_distance = @euclidean_distance;
    h.percentile_value   = @percentile_value;
    h.normalize01        = @normalize01;
    h.ensure_folder      = @ensure_folder;
    h.t_pvalue           = @t_pvalue;
    h.make_colormap      = @make_colormap;
end

% -------------------------------------------------------------------------
function u = white_rand(count, seed)
%WHITE_RAND  Deterministic uniform (0,1) numbers, column vector.
%   Park-Miller "minimal standard" generator, x <- 16807*x mod (2^31-1).
%   Implemented with exact double arithmetic, so the output is identical on
%   every MATLAB version, operating system and CPU. This is what makes the
%   whole project bit-for-bit reproducible.
    modulus = 2147483647;
    x = mod(1103515245 * (seed + 977) + 12345, modulus - 1) + 1;   % scramble seed
    for k = 1:20                       % burn-in
        x = mod(16807 * x, modulus);
    end
    u = zeros(count, 1);
    for k = 1:count
        x = mod(16807 * x, modulus);
        u(k) = x / modulus;
    end
end

% -------------------------------------------------------------------------
function z = white_randn(rows, cols, seed)
%WHITE_RANDN  Deterministic standard-normal matrix (Box-Muller transform).
    count = rows * cols;
    m = ceil(count / 2);
    u = white_rand(2 * m, seed);
    u1 = u(1:m);
    u2 = u(m+1:2*m);
    r  = sqrt(-2 * log(u1));
    zz = [r .* cos(2*pi*u2); r .* sin(2*pi*u2)];
    z  = reshape(zz(1:count), rows, cols);
end

% -------------------------------------------------------------------------
function field = fractal_noise(n, beta, seed)
%FRACTAL_NOISE  Spatially correlated noise with power spectrum ~ 1/f^beta.
%   beta ~ 1 : rough / small-scale texture
%   beta ~ 2 : natural-looking terrain-like fields
%   beta ~ 3 : very smooth, large blobs
%   Output is standardised to zero mean and unit standard deviation.
    k = ifftshift(-floor(n/2):ceil(n/2)-1);
    [kx, ky] = meshgrid(k, k);
    f = sqrt(kx.^2 + ky.^2);
    f(1,1) = 1;                         % avoid division by zero at DC
    amplitude = f .^ (-beta / 2);
    amplitude(1,1) = 0;                 % remove the mean
    spectrum = fft2(white_randn(n, n, seed)) .* amplitude;
    field = real(ifft2(spectrum));
    field = (field - mean(field(:))) / std(field(:));
end

% -------------------------------------------------------------------------
function B = gaussian_smooth(A, sigma)
%GAUSSIAN_SMOOTH  Separable Gaussian filter with edge replication.
    if sigma <= 0
        B = A;
        return;
    end
    r = ceil(3 * sigma);
    x = -r:r;
    kernel = exp(-x.^2 / (2 * sigma^2));
    kernel = kernel / sum(kernel);
    [nr, nc] = size(A);
    ri = [ones(1, r), 1:nr, nr * ones(1, r)];
    ci = [ones(1, r), 1:nc, nc * ones(1, r)];
    P = A(ri, ci);
    P = conv2(P, kernel(:), 'valid');   % filter down the columns
    B = conv2(P, kernel(:).', 'valid'); % filter along the rows
end

% -------------------------------------------------------------------------
function d = euclidean_distance(mask, cell_size)
%EUCLIDEAN_DISTANCE  Distance (metres) from every cell to the nearest true cell.
%   Exact Euclidean distance transform (Felzenszwalb & Huttenlocher, 2012),
%   implemented in plain MATLAB so that bwdist (Image Processing Toolbox)
%   is not required.
    if ~any(mask(:))
        d = inf(size(mask));
        return;
    end
    [nr, nc] = size(mask);
    big = 1e10;
    f = big * ones(nr, nc);
    f(mask) = 0;
    g = zeros(nr, nc);
    for j = 1:nc
        g(:, j) = squared_dt_1d(f(:, j));
    end
    q = zeros(nr, nc);
    for i = 1:nr
        q(i, :) = squared_dt_1d(g(i, :).').';
    end
    d = sqrt(q) * cell_size;
end

function d = squared_dt_1d(f)
    n = numel(f);
    d = zeros(n, 1);
    v = zeros(n, 1);          % locations of parabolas in the lower envelope
    z = zeros(n + 1, 1);      % boundaries between parabolas
    k = 1;
    v(1) = 1;
    z(1) = -Inf;
    z(2) = Inf;
    for q = 2:n
        s = ((f(q) + q^2) - (f(v(k)) + v(k)^2)) / (2*q - 2*v(k));
        while s <= z(k)
            k = k - 1;
            s = ((f(q) + q^2) - (f(v(k)) + v(k)^2)) / (2*q - 2*v(k));
        end
        k = k + 1;
        v(k) = q;
        z(k) = s;
        z(k + 1) = Inf;
    end
    k = 1;
    for q = 1:n
        while z(k + 1) < q
            k = k + 1;
        end
        d(q) = (q - v(k))^2 + f(v(k));
    end
end

% -------------------------------------------------------------------------
function v = percentile_value(x, p)
%PERCENTILE_VALUE  p-th percentile (0-100) of the values in x.
    xs = sort(x(:));
    idx = min(max(ceil(p / 100 * numel(xs)), 1), numel(xs));
    v = xs(idx);
end

% -------------------------------------------------------------------------
function B = normalize01(A)
%NORMALIZE01  Linearly rescale to the range [0, 1].
    lo = min(A(:));
    hi = max(A(:));
    B = (A - lo) / (hi - lo);
end

% -------------------------------------------------------------------------
function ensure_folder(folder)
%ENSURE_FOLDER  Create a folder if it does not already exist.
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end

% -------------------------------------------------------------------------
function p = t_pvalue(t, dof)
%T_PVALUE  Two-sided p-value for a t statistic using the incomplete beta
%   function (base MATLAB; avoids tcdf from the Statistics Toolbox).
    p = betainc(dof ./ (dof + t.^2), dof / 2, 0.5);
end

% -------------------------------------------------------------------------
function cmap = make_colormap(colors, n)
%MAKE_COLORMAP  Build an n-row colormap by linear interpolation of colors.
    m = size(colors, 1);
    xi = linspace(1, m, n);
    cmap = zeros(n, 3);
    for c = 1:3
        cmap(:, c) = interp1(1:m, colors(:, c), xi, 'linear');
    end
    cmap = min(max(cmap, 0), 1);
end
