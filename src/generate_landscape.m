function L = generate_landscape(cfg, h)
%GENERATE_LANDSCAPE  Build a synthetic mixed urban landscape (land-cover map).
%
%   L = generate_landscape(cfg, h)
%
%   Land-cover classes
%       1 = Water              (river + lake)
%       2 = Dense vegetation   (parks, forest patches)
%       3 = Sparse vegetation  (grass, scattered trees)
%       4 = Built-up           (urban cores + transport corridors)
%       5 = Bare soil          (fallow land, construction, urban fringe)
%
%   Strategy: every class is derived from smooth, spatially coherent
%   "driver" fields (Gaussian urban cores, line corridors, a sinuous river,
%   1/f^beta noise), so the map contains connected, irregular shapes and
%   not salt-and-pepper pixels.

    n  = cfg.grid_size;
    cs = cfg.cell_size_m;

    % Normalised coordinates in [0,1]. Row 1 is the southern edge.
    [col, row] = meshgrid(1:n, 1:n);
    Xn = (col - 0.5) / n;
    Yn = (row - 0.5) / n;

    % ---------------------------------------------------------------
    % 1) Urban intensity field: cores + corridors + texture
    % ---------------------------------------------------------------
    urban = zeros(n);
    for k = 1:size(cfg.urban_cores, 1)              % [cx cy sigma weight]
        c = cfg.urban_cores(k, :);
        d2 = (Xn - c(1)).^2 + (Yn - c(2)).^2;
        urban = urban + c(4) * exp(-d2 / (2 * c(3)^2));
    end
    for k = 1:size(cfg.corridors, 1)                % [x1 y1 x2 y2 halfwidth weight]
        c = cfg.corridors(k, :);
        d = point_segment_distance(Xn, Yn, c(1), c(2), c(3), c(4));
        urban = urban + c(6) * exp(-d.^2 / (2 * c(5)^2));
    end
    urban = urban + cfg.urban_noise_weight * h.fractal_noise(n, 2.6, cfg.seed + 11);
    urban = h.gaussian_smooth(urban, 1.5);
    urban_norm = h.normalize01(urban);

    % ---------------------------------------------------------------
    % 2) Water: a meandering river plus an irregular lake
    % ---------------------------------------------------------------
    centre_y = cfg.river_y0 ...
             + 0.10 * sin(2*pi*1.3*Xn + 1.0) ...
             + 0.04 * sin(2*pi*3.1*Xn + 2.0);
    half_width = cfg.river_halfwidth * (1 + 0.5 * sin(2*pi*2.2*Xn + 0.5));
    river = abs(Yn - centre_y) < half_width;

    lake_d2 = (Xn - cfg.lake_center(1)).^2 + (Yn - cfg.lake_center(2)).^2;
    lake_field = exp(-lake_d2 / (2 * cfg.lake_sigma^2)) ...
               + cfg.lake_roughness * h.fractal_noise(n, 3.0, cfg.seed + 21);
    lake = lake_field > cfg.lake_threshold;

    water = river | lake;

    % ---------------------------------------------------------------
    % 3) Built-up: highest urban intensity among non-water cells
    % ---------------------------------------------------------------
    thr_urban = h.percentile_value(urban_norm(~water), 100 * (1 - cfg.builtup_fraction));
    builtup = (urban_norm > thr_urban) & ~water;

    % ---------------------------------------------------------------
    % 4) Vegetation and bare soil in the remaining land
    % ---------------------------------------------------------------
    remaining = ~water & ~builtup;

    % Vegetation prefers areas away from the urban core.
    veg_field = 0.75 * h.fractal_noise(n, 3.0, cfg.seed + 31) ...
              - 0.90 * urban_norm;
    % Bare soil: patchy, more common on the urban fringe.
    bare_field = h.fractal_noise(n, 2.6, cfg.seed + 41) + 0.8 * urban_norm;

    thr_bare = h.percentile_value(bare_field(remaining), 100 * (1 - cfg.bare_fraction_of_land));
    bare = remaining & (bare_field > thr_bare);

    vegetated = remaining & ~bare;
    thr_dense = h.percentile_value(veg_field(vegetated), 100 * (1 - cfg.dense_share_of_vegetation));
    dense  = vegetated & (veg_field > thr_dense);
    sparse = vegetated & ~dense;

    % ---------------------------------------------------------------
    % 5) Assemble the classified map
    % ---------------------------------------------------------------
    landcover = zeros(n, 'uint8');
    landcover(water)   = 1;
    landcover(dense)   = 2;
    landcover(sparse)  = 3;
    landcover(builtup) = 4;
    landcover(bare)    = 5;

    L.landcover     = landcover;
    L.urban_norm    = urban_norm;     % continuous urban intensity (0-1)
    L.veg_field     = veg_field;
    L.x_km          = ((1:n) - 0.5) * cs / 1000;
    L.y_km          = ((1:n) - 0.5) * cs / 1000;
    L.Xn            = Xn;
    L.Yn            = Yn;
    L.class_names   = {'Water', 'Dense vegetation', 'Sparse vegetation', 'Built-up', 'Bare soil'};
end

% -------------------------------------------------------------------------
function d = point_segment_distance(X, Y, x1, y1, x2, y2)
%POINT_SEGMENT_DISTANCE  Distance from each grid point to a line segment.
    dx = x2 - x1;
    dy = y2 - y1;
    t = ((X - x1) * dx + (Y - y1) * dy) / (dx^2 + dy^2);
    t = min(max(t, 0), 1);
    d = sqrt((X - (x1 + t * dx)).^2 + (Y - (y1 + t * dy)).^2);
end
