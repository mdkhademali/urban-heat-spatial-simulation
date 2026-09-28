function S = simulate_lst(cfg, L, h)
%SIMULATE_LST  Simulate spectral indices (NDVI, NDBI, NDWI) and LST.
%
%   S = simulate_lst(cfg, L, h)
%
%   Design principle
%   ----------------
%   LST is NOT assigned from the land-cover class. It is built from physical
%   drivers (sub-pixel fractions of each cover type, distance to water and
%   vegetation, an urban heat-island term, a regional gradient and spatially
%   correlated noise). The spectral indices are generated from the same
%   landscape but with their OWN independent noise, so NDVI/NDBI/NDWI are
%   correlated with LST only indirectly, through the shared landscape --
%   just as in real remote-sensing data.
%
%   LST equation (degrees Celsius)
%     LST = T0 + a_b*f_built + a_s*f_soil - a_dv*f_dense - a_sv*f_sparse
%              - a_w*f_water - c_w*exp(-d_w/Lw)*(1-f_water)
%              - c_v*exp(-d_v/Lv)*(1-f_dense) + UHI + gradient + noise

    lc = L.landcover;
    n  = size(lc, 1);
    cs = cfg.cell_size_m;

    % ---------------------------------------------------------------
    % 1) Sub-pixel cover fractions (Gaussian-blurred class indicators)
    % ---------------------------------------------------------------
    %   'optical' blur  -> used for NDVI / NDBI / NDWI (sensor mixing)
    %   'thermal' blur  -> used for LST (heat diffuses over a wider area)
    frac_opt = zeros(n, n, 5);
    frac_th  = zeros(n, n, 5);
    for k = 1:5
        indicator = double(lc == k);
        frac_opt(:, :, k) = h.gaussian_smooth(indicator, cfg.optical_blur_cells);
        frac_th(:, :, k)  = h.gaussian_smooth(indicator, cfg.thermal_blur_cells);
    end

    % ---------------------------------------------------------------
    % 2) Spectral indices: mixing of class-typical values + noise
    % ---------------------------------------------------------------
    ndvi = zeros(n);  ndbi = zeros(n);  ndwi = zeros(n);
    for k = 1:5
        ndvi = ndvi + frac_opt(:, :, k) * cfg.class_ndvi(k);
        ndbi = ndbi + frac_opt(:, :, k) * cfg.class_ndbi(k);
        ndwi = ndwi + frac_opt(:, :, k) * cfg.class_ndwi(k);
    end
    ndvi = ndvi + cfg.ndvi_noise_sd * h.fractal_noise(n, 1.4, cfg.seed + 101);
    ndbi = ndbi + cfg.ndbi_noise_sd * h.fractal_noise(n, 1.4, cfg.seed + 102);
    ndwi = ndwi + cfg.ndwi_noise_sd * h.fractal_noise(n, 1.4, cfg.seed + 103);
    ndvi = min(max(ndvi, -1), 1);
    ndbi = min(max(ndbi, -1), 1);
    ndwi = min(max(ndwi, -1), 1);

    % ---------------------------------------------------------------
    % 3) Distances to cooling features (metres)
    % ---------------------------------------------------------------
    dist_water = h.euclidean_distance(lc == 1, cs);
    dist_dense = h.euclidean_distance(lc == 2, cs);

    % ---------------------------------------------------------------
    % 4) LST components (each one is stored so the model is transparent)
    % ---------------------------------------------------------------
    f_water  = frac_th(:, :, 1);
    f_dense  = frac_th(:, :, 2);
    f_sparse = frac_th(:, :, 3);
    f_built  = frac_th(:, :, 4);
    f_soil   = frac_th(:, :, 5);

    comp.base       = cfg.base_temp_C * ones(n);
    comp.built_up   =  cfg.a_built  * f_built;
    comp.soil       =  cfg.a_soil   * f_soil;
    comp.dense_veg  = -cfg.a_dense  * f_dense;
    comp.sparse_veg = -cfg.a_sparse * f_sparse;
    comp.water_local = -cfg.a_water * f_water;
    % Cooling that spreads from water into neighbouring cells (exponential decay)
    comp.water_advect = -cfg.c_water * exp(-dist_water / cfg.decay_water_m) .* (1 - f_water);
    % Shading / evapotranspiration halo around dense vegetation
    comp.veg_halo = -cfg.c_veg * exp(-dist_dense / cfg.decay_veg_m) .* (1 - f_dense);
    % Urban heat island: regional warming that follows the smoothed urban intensity
    comp.uhi = cfg.uhi_amplitude_C * h.gaussian_smooth(L.urban_norm, cfg.uhi_blur_cells);
    % Large-scale gradient (e.g. elevation / sea-breeze / regional advection)
    comp.gradient = cfg.gradient_x_C * (L.Xn - 0.5) + cfg.gradient_y_C * (L.Yn - 0.5);
    % Spatially correlated random variability (unresolved processes, sensor noise)
    comp.noise = cfg.lst_noise_sd * h.fractal_noise(n, 1.6, cfg.seed + 201);

    lst = comp.base + comp.built_up + comp.soil + comp.dense_veg + comp.sparse_veg ...
        + comp.water_local + comp.water_advect + comp.veg_halo ...
        + comp.uhi + comp.gradient + comp.noise;

    % ---------------------------------------------------------------
    % 5) Output
    % ---------------------------------------------------------------
    S.ndvi = ndvi;
    S.ndbi = ndbi;
    S.ndwi = ndwi;
    S.lst  = lst;
    S.lst_anomaly = lst - mean(lst(:));
    S.dist_water_m = dist_water;
    S.dist_dense_m = dist_dense;
    S.fractions_thermal = frac_th;
    S.components = comp;
end
