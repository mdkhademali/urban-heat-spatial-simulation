function A = analyze_spatial_relationships(cfg, L, S, h)
%ANALYZE_SPATIAL_RELATIONSHIPS  Statistics, regression and cooling profiles.
%
%   A = analyze_spatial_relationships(cfg, L, S, h)
%
%   Outputs (fields of A)
%     summary      global mean/SD/min/max of LST and mean NDVI/NDBI/NDWI
%     corr         Pearson r (and p) of NDVI, NDBI, NDWI versus LST
%     regression   LST ~ NDVI + NDBI + NDWI (coefficients, SE, t, p, R2, RMSE)
%     class_stats  per land-cover class statistics
%     water_profile / veg_profile   distance-zone cooling profiles
%     scatter_idx  reproducible subsample used for scatter plots

    lc   = L.landcover;
    lst  = S.lst(:);
    ndvi = S.ndvi(:);
    ndbi = S.ndbi(:);
    ndwi = S.ndwi(:);
    N    = numel(lst);

    % ---------------------------------------------------------------
    % 1) Global descriptive statistics
    % ---------------------------------------------------------------
    A.summary.lst_mean  = mean(lst);
    A.summary.lst_std   = std(lst);
    A.summary.lst_min   = min(lst);
    A.summary.lst_max   = max(lst);
    A.summary.ndvi_mean = mean(ndvi);
    A.summary.ndbi_mean = mean(ndbi);
    A.summary.ndwi_mean = mean(ndwi);
    A.summary.n_cells   = N;

    % ---------------------------------------------------------------
    % 2) Pearson correlations with LST
    % ---------------------------------------------------------------
    names = {'NDVI', 'NDBI', 'NDWI'};
    preds = [ndvi, ndbi, ndwi];
    A.corr.names = names;
    A.corr.r = zeros(1, 3);
    A.corr.p = zeros(1, 3);
    for k = 1:3
        R = corrcoef(preds(:, k), lst);
        r = R(1, 2);
        t = r * sqrt((N - 2) / (1 - r^2));
        A.corr.r(k) = r;
        A.corr.p(k) = h.t_pvalue(t, N - 2);
    end
    Rall = corrcoef([ndvi, ndbi, ndwi, lst]);
    A.corr.matrix = Rall;
    A.corr.matrix_names = {'NDVI', 'NDBI', 'NDWI', 'LST'};

    % ---------------------------------------------------------------
    % 3) Multiple linear regression: LST ~ NDVI + NDBI + NDWI
    %    Ordinary least squares in base MATLAB (no Statistics Toolbox).
    % ---------------------------------------------------------------
    X    = [ones(N, 1), preds];
    beta = X \ lst;
    fitted = X * beta;
    resid  = lst - fitted;
    dof    = N - size(X, 2);
    sse    = sum(resid.^2);
    sst    = sum((lst - mean(lst)).^2);
    sigma2 = sse / dof;
    covb   = sigma2 * inv(X' * X);
    se     = sqrt(diag(covb));
    tval   = beta ./ se;
    pval   = h.t_pvalue(tval, dof);

    A.regression.terms    = {'Intercept', 'NDVI', 'NDBI', 'NDWI'};
    A.regression.beta     = beta(:).';
    A.regression.se       = se(:).';
    A.regression.t        = tval(:).';
    A.regression.p        = pval(:).';
    A.regression.r2       = 1 - sse / sst;
    A.regression.adj_r2   = 1 - (sse / dof) / (sst / (N - 1));
    A.regression.rmse     = sqrt(sse / N);
    A.regression.dof      = dof;
    % Standardised coefficients (effect of +1 SD in the predictor, in SD of LST)
    A.regression.beta_std = beta(2:4).' .* std(preds, 0, 1) / std(lst);

    % Simple regression lines used in the scatter plots
    A.fit_ndvi = polyfit(ndvi, lst, 1);
    A.fit_ndbi = polyfit(ndbi, lst, 1);

    % ---------------------------------------------------------------
    % 4) Land-cover class statistics
    % ---------------------------------------------------------------
    cell_area_km2 = (cfg.cell_size_m / 1000)^2;
    K = 5;
    cs_stats = struct('name', {}, 'n', {}, 'area_km2', {}, 'percent', {}, ...
        'lst_mean', {}, 'lst_std', {}, 'lst_min', {}, 'lst_max', {}, ...
        'ndvi_mean', {}, 'ndbi_mean', {}, 'ndwi_mean', {});
    for k = 1:K
        m = (lc(:) == k);
        cs_stats(k).name      = L.class_names{k};
        cs_stats(k).n         = sum(m);
        cs_stats(k).area_km2  = sum(m) * cell_area_km2;
        cs_stats(k).percent   = 100 * sum(m) / N;
        cs_stats(k).lst_mean  = mean(lst(m));
        cs_stats(k).lst_std   = std(lst(m));
        cs_stats(k).lst_min   = min(lst(m));
        cs_stats(k).lst_max   = max(lst(m));
        cs_stats(k).ndvi_mean = mean(ndvi(m));
        cs_stats(k).ndbi_mean = mean(ndbi(m));
        cs_stats(k).ndwi_mean = mean(ndwi(m));
    end
    A.class_stats = cs_stats;

    % ---------------------------------------------------------------
    % 5) Distance-zone cooling profiles
    %    Water profile: land cells only (water cells excluded, otherwise the
    %    result would be trivially cold). Vegetation profile: cells that are
    %    neither water nor dense vegetation.
    % ---------------------------------------------------------------
    A.water_profile = zone_profile(S.dist_water_m, S.lst, lc ~= 1, cfg.zone_edges_m);
    A.veg_profile   = zone_profile(S.dist_dense_m, S.lst, (lc ~= 1) & (lc ~= 2), cfg.zone_edges_m);

    % Correlation between distance and LST inside the analysed range
    A.water_profile.r_distance_lst = distance_correlation(S.dist_water_m, S.lst, lc ~= 1, cfg.zone_edges_m(end));
    A.veg_profile.r_distance_lst   = distance_correlation(S.dist_dense_m, S.lst, (lc ~= 1) & (lc ~= 2), cfg.zone_edges_m(end));

    % ---------------------------------------------------------------
    % 6) Reproducible random subsample for scatter plots
    % ---------------------------------------------------------------
    [~, order] = sort(h.white_rand(N, cfg.seed + 301));
    A.scatter_idx = order(1:min(cfg.scatter_points, N));
end

% -------------------------------------------------------------------------
function P = zone_profile(dist_m, values, valid, edges)
%ZONE_PROFILE  Mean of values within distance zones (lo < d <= hi).
    nz = numel(edges) - 1;
    P.edges  = edges;
    P.labels = cell(1, nz);
    P.mid_m  = zeros(1, nz);
    P.mean   = nan(1, nz);
    P.std    = nan(1, nz);
    P.sem    = nan(1, nz);
    P.n      = zeros(1, nz);
    for z = 1:nz
        lo = edges(z);
        hi = edges(z + 1);
        m = valid & (dist_m > lo) & (dist_m <= hi);
        v = values(m);
        P.labels{z} = sprintf('%d-%d m', lo, hi);
        P.mid_m(z)  = (lo + hi) / 2;
        P.n(z)      = numel(v);
        if ~isempty(v)
            P.mean(z) = mean(v);
            P.std(z)  = std(v);
            P.sem(z)  = std(v) / sqrt(numel(v));
        end
    end
    P.cooling_first_to_last_C = P.mean(end) - P.mean(1);   % positive = warmer far away
end

% -------------------------------------------------------------------------
function r = distance_correlation(dist_m, lst, valid, max_dist)
    m = valid & (dist_m > 0) & (dist_m <= max_dist);
    R = corrcoef(dist_m(m), lst(m));
    r = R(1, 2);
end
