%% MAIN_SIMULATION  Spatial Simulation of Urban Heat Patterns (synthetic data)
%
%  Run this script (the ONLY one you need):   >> main_simulation
%
%  Workflow
%    1. create the synthetic landscape       (generate_landscape.m)
%    2. simulate NDVI / NDBI / NDWI and LST  (simulate_lst.m)
%    3. spatial + statistical analysis       (analyze_spatial_relationships.m)
%    4. create and save all figures          (create_figures.m)
%    5. save results (.mat) and summary CSV
%    6. print a concise summary in the Command Window
%
%  NOTE: everything here is SYNTHETIC. The outputs are NOT satellite-derived.
%  Requirements: base MATLAB R2016b or newer. No toolboxes are needed.

clear; clc; close all;
t_start = tic;

%% ------------------------------------------------------------------------
%  CONFIGURATION  (all tunable parameters live here)
%  ------------------------------------------------------------------------
cfg.seed        = 42;      % master random seed (whole project is deterministic)
cfg.grid_size   = 500;     % grid is grid_size x grid_size cells
cfg.cell_size_m = 10;      % cell size in metres (500 cells -> 5 km x 5 km)
cfg.figure_dpi  = 300;     % resolution of saved PNG figures
cfg.scatter_points = 8000; % number of cells drawn in scatter plots

% ---- Landscape (see generate_landscape.m) ---------------------------------
% Urban cores: [centre_x  centre_y  sigma  weight]   (normalised 0-1 units)
cfg.urban_cores = [0.50 0.45 0.15 1.00;
                   0.20 0.25 0.08 0.55;
                   0.80 0.70 0.09 0.50];
% Urban corridors (roads / ribbon development):
%                  [x1   y1   x2   y2   half_width  weight]
cfg.corridors   = [0.50 0.45 0.20 0.25 0.012 0.55;
                   0.50 0.45 0.80 0.70 0.012 0.55;
                   0.05 0.45 0.95 0.45 0.010 0.40];
cfg.urban_noise_weight     = 0.12;   % texture added to the urban field
cfg.builtup_fraction       = 0.26;   % share of non-water cells that are built-up
cfg.bare_fraction_of_land  = 0.13;   % share of remaining land that is bare soil
cfg.dense_share_of_vegetation = 0.50;% share of vegetated land that is dense

cfg.river_y0        = 0.78;          % mean river position (normalised)
cfg.river_halfwidth = 0.011;         % river half-width (normalised units)
cfg.lake_center     = [0.78 0.20];   % lake centre (normalised)
cfg.lake_sigma      = 0.085;
cfg.lake_roughness  = 0.30;          % irregularity of the lake shoreline
cfg.lake_threshold  = 0.55;

% ---- Spectral indices (class-typical values, see simulate_lst.m) -----------
%                      Water  Dense  Sparse  Built  Bare
cfg.class_ndvi = [     -0.25   0.78    0.42  0.10  0.15];
cfg.class_ndbi = [     -0.30  -0.38   -0.16  0.28  0.12];
cfg.class_ndwi = [      0.55  -0.50   -0.30 -0.20 -0.25];
cfg.ndvi_noise_sd = 0.05;
cfg.ndbi_noise_sd = 0.05;
cfg.ndwi_noise_sd = 0.05;
cfg.optical_blur_cells = 1.2;        % sensor mixing (sigma in cells)

% ---- LST model parameters (degrees C) -----------------------------------
cfg.base_temp_C     = 31.0;          % T0: reference temperature of open land
cfg.a_built         = 7.0;           % warming by built-up fraction
cfg.a_soil          = 4.0;           % warming by bare-soil fraction
cfg.a_dense         = 6.0;           % cooling by dense vegetation fraction
cfg.a_sparse        = 2.5;           % cooling by sparse vegetation fraction
cfg.a_water         = 8.0;           % local cooling by water fraction
cfg.c_water         = 3.0;           % amplitude of water cooling halo
cfg.decay_water_m   = 250;           % e-folding distance of water halo (m)
cfg.c_veg           = 1.5;           % amplitude of vegetation cooling halo
cfg.decay_veg_m     = 150;           % e-folding distance of vegetation halo (m)
cfg.uhi_amplitude_C = 2.5;           % urban heat island amplitude
cfg.uhi_blur_cells  = 15;            % smoothing of the urban-intensity field
cfg.gradient_x_C    = 1.2;           % west-to-east warming across the scene
cfg.gradient_y_C    = -0.8;          % south-to-north change across the scene
cfg.lst_noise_sd    = 0.7;           % spatially correlated noise (SD)
cfg.thermal_blur_cells = 3;          % thermal mixing (sigma in cells)

% ---- Distance zones for the cooling analysis (metres) -----------------------
cfg.zone_edges_m = [0 50 100 200 300 500 1000];

%% ------------------------------------------------------------------------
%  PATHS  (folders are created automatically)
%  ------------------------------------------------------------------------
src_dir  = fileparts(mfilename('fullpath'));
root_dir = fileparts(src_dir);
addpath(src_dir);
fig_dir  = fullfile(root_dir, 'figures');
res_dir  = fullfile(root_dir, 'results');

h = helper_functions();
h.ensure_folder(fig_dir);
h.ensure_folder(res_dir);
rng(cfg.seed);   % fixed seed (all fields also use a deterministic built-in generator)

fprintf('\n=== Urban Heat Spatial Simulation (synthetic data) ===\n');

%% 1) Landscape --------------------------------------------------------------
fprintf('[1/5] Generating synthetic landscape ...\n');
L = generate_landscape(cfg, h);

%% 2) NDVI / NDBI / NDWI / LST ----------------------------------------------
fprintf('[2/5] Simulating NDVI, NDBI, NDWI and LST ...\n');
S = simulate_lst(cfg, L, h);

%% 3) Analysis -------------------------------------------------------------
fprintf('[3/5] Analysing spatial relationships ...\n');
A = analyze_spatial_relationships(cfg, L, S, h);

%% 4) Figures --------------------------------------------------------------
fprintf('[4/5] Creating figures ...\n');
fig_files = create_figures(cfg, L, S, A, h, fig_dir);

%% 5) Save results ---------------------------------------------------------
fprintf('[5/5] Saving results ...\n');

% --- MAT file (large grids stored as single precision to keep the file small)
landcover        = L.landcover;                       %#ok<NASGU>
ndvi             = single(S.ndvi);                    %#ok<NASGU>
ndbi             = single(S.ndbi);                    %#ok<NASGU>
ndwi             = single(S.ndwi);                    %#ok<NASGU>
lst              = single(S.lst);                     %#ok<NASGU>
lst_anomaly      = single(S.lst_anomaly);             %#ok<NASGU>
dist_water_m     = single(S.dist_water_m);            %#ok<NASGU>
dist_dense_veg_m = single(S.dist_dense_m);            %#ok<NASGU>
lst_components   = structfun(@single, S.components, 'UniformOutput', false); %#ok<NASGU>
class_names      = L.class_names;                     %#ok<NASGU>
x_km             = L.x_km;                            %#ok<NASGU>
y_km             = L.y_km;                            %#ok<NASGU>
analysis         = rmfield(A, 'scatter_idx');         %#ok<NASGU>
mat_file = fullfile(res_dir, 'simulation_results.mat');
save(mat_file, 'cfg', 'landcover', 'ndvi', 'ndbi', 'ndwi', 'lst', 'lst_anomaly', ...
    'dist_water_m', 'dist_dense_veg_m', 'lst_components', 'class_names', ...
    'x_km', 'y_km', 'analysis', '-v7');

% --- CSV file --------------------------------------------------------------
csv_file = fullfile(res_dir, 'summary_statistics.csv');
fid = fopen(csv_file, 'w');
fprintf(fid, 'category,metric,value,unit\n');
w = @(cat, met, val, unit) fprintf(fid, '%s,%s,%.6g,%s\n', cat, met, val, unit);

w('global', 'grid_size_cells',   cfg.grid_size,          'cells');
w('global', 'cell_size',         cfg.cell_size_m,        'm');
w('global', 'lst_mean',          A.summary.lst_mean,     'degC');
w('global', 'lst_std',           A.summary.lst_std,      'degC');
w('global', 'lst_min',           A.summary.lst_min,      'degC');
w('global', 'lst_max',           A.summary.lst_max,      'degC');
w('global', 'ndvi_mean',         A.summary.ndvi_mean,    'unitless');
w('global', 'ndbi_mean',         A.summary.ndbi_mean,    'unitless');
w('global', 'ndwi_mean',         A.summary.ndwi_mean,    'unitless');

for k = 1:3
    w('correlation', ['pearson_r_' A.corr.names{k} '_vs_LST'], A.corr.r(k), 'unitless');
    w('correlation', ['p_value_'   A.corr.names{k} '_vs_LST'], A.corr.p(k), 'unitless');
end

for k = 1:4
    term = A.regression.terms{k};
    w('regression', ['coef_'    term], A.regression.beta(k), 'degC per unit');
    w('regression', ['std_err_' term], A.regression.se(k),   'degC per unit');
    w('regression', ['t_stat_'  term], A.regression.t(k),    'unitless');
    w('regression', ['p_value_' term], A.regression.p(k),    'unitless');
end
w('regression', 'R2',      A.regression.r2,     'unitless');
w('regression', 'adj_R2',  A.regression.adj_r2, 'unitless');
w('regression', 'RMSE',    A.regression.rmse,   'degC');

for k = 1:5
    c = A.class_stats(k);
    tag = strrep(c.name, ' ', '_');
    w(['class_' tag], 'area',      c.area_km2,  'km2');
    w(['class_' tag], 'percent',   c.percent,   'percent');
    w(['class_' tag], 'lst_mean',  c.lst_mean,  'degC');
    w(['class_' tag], 'lst_std',   c.lst_std,   'degC');
    w(['class_' tag], 'ndvi_mean', c.ndvi_mean, 'unitless');
    w(['class_' tag], 'ndbi_mean', c.ndbi_mean, 'unitless');
    w(['class_' tag], 'ndwi_mean', c.ndwi_mean, 'unitless');
end

profiles = {A.water_profile, 'water_distance'; A.veg_profile, 'dense_veg_distance'};
for p = 1:2
    P = profiles{p, 1};
    for z = 1:numel(P.mean)
        tag = strrep(P.labels{z}, ' ', '');
        w(['zone_' profiles{p, 2}], ['lst_mean_' tag], P.mean(z), 'degC');
        w(['zone_' profiles{p, 2}], ['n_cells_'  tag], P.n(z),    'cells');
    end
    w(['zone_' profiles{p, 2}], 'warming_first_to_last_zone', P.cooling_first_to_last_C, 'degC');
    w(['zone_' profiles{p, 2}], 'pearson_r_distance_vs_LST',  P.r_distance_lst,          'unitless');
end
fclose(fid);

%% 6) Command-window summary -------------------------------------------------
deg = char(176);
fprintf('\n---------------- SUMMARY ----------------\n');
fprintf('Grid            : %d x %d cells @ %d m  (%.1f x %.1f km)\n', cfg.grid_size, ...
    cfg.grid_size, cfg.cell_size_m, cfg.grid_size*cfg.cell_size_m/1000, cfg.grid_size*cfg.cell_size_m/1000);
fprintf('LST             : mean %.2f %sC | SD %.2f | min %.2f | max %.2f\n', ...
    A.summary.lst_mean, deg, A.summary.lst_std, A.summary.lst_min, A.summary.lst_max);
fprintf('Mean indices    : NDVI %.3f | NDBI %.3f | NDWI %.3f\n', ...
    A.summary.ndvi_mean, A.summary.ndbi_mean, A.summary.ndwi_mean);
fprintf('Pearson r (LST): NDVI %+.3f | NDBI %+.3f | NDWI %+.3f\n', A.corr.r);
fprintf('Regression      : LST = %.2f %+.2f*NDVI %+.2f*NDBI %+.2f*NDWI\n', A.regression.beta);
fprintf('                  R2 = %.3f | RMSE = %.2f %sC | n = %d\n', ...
    A.regression.r2, A.regression.rmse, deg, A.summary.n_cells);
fprintf('Mean LST by class:\n');
for k = 1:5
    fprintf('   %-18s %6.2f %sC  (%4.1f %% of area)\n', A.class_stats(k).name, ...
        A.class_stats(k).lst_mean, deg, A.class_stats(k).percent);
end
fprintf('Water cooling profile (mean LST by distance zone):\n');
for z = 1:numel(A.water_profile.mean)
    fprintf('   %-10s %6.2f %sC  (n = %d)\n', A.water_profile.labels{z}, ...
        A.water_profile.mean(z), deg, A.water_profile.n(z));
end
fprintf('Outputs         : %d figures -> %s\n', numel(fig_files), fig_dir);
fprintf('                  %s\n                  %s\n', mat_file, csv_file);
fprintf('Elapsed time    : %.1f s\n', toc(t_start));
fprintf('------------------------------------------\n');
