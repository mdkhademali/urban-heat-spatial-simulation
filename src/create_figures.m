function fig_files = create_figures(cfg, L, S, A, h, fig_dir)
%CREATE_FIGURES  Produce all publication-style figures and save them as PNG.
%
%   fig_files = create_figures(cfg, L, S, A, h, fig_dir)
%
%   Figures are rendered off-screen (Visible = 'off') with a white
%   background, consistent fonts and 300 dpi (summary figure: 250 dpi).

    h.ensure_folder(fig_dir);
    fig_files = {};

    % ---- shared style --------------------------------------------------
    st.font   = 'Helvetica';
    st.fs     = 11;
    st.dpi    = cfg.figure_dpi;
    st.x_km   = L.x_km;
    st.y_km   = L.y_km;

    % ---- colormaps -----------------------------------------------------
    cm.landcover = [0.16 0.42 0.80;    % 1 water
                    0.05 0.36 0.14;    % 2 dense vegetation
                    0.58 0.77 0.36;    % 3 sparse vegetation
                    0.52 0.52 0.56;    % 4 built-up
                    0.86 0.72 0.50];   % 5 bare soil
    cm.ndvi = h.make_colormap([0.35 0.24 0.15; 0.80 0.68 0.42; 0.96 0.93 0.62; ...
                               0.55 0.78 0.32; 0.10 0.50 0.18; 0.02 0.30 0.10], 256);
    cm.ndbi = h.make_colormap([0.13 0.36 0.62; 0.62 0.80 0.90; 0.97 0.97 0.95; ...
                               0.98 0.73 0.40; 0.80 0.20 0.10], 256);
    cm.lst  = h.make_colormap([0.10 0.20 0.55; 0.20 0.55 0.80; 0.60 0.85 0.75; ...
                               0.98 0.93 0.55; 0.98 0.62 0.22; 0.80 0.12 0.12; 0.40 0.02 0.10], 256);
    cm.anom = h.make_colormap([0.10 0.25 0.60; 0.45 0.68 0.85; 0.97 0.97 0.97; ...
                               0.98 0.65 0.45; 0.72 0.10 0.12], 256);

    % =====================================================================
    % Figure 1 - land cover
    % =====================================================================
    [fig, ax] = new_figure(st, 7.4, 6.2);
    imagesc(ax, st.x_km, st.y_km, double(L.landcover));
    style_map(ax, st);
    colormap(ax, cm.landcover);
    caxis(ax, [0.5 5.5]);
    cb = colorbar(ax);
    set(cb, 'Ticks', 1:5, 'TickLabels', L.class_names);
    label_colorbar(cb, 'Land-cover class');
    title(ax, 'Figure 1. Synthetic Land-Cover Map', 'FontWeight', 'bold');
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '01_synthetic_landcover.png'), st.dpi);

    % =====================================================================
    % Figure 2 - NDVI
    % =====================================================================
    [fig, ax] = new_figure(st, 7.4, 6.2);
    imagesc(ax, st.x_km, st.y_km, S.ndvi);
    style_map(ax, st);
    colormap(ax, cm.ndvi);
    caxis(ax, [-0.4 0.9]);
    cb = colorbar(ax);
    label_colorbar(cb, 'NDVI (unitless)');
    title(ax, 'Figure 2. Simulated NDVI', 'FontWeight', 'bold');
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '02_ndvi.png'), st.dpi);

    % =====================================================================
    % Figure 3 - NDBI
    % =====================================================================
    [fig, ax] = new_figure(st, 7.4, 6.2);
    imagesc(ax, st.x_km, st.y_km, S.ndbi);
    style_map(ax, st);
    colormap(ax, cm.ndbi);
    caxis(ax, [-0.5 0.5]);
    cb = colorbar(ax);
    label_colorbar(cb, 'NDBI (unitless)');
    title(ax, 'Figure 3. Simulated NDBI', 'FontWeight', 'bold');
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '03_ndbi.png'), st.dpi);

    % =====================================================================
    % Figure 4 - LST
    % =====================================================================
    [fig, ax] = new_figure(st, 7.4, 6.2);
    imagesc(ax, st.x_km, st.y_km, S.lst);
    style_map(ax, st);
    colormap(ax, cm.lst);
    lst_lim = [h.percentile_value(S.lst, 0.5), h.percentile_value(S.lst, 99.5)];
    caxis(ax, lst_lim);
    cb = colorbar(ax);
    label_colorbar(cb, ['Land Surface Temperature (' char(176) 'C)']);
    title(ax, 'Figure 4. Simulated Land Surface Temperature (LST)', 'FontWeight', 'bold');
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '04_lst.png'), st.dpi);

    % =====================================================================
    % Figure 5 - LST anomaly (heat-island pattern)
    % =====================================================================
    [fig, ax] = new_figure(st, 7.4, 6.2);
    imagesc(ax, st.x_km, st.y_km, S.lst_anomaly);
    style_map(ax, st);
    colormap(ax, cm.anom);
    amax = max(abs([h.percentile_value(S.lst_anomaly, 0.5), h.percentile_value(S.lst_anomaly, 99.5)]));
    caxis(ax, [-amax amax]);
    cb = colorbar(ax);
    label_colorbar(cb, ['LST anomaly relative to scene mean (' char(176) 'C)']);
    title(ax, 'Figure 5. LST Anomaly / Heat-Island Pattern', 'FontWeight', 'bold');
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '05_lst_anomaly.png'), st.dpi);

    % =====================================================================
    % Figure 6 - NDVI vs LST
    % =====================================================================
    idx = A.scatter_idx;
    [fig, ax] = new_figure(st, 6.8, 5.6);
    plot(ax, S.ndvi(idx), S.lst(idx), '.', 'MarkerSize', 4, 'Color', [0.20 0.45 0.75]);
    hold(ax, 'on');
    xx = [min(S.ndvi(:)), max(S.ndvi(:))];
    plot(ax, xx, polyval(A.fit_ndvi, xx), '-', 'Color', [0.80 0.10 0.10], 'LineWidth', 2.2);
    style_axes(ax, st);
    xlabel(ax, 'NDVI (unitless)');
    ylabel(ax, ['LST (' char(176) 'C)']);
    title(ax, 'Figure 6. NDVI vs LST', 'FontWeight', 'bold');
    txt = sprintf('r = %.3f\nslope = %.2f %sC per NDVI unit\nn = %d cells (%d plotted)', ...
        A.corr.r(1), A.fit_ndvi(1), char(176), A.summary.n_cells, numel(idx));
    text(ax, 0.03, 0.05, txt, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
        'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'FontSize', st.fs - 1);
    legend(ax, {'Sampled cells', 'OLS regression line'}, 'Location', 'northeast', 'Box', 'off');
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '06_ndvi_lst_scatter.png'), st.dpi);

    % =====================================================================
    % Figure 7 - NDBI vs LST
    % =====================================================================
    [fig, ax] = new_figure(st, 6.8, 5.6);
    plot(ax, S.ndbi(idx), S.lst(idx), '.', 'MarkerSize', 4, 'Color', [0.85 0.45 0.15]);
    hold(ax, 'on');
    xx = [min(S.ndbi(:)), max(S.ndbi(:))];
    plot(ax, xx, polyval(A.fit_ndbi, xx), '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 2.2);
    style_axes(ax, st);
    xlabel(ax, 'NDBI (unitless)');
    ylabel(ax, ['LST (' char(176) 'C)']);
    title(ax, 'Figure 7. NDBI vs LST', 'FontWeight', 'bold');
    txt = sprintf('r = %.3f\nslope = %.2f %sC per NDBI unit\nn = %d cells (%d plotted)', ...
        A.corr.r(2), A.fit_ndbi(1), char(176), A.summary.n_cells, numel(idx));
    text(ax, 0.03, 0.95, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'FontSize', st.fs - 1);
    legend(ax, {'Sampled cells', 'OLS regression line'}, 'Location', 'southeast', 'Box', 'off');
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '07_ndbi_lst_scatter.png'), st.dpi);

    % =====================================================================
    % Figure 8 - mean LST by land-cover class
    % =====================================================================
    cs_stats = A.class_stats;
    means = [cs_stats.lst_mean];
    sds   = [cs_stats.lst_std];
    [fig, ax] = new_figure(st, 7.6, 5.6);
    hold(ax, 'on');
    for k = 1:5
        bar(ax, k, means(k), 0.65, 'FaceColor', cm.landcover(k, :), 'EdgeColor', [0.2 0.2 0.2]);
    end
    errorbar(ax, 1:5, means, sds, 'k.', 'LineWidth', 1.4, 'CapSize', 9);
    for k = 1:5
        text(ax, k, means(k) + sds(k) + 0.35, sprintf('%.1f', means(k)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', st.fs);
    end
    style_axes(ax, st);
    set(ax, 'XTick', 1:5, 'XTickLabel', L.class_names);
    ylabel(ax, ['Mean LST (' char(176) 'C)']);
    xlabel(ax, 'Land-cover class');
    ylim(ax, [min(means - sds) - 1, max(means + sds) + 2]);
    title(ax, 'Figure 8. Mean LST by Land-Cover Class (error bars: \pm1 SD)', 'FontWeight', 'bold');
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '08_landcover_mean_lst.png'), st.dpi);

    % =====================================================================
    % Figure 9 - water distance cooling profile
    % =====================================================================
    [fig, ax] = new_figure(st, 7.4, 5.6);
    profile_plot(ax, A.water_profile, [0.16 0.42 0.80], st);
    title(ax, 'Figure 9. Cooling Profile: Distance from Water Bodies', 'FontWeight', 'bold');
    xlabel(ax, 'Distance from nearest water body (m)');
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '09_water_cooling_profile.png'), st.dpi);

    % =====================================================================
    % Figure 10 - research summary (multi-panel)
    % =====================================================================
    fig = figure('Visible', 'off', 'Color', 'w', 'InvertHardcopy', 'off', ...
        'Units', 'inches', 'Position', [0.5 0.5 13.5 10.5], ...
        'PaperUnits', 'inches', 'PaperPosition', [0 0 13.5 10.5], 'PaperSize', [13.5 10.5]);

    ax1 = subplot(2, 2, 1);
    imagesc(ax1, st.x_km, st.y_km, double(L.landcover));
    style_map(ax1, st);  colormap(ax1, cm.landcover);  caxis(ax1, [0.5 5.5]);
    cb = colorbar(ax1);  set(cb, 'Ticks', 1:5, 'TickLabels', L.class_names);
    title(ax1, '(a) Land cover', 'FontWeight', 'bold');

    ax2 = subplot(2, 2, 2);
    imagesc(ax2, st.x_km, st.y_km, S.ndvi);
    style_map(ax2, st);  colormap(ax2, cm.ndvi);  caxis(ax2, [-0.4 0.9]);
    cb = colorbar(ax2);  label_colorbar(cb, 'NDVI');
    title(ax2, '(b) Simulated NDVI', 'FontWeight', 'bold');

    ax3 = subplot(2, 2, 3);
    imagesc(ax3, st.x_km, st.y_km, S.lst);
    style_map(ax3, st);  colormap(ax3, cm.lst);  caxis(ax3, lst_lim);
    cb = colorbar(ax3);  label_colorbar(cb, ['LST (' char(176) 'C)']);
    title(ax3, '(c) Simulated LST', 'FontWeight', 'bold');

    ax4 = subplot(2, 2, 4);
    profile_plot(ax4, A.water_profile, [0.16 0.42 0.80], st);
    title(ax4, '(d) Water cooling profile', 'FontWeight', 'bold');
    xlabel(ax4, 'Distance from water (m)');

    annotation(fig, 'textbox', [0 0.965 1 0.035], 'String', ...
        'Spatial Simulation of Urban Heat Patterns: Research Summary (synthetic data)', ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 15, ...
        'FontWeight', 'bold', 'FontName', st.font);
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '10_research_summary.png'), 250);

    % =====================================================================
    % Figure 11 - dense vegetation cooling profile
    % =====================================================================
    [fig, ax] = new_figure(st, 7.4, 5.6);
    profile_plot(ax, A.veg_profile, [0.05 0.45 0.15], st);
    title(ax, 'Figure 11. Cooling Profile: Distance from Dense Vegetation', 'FontWeight', 'bold');
    xlabel(ax, 'Distance from nearest dense vegetation (m)');
    fig_files{end+1} = save_figure(fig, fullfile(fig_dir, '11_vegetation_cooling_profile.png'), st.dpi);
end

% =========================================================================
%  Local plotting helpers
% =========================================================================
function [fig, ax] = new_figure(st, w_in, h_in)
    fig = figure('Visible', 'off', 'Color', 'w', 'InvertHardcopy', 'off', ...
        'Units', 'inches', 'Position', [0.5 0.5 w_in h_in], ...
        'PaperUnits', 'inches', 'PaperPosition', [0 0 w_in h_in], ...
        'PaperSize', [w_in h_in]);
    ax = axes('Parent', fig);
end

function style_axes(ax, st)
    set(ax, 'FontName', st.font, 'FontSize', st.fs, 'LineWidth', 0.8, ...
        'Box', 'on', 'TickDir', 'out', 'Color', 'w', ...
        'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15]);
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.25);
end

function style_map(ax, st)
    set(ax, 'YDir', 'normal', 'FontName', st.font, 'FontSize', st.fs, ...
        'LineWidth', 0.8, 'Box', 'on', 'TickDir', 'out');
    axis(ax, 'image');
    xlabel(ax, 'Easting (km)');
    ylabel(ax, 'Northing (km)');
end

function label_colorbar(cb, txt)
    ylabel(cb, txt);
end

function profile_plot(ax, P, color, st)
%PROFILE_PLOT  Mean LST versus distance zone with +/-1 SD bars.
    x = 1:numel(P.mean);
    hold(ax, 'on');
    plot(ax, x, P.mean, '-o', 'Color', color, 'LineWidth', 2.4, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
    errorbar(ax, x, P.mean, P.std, 'LineStyle', 'none', 'Color', color * 0.7 + 0.3, ...
        'LineWidth', 1.2, 'CapSize', 7);
    style_axes(ax, st);
    set(ax, 'XTick', x, 'XTickLabel', P.labels, 'XLim', [0.5 numel(x) + 0.5]);
    ylabel(ax, ['Mean LST (' char(176) 'C)']);
    for k = x
        text(ax, k, P.mean(k) - P.std(k) - 0.25, sprintf('%.2f', P.mean(k)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', st.fs - 2, 'Color', [0.2 0.2 0.2]);
    end
    lo = min(P.mean - P.std) - 1.0;
    hi = max(P.mean + P.std) + 0.6;
    ylim(ax, [lo hi]);
    txt = sprintf('Warming over profile: %+.2f %sC\nr(distance, LST) = %.2f', ...
        P.cooling_first_to_last_C, char(176), P.r_distance_lst);
    text(ax, 0.97, 0.06, txt, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], ...
        'FontSize', st.fs - 1);
    legend(ax, {'Zone mean', '\pm1 SD'}, 'Location', 'northwest', 'Box', 'off');
end

function fname = save_figure(fig, filename, dpi)
    print(fig, filename, '-dpng', sprintf('-r%d', dpi));
    close(fig);
    fname = filename;
end
