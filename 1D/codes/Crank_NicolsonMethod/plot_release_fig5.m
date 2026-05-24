%% Regenerate 1D_case5_computed.png and 1D_case5_comparison.png
% Matches the style of BloodFlow_release/docs/figures/1D_case5_computed.png:
%   blue lines, artery-number badge top-left, grid on, headers "Flux"/"Pressure".

RESULTS_FILE = 'P6m_results.mat';
if ~isfile(RESULTS_FILE)
    P6_m;  % NB: clears workspace, so RESULTS_FILE must be re-assigned below
    save('P6m_results.mat', ...
        'captureQ1','captureQ2','captureQ3','captureQ4', ...
        'captureA1','captureA2','captureA3','captureA4', ...
        'A','Q','beta','A0','dt','number_of_time_step','-v7.3');
end
results_file = 'P6m_results.mat';
S = load(results_file);
dt = S.dt; number_of_time_step = S.number_of_time_step;
beta = S.beta; A0 = S.A0;
captureQ1 = S.captureQ1; captureQ2 = S.captureQ2;
captureQ3 = S.captureQ3; captureQ4 = S.captureQ4;
captureA1 = S.captureA1; captureA2 = S.captureA2;
captureA3 = S.captureA3; captureA4 = S.captureA4;

% Tc = 0.8 s; 10 beats total in 8 s; last beat (steady state) is 7.2–8.0 s.
Tc = 0.8;
steps_per_beat = round(Tc/dt);
n_beats        = floor(number_of_time_step/steps_per_beat);
beat_idx = (n_beats-1)*steps_per_beat + 1 : n_beats*steps_per_beat;
tt = beat_idx*dt;

P1 = beta(1) *(sqrt(captureA1)-sqrt(A0(1))) /1000;
P2 = beta(8) *(sqrt(captureA2)-sqrt(A0(8))) /1000;
P3 = beta(37)*(sqrt(captureA3)-sqrt(A0(37)))/1000;
P4 = beta(54)*(sqrt(captureA4)-sqrt(A0(54)))/1000;

art_ids = [1, 8, 37, 54];
Q_traces = {captureQ1, captureQ2, captureQ3, captureQ4};
P_traces = {P1, P2, P3, P4};

fig = figure('Visible','off','Color','w','Position',[100 100 600 720]);

for r = 1:4
    % Flux subplot
    axQ = subplot(4,2,2*r-1);
    plot(tt, Q_traces{r}(beat_idx), 'b', 'LineWidth', 1.6); hold on; grid on;
    set(axQ,'XLim',[tt(1) tt(end)],'GridAlpha',0.25,'FontSize',9);
    ylabel('ml/s','FontSize',10);
    if r < 4, set(axQ,'XTickLabel',[]); else, xlabel('seconds','FontSize',10); end
    if r == 1, title('Flux','FontSize',12,'FontWeight','bold'); end
    yl = ylim; xl = xlim;
    text(xl(1)+(xl(2)-xl(1))*0.04, yl(2)-(yl(2)-yl(1))*0.10, num2str(art_ids(r)), ...
        'FontWeight','bold','FontSize',11, 'BackgroundColor','w','Margin',1, 'EdgeColor','none');

    % Pressure subplot
    axP = subplot(4,2,2*r);
    plot(tt, P_traces{r}(beat_idx), 'b', 'LineWidth', 1.6); hold on; grid on;
    set(axP,'XLim',[tt(1) tt(end)],'GridAlpha',0.25,'FontSize',9);
    ylabel('kPa','FontSize',10);
    if r < 4, set(axP,'XTickLabel',[]); else, xlabel('seconds','FontSize',10); end
    if r == 1, title('Pressure','FontSize',12,'FontWeight','bold'); end
    yl = ylim; xl = xlim;
    text(xl(1)+(xl(2)-xl(1))*0.04, yl(2)-(yl(2)-yl(1))*0.10, num2str(art_ids(r)), ...
        'FontWeight','bold','FontSize',11, 'BackgroundColor','w','Margin',1, 'EdgeColor','none');
end

release_figs = '/data/dsq1/BloodFlow/BloodFlow_release/docs/figures';
computed_png = fullfile(release_figs, '1D_case5_computed.png');
exportgraphics(fig, computed_png, 'Resolution', 150);
fprintf('saved: %s\n', computed_png);

% --- side-by-side with literature panel via tiledlayout ---
lit_png = fullfile(release_figs, '1D_case5_literature.png');
if isfile(lit_png) && isfile(computed_png)
    lit  = imread(lit_png);
    comp = imread(computed_png);

    figS = figure('Visible','off','Color','w','Position',[50 50 1600 1100]);
    tl = tiledlayout(figS, 1, 2, 'Padding','compact', 'TileSpacing','compact');

    nexttile(tl);  imshow(lit);
    title('Literature — Wang 2014, Fig. 3.10 (4-method overlay)', 'FontSize',12);
    nexttile(tl);  imshow(comp);
    title(sprintf('Computed — 1D/codes/Crank\\_NicolsonMethod/P6\\_m.m  (data: P6m\\_results.mat, dt=%g)', dt), ...
        'FontSize',12);

    title(tl, '1D case 5 — 55-vessel systemic arterial network', ...
        'FontSize',16, 'FontWeight','bold');
    subtitle(tl, 'Flux (ml/s) and Pressure (kPa) at arteries 1, 8, 37, 54 over one cardiac cycle (10th heartbeat shown)', ...
        'FontSize',11, 'Color',[0.25 0.25 0.25]);

    out_png = fullfile(release_figs, '1D_case5_comparison.png');
    exportgraphics(figS, out_png, 'Resolution', 130);
    fprintf('saved: %s\n', out_png);
end
