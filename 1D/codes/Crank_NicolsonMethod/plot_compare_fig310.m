%% Plot the last heart beat in the Figure 3.10 layout.
% Re-uses P6m_results.mat if available; otherwise runs P6_m and saves.
results_file = 'P6m_results.mat';
if isfile(results_file)
    fprintf('loading: %s\n', results_file);
    load(results_file);
else
    P6_m;
    save(results_file, ...
        'captureQ1','captureQ2','captureQ3','captureQ4', ...
        'captureA1','captureA2','captureA3','captureA4', ...
        'A','Q','beta','A0','dt','number_of_time_step','-v7.3');
    fprintf('saved: %s\n', results_file);
end

% --- last heart beat indices (steady state) ---
% InflowFuncation uses Tc = 0.8 s (half-sine over 0..0.4 s, zero 0.4..0.8 s).
% At 8 s of physical time we have 10 beats; plot the last (7.2 - 8.0 s) to
% match Figure 3.10 of the reference.
Tc = 0.8;
steps_per_beat = round(Tc/dt);
n_beats        = floor(number_of_time_step / steps_per_beat);
beat_start     = (n_beats-1)*steps_per_beat + 1;
beat_end       = n_beats*steps_per_beat;
idx = beat_start:beat_end;
tt  = idx*dt;
t_xlim = [tt(1) tt(end)];

% --- pressure (kPa) ---
P1 = beta(1) *(sqrt(captureA1)-sqrt(A0(1))) /1000;
P2 = beta(8) *(sqrt(captureA2)-sqrt(A0(8))) /1000;
P3 = beta(37)*(sqrt(captureA3)-sqrt(A0(37)))/1000;
P4 = beta(54)*(sqrt(captureA4)-sqrt(A0(54)))/1000;

fig = figure('Position',[100 100 900 1000],'Visible','off','Color','w');

% Reference Figure 3.10 y-axis ranges (Wang & Parker / Sherwin 55-vessel set)
%   artery 1:  Q [0  500],  P [10  16]
%   artery 8:  Q [9   14],  P [11.5 15.5]
%   artery 37: Q [-50 125], P [10   15]
%   artery 54: Q [8   16],  P [10   16]
panels = {
    1, captureQ1, 'Artery 1  (aortic root)',  'Q [ml/s]',  [0 500]; ...
    2, P1,        'Artery 1',                  'P [kPa]',  [10 16]; ...
    3, captureQ2, 'Artery 8',                  'Q [ml/s]', [9 14]; ...
    4, P2,        'Artery 8',                  'P [kPa]',  [10 18]; ...
    5, captureQ3, 'Artery 37  (reverse flow)', 'Q [ml/s]', [-50 130]; ...
    6, P3,        'Artery 37',                 'P [kPa]',  [9 16]; ...
    7, captureQ4, 'Artery 54',                 'Q [ml/s]', [8 16]; ...
    8, P4,        'Artery 54',                 'P [kPa]',  [9 18]; ...
};
for k = 1:size(panels,1)
    subplot(4,2,panels{k,1});
    plot(tt, panels{k,2}(idx),'b','LineWidth',1.4); grid on; hold on;
    ax = gca; ax.XLim = t_xlim; ax.YLim = panels{k,5};
    xlabel('t [s]'); ylabel(panels{k,4}); title(panels{k,3});
end
sgtitle(sprintf('Last heart beat (beat %d, t = %.1f–%.1f s),  dt=%.0e,  MacCormack+CN hybrid', ...
    n_beats, t_xlim(1), t_xlim(2), dt));

our_png = 'P6m_fig310_ours.png';
exportgraphics(fig, our_png, 'Resolution', 130);
fprintf('saved: %s\n', our_png);

% --- side-by-side with reference Figure 3.10 ---
if isfile('reference_fig310.png')
    ours = imread(our_png);
    ref  = imread('reference_fig310.png');
    % pad/crop to same height
    h = max(size(ours,1), size(ref,1));
    ours = padarray(ours, [h-size(ours,1) 0 0], 255, 'post');
    ref  = padarray(ref,  [h-size(ref,1)  0 0], 255, 'post');
    side = [ref, 255*ones(h, 30, 3, 'uint8'), ours];

    figS = figure('Visible','off','Color','w','Position',[100 100 2400 1000]);
    imshow(side);
    title('Reference: Figure 3.10 (PDF page 50)        |        Ours: MacCormack+CN, dt=1e-4, optimized', ...
        'FontSize',13);
    exportgraphics(figS, 'P6m_fig310_compare.png', 'Resolution', 110);
    fprintf('saved: P6m_fig310_compare.png\n');
end
