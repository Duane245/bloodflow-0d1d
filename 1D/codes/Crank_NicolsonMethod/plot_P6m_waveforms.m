%% Run P6_m and plot the 4 capture locations
P6_m;

time = (1:number_of_time_step)*dt;

fig = figure('Position',[100 100 1200 800], 'Visible','off');

subplot(4,2,1); plot(time, captureQ1); grid on;
xlabel('t (s)'); ylabel('Q [cm^3/s]');
title(sprintf('Q at vessel 1, node 20  (max=%.1f)', max(captureQ1)));

subplot(4,2,2); plot(time, captureA1); grid on;
xlabel('t (s)'); ylabel('A [cm^2]');
title(sprintf('A at vessel 1, node 20  (max=%.2f)', max(captureA1)));

subplot(4,2,3); plot(time, captureQ2); grid on;
xlabel('t (s)'); ylabel('Q [cm^3/s]');
title(sprintf('Q at vessel 8, node 117  (max=%.2f)', max(captureQ2)));

subplot(4,2,4); plot(time, captureA2); grid on;
xlabel('t (s)'); ylabel('A [cm^2]');
title(sprintf('A at vessel 8, node 117  (max=%.3f)', max(captureA2)));

subplot(4,2,5); plot(time, captureQ3); grid on;
xlabel('t (s)'); ylabel('Q [cm^3/s]');
title(sprintf('Q at vessel 37, node 5  (max=%.2f)', max(captureQ3)));

subplot(4,2,6); plot(time, captureA3); grid on;
xlabel('t (s)'); ylabel('A [cm^2]');
title(sprintf('A at vessel 37, node 5  (max=%.3f)', max(captureA3)));

subplot(4,2,7); plot(time, captureQ4); grid on;
xlabel('t (s)'); ylabel('Q [cm^3/s]');
title(sprintf('Q at vessel 54, node 161  (max=%.2f)', max(captureQ4)));

subplot(4,2,8); plot(time, captureA4); grid on;
xlabel('t (s)'); ylabel('A [cm^2]');
title(sprintf('A at vessel 54, node 161  (max=%.3f)', max(captureA4)));

sgtitle(sprintf('P6\\_m  55-vessel CN+MacCormack,  dt=%.0e,  %d steps,  %.1f s physical', ...
    dt, number_of_time_step, number_of_time_step*dt));

outpng = 'P6m_waveforms.png';
exportgraphics(fig, outpng, 'Resolution', 120);
fprintf('saved: %s\n', outpng);
