%% EEG signal simulation for teaching
% Goal:
% 1) Start from pure beta.
% 2) Add one component at each stage.
% 3) Clearly mark where the new change appears.

clear; clc; close all;

%% 1. Basic settings
fs = 1000;                  % sampling rate (Hz)
duration = 3;               % total duration (s)
t = 0:1/fs:duration-1/fs;
N = numel(t);
rng(1);

%% 2. Build components
center = 1.50;
short_win = [1.35 1.65];
mid_win = [1.20 1.80];
long_win = [1.05 1.95];

beta = 8 * sin(2*pi*20*t);

alpha = 6 * sin(2*pi*10*t) .* smooth_gate(t, mid_win(1), mid_win(2), 0.03);

delta = 14 * sin(2*pi*2.2*t) .* smooth_gate(t, long_win(1), long_win(2), 0.06);

gamma = 2.8 * sin(2*pi*40*t) .* smooth_gate(t, short_win(1), short_win(2), 0.02);

[b_emg, a_emg] = butter(4, [35 100]/(fs/2), 'bandpass');
emg = filtfilt(b_emg, a_emg, randn(size(t)));
emg = emg / std(emg);
emg = 14 * emg .* smooth_gate(t, short_win(1), short_win(2), 0.015);

blink = 90 * exp(-((t - center)/0.13).^2) ...
      - 28 * exp(-((t - (center + 0.12))/0.06).^2);

motion_rise = smooth_gate(t, long_win(1), center, 0.04);
motion_fall = smooth_gate(t, center, long_win(2), 0.06);
motion = 20 * motion_rise + 10 * motion_fall ...
       + 4 * pink_noise_simple(N, 1.0) .* smooth_gate(t, long_win(1), long_win(2), 0.05);

line_noise = 7 * sin(2*pi*50*t) .* smooth_gate(t, short_win(1), short_win(2), 0.01);

%% 3. Cumulative stages: add one component each time
stage_added = {
    zeros(size(t))
    alpha
    delta
    gamma
    emg
    blink
    motion
    line_noise
    };

stage_signals = cell(1, 8);
stage_signals{1} = beta;
for k = 2:8
    stage_signals{k} = stage_signals{k-1} + stage_added{k};
end

stage_titles = {
    '1. Pure Beta'
    '2. + Alpha'
    '3. + Delta'
    '4. + Gamma'
    '5. + EMG'
    '6. + EOG / Blink'
    '7. + Motion artifact'
    '8. + 50 Hz line noise'
    };

change_windows = {
    []
    mid_win
    long_win
    short_win
    short_win
    [1.25 1.75]
    long_win
    short_win
    };

added_titles = {
    ''
    'Added component: Alpha'
    'Added component: Delta'
    'Added component: Gamma'
    'Added component: EMG'
    'Added component: EOG / Blink'
    'Added component: Motion artifact'
    'Added component: 50 Hz line noise'
    };

%% 4. Main overview: all stages
figure('Color', 'w', 'Position', [70 30 1300 980]);
for k = 1:8
    subplot(8,1,k);
    plot(t, stage_signals{k}, 'k', 'LineWidth', 0.95); hold on;
    xlim([0 duration]);
    yk = max(abs(stage_signals{k})) * 1.15 + eps;
    ylim([-yk yk]);
    grid on;
    ylabel('\muV');
    title(stage_titles{k});
    if ~isempty(change_windows{k})
        mark_change_region(gca, change_windows{k}, yk);
    end
end
xlabel('Time (s)');

%% 5. What was added at each stage
figure('Color', 'w', 'Position', [90 40 1200 920]);
for k = 2:8
    subplot(7,1,k-1);
    plot(t, stage_added{k}, 'Color', [0.85 0.2 0.2], 'LineWidth', 0.95); hold on;
    xlim([0 duration]);
    yk = max(abs(stage_added{k})) * 1.15 + eps;
    ylim([-yk yk]);
    grid on;
    ylabel('\muV');
    title(added_titles{k});
    mark_change_region(gca, change_windows{k}, yk);
end
xlabel('Time (s)');

%% 6. Before/after comparisons for each added component
for k = 2:8
    make_comparison_figure( ...
        t, ...
        stage_signals{k-1}, ...
        stage_signals{k}, ...
        stage_added{k}, ...
        change_windows{k}, ...
        stage_titles{k-1}, ...
        stage_titles{k}, ...
        added_titles{k});
end

%% 7. PSD of first and last stage
figure('Color', 'w', 'Position', [120 100 900 500]);
[p_first, f_first] = pwelch(stage_signals{1}, hanning(2*fs), fs, 2048, fs);
[p_last, f_last] = pwelch(stage_signals{8}, hanning(2*fs), fs, 2048, fs);
loglog(f_first, p_first, 'LineWidth', 1.4); hold on;
loglog(f_last, p_last, 'LineWidth', 1.4);
grid on;
xlim([1 120]);
xlabel('Frequency (Hz)');
ylabel('Power spectral density');
title('PSD: Pure Beta vs Final contaminated signal');
legend('Pure Beta', 'Final stage', 'Location', 'southwest');
xline(10, '--', 'Alpha');
xline(20, '--', 'Beta');
xline(40, '--', 'Gamma');
xline(50, '--r', '50 Hz');

%% ===== Local functions =====
function make_comparison_figure(t, sig_before, sig_after, added_comp, change_win, left_title, right_title, mid_title)
    figure('Color', 'w', 'Position', [100 80 1250 760]);

    subplot(3,1,1);
    plot(t, sig_before, 'k', 'LineWidth', 0.95);
    xlim([0 t(end)]);
    y1 = max(abs(sig_before)) * 1.15 + eps;
    ylim([-y1 y1]);
    grid on;
    ylabel('\muV');
    title(left_title);
    if ~isempty(change_win)
        mark_change_region(gca, change_win, y1);
    end

    subplot(3,1,2);
    plot(t, added_comp, 'Color', [0.85 0.2 0.2], 'LineWidth', 0.95); hold on;
    xlim([0 t(end)]);
    y2 = max(abs(added_comp)) * 1.15 + eps;
    ylim([-y2 y2]);
    grid on;
    ylabel('\muV');
    title(mid_title);
    mark_change_region(gca, change_win, y2);

    subplot(3,1,3);
    plot(t, sig_before, 'Color', [0.65 0.65 0.65], 'LineWidth', 0.9); hold on;
    plot(t, sig_after, 'k', 'LineWidth', 1.0);
    xlim([0 t(end)]);
    y3 = max(abs([sig_before sig_after])) * 1.15 + eps;
    ylim([-y3 y3]);
    grid on;
    ylabel('\muV');
    xlabel('Time (s)');
    title(right_title);
    mark_change_region(gca, change_win, y3);
    legend('Before', 'After', 'Location', 'northwest');
end

function mark_change_region(ax, win, ymax)
    axes(ax);
    patch([win(1) win(2) win(2) win(1)], ...
          [-ymax -ymax ymax ymax], ...
          [1.0 0.92 0.92], ...
          'EdgeColor', 'none', ...
          'FaceAlpha', 0.35);
    xline(win(1), '--r', 'start', 'LineWidth', 1);
    xline(win(2), '--r', 'end', 'LineWidth', 1);
    uistack(findobj(ax, 'Type', 'Line'), 'top');
end

function g = smooth_gate(t, t1, t2, tau)
    g = 0.5 * (tanh((t - t1)/tau) - tanh((t - t2)/tau));
end

function x = pink_noise_simple(N, alpha)
    if nargin < 2
        alpha = 1.0;
    end

    w = randn(1, N);
    W = fft(w);
    f = [0:floor(N/2), -floor((N-1)/2):-1];
    scale = 1 ./ (abs(f) + 1).^(alpha/2);
    X = W .* scale;
    x = real(ifft(X));
    x = x - mean(x);
    x = x / std(x);
end
