% ANIMATE_3D_HEATMAPS  Generate animated GIFs of p-y and t-z 3D surfaces
%
% Run this script AFTER running THREE_POD_KIM2014_VALIDATION_LD123 or
% FOUR_POD_BARARI_VALIDATION. The simulation stores per-step data in the
% workspace variables: anim_s_fixed, anim_KKx, anim_KKy, anim_KKz, etc.
%
% Outputs:
%   py_3d_animation.gif  - p-y force surface evolving with rotation angle
%   tz_3d_animation.gif  - t-z force surface evolving with rotation angle

%% ---- Configuration ----
delay_time = 0.33;       % seconds per frame in GIF
view_angle = [0, 0]; % 3D view azimuth & elevation
out_dir = fileparts(mfilename('fullpath'));  % save GIFs next to this script

%% ---- Detect configuration from workspace ----
n_steps = length(anim_KKx);
n_buckets = length(anim_KKx{1});  % 3 for tripod, 4 for tetrapod

if n_buckets == 3
    config_name = 'Tripod';
    bucket_labels = {'B1 (60)', 'B2 (180)', 'B3 (300)'};
    th_angles = [60, 180, 300];
elseif n_buckets == 4
    config_name = 'Tetrapod';
    bucket_labels = {'B1 (45)', 'B2 (135)', 'B3 (225)', 'B4 (315)'};
    th_angles = [45, 135, 225, 315];
else
    error('Unexpected number of buckets: %d', n_buckets);
end

fprintf('Detected %s (%d buckets, %d rotation steps)\n', config_name, n_buckets, n_steps);

%% ---- Compute global color limits across all steps ----
global_max_py = 0;
global_max_tz = 0;
for s = 1:n_steps
    for b = 1:n_buckets
        fpy = sqrt(anim_KKx{s}{b}.^2 + anim_KKy{s}{b}.^2);
        global_max_py = max(global_max_py, max(fpy));
        global_max_tz = max(global_max_tz, max(abs(anim_KKz{s}{b})));
    end
end

%% ---- Blue-white-red diverging colormap ----
n_cmap = 256;
blue_white_red = [linspace(0,1,n_cmap/2)', linspace(0,1,n_cmap/2)', ones(n_cmap/2,1); ...
                  ones(n_cmap/2,1), linspace(1,0,n_cmap/2)', linspace(1,0,n_cmap/2)'];

%% ========== Animate p-y 3D surface ==========
gif_py = fullfile(out_dir, 'py_3d_animation.gif');
vid_py = fullfile(out_dir, 'py_3d_animation.mp4');
fig_py = figure('Position', [100 100 1200 800], 'Color', 'w');
vw_py = VideoWriter(vid_py, 'MPEG-4');
vw_py.FrameRate = round(1 / delay_time);
vw_py.Quality = 95;
open(vw_py);
fprintf('Generating p-y animation...\n');

for s = 1:n_steps
    clf(fig_py);
    ax = axes(fig_py);
    hold(ax, 'on');

    for b = 1:n_buckets
        pos = anim_s_fixed{s}{b};
        Xg = reshape(pos(1,:), n, n_l);
        Yg = reshape(pos(2,:), n, n_l);
        Zg = reshape(pos(3,:), n, n_l);
        Cg = reshape(sqrt(anim_KKx{s}{b}.^2 + anim_KKy{s}{b}.^2), n, n_l);
        % Close cylinder
        Xg = [Xg; Xg(1,:)]; Yg = [Yg; Yg(1,:)];
        Zg = [Zg; Zg(1,:)]; Cg = [Cg; Cg(1,:)];
        surf(ax, Xg, Yg, Zg, Cg, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.95);
    end

    % Rotation center trajectory (up to current step)
    plot3(ax, CENTE_OF_ROTATION(1,1:s), CENTE_OF_ROTATION(2,1:s), CENTE_OF_ROTATION(3,1:s), ...
        '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.5, 'MarkerSize', 7, ...
        'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k');
    scatter3(ax, CENTE_OF_ROTATION(1,s), CENTE_OF_ROTATION(2,s), CENTE_OF_ROTATION(3,s), ...
        200, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.2 0.9 0.2], 'LineWidth', 2);

    % Bucket labels
    for b = 1:n_buckets
        text(ax, l*cosd(th_angles(b)), l*sind(th_angles(b)), 1.5, bucket_labels{b}, ...
            'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end

    colormap(ax, jet);
    clim(ax, [0 global_max_py]);
    cb = colorbar(ax);
    cb.Label.String = '|F_{py}| (kN)'; cb.Label.FontSize = 14; cb.Label.FontWeight = 'bold';
    title(ax, sprintf('%s p-y Force Distribution  |  \\theta = %.2f deg', config_name, ddd(s)), ...
        'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    xlabel(ax, 'X (m)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    ylabel(ax, 'Y (m)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    zlabel(ax, 'Z (m)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    set(ax, 'FontSize', 12, 'FontWeight', 'bold');
    grid(ax, 'on'); axis(ax, 'equal'); view(ax, view_angle);
    hold(ax, 'off');
    drawnow;

    % Capture frame
    frame = getframe(fig_py);
    writeVideo(vw_py, frame);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);
    if s == 1
        imwrite(A, map, gif_py, 'gif', 'LoopCount', Inf, 'DelayTime', delay_time);
    else
        imwrite(A, map, gif_py, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
    fprintf('\r  p-y frame %d / %d', s, n_steps);
end
close(vw_py);
fprintf('\n  Saved: %s\n', gif_py);
fprintf('  Saved: %s\n', vid_py);

%% ========== Animate t-z 3D surface ==========
gif_tz = fullfile(out_dir, 'tz_3d_animation.gif');
vid_tz = fullfile(out_dir, 'tz_3d_animation.mp4');
fig_tz = figure('Position', [100 100 1200 800], 'Color', 'w');
vw_tz = VideoWriter(vid_tz, 'MPEG-4');
vw_tz.FrameRate = round(1 / delay_time);
vw_tz.Quality = 98;
open(vw_tz);
fprintf('Generating t-z animation...\n');

for s = 1:n_steps
    clf(fig_tz);
    ax = axes(fig_tz);
    hold(ax, 'on');

    for b = 1:n_buckets
        pos = anim_s_fixed{s}{b};
        Xg = reshape(pos(1,:), n, n_l);
        Yg = reshape(pos(2,:), n, n_l);
        Zg = reshape(pos(3,:), n, n_l);
        Cg = reshape(anim_KKz{s}{b}, n, n_l);
        % Close cylinder
        Xg = [Xg; Xg(1,:)]; Yg = [Yg; Yg(1,:)];
        Zg = [Zg; Zg(1,:)]; Cg = [Cg; Cg(1,:)];
        surf(ax, Xg, Yg, Zg, Cg, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.95);
    end

    % Rotation center trajectory
    plot3(ax, CENTE_OF_ROTATION(1,1:s), CENTE_OF_ROTATION(2,1:s), CENTE_OF_ROTATION(3,1:s), ...
        '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.5, 'MarkerSize', 7, ...
        'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k');
    scatter3(ax, CENTE_OF_ROTATION(1,s), CENTE_OF_ROTATION(2,s), CENTE_OF_ROTATION(3,s), ...
        200, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.2 0.9 0.2], 'LineWidth', 2);

    % Bucket labels
    for b = 1:n_buckets
        text(ax, l*cosd(th_angles(b)), l*sind(th_angles(b)), 1.5, bucket_labels{b}, ...
            'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end

    colormap(ax, blue_white_red);
    clim(ax, [-global_max_tz, global_max_tz]);
    cb = colorbar(ax);
    cb.Label.String = 'F_{tz} (kN)'; cb.Label.FontSize = 14; cb.Label.FontWeight = 'bold';
    title(ax, sprintf('%s t-z Force Distribution  |  \\theta = %.2f deg', config_name, ddd(s)), ...
        'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    xlabel(ax, 'X (m)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    ylabel(ax, 'Y (m)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    zlabel(ax, 'Z (m)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    set(ax, 'FontSize', 12, 'FontWeight', 'bold');
    grid(ax, 'on'); axis(ax, 'equal'); view(ax, view_angle);
    hold(ax, 'off');
    drawnow;

    % Capture frame
    frame = getframe(fig_tz);
    writeVideo(vw_tz, frame);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);
    if s == 1
        imwrite(A, map, gif_tz, 'gif', 'LoopCount', Inf, 'DelayTime', delay_time);
    else
        imwrite(A, map, gif_tz, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
    fprintf('\r  t-z frame %d / %d', s, n_steps);
end
close(vw_tz);
fprintf('\n  Saved: %s\n', gif_tz);
fprintf('  Saved: %s\n', vid_tz);

fprintf('\nDone! Generated:\n  %s\n  %s\n  %s\n  %s\n', gif_py, vid_py, gif_tz, vid_tz);
