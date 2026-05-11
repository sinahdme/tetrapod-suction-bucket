function fg_vec = detect_fg_per_bucket(s_end_cells, s_fixed_cells, adj_matrix, SD_ratio, delta_ref, LD_ratio, config)
% DETECT_FG_PER_BUCKET  Displacement-dependent group p-multiplier
%
% f_g ramps linearly from 1.0 (no reduction) to f_g_min (full reduction)
% as the differential compression displacement increases:
%
%   f_g_i = 1 - (1 - f_g_min) * min(delta_i / delta_ref, 1)
%
% where delta_i = how much more bucket i settles than the group mean (m)
%
% f_g_min is computed from the enriched Barari et al. (2023) database using
% 2D interpolation of E_cM (moment group efficiency) across all L/D and S/D.
% Falls back to the original single-equation method if barari_database.m
% is not available.
%
% Inputs:
%   s_end_cells    - cell array of body-frame positions (3xN per bucket)
%   s_fixed_cells  - cell array of fixed-frame positions (3xN per bucket)
%   adj_matrix     - n x n logical adjacency matrix
%   SD_ratio       - spacing-to-diameter ratio S/D
%   delta_ref      - reference displacement (m) at which full reduction is reached
%   LD_ratio       - (optional) embedment ratio L/D, default 1.0
%   config         - (optional) 'tri' or 'tetra', default 'tetra'
%
% Returns fg_vec: per-bucket f_g values (between f_g_min and 1.0)

% Handle optional arguments for backward compatibility
if nargin < 6 || isempty(LD_ratio)
    LD_ratio = 1.0;
end
if nargin < 7 || isempty(config)
    config = 'tetra';
end

n_buckets = length(s_end_cells);

% Compute f_g_min from enriched database or fallback
f_g_min = compute_fg_min(SD_ratio, LD_ratio, config);

% Step 1: Compute vertical displacement per bucket
disp_z = zeros(1, n_buckets);
for i = 1:n_buckets
    disp_z(i) = mean(abs(s_fixed_cells{i}(3,:)) - abs(s_end_cells{i}(3,:)));
end

% Step 2: Relative compression detection (above group mean = compression)
mean_disp = mean(disp_z);
is_compression = disp_z > mean_disp;

% Step 3: Displacement-dependent f_g for compression buckets with adjacent neighbors
fg_vec = ones(1, n_buckets);
for i = 1:n_buckets
    if is_compression(i)
        has_adj_compression = false;
        for j = 1:n_buckets
            if j ~= i && is_compression(j) && adj_matrix(i,j)
                has_adj_compression = true;
                break;
            end
        end
        if has_adj_compression
            delta_i = disp_z(i) - mean_disp;  % differential compression (m)
            ratio = min(delta_i / delta_ref, 1.0);
            fg_vec(i) = 1 - (1 - f_g_min) * ratio;
        end
    end
end

end


function f_g_min = compute_fg_min(SD_ratio, LD_ratio, config)
% COMPUTE_FG_MIN  Compute f_g_min from enriched database with fallback
%
% Uses barari_database() for 2D interpolation of E_cM across L/D and S/D.
% The group penalty is derived from the ratio of E_cM at the given S/D
% to E_cM at S/D=5.0 (where group interaction is negligible).
%
% Scale factor 0.723 calibrated at S/D=3.0, L/D=1.0 where f_g_min=0.70.

scale = 0.723;

try
    db = barari_database();

    % Get E_cM at requested S/D and at reference S/D=5.0
    EcM_current = db.lookup.EcM(config, LD_ratio, SD_ratio);
    EcM_ref     = db.lookup.EcM(config, LD_ratio, 5.0);

    % Guard against division by zero or negative reference
    if EcM_ref <= 0
        EcM_ref = 1.0;
    end

    penalty = max(0, 1 - EcM_current / EcM_ref);
    f_g_min = max(0.30, 1 - scale * penalty);
    f_g_min = min(1.0, f_g_min);

catch
    % Fallback: original single-equation method (tetrapod, L/D=1.0)
    eta_M = 0.05 * SD_ratio^2 + 1.00 * SD_ratio + 0.50;
    eta_M_ref = 6.75;  % eta_M at S/D=5.0
    penalty = max(0, 1 - eta_M / eta_M_ref);
    f_g_min = max(0.30, 1 - scale * penalty);
    f_g_min = min(1.0, f_g_min);
end

end
