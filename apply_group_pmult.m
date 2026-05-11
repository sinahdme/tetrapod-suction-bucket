function [KKx_out, KKy_out, fg_applied] = ...
    apply_group_pmult(KKx_cells, KKy_cells, s_end_cells, s_fixed_cells, ...
    adj_matrix, f_g)
% APPLY_GROUP_PMULT  Apply group interaction p-multiplier to p-y forces
%
% Compression detection: displacement = abs(z_fixed) - abs(z_body)
%   positive = bucket pushed deeper = compression
%
% f_g applied only to compression buckets with adjacent compression neighbors
%
% Inputs:
%   KKx_cells, KKy_cells  - cell arrays of p-y force vectors per bucket
%   s_end_cells            - cell array of body node positions (3xN)
%   s_fixed_cells          - cell array of fixed node positions (3xN)
%   adj_matrix             - n_buckets x n_buckets logical adjacency matrix
%   f_g                    - p-multiplier (0 < f_g <= 1)
%
% Outputs:
%   KKx_out, KKy_out      - cell arrays of modified p-y forces
%   fg_applied             - logical array indicating which buckets got f_g

n_buckets = length(KKx_cells);

% Step 1: Detect compression per bucket
is_compression = false(1, n_buckets);
for i = 1:n_buckets
    disp_z = mean(abs(s_fixed_cells{i}(3,:)) - abs(s_end_cells{i}(3,:)));
    is_compression(i) = disp_z > 0;  % positive = pushed deeper
end

% Step 2: Check for adjacent compression neighbors
fg_applied = false(1, n_buckets);
for i = 1:n_buckets
    if is_compression(i)
        for j = 1:n_buckets
            if j ~= i && is_compression(j) && adj_matrix(i,j)
                fg_applied(i) = true;
                break;
            end
        end
    end
end

% Step 3: Apply f_g
KKx_out = KKx_cells;
KKy_out = KKy_cells;
for i = 1:n_buckets
    if fg_applied(i)
        KKx_out{i} = f_g * KKx_out{i};
        KKy_out{i} = f_g * KKy_out{i};
    end
end

end
