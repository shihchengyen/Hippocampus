function map_smooth = laplace_smoothing(map, x_bins, y_bins, pillar_centers, pillar_half_length, n_iter)
%LAPLACE_SMOOTHING Smooth a place map using graph Laplacian diffusion
%
%   map_smooth = laplace_smoothing(map, x_bins, y_bins, pillar_centers, pillar_half_length, n_iter)
%
%   INPUTS:
%       map               - 1x1600 raw firing rate map (NaN for unvisited/pillar bins)
%       x_bins            - 1x40 bin centers in x direction
%       y_bins            - 1x40 bin centers in y direction
%       pillar_centers    - Nx2 matrix of [x, y] pillar centers
%                           e.g. [-5,-5; -5,5; 5,5; 5,-5]
%       pillar_half_length- half side length of each pillar (e.g. 2.5)
%       n_iter            - number of Laplacian smoothing iterations
%
%   OUTPUT:
%       map_smooth        - 1x1600 smoothed firing rate map
%
%   NOTES:
%       - Bins inside pillars are excluded from smoothing
%       - Boundary bins (edges and pillar edges) only connect to valid neighbors
%       - The Laplacian matrix L = D^-1 * A (row-normalized)
%         so each bin gets the average of its neighbors at each iteration
%
%   EXAMPLE:
%       x_bins = linspace(-12.1875, 12.1875, 40);
%       y_bins = linspace(-12.1875, 12.1875, 40);
%       pillar_centers = [-5,-5; -5,5; 5,5; 5,-5];
%       map_smooth = laplace_smoothing(map, x_bins, y_bins, pillar_centers, 2.5, 10);

    n_x = length(x_bins);
    n_y = length(y_bins);
    n_bins = n_x * n_y; % 1600

    % ------------------------------------------------------------------ %
    %  1. Identify pillar bins
    %     Bins are indexed as: bin = (col-1)*n_y + row
    %     where col = x index, row = y index
    %     i.e. column-major order matching MATLAB's default
    % ------------------------------------------------------------------ %
    is_pillar = false(n_y, n_x); % (row=y, col=x)

    for p = 1:size(pillar_centers, 1)
        px = pillar_centers(p, 1);
        py = pillar_centers(p, 2);
        % Find bins whose centers fall inside pillar boundaries
        in_x = (x_bins >= px - pillar_half_length) & (x_bins <= px + pillar_half_length);
        in_y = (y_bins >= py - pillar_half_length) & (y_bins <= py + pillar_half_length);
        is_pillar(in_y, in_x) = true;
    end

    % Flatten to 1D (column-major)
    is_pillar_flat = is_pillar(:); % 1600x1

    % ------------------------------------------------------------------ %
    %  2. Build adjacency matrix A (1600 x 1600)
    %     4-connectivity: up, down, left, right
    %     No connections to/from pillar bins
    %     No connections across edges
    % ------------------------------------------------------------------ %
    fprintf('Building adjacency matrix...\n');

    % Pre-allocate sparse adjacency matrix
    % Max 4 neighbors per bin, so max 4*1600 = 6400 entries
    I = zeros(4 * n_bins, 1);
    J = zeros(4 * n_bins, 1);
    idx = 0;

    for col = 1:n_x       % x direction
        for row = 1:n_y   % y direction

            bin = (col-1)*n_y + row;

            % Skip pillar bins entirely
            if is_pillar_flat(bin)
                continue;
            end

            % Check 4 neighbors: up, down, left, right
            neighbors = [...
                row-1, col;   % down  (y decreases)
                row+1, col;   % up    (y increases)
                row,   col-1; % left  (x decreases)
                row,   col+1; % right (x increases)
            ];

            for n = 1:size(neighbors, 1)
                nr = neighbors(n, 1);
                nc = neighbors(n, 2);

                % Check within grid bounds
                if nr < 1 || nr > n_y || nc < 1 || nc > n_x
                    continue;
                end

                neighbor_bin = (nc-1)*n_y + nr;

                % Skip if neighbor is a pillar bin
                if is_pillar_flat(neighbor_bin)
                    continue;
                end

                % Add edge
                idx = idx + 1;
                I(idx) = bin;
                J(idx) = neighbor_bin;
            end
        end
    end

    % Trim pre-allocated zeros
    I = I(1:idx);
    J = J(1:idx);
    A = sparse(I, J, ones(idx,1), n_bins, n_bins);

    % ------------------------------------------------------------------ %
    %  3. Build symmetric normalized Laplacian: Ls = I - Dp * A * Dp
    %     where Dp_i = 1/sqrt(sum(A(:,i)))  (column sum = degree)
    % ------------------------------------------------------------------ %
    degree = full(sum(A, 1))'; % column sums, 1600x1
    degree(degree == 0) = 1;   % avoid division by zero for isolated/pillar bins
    Dp = spdiags(1./sqrt(degree), 0, n_bins, n_bins);
    Ls = speye(n_bins) - Dp * A * Dp;

    % ------------------------------------------------------------------ %
    %  4. Apply smoothing iteratively
    %     map_smooth = L^n_iter * map
    %     NaN bins (pillar + unvisited) handled by replacing with 0,
    %     smoothing, then restoring NaNs for pillar bins
    % ------------------------------------------------------------------ %
    map_col = map(:); % ensure column vector 1600x1

    % Track which bins are NaN in original map (unvisited)
    nan_mask = isnan(map_col);

    % Replace NaN with 0 for matrix multiplication
    map_col(nan_mask) = 0;

    % Zero out pillar bins so they don't contribute
    map_col(is_pillar_flat) = 0;

    % Iterative smoothing: repeatedly apply (I - Ls) = Dp*A*Dp
    % which is the diffusion operator derived from Ls
    diffusion = speye(n_bins) - Ls; % = Dp * A * Dp

    for iter = 1:n_iter
        map_col = diffusion * map_col;
        % Re-zero pillar bins after each iteration
        map_col(is_pillar_flat) = 0;
    end

    % Restore NaN for pillar bins
    map_col(is_pillar_flat) = NaN;

    % Restore NaN for bins that were unvisited in original map
    % (optional: comment out if you want smoothing to fill unvisited bins)
    map_col(nan_mask) = NaN;

    % Return as row vector to match input format
    map_smooth = map_col';

end