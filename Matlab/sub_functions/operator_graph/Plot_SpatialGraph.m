function Plot_SpatialGraph(n1, n2, G_or_ijw, varargin)
% PlotSpatialGraph
% -------------------------------------------------------------------------
% Visualize a spatial graph over an n1-by-n2 grid.
% Accepts either:
%   - G_sp: [N x E] weighted incidence matrix (sparse), or
%   - ijw : [E x 3] edge list, columns = [i, j, w], 1-based nodes
%
% Nodes are placed on a regular grid:
%   Node index p maps to (x,y) with:
%       y = mod(p-1, n1) + 1
%       x = floor((p-1)/n1) + 1
%
% OPTIONS (name-value):
%   'MaxEdges'     (default: 20000)  % draw at most this many strongest edges
%   'WeightThresh' (default: 0)      % draw edges with w >= this threshold
%   'EdgeAlpha'    (default: 0.25)   % transparency for edges
%   'EdgeColor'    (default: [0 0.447 0.741]) % MATLAB default blue
%   'LineWidth'    (default: [])     % if empty, auto scales by weight
%   'ShowNodes'    (default: true)
%   'NodeMarker'   (default: '.')
%   'NodeAlpha'    (default: 0.6)
%   'FlipY'        (default: true)   % image-style y axis (top=1)
%
% Example:
%   PlotSpatialGraph(n1, n2, G_sp, 'MaxEdges', 5000, 'WeightThresh', 0.1);
%   PlotSpatialGraph(n1, n2, info.ijw, 'EdgeAlpha', 0.15, 'ShowNodes', false);
% -------------------------------------------------------------------------

% ----- parse options -----
p = inputParser;
addParameter(p, 'MaxEdges', 20000);
addParameter(p, 'WeightThresh', 0);
addParameter(p, 'EdgeAlpha', 0.25);
addParameter(p, 'EdgeColor', [0 0.447 0.741]);
addParameter(p, 'LineWidth', []);
addParameter(p, 'ShowNodes', true);
addParameter(p, 'NodeMarker', '.');
addParameter(p, 'NodeAlpha', 0.6);
addParameter(p, 'FlipY', true);
parse(p, varargin{:});
opt = p.Results;

N = n1*n2;

% ----- build edge list [i j w] -----
if ismatrix(G_or_ijw) && size(G_or_ijw,2)==3
    ijw = G_or_ijw;
else
    G = G_or_ijw;
    if ~issparse(G), error('G must be sparse or pass ijw as [E x 3].'); end
    % Each column e has +w at node i and -w at node j
    [rows, cols, vals] = find(G);
    % group by column (edge)
    [e_unique, ~, grp] = unique(cols, 'stable');
    if numel(e_unique) == 0
        warning('No edges to plot.'); return;
    end
    % Per edge, find positive row (tail) and negative row (head); weight = abs(val)
    % Accumulate: first positive then negative (2 entries per edge)
    % Safer approach: split by sign then join by columns
    posMask = vals > 0;
    negMask = vals < 0;
    % map col -> row for + and -
    i_pos = accumarray(cols(posMask), rows(posMask), [max(cols) 1], @(x) x(1), 0);
    i_neg = accumarray(cols(negMask), rows(negMask), [max(cols) 1], @(x) x(1), 0);
    w_col = accumarray(cols, abs(vals), [max(cols) 1], @(x) max(x), 0);
    keep  = (i_pos>0) & (i_neg>0) & (w_col>0);
    ijw   = [i_pos(keep), i_neg(keep), w_col(keep)];
end

% ----- filter by threshold -----
if opt.WeightThresh > 0
    ijw = ijw(ijw(:,3) >= opt.WeightThresh, :);
end
if isempty(ijw)
    warning('All edges filtered out by threshold.'); return;
end

% ----- keep strongest edges (by weight) -----
E = size(ijw,1);
if E > opt.MaxEdges
    [~,ord] = maxk(ijw(:,3), opt.MaxEdges);
    ijw = ijw(ord,:);
    E = size(ijw,1);
end

% ----- node coordinates -----
% p -> (x,y)
toXY = @(p) deal( floor((p-1)/n1) + 1 , mod(p-1, n1) + 1 );
[i_x, i_y] = toXY(ijw(:,1));
[j_x, j_y] = toXY(ijw(:,2));

% ----- build line segments for vectorized plotting -----
% Make arrays [x1 x2 NaN] and [y1 y2 NaN] per edge and plot in one call
Xseg = [i_x, j_x, nan(E,1)].';
Yseg = [i_y, j_y, nan(E,1)].';

% ----- draw -----
holdState = ishold;
hold on;

% auto linewidth scaling by weight (optional)
if isempty(opt.LineWidth)
    % normalize weights to [0.5, 2.5] for visibility
    w = ijw(:,3);
    w = (w - min(w)) / max(eps, (max(w)-min(w)));
    lw = 0.5 + 2.0 * w;
else
    lw = opt.LineWidth * ones(E,1);
end

% draw edges (one batch per unique LW for performance, or simple loop)
% Simple and reasonably fast:
h = plot(Xseg(:), Yseg(:), '-', 'Color', opt.EdgeColor, 'LineWidth', 1);
% apply alpha by patch object (line alpha via patch workaround)
if verLessThan('matlab','9.6')
    % older MATLAB: alpha on lines is limited; ignore
else
    set(h, 'Color', [opt.EdgeColor, opt.EdgeAlpha]);
end

% If per-edge linewidth is desired, redraw in chunks (slower):
% for b = 1:E
%     plot([i_x(b) j_x(b)], [i_y(b) j_y(b)], '-', ...
%         'Color', [opt.EdgeColor opt.EdgeAlpha], 'LineWidth', lw(b));
% end

% draw nodes
if opt.ShowNodes
    [XX,YY] = meshgrid(1:n2, 1:n1);
    scatter(XX(:), YY(:), 6, 'k', opt.NodeMarker, 'MarkerEdgeAlpha', opt.NodeAlpha);
end

axis equal tight;
xlim([0.5 n2+0.5]); ylim([0.5 n1+0.5]);

if opt.FlipY
    set(gca,'YDir','reverse'); % image-style coordinates (row 1 at top)
end
box on; grid on;
xlabel('x (column)'); ylabel('y (row)');
title(sprintf('Spatial Graph (|E| = %d)', E));

if ~holdState, hold off; end
end
