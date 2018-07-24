% Plot the prctiles around a mean or median
function plotErrBnd(currAx, x, y, yL, yU, col)

% Assumptions and notes
% - uses boundedline package and plots on currAx
% - upper and lower bounds yL, yU are absolute (not rel to y)
% - y can be mean, median or any central measure
% - column vectors expected

% Check input structure
if ~iscolumn(x) || ~iscolumn(y) || ~iscolumn(yL) || ~iscolumn(yU)
    error('Expect column vectors - transposing');
end

% Bounded line with formatting and outline
axes(currAx);
[l,p] = boundedline(x, y, [y-yL yU-y], '-');
p.FaceAlpha = 0.08; p.FaceColor = col;
l.Color = col; l.LineWidth = 1.5;
outlinebounds(l,p);