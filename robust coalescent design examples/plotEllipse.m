% Plots an ellipse in parametric form given mean and covariance matrix
function [x, y] = plotEllipse(M, C, s)

% Mean gives ellipse centre
xCen = M(1); yCen = M(2);

% Covariance and confidence value (s) gives radii
rad = sqrt(s*C); rad = diag(rad);
xRad = rad(1); yRad = rad(2);

% Parametric space
t = 0:0.01:2*pi;
% Ellipse points ((x - xCen)/xRad)^2 + ((y - yCen)/yRad)^2 = 1
x = xRad*cos(t) + xCen;
y = yRad*sin(t) + yCen;