% Draw a multivariate Gaussian ellipse to illustrate optimality designs
clearvars; clc; close all

% Assumptions and notes
% - a 2 segment population and coalescents divide between them
% - D and E optimality of interest
% - assume Fisher confidence intervals from Cramer Rao bound

% Error ellipse package
addpath(genpath('/Users/kris/Documents/MATLAB'));

% Mean and vector for 2 populations
N1 = 100; N2 = 200;
%M = [N1 N2]; Mlog = log(M);

M = [0 0]; Mlog = [0 0];

% Total no. coalescents, n and split variable m
n = 100;
mset = [round(linspace(1, n-1, 10)), n/2];
mset = unique(sort(mset)); % ensure have optima
idmid = find(mset == n/2);
lenm = length(mset);

% Want 95% confidence so use relative chi-squared
s = 9.210; % 99%
%s = 5.991; % 95%

% Ellipse coordinates with m values
xset = cell(1, lenm); yset = cell(1, lenm);
xset0 = cell(1, lenm); yset0 = cell(1, lenm);

% Get set of log ellipses at various m
for i = 1:lenm
    % Covariance at choice of m (based on inverse Fisher)
    Clog = diag([1/mset(i), 1/(n-mset(i))]);
    C = diag([N1^2/mset(i), N2^2/(n - mset(i))]);
    % Ellipses for normal and log parametrisations
    [xset{i}, yset{i}] = plotEllipse(Mlog, Clog, s);
    [xset0{i}, yset0{i}] = plotEllipse(M, C, s);
end

% Concatenate cells for plotting
xcat = cat(1, xset{:}); xcat = xcat';
ycat = cat(1, yset{:}); ycat = ycat';
xcat0 = cat(1, xset0{:}); xcat0 = xcat0';
ycat0 = cat(1, yset0{:}); ycat0 = ycat0';

% Separate optimal ellipse
xrest = xcat(:, setdiff(1:lenm, idmid));
yrest = ycat(:, setdiff(1:lenm, idmid));
xopt = xcat(:, idmid); yopt = ycat(:, idmid);
xrest0 = xcat0(:, setdiff(1:lenm, idmid));
yrest0 = ycat0(:, setdiff(1:lenm, idmid));
xopt0 = xcat0(:, idmid); yopt0 = ycat0(:, idmid);

% Also get an e optimal ellipse for normal case
me = round(n/((N2/N1)^2 + 1));
C0e = diag([N1^2/me, N2^2/(n - me)]);
[x0e, y0e] = plotEllipse(M, C0e, s);


% Plot log ellipses
figure;
plot(xrest, yrest, '-', 'color', [0.8 0.8 0.8], 'linewidth', 2);
hold on;
plot(xopt, yopt, 'r', 'linewidth', 2);
hold off; grid off; box off; 
axis square;
h = gca;
h.XLim = [min(min(xcat))-0.2, max(max(xcat))+0.2];
h.YLim = h.XLim;
xlabel('log(N_1)'); ylabel('log(N_2)'); 

% Plot normal ellipses
figure;
plot(xrest0, yrest0, '-', 'color', [0.8 0.8 0.8], 'linewidth', 2);
hold on;
plot(xopt0, yopt0, 'r', 'linewidth', 2);
plot(x0e, y0e, 'color', [0.31 0.31 0.31], 'linewidth', 2);
hold off; grid off; box off; 
axis square;
h = gca;
h.XLim = [min(min(xcat0))-20, max(max(xcat0))+20];
h.YLim = h.XLim;
xlabel('N_1'); ylabel('N_2'); 









