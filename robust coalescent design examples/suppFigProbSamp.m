% Extra plots for supplement
clearvars; clc; close all;

% Assumptions and notes
% - plots based around square wave robust design figure
% - uses settings from probSampDetDiffN

% Period (1/1000 months), amp = Nmax/Nmin
T = 12000;
amp = 2;

% High and low populations
N1 = [T/8 T/4 T/2 T]; 
N2 = amp*N1; lenN = length(N1);

% Define square waves for t
t = linspace(0, 2.5*T, 1000);
lent = length(t); Nt = zeros(4, lent);
for i = 1:lenN
    Nt(i, :) = getSqWave(t, T, 0.5, N1(i), N2(i));
end

% Colours of curves
col = {'k', 'r', 'g', 'm'};

% Plot all 4 trajectories
figure;
hold on
for i = 1:lenN
    plot(t/T, Nt(i, :)/T, 'color', col{i}, 'linewidth', 2);
end
hold off
xlabel('relative time, t/T'); ylabel('scaled population size, N(t)/T');
grid off; box off;
legend('T/8', 'T/4', 'T/2', 'T', 'location', 'best');
ylim([0.1 2.1]);

% Conceptual trajectory
figure;
plot(t/T, Nt(2, :)/T, 'color', col{2}, 'linewidth', 2);
xlabel('relative time, t/T'); ylabel('scaled population size, N(t)/T');
grid off; box off;
ylim([0.1 2.1]);


% Plot sampling strategy for specific n case
idplt = 1;
% Strip meanSampIntro of small values and round
figRan = round(linspace(1, length(q1), 3));
meanSampRnd = round(D{idplt}.meanSampIntro(figRan, :), 2);

% Plot how the sample distributions look in terms of half periods
figure;
% Control figure size
smallFig = length(figRan); 
for jj = 1:3
    subplot(3, 1, jj);
    h = stem(meanSampRnd(jj, :), 'linewidth', 2);
    % Accessorise colours
    h.Color = [0.8 0.8 0.8];
    h.MarkerEdgeColor = 'm';
    h.Marker = '.';
    h.MarkerSize = 22;
   
    box off; grid off;
    %set(gca,'ytick',[]);
    if jj ~= 3
        set(gca,'xtick',[]);
    else
        xlabel('half-period number');
    end
    ylabel(['p_1 = ' num2str(round(D{idplt}.p1(figRan(jj)), 2))]);
end
