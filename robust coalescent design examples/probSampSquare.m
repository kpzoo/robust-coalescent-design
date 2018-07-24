% Deterministically place groups or single samples along half-periods
clearvars
clc; close all;

% Assumptions and notes
% - varies N1 as a fraction of T instead of samples introduced
% - all populations are in log form so aim for even distribution
% - period of 1 year, allows change of population ratio
% - sampling schemes are deterministic, and uniform
% - p1 = 0 means all samples in N1, p1 = 1 in N2
% - p1 = a means a*(n-1) in N1 and 1-p1 in N2

%% Square wave and sample prob setup

% Decide if saving data and if using mean or bin noise with max constraint
wantSave = 0;
useMean = 0;

% Clock the time when started code and path
clk = datetime('now'); tic;
filePath = [pwd,'/',mfilename];

% Partition and notBoxPlot package <---- choose directory
addpath(genpath('/Users/kris/Documents/MATLAB'));

% Period (1/1000 months), amp = Nmax/Nmin
T = 12000;
amp = 2;

% High and low populations
%N1 = [T/8 T/4 T/2 T];
N1 = T/8;
N2 = amp*N1; lenN = length(N1);
% Want log populations
N1 = log(N1); N2 = log(N2);

% Sample probs, p1 = tip in N1 region, p2 = 1-p1
p1 = 0:0.01:1;
p2 = 1 - p1;
nPrs = length(p1);

% Size of tree and repetitions
M = 5000; n = 100;

% Max sample no. samples introduced each T/2
nMaxSet = 1;
% Get coalescent binomial factors for lineages and make no. lineages an
% index so fac(1) = 0, fac(k) = k choose 2
nset = n:-1:2;
fac = nset.*(nset-1)/2;
fac = sort([0 fac]);

for kk = 1:lenN
    
    % Store sample schemes for each p1, and no. times tintro
    sampCell = cell(nPrs, M);
    tintro = zeros(nPrs, M);
    % Fractions of tips actually introduced in each segment
    sampFrac = zeros(nPrs, M);
    
    % M sample schemes at each (p1, p2), each introduced at a T/2
    for i = 1:nPrs
        % Sample probs in this case
        p0 = p1(i); q0 = p2(i);
        % Sample number
        for j = 1:M
            % Alternates sampling prob between p1 and 1-p1 as diff intervals
            [sampCell{i, j}, tintro(i, j)] = sampDetSeq(n, p0, nMaxSet);
            % Get fraction of tips actually introduced in N1 vs N2
            sampFrac(i, j) = sum(sampCell{i, j}(1:2:end))/n;
        end
    end
    
    % Pad all sample schemes with extra 0s so same (max) length
    nElemMax = max(max(tintro));
    sampPad = cell(1, 1);
    for i = 1:nPrs
        for j = 1:M
            tempPad = sampCell{i, j};
            % Zero padding is on future T/2s
            sampPad{i, j} = padarray(tempPad, [0 nElemMax - length(tempPad)], 'post');
        end
    end
    % Convert to matrix and reshape
    sampPad = cell2mat(sampPad);
    sampPad = reshape(sampPad, [nPrs nElemMax M]);
    
    % On every p1 get a mean sample introduction list across nElemMax
    meanSampIntro = zeros(nPrs, nElemMax);
    for i = 1:nPrs
        tempPad = squeeze(sampPad(i, :, :));
        meanSampIntro(i, :) = mean(tempPad, 2);
    end
    
    %% Coalescent and sampling simulation
    
    % Coalescent event times and lineage counts
    tcoal = zeros(nPrs, M, n-1); % n-1 coalescences on each tree
    nLin = cell(nPrs, M); tLin = nLin; tcoalRef = tcoal;
    
    % No. coalescences within an odd or even half period
    m1 = zeros(nPrs, M); m2 = m1;
    
    % Construct a coalescent tree for each sampled scheme
    parfor i = 1:nPrs
        for j = 1:M
            % Extract sample numbers and times (in T/2)
            svec = 1:tintro(i, j);
            nvec = sampCell{i, j};
            svec = svec - 1; % (0 is first time)
            
            % Set nData as distinct coalescent events and sample times
            nSampTimes = tintro(i, j);
            nData = n + nSampTimes - 2;
            
            % Simulate heterochronous coalescent
            [~, ~, tcoalTemp] = sampHetSqHalfPeriod(T, N1(kk),...
                N2(kk), fac, svec, nvec, nData, nSampTimes, n);
            % Remove the first 0
            tcoal(i, j, :) = tcoalTemp(2:end);
            
            % Find all events within high and low populations
            tcoalRef(i, j, :) = mod(tcoalTemp(2:end), T)/T;
            m1(i, j) = sum(tcoalRef(i, j, :) <= 0.5);
            m2(i, j) = sum(tcoalRef(i, j, :) > 0.5);
        end
    end
    
    % Convert m1 to a raw fraction and one around 0.5
    pm1.abs = m1/(n-1);
    pm1.rel = abs(pm1.abs - 0.5);
    
    % Check that m1 and m2 consistent
    mdiff = m1 + m2 - (n-1);
    if any(mdiff)
        error('m_i do not sum to n-1');
    end
    % Check that maxima and minima are within bounds
    if max(max(pm1.abs)) > 1 || min(min(pm1.abs)) < 0
        error('pm1.abs out of bounds');
    end
    
    % Statistics on pm1 and Rpm1 - mean and 95% range
    pm1M.abs = mean(pm1.abs, 2)';
    pm1M.rel = mean(pm1.rel, 2)';
    pm1L.abs = prctile(pm1.abs', 2.5);
    pm1L.rel = prctile(pm1.rel', 2.5);
    pm1U.abs = prctile(pm1.abs', 97.5);
    pm1U.rel = prctile(pm1.rel', 97.5);
    
    % Correct for 1st half period as 0 id
    tintro = tintro - 1;
    tintroM = mean(tintro, 2);
    % Correlation betw samples and coalescents
    corrs.emp_samp_coal = corr2(sampFrac', pm1.abs');
    % Mean correlations
    sampFracM = mean(sampFrac, 2);
    corrs.samp_coal = corr2(sampFracM, pm1M.abs');
    corrs.prob_coal = corr2(p1', pm1M.abs');
    corrs.prob_samp = corr2(p1', sampFracM);
    disp(corrs);
    % Store sample struc means
    sampStruc.tintro = tintroM;
    sampStruc.sampFrac = sampFracM;
    
    % Find point on mean curve where achieve fraction of 0.5 for m1
    [m1best.minAbsDiff, m1best.idMin] = min(pm1M.rel);
    m1best.frac = pm1M.abs(m1best.idMin);
    m1best.p1 = p1(m1best.idMin);
    m1best.samp = sampFracM(m1best.idMin);
    disp(['m_1 = 0.5 is at p_1 = ' num2str(m1best.p1)]);
    disp(['m_1 = 0.5 is at samp frac = ' num2str(m1best.samp)]);
    
    % Save data in dedicated folder
    varList = {'N1', 'N2', 'amp', 'M', 'n', 'T', 'p1', 'pm1', 'tintroM',...
        'pm1M', 'pm1L', 'pm1U', 'clk', 'filePath', 'sampStruc', 'corrs',...
        'm1best', 'useMean', 'meanSampIntro', 'sampFracM', 'nPrs', 'lenN'};
    % Create folder to save
    homeName = cd;
    dirName = 'ntipsExample'; % <---- change name as needed
    mkdir(dirName);
    cd(dirName);
    save(['detSamples_' num2str(M) '_' num2str(n) '_' num2str(kk)], varList{:});
    cd(homeName);
    
    disp(['Progress: ' num2str(kk) ' of ' num2str(lenN)]);
end

% Sim time
trun = toc/60;
disp(['Sim time = ' num2str(trun) ' mins']);


%% Data visualisation

% Grab data file names from folder
cd(dirName);
files = dir('*.mat');
nFile = length(files);
D = cell(1, nFile);
% Load data from files
for i = 1:nFile
    D{i} = load(files(i).name, 'pm1M', 'pm1L', 'pm1U', 'sampFracM',...
        'tintroM', 'p1', 'm1best', 'meanSampIntro', 'n', 'N1', 'N2');
end
cd(homeName);

% Colours of curves
col = {'k', 'r', 'g', 'm'};

% Plot the d(m1) statistics against p1
figure;
hold on;
for i = 1:nFile
    plotErrBnd(gca, D{i}.p1', D{i}.pm1M.rel', D{i}.pm1L.rel', D{i}.pm1U.rel', col{i});
    % The optimal p1 value
    mID = D{i}.m1best.idMin;
    plot(D{i}.p1(mID), D{i}.pm1M.rel(mID), '.', 'color', [0.5 0.5 0.5], 'markersize', 30);
end
hold off; grid off; box off;
xlim([0 1]);
xlabel('p_1'); ylabel('d(m_1)');

% Absolute m1 fractions against p1
figure;
hold on;
for i = 1:nFile
    plotErrBnd(gca, D{i}.p1', D{i}.pm1M.abs', D{i}.pm1L.abs', D{i}.pm1U.abs', col{i});
    % The optimal p1 value
    mID = D{i}.m1best.idMin;
    plot(D{i}.p1(mID), D{i}.pm1M.abs(mID), '.', 'color', [0.5 0.5 0.5], 'markersize', 30);
end
hold off; grid off; box off;
xlim([0 1]);
xlabel('p_1'); ylabel('p_{m1}');


% Plot the d(m1) statistics of a subset of the data against mean q1
figure;
hold on;
for i = 1:nFile
    q1 = D{i}.sampFracM;
    plotErrBnd(gca, q1, D{i}.pm1M.rel', D{i}.pm1L.rel', D{i}.pm1U.rel', col{i});
    % The optimal p1 value
    mID = D{i}.m1best.idMin;
    plot(q1(mID), D{i}.pm1M.rel(mID), '.', 'color', [0.5 0.5 0.5], 'markersize', 30);
end
hold off; grid off; box off;
xlim([0 1]);
xlabel('q_1'); ylabel('d(m_1)');

% Plot sampling strategy for specific n case
idplt = 1;
% Strip meanSampIntro of small values and round
figRan = round(linspace(1, length(q1), 10));
meanSampRnd = round(D{idplt}.meanSampIntro(figRan, :), 2);

% Plot how the sample distributions look in terms of half periods
figure;
% Control figure size
smallFig = length(figRan);
for jj = 1:10
    % Assumes 10 separate schemes
    if jj <= 5
        subplot(5, 2, 2*jj-1);
    else
        subplot(5, 2, 2*(jj-5));
    end
    h = stem(meanSampRnd(jj, :), 'linewidth', 2);
    % Accessorise colours
    h.Color = [0.8 0.8 0.8];
    h.MarkerEdgeColor = 'm';
    h.Marker = '.';
    h.MarkerSize = 22;
    
    box off; grid off;
    %set(gca,'ytick',[]);
    if jj ~= 5 && jj ~= 10
        set(gca,'xtick',[]);
    else
        xlabel('no. T/2');
    end
    ylabel(['p_1 = ' num2str(round(D{idplt}.p1(figRan(jj)), 2))]);
end

% Cumsum across half periods with p1
figure;
hold all
%jj = 1:10:length(D{idplt}.p1)
for jj = round(linspace(1, length(D{idplt}.p1), 5))
    h = plot(1:size(D{idplt}.meanSampIntro, 2), cumsum(D{idplt}.meanSampIntro(jj, 1:1:end)), 'linewidth', 2);
    if D{idplt}.p1(jj) <= 0.5
        set(h, 'color', 'b', 'marker', 'o');
    else
        set(h, 'color', 'r', 'marker', 'x');
    end
end
hold off
xlabel('no. T/2');
ylabel('sum of sample introductions');



















