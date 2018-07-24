% Deterministic sequential sampling
function [nIntro, tRes] = sampDetSeq(nS, p, sampNo)

% Assumptions and notes
% - fixed no. samples introduced at any time
% - sequentially introduce based on p and 1-p
% - alternates between p and 1-p as other samples for N2 vs N1
% - p is fraction of samples for N1 and nS is total samples

% Get fractions of N1 and N2 samples and remaining
f1 = round(p*nS); f2 = nS-f1;
fRem = [f1 f2];
% No. half periods for each case with full sample introductions
h1 = floor(f1/sampNo); h2 = floor(f2/sampNo);
hset = [h1 h2];

% Cells to store times and resource usuage (successes in binomials)
nIntro = cell(1, 1);  tRes = 0; 

while any(fRem > 0)
    % Use resources, fRem is what left
    fRemOld = fRem;
    % Assign samples according to if even or odd T/2
    if ~rem(tRes, 2)
        % Assign N1 samples, account for h1+1 introductions of < f1    
        fRem(1) = fRemOld(1) - min(sampNo, fRemOld(1));
        % Increment time and get introduction
        tRes = tRes + 1;
        nIntro{tRes} = fRemOld(1) - fRem(1);
    else
        % Assign N2 samples, account for h2+1 introductions of < f2
        fRem(2) = fRemOld(2) - min(sampNo, fRemOld(2));
        % Increment time and get introduction
        tRes = tRes + 1;
        nIntro{tRes} = fRemOld(2) - fRem(2);
    end
end

% Get sample vector
nIntro = cell2mat(nIntro);
% Check used all resources
if sum(sum(nIntro) ~= nS)
    assignin('base', 'nIntro', nIntro); 
    error('Not all sampling resources were used');
end