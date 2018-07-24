% Heterochronously sampled coalescent for square waves
function [nLin, tLin, tcoal] = sampHetSqHalfPeriod(T, N1, N2, fac, svec, nvec,...
    nData, nSampTimes, n)

% Assumptions and notes
% - condition on nLinCurr changed to account for no samples at 0
% - duty is 0.5 and svec is in units of T/2
% - uses rejection sampling to better account for switch
% - the square wave switch times are known (T & duty)
% - only 2 possible levels N1 and N2, given in log of population

% Maximim non-lineage rate (log populations)
lmax = max(-N1, -N2);
emax = exp(lmax);
eN1 = exp(-N1); eN2 = exp(-N2);
% Convert svec into real time
svec = 0.5*svec*T;

% Initialise coalescent times and LTT plot variables
tcoal = zeros(1, n); % defined so tcoal(1) = 0 = svec(1)
nLin = zeros(1, nData + 1);
tLin = zeros(1, nData + 1); % all lineage values
nLin(1) = nvec(1); tLin(1) = svec(1);

% Count variables for simulation and booleans
iev = 1; % total no. loops
nLinCurr = nLin(1); % current no. lineages
tLinCurr = tLin(1); % start time
isamp = 1; % count no. samples <= nSampTimes
icoal = 0; % count no. coalescents <= n-1
sampVec = [svec inf]; % like svec but with inf for when exhaust samples in condition

% Main loop to simulate heterochronous coalescent
while(iev < nData+1)

    % Lineages fall to 1 before end then use next sample time and update
    if nLinCurr <= 1 && isamp < nSampTimes + 1
        % Update to next sample, no coalescents can occur
        sampTrue = 1;
    else
        % Reject boolean, change when tnext accepted
        reject = 1; tnext = 0;
        while(reject)
            % Get next time according to bounding rate
            r = rand(1, 2); lammax = fac(nLinCurr)*emax;
            tnext = tLinCurr -log(r(1))/lammax;
            % True rate at this tnext (modulo period to determine location)
            tloc = mod(tnext, T)/T;
            if tloc <= 0.5
                % First linear segment
                lamtrue = fac(nLinCurr)*eN1;
            else
                % Second linear segment
                lamtrue = fac(nLinCurr)*eN2;
            end 
            % Test if reject or accept next time
            if r(2) <= lamtrue/lammax
                % Event can be accepted
                reject = 0;
            end
        end
        
        % If accepted coalescent time is after next sample time then sample
        if tnext > sampVec(isamp+1)
            % Update to next sample, happens before next coalescent time
            sampTrue = 1;
        else
            % A coalescent event has occurred
            nLinCurr = nLinCurr - 1;
            tLinCurr = tnext;
            icoal = icoal + 1;
            sampTrue = 0;
            tcoal(icoal+1) = tLinCurr;
            %disp(['Completed coalescent ' num2str(icoal+1) ' of ' num2str(n-1)]);
        end        
    end
    
    % Sample update involves lineages and times from svec and nvec
    if sampTrue
        isamp = isamp + 1;
        nLinCurr = nLinCurr + nvec(isamp);
        tLinCurr = svec(isamp); % svec should never be exceeded
    end
    
    % Store data
    iev = iev + 1;
    nLin(iev) = nLinCurr;
    tLin(iev) = tLinCurr;
end

% Check simulation made use of all data
if(iev ~= isamp + icoal)
    error('Not all time data used in simulation');
end
% Check for consistency in event times
tcheck = sort([tcoal svec(2:end)]);
if ~all(tcheck == tLin)
    disp('The times are inconsistent');
end
if any(isnan(tLin))
    error('NaN values in the simulated times');
end