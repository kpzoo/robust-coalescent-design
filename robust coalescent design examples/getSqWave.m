% Get true square wave demographics
function Nt = getSqWave(t, T, duty, N1, N2)

% Assumption and notes
% - T is period and duty is high time fraction
% - N1 is starting level and N2 the next switched level
% - works for vector t

% Get location modulo period
tloc = mod(t, T)/T;

% If in duty part of cycle then it is N1 else N2
Nt = N1*ones(size(t));
Nt(tloc > duty) = N2;
    