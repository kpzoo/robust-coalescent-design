function cellPart = intpartitions(varargin)
%INTPARTITION performs integer partition, i.e. the partition of of a set 
%containing homogenous elements. The function generates a cell array 
%containing a list of vectors representing all possible ways 
%of partitioning a set containing intIn number of identical elements 
%without order awareness. The numerical representation in the 
%output describes the partitions as: {[3 1]} = [1 1 1 | 1].
%   
%   Syntax:
%   intpartition(n)
%   intpartition(n,s)
%
%   The resulting partions can also be seen as writing the input n as a sum
%   of positive integers. An optional argument s can be supplied to output
%   a subset of partitions with number of parts less than or equal to s.
%
%   Number of ways of partitioning is according to sequence:
%   http://oeis.org/A000041 - (1), 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, ...
%
%   Example 1: intpartition(4) gives {[1 1 1 1],[1 1 2],[1 3],[2 2],4}
%   Example 2: intpartition(10,2) gives {[3,7],[4,6],[5,5],10}
%
    n=varargin{1};
    if nargin >= 2
        s=varargin{2};
    else
        s=inf;
    end
    
    n = round(n);
    if ~isscalar(n)
        error('Invalid input. Input must be a scalar integer')
    end
    
    cellPart={};
	a = zeros(n,1);
    k = 1;
	y = n - 1;
    while k ~= 0
        x = a(k) + 1;
        k = k-1;
        while 2*x <= y
            a(k+1) = x;
            y = y - x;
            k = k + 1;
        end
        l = k + 1;
        while x <= y
            a(k+1) = x;
            a(l+1) = y;
            if k+2<=s
                cellPart(end+1) = {a(1:k + 2)};
            end
            x = x + 1;
            y = y - 1;
        end
        a(k+1) = x + y;
        y = x + y - 1;
        if k+1<=s
            cellPart(end+1) = {a(1:k + 1)};
        end
   end
end