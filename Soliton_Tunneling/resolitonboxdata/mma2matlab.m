function [ f ] = mma2matlab( f )
% Input:
%   f: expression from mathematica (as a string)
% Output:
%   f: expression converted to Mathlab syntax (as a string)

f = strrep(f, '*', '.*');
f = strrep(f, '/', './');
f = strrep(f, '^', '.^');
f = strrep(f, '[', '(');
f = strrep(f, ']', ')');
f = lower(f);
disp(f);



end

