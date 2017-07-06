function [ w ] = numunwrap( v, zmax )
%numunwrap unwraps a vector v, assuming its period is zmax
% assumes v is supposed to be monotonic increasingw
w = v;
d = diff(w);
badinds = find(d<0)+1;

for ii = 1:length(badinds)
    % Ensures actually a wraparound problem
    % (sometimes wrong soliton picked up)
        prevind = badinds(ii) - 1;
        while ismember(prevind,badinds)
            prevind = prevind - 1;
        end
        if abs(w(prevind)-w(badinds(ii)))>0.87*zmax %0.87 because lazy and it works
            w(badinds(ii):end) = w(badinds(ii):end) + zmax;
        end
end

end

