% Define Heaviside function such that:
% H(x) = 0 if x <= 0
% H(x) = 1 if x > 0

function out = heaviside(input)
assert(isnumeric(input))
if input <= 0
    out = 0;
else
    out = 1;
end
end
