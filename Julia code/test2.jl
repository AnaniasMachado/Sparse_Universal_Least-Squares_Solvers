using MATLAB

M = 2
m = 10
n = 5
r = 2
density = 0.1

rc = zeros(r);

for i = 1:r
    rc(i) = M * ((1/M) ^ (2/(i+1)))
end

rc = mxarray(rc)

# Calls Matlab function
mat"$A = sprand($m, $n, $density, $rc)"

@show A
