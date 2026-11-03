function x = tdma_batch(a, b, c, d)
%TDMA_BATCH  Thomas algorithm for B independent tridiagonal systems.
%   x = tdma_batch(a, b, c, d)   all inputs B x N
%   a : sub-diagonal   (a(:,1) unused)
%   b : main diagonal
%   c : super-diagonal (c(:,N) unused)
%   d : right-hand side
%   Loops over the N cells and vectorises over the B systems (column slices
%   are contiguous in MATLAB, so this is fast).

[B, N] = size(b);
cprime = zeros(B, N);
dprime = zeros(B, N);
cprime(:, 1) = c(:, 1) ./ b(:, 1);
dprime(:, 1) = d(:, 1) ./ b(:, 1);
for i = 2:N
    m = b(:, i) - a(:, i) .* cprime(:, i-1);
    cprime(:, i) = c(:, i) ./ m;
    dprime(:, i) = (d(:, i) - a(:, i) .* dprime(:, i-1)) ./ m;
end
x = zeros(B, N);
x(:, N) = dprime(:, N);
for i = N-1:-1:1
    x(:, i) = dprime(:, i) - cprime(:, i) .* x(:, i+1);
end
end
