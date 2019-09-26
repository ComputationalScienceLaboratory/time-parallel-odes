function J = lorenz96_jacobian(~, y)

N = numel(y);

diagmat = [[y(end - 1); zeros(N - 1, 1)], [-y(2:(end - 1)) ; 0 ; 0], ...
    [y([3:end , 1]) - y([end, 1:(end - 2)]) ; 0], -1*ones(N, 1), [0; y([end, 1:(end - 2)])], ...
    [zeros(N - 2, 1); -y([end, 1])], [zeros(N - 1, 1); y(2)-y(end - 1)]];
dnum = [-(N - 1), -2, -1, 0, 1, N - 2, N - 1];

J = spdiags(diagmat, dnum, N, N);

end
