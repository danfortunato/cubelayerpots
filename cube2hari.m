function [HfC, CfH] = cube2hari(ncheb, nside)

% - Reorder indices to match Hari
% - Shift from [-0.5, 0.5]^3 to [0, 1]^3
% - Interpolate from Chebyshev type-1 grid to Hari's grid (Chebyshev type-2 grid
%   with the endpoints dropped)

idx = reorder(ncheb, nside);
[x, ~, w] = chebpts(ncheb, 1);
xhari = chebpts(ncheb+2, 2);
xhari([1 end]) = [];
whari = 1./prod(xhari-xhari.'+eye(ncheb)).';
whari = whari / max(abs(whari));
P = barymat(xhari, x, w);
Q = barymat(x, xhari, whari);
PP = kron(P, P);
QQ = kron(Q, Q);
Pblock = repmat({PP}, 6*nside^2, 1);
Qblock = repmat({QQ}, 6*nside^2, 1);
PPP = matlab.internal.math.blkdiag(Pblock{:});
QQQ = matlab.internal.math.blkdiag(Qblock{:});
[pi, pj, pv] = find(PPP);
[qi, qj, qv] = find(QQQ);
n = 6*ncheb^2*nside^2;
HfC = sparse(pi, idx(pj), pv, n, n);
CfH = sparse(idx(qi), qj, qv, n, n);

end
