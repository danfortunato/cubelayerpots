ncheb = 5;
nref = 3;
nq = 16;
rlam = 0;
nside = 2^nref;
ndisc = ncheb^2;
nsq_per_face = nside^2;
nsys = 6*nsq_per_face*ndisc;
[HfC, CfH] = cube2hari(ncheb, nside);

d = 1/800;
x0 =  0.1;
y0 = -0.1;
z0 =  0;
f = @(x,y,z) (4*pi*d)^(-3/2)*exp(-((x-x0).^2+(y-y0).^2+(z-z0).^2)/(4*d));
sol = @(x,y,z) pfreegauss(d, x0, y0, z0, x, y, z);

t = tic;
[xx, yy, zz] = cubepts(ncheb, nside);
[S, D] = cubelayermats(rlam, ncheb, nside, nq);

% Convert layer potentials from our cube to Hari's
S = HfC * S * CfH;
D = HfC * D * CfH;
xx = reshape(HfC * xx(:), size(xx)) + 0.5;
yy = reshape(HfC * yy(:), size(yy)) + 0.5;
zz = reshape(HfC * zz(:), size(zz)) + 0.5;
fprintf('Total time for layer potentials: %gs.\n', toc(t));

% slpval = reshape(slpmat*slpdens(:), [ndisc nsq_per_face 6]);
% dlpval = reshape(dlpmat*dlpdens(:), [ndisc nsq_per_face 6]);

% Pass f to Hari's code...
% addpath ../hpstr

pdo = struct();
pdo.dxx = -1; pdo.dyy = -1; pdo.dzz = -1;
pdo.dxy =  0; pdo.dyz =  0; pdo.dxz =  0;
pdo.dx  =  0; pdo.dy  =  0; pdo.dz  =  0;
pdo.b   =  rlam;
dim = 3;
order = ncheb+1;
rr = refel(dim, order);
anchor = [-0.5 -0.5 -0.5];
sz = 1;
level = nref;

root = block(rr, anchor, sz);
root.split(level);
root.update_idx(1);
tstart = tic;
root.build(pdo, f);
DtN_int = full(root.D2N);
neu_part = DtN_int(:,end);
DtN_int = DtN_int(:,1:end-1);
t = toc(tstart);

fprintf('Problem size: %d\n', ((2^level)*order + 1)^dim);
fprintf('Total time to build: %gs\n', t);

I = eye(size(S));
A = I/2 - D + S*DtN_int;

rhs = S * (-neu_part);
dir_tot = A \ rhs;

% Evaluate the computed solution at a target
xt =  0.3;
yt =  0;
zt = -0.5;

% Evaluate on the interior using Hari's code:
tstart = tic;
u = -root.solve(dir_tot);
utarg = -root.solve_at(dir_tot, [xt yt zt]);
t = toc(tstart);
fprintf('Total time to solve: %g\n', t);

% Evaluate at exterior targets:
% dlpdens = dir_scat;
% slpdens = neu_scat;
% D_u  = evaldlp(rlam, ncheb, nside, nq, dlpdens, xt, yt, zt);
% S_du = evalslp(rlam, ncheb, nside, nq, slpdens, xt, yt, zt);
% u = D_u - S_du;

uexact = sol(xt,yt,zt);
err = norm(utarg - uexact, inf) / norm(uexact, inf);
fprintf('Maximum error at target = %e\n', norm(err(:), inf));

usol = root.eval(sol);
err = u - usol;
fprintf('Maximum error at nodes  = %e\n', norm(err(:), inf));

figure(1), clf
set(gcf, 'Position', [100 100 1200 400])

subplot(131)
root.slice(u, 'x', 0.2);
title('HPS solution')

subplot(132);
root.slice(usol, 'x', 0.2);
title('Exact solution')

subplot(133)
root.slice(err, 'x', 0.2);
colorbar
title('Error')
