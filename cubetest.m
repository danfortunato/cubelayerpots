%
%     test Green's identity on constant density on cube:
%     4 squares with 5x5 grids on each face.
%
%
%             ____
%            /| 6 /|      6 = top                  y
%           /_|_ / |      5 = bottom           z| /
%         1|  |_ |_|2     4 = (y=1 face)        |/___ x
%          | /   | /      3 = (y=-1 face)
%          |/__5_|/       2 = (x=1 face)
%                         1 = (x=-1 face)
%     for 6: x=fast,y=slow
%          centers (-hh,-hh,-1),(hh,-hh,-1),(-hh,hh,-1),(hh,hh,-1),
%           e1 = (1,0,0),e2 = (0,1,0),e3 = (0,0,1),
%     for 5: y=fast,x=slow
%          centers (-hh,-hh,1),(hh,-hh,1),(-hh,hh,1),(hh,hh,1),
%           e1 = (0,1,0),e2 = (1,0,0),e3 = (0,0,-1),
%     for 4: z=fast,x=slow
%          centers (-hh,1,-hh),(hh,1,-hh),(-hh,1,hh),(hh,1,hh),
%           e1 = (0,0,1),e2 = (1,0,0),e3 = (0,1,0),
%     for 3: x=fast,z=slow
%          centers (-hh,-1,-hh),(hh,-1,-hh),(-hh,-1,hh),(hh,1,hh),
%           e1 = (1,0,0),e2 = (0,0,1),e3 = (0,-1,0),
%     for 2: y=fast,z=slow
%          centers (1,-hh,-hh),(1,hh,-hh),(1,-hh,hh),(1,hh,hh),
%           e1 = (0,1,0),e2 = (0,0,1),e3 = (1,0,0),
%     for 1: z=fast,y=slow
%          centers (-1,-hh,-hh),(-1,hh,-hh),(-1,-hh,hh),(1,hh,hh),
%           e1 = (0,0,1),e2 = (0,0,1),e3 = (-1,0,0),

tic

ncheb = 5;
nref = 2;
nside = 2^nref;
nq = 16;
rlam = 0;
xpp = -1.2;
ypp =  2.2;
zpp =  3.2;
[xx, yy, zz] = cubepts(ncheb, nside);

rx = xx - xpp;
ry = yy - ypp;
rz = zz - zpp;
rr2 = rx.^2 + ry.^2 + rz.^2;
rr = sqrt(rr2);
uu = exp(-rlam*rr)./rr;
ux = -rx.*exp(-rlam*rr)./(rr.*rr2) - (rlam*rx.*exp(-rlam*rr))./rr2;
uy = -ry.*exp(-rlam*rr)./(rr.*rr2) - (rlam*ry.*exp(-rlam*rr))./rr2;
uz = -rz.*exp(-rlam*rr)./(rr.*rr2) - (rlam*rz.*exp(-rlam*rr))./rr2;

slpdens = zeros(ncheb^2, nside^2, 6);
slpdens(:,:,1) = -ux(:,:,1);
slpdens(:,:,2) =  ux(:,:,2);
slpdens(:,:,3) = -uy(:,:,3);
slpdens(:,:,4) =  uy(:,:,4);
slpdens(:,:,5) = -uz(:,:,5);
slpdens(:,:,6) =  uz(:,:,6);
dlpdens = uu;

[slpmat, dlpmat] = cubelayermats(rlam, ncheb, nside, nq);

% Now de-interleave target and source data and reshape into matrices
% slpmat = permute(slpmat, [1 3 5 2 4 6]);
% slpmat = reshape(slpmat, nsys, nsys);
% dlpmat = permute(dlpmat, [1 3 5 2 4 6]);
% dlpmat = reshape(dlpmat, nsys, nsys);

slpval = reshape(slpmat*slpdens(:), [ncheb^2 nside^2 6]);
dlpval = reshape(dlpmat*dlpdens(:), [ncheb^2 nside^2 6]);

ucomp = slpval - dlpval;
err = sum((0.5*uu(:) - ucomp(:)).^2);

fprintf('L2 error = %g\n', sqrt(err));
 
toc
