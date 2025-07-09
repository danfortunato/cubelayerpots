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

ndisc = ncheb^2;
nsq_per_face = nside^2;
nsys = 6 * nsq_per_face * ndisc;

patchsize = 1/nside;
A1 = 0;
B1 = patchsize;
U1 = (B1-A1)/2;
V1 = (B1+A1)/2;
xcheb = chebpts(ncheb, 1);
x1 = U1 * xcheb;
[xs, xf] = meshgrid(x1);
xs = xs(:);
xf = xf(:);

% Set point locations and density values
cent = -0.5 + ((1:nside)-0.5)*patchsize;
centers = zeros(3, nsq_per_face, 6);
sources = zeros(3, ncheb^2, nsq_per_face, 6);

n = 0;
for iface = 1:6
    for nslow = 1:nside
        for nfast = 1:nside
            jsq = (nslow-1)*nside + nfast;
            if (iface == 1)
                xcent = -0.5;
                ycent = cent(nslow);
                zcent = cent(nfast);
                sources(1,:,jsq,iface) = -0.5;
                sources(2,:,jsq,iface) = xs + ycent;
                sources(3,:,jsq,iface) = xf + zcent;
            elseif (iface == 2)
                xcent = 0.5;
                ycent = cent(nfast);
                zcent = cent(nslow);
                sources(1,:,jsq,iface) = 0.5;
                sources(2,:,jsq,iface) = xf + ycent;
                sources(3,:,jsq,iface) = xs + zcent;
            elseif (iface == 3)
                xcent = cent(nfast);
                ycent = -0.5;
                zcent = cent(nslow);
                sources(1,:,jsq,iface) = xf + xcent;
                sources(2,:,jsq,iface) = -0.5;
                sources(3,:,jsq,iface) = xs + zcent;
            elseif (iface == 4)
                xcent = cent(nslow);
                ycent = 0.5;
                zcent = cent(nfast);
                sources(1,:,jsq,iface) = xs + xcent;
                sources(2,:,jsq,iface) = 0.5;
                sources(3,:,jsq,iface) = xf + zcent;
            elseif (iface == 5)
                xcent = cent(nslow);
                ycent = cent(nfast);
                zcent = -0.5;
                sources(1,:,jsq,iface) = xs + xcent;
                sources(2,:,jsq,iface) = xf + ycent;
                sources(3,:,jsq,iface) = -0.5;
            elseif (iface == 6)
                xcent = cent(nfast);
                ycent = cent(nslow);
                zcent = 0.5;
                sources(1,:,jsq,iface) = xf + xcent;
                sources(2,:,jsq,iface) = xs + ycent;
                sources(3,:,jsq,iface) = 0.5;
            end
            centers(1,jsq,iface) = xcent;
            centers(2,jsq,iface) = ycent;
            centers(3,jsq,iface) = zcent;
        end
    end
end

xpp = -1.2;
ypp =  2.2;
zpp =  3.2;

xx = squeeze(sources(1,:,:,:));
yy = squeeze(sources(2,:,:,:));
zz = squeeze(sources(3,:,:,:));
rx = xx - xpp;
ry = yy - ypp;
rz = zz - zpp;
rr = rx.^2 + ry.^2 + rz.^2;
uu = 1./sqrt(rr);
ux = -rx.*uu./rr;
uy = -ry.*uu./rr;
uz = -rz.*uu./rr;

slpdens = zeros(ncheb^2, nsq_per_face, 6);
slpdens(:,:,1) = -ux(:,:,1);
slpdens(:,:,2) =  ux(:,:,2);
slpdens(:,:,3) = -uy(:,:,3);
slpdens(:,:,4) =  uy(:,:,4);
slpdens(:,:,5) = -uz(:,:,5);
slpdens(:,:,6) =  uz(:,:,6);
dlpdens = uu;

[slpmat, dlpmat] = cubelayermats(ncheb, nside, nq);

% Now de-interleave target and source data and reshape into matrices
slpmat = permute(slpmat, [1 3 5 2 4 6]);
slpmat = reshape(slpmat, nsys, nsys);
dlpmat = permute(dlpmat, [1 3 5 2 4 6]);
dlpmat = reshape(dlpmat, nsys, nsys);

slpval = reshape(slpmat*slpdens(:), [ndisc nsq_per_face 6]);
dlpval = reshape(dlpmat*dlpdens(:), [ndisc nsq_per_face 6]);

ucomp = slpval - dlpval;
err = sum((0.5*uu(:) - ucomp(:)).^2);

fprintf('L2 error = %g\n', sqrt(err));
 
toc
