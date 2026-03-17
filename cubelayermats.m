function [slpmat, dlpmat] = cubelayermats(rlam, ncheb, nside, nq)

ibflag = mkibflag(nside);

boxsize = 1/nside;
rlamuse = rlam*boxsize/2;
tab_file = sprintf('tab/tab_rlam=%.16g_ncheb=%d.mat', rlamuse, ncheb);
if isfile(tab_file)
    load(tab_file, 'slpsame', 'slpside', 'dlpside');
else
    [slpsame, slpside, dlpside] = gentab(rlamuse, ncheb);
    save(tab_file, 'slpsame', 'slpside', 'dlpside');
end

nsq_per_face = nside^2;
slpmat = zeros(ncheb^2, ncheb^2, nsq_per_face, nsq_per_face, 6, 6);
dlpmat = zeros(ncheb^2, ncheb^2, nsq_per_face, nsq_per_face, 6, 6);

% Evaluate Chebyshev polynomials at the quadrature nodes
xcheb = chebpts(ncheb, 1);
[xq, wq] = legpts(nq);
tx = zeros(nq, ncheb);
tx(:,1) = 1;
tx(:,2) = xq;
for jj = 3:ncheb
    tx(:,jj) = 2*xq.*tx(:,jj-1) - tx(:,jj-2);
end

% Convert the quadrature weights from coefficient space to grid space
cfv = chebtech1.vals2coeffs(eye(ncheb));
txw = wq(:) .* tx * cfv;
txyw = reshape(txw, [nq 1 ncheb 1]) .* reshape(txw, [1 nq 1 ncheb]);
txyw = reshape(txyw, [nq^2 ncheb^2]);

for ifacet = 1:6
    for ifaces = 1:6
        [slpfar, dlpfar] = genlpfar(rlam, ifaces, ifacet, nside, xcheb, xq, txyw);
        for jsq = 1:nsq_per_face
            for jtq = 1:nsq_per_face
                ifl = ibflag(jtq, jsq, ifacet, ifaces);
                if (ifl == -2)
                    slpmat(:,:,jtq,jsq,ifacet,ifaces) = slpfar(:,:,jtq,jsq);
                    dlpmat(:,:,jtq,jsq,ifacet,ifaces) = dlpfar(:,:,jtq,jsq);
                elseif (ifacet ~= ifaces)
                    slpmat(:,:,jtq,jsq,ifacet,ifaces) = slpside(:,:,ifacet,ifaces,ifl+2);
                    dlpmat(:,:,jtq,jsq,ifacet,ifaces) = dlpside(:,:,ifacet,ifaces,ifl+2);
                else
                    slpmat(:,:,jtq,jsq,ifacet,ifaces) = slpsame(:,:,ifl);
                end
            end
        end
    end
end

boxsize = 1/nside;
slpmat = slpmat * boxsize/2;

% Now de-interleave target and source data and reshape into matrices
ndof = ncheb^2 * nsq_per_face * 6;
slpmat = permute(slpmat, [1 3 5 2 4 6]);
slpmat = reshape(slpmat, ndof, ndof);
dlpmat = permute(dlpmat, [1 3 5 2 4 6]);
dlpmat = reshape(dlpmat, ndof, ndof);

% reshape(slpmat, [ncheb^2 nsq_per_face 6 ncheb^2 nsq_per_face 6]);

end

function [tabslp, tabdlp] = genlpfar(rlam, ifaces, ifacet, nside, xcheb, xq, txyw)
%
%     generates SLP, DLP far field NxN interaction matrix for all sources
%     on the IFACES side of the unit cube and all targets on the IFACET
%     side of the unit cube. We assume the discretization has NSIDE*NSIDE
%     squares on each face with tensor-product Chebyshev nodes on each square.
%
%     INPUT:
%     ifaces        cube face of source
%     ifacet        cube face of target
%     nside         nside*nside is number of squares on each face
%     xcheb         Chebyshev nodes (type-1, i.e. no endpoints) on [-1, 1]
%     xleg          Legendre nodes on [-1, 1]
%     
%
%     OUTPUT:
%     tabslp        map from densities at source to potentials at
%                   targets of dimension (nnodes*nnodes,nnodes*nnodes)
%     tabdlp        map from densities at source to potentials at
%                   targets of dimension (nnodes*nnodes,nnodes*nnodes)

%   create grid pts in 1d
%
%     |---|---|---|---|
%     boxsize = boxs/nside
%     ns = ncheb*(nslows-1)+ nfast
%
%   finest(TARG,j,k) =
%     int_Box T_j(x1)T_k(x2) R(TARG,x1,x2,zz) dx1 dx2
%    = sum_{m1,m2} Tj( xx(m1) ) T_k( yy(m2))
%           R(TARG,xx(m1),yy(m2),zz) w(m1) w(m2)
%    = sum_{m1} Tj( xx(m1)) w(m1) S1(m1,TARG,k)
%    where  S1(m1,TARG,k) =
%           sum_{m2}  T_k( yy(m2)) R(TARG,xx(m1),xx(m2),zz) w(m2)

ncheb = length(xcheb);
nq    = length(xq);

nsq_per_face = nside^2;
ns = 1:nsq_per_face;
nt = 1:nsq_per_face;

nfasts = mod(ns, nside);
nfasts(nfasts == 0) = nside;
nslows = (ns-nfasts)/nside + 1;

nfastt = mod(nt, nside);
nfastt(nfastt == 0) = nside;
nslowt = (nt-nfastt)/nside + 1;

a1 = -nside + (nfasts-1)*2;
b1 = a1 + 2;
u1 = (b1-a1)/2;
v1 = (b1+a1)/2;
xfs = u1.*xq + v1;

a2 = -nside + (nslows-1)*2;
b2 = a2 + 2;
u2 = (b2-a2)/2;
v2 = (b2+a2)/2;
xss = u2.*xq + v2;

at1 = -nside + (nfastt-1)*2;
bt1 = at1 + 2;
ut1 = (bt1-at1)/2;
vt1 = (bt1+at1)/2;
xft = ut1.*xcheb + vt1;

at2 = -nside + (nslowt-1)*2;
bt2 = at2 + 2;
ut2 = (bt2-at2)/2;
vt2 = (bt2+at2)/2;
xst = ut2.*xcheb + vt2;

xfs = reshape(xfs, [nq 1 1 1 1 nsq_per_face]);
xss = reshape(xss, [1 nq 1 1 1 nsq_per_face]);
xft = reshape(xft, [1 1 ncheb 1 nsq_per_face 1]);
xst = reshape(xst, [1 1 1 ncheb nsq_per_face 1]);

switch ifacet
    case 1
        xtarg = -nside;
        ytarg = xst;
        ztarg = xft;
    case 2
        xtarg = nside;
        ytarg = xft;
        ztarg = xst;
    case 3
        ytarg = -nside;
        ztarg = xst;
        xtarg = xft;
    case 4
        ytarg = nside;
        ztarg = xft;
        xtarg = xst;
    case 5
        ztarg = -nside;
        ytarg = xft;
        xtarg = xst;
    case 6
        ztarg = nside;
        ytarg = xst;
        xtarg = xft;
end

switch ifaces
    case 1
        rx = xtarg + nside;
        ry = ytarg - xss;
        rz = ztarg - xfs;
        dr = -rx;
    case 2
        rx = xtarg - nside;
        ry = ytarg - xfs;
        rz = ztarg - xss;
        dr = rx;
    case 3
        rx = xtarg - xfs;
        ry = ytarg + nside;
        rz = ztarg - xss;
        dr = -ry;
    case 4
        rx = xtarg - xss;
        ry = ytarg - nside;
        rz = ztarg - xfs;
        dr = ry;
    case 5
        rx = xtarg - xss;
        ry = ytarg - xfs;
        rz = ztarg + nside;
        dr = -rz;
    case 6
        rx = xtarg - xfs;
        ry = ytarg - xss;
        rz = ztarg - nside;
        dr = rz;
end

rr2 = rx.^2 + ry.^2 + rz.^2;
rr = sqrt(rr2);

boxsize = 1/nside;
rlamuse = rlam*boxsize/2;

kern = exp(-rlamuse*rr)./(4*pi*rr);
kern = reshape(kern, [nq^2 ncheb^2 nsq_per_face nsq_per_face]);
f = tensorprod(txyw, kern, 1, 1);
tabslp = permute(f, [2 1 3 4]);

kern = dr.*(1./rr+rlamuse).*exp(-rlamuse*rr)./(4*pi*rr2);
kern = reshape(kern, [nq^2 ncheb^2 nsq_per_face nsq_per_face]);
f = tensorprod(txyw, kern, 1, 1);
tabdlp = permute(f, [2 1 3 4]);

end
