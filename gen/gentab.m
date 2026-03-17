function [slpsame, slpside, dlpside] = gentab(rlam, ncheb, opts)
% GENTAB  Generate Laplace/Yukawa layer potential quadrature tables on a cube.
%
%   [SLPSAME, SLPSIDE, DLPSIDE] = GENTAB(RLAM, NCHEB) generates
%   near-interaction quadrature tables for a cube mesh with Chebyshev
%   order NCHEB and Yukawa parameter RLAM (use RLAM=0 for Laplace).
%
%   [SLPSAME, SLPSIDE, DLPSIDE] = GENTAB(RLAM, NCHEB, NTHREADS) also
%   sets the number of OpenMP threads (default: use OMP_NUM_THREADS).
%
%   Outputs:
%     SLPSAME  (nf,nf,9)     - SLP weights, same-face neighbors
%     SLPSIDE  (nf,nf,6,6,3) - SLP weights, adjacent-face neighbors
%     DLPSIDE  (nf,nf,6,6,3) - DLP weights, adjacent-face neighbors
%   where nf = NCHEB^2. The last index of SLPSIDE/DLPSIDE corresponds
%   to istep = -1, 0, +1 (indices 1, 2, 3).

arguments
    rlam          (1,1) {mustBeReal}
    ncheb         (1,1) {mustBePositive, mustBeInteger}
    opts.nthreads (1,1) {mustBeNonnegative, mustBeInteger} = maxNumCompThreads
end

if ~isscalar(rlam) || ~isreal(rlam)
    error('rlam must be a real scalar.');
end

if ~isscalar(ncheb) || ncheb < 1 || ncheb ~= floor(ncheb)
    error('ncheb must be a positive integer scalar.');
end

nf = ncheb^2;
sz_same = [nf nf 9];
sz_side = [nf nf 6 6 3];
ntab_same = prod(sz_same);
ntab_side = prod(sz_side);
tab_slpsame = zeros(ntab_same, 1);
tab_slpside = zeros(ntab_side, 1);
tab_dlpside = zeros(ntab_side, 1);
nthreads = opts.nthreads;

mex_id_ = 'gentab_omp(i double, i int, i int, io double[x], io double[x], io double[x])';
[tab_slpsame, tab_slpside, tab_dlpside] = gentab_mex(mex_id_, rlam, ncheb, nthreads, tab_slpsame, tab_slpside, tab_dlpside, ntab_same, ntab_side, ntab_side);

slpsame = reshape(tab_slpsame, sz_same);
slpside = reshape(tab_slpside, sz_side);
dlpside = reshape(tab_dlpside, sz_side);

end
