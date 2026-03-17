function u = pfreegauss(d, x0, y0, z0, x, y, z)
%PFREEGAUSS
%   Computes exact solution to Poisson equation in free space due to Gaussian
%   source distribution centered at (x0, y0, z0) with variance d, i.e.:
%
%      f(x,y,z) = (4*pi*d)^(-3/2) * exp( -((x-x0)^2+(y-y0)^2+(z-z0)^2) / (4d) )
%
%   Input:
%
%      d          = Gaussian variance
%      x0, y0, z0 = Gaussian center
%      x, y, z    = Targets
%
%   Output:
%
%      u = Potential at targets

r = sqrt((x-x0).^2 + (y-y0).^2 + (z-z0).^2);
arg = r/(2*sqrt(d));
ddd = (4*pi*d)^(3/2);

u = zeros(size(arg));
idx = (arg > 1e-3);
argi  = arg( idx);
argni = arg(~idx);
u( idx) = erf(argi)./(4*pi*r(idx)) - 2*d*exp(-argi.*argi)/ddd;
u(~idx) = 2*d*(2*argni.*argni/3 - 2*argni.*argni.*argni.*argni/5 + (argni.^6)/7)/ddd;
u = u + 2*d*exp(-arg.*arg)/ddd;

end
