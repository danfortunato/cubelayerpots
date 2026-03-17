function [xx, yy, zz] = cubepts(ncheb, nside)

%             ____
%            /| 6 /|      6 = top                  y
%           /_|_ / |      5 = bottom           z| /
%         1|  |_ |_|2     4 = (y= 1/2 face)     |/___ x
%          | /   | /      3 = (y=-1/2 face)
%          |/__5_|/       2 = (x= 1/2 face)
%                         1 = (x=-1/2 face)
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

nsq_per_face = nside^2;
patchsize = 1/nside;
dom = [-patchsize/2 patchsize/2];
x1 = chebpts(ncheb, dom, 1);
[xs, xf] = meshgrid(x1);
xs = xs(:);
xf = xf(:);

% Set point locations and density values
cent = -0.5 + ((1:nside)-0.5)*patchsize;
centers = zeros(3, nsq_per_face, 6);
sources = zeros(3, ncheb^2, nsq_per_face, 6);

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

xx = squeeze(sources(1,:,:,:));
yy = squeeze(sources(2,:,:,:));
zz = squeeze(sources(3,:,:,:));

end
