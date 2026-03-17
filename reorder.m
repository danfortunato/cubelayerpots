function [idx, chebidx, sqidx, faceidx] = reorder(ncheb, nside)

if ( nargin == 0 )
    test();
    return
end

cheb = (1:ncheb^2).';
chebt = reshape(cheb, ncheb, ncheb).';
chebt = chebt(:);

x = uint64(0:nside-1);
[xx, yy] = meshgrid(x);
m = cartesian2morton(xx, yy);
mt = m.';
[~, sq ] = sort(m(:));
[~, sqt] = sort(mt(:));

chebidx = [ chebt cheb  cheb  chebt chebt cheb ];
sqidx   = [ sq    sqt   sqt   sq    sq    sqt  ];
faceidx = [ 1     3     2     5     4     6    ]; % Hari's face ordering

chebidx = reshape(chebidx, [ncheb^2 1 6]);
sqidx   = reshape(sqidx,   [1 nside^2 6]);
faceidx = reshape(faceidx, [1 1 6]);
idx = chebidx + (sqidx-1)*ncheb^2 + (faceidx-1)*nside^2*ncheb^2;
idx = idx(:);

chebidx = reshape(chebidx, [ncheb^2 6]);
sqidx   = reshape(sqidx,   [nside^2 6]);
faceidx = reshape(faceidx,   [1 6]);

end

function test()

nref = 3;
nside = 2^nref;
ncheb = 5;
[~, ~, sqidx, ~] = reorder(ncheb, nside);
idx = sqidx(:,1);

x = uint64(0:nside-1);
[xx, yy] = meshgrid(x);
xm = xx(idx);
ym = yy(idx);

cmap = colormap('redblue');
P = length(cmap);
cmap = interp1(1:P, cmap, linspace(1, P, nside^2), 'linear');

clf
axis equal
axis([0 nside 0 nside])
hold on
drawnow
shg
for k = 1:length(xm)
    p = polyshape([xm(k) xm(k) xm(k)+1 xm(k)+1], [ym(k) ym(k)+1 ym(k)+1 ym(k)]);
    plot(p, facecolor=cmap(k,:))
    text(double(xm(k))+0.5, double(ym(k))+0.5, 0, sprintf('%d', idx(k)))
    drawnow
    hold on
    pause(0.02)
end
shg

end

function [x, y] = morton2cartesian(morton)
    x = Compact1By1(morton);
    y = Compact1By1(bitshift(morton, -1));
end

function morton = cartesian2morton(x, y)
    morton = bitshift(Part1By1(y), 1) + Part1By1(x);
end

function x = Part1By1(x)
% "Insert" a 0 bit after each of the 16 low bits of x
    x = bitand(x, 0x00000000ffffffff);
    x = bitand(bitxor(x, bitshift(x, 16)), 0x0000ffff0000ffff);
    x = bitand(bitxor(x, bitshift(x, 8)),  0x00ff00ff00ff00ff);
    x = bitand(bitxor(x, bitshift(x, 4)),  0x0f0f0f0f0f0f0f0f);
    x = bitand(bitxor(x, bitshift(x, 2)),  0x3333333333333333);
    x = bitand(bitxor(x, bitshift(x, 1)),  0x5555555555555555);
end

function x = Compact1By1(x)
% Inverse of Part1By1 - "delete" all odd-indexed bits
    x = bitand(x, 0x5555555555555555);
    x = bitand(bitxor(x, bitshift(x, -1)),  0x3333333333333333);
    x = bitand(bitxor(x, bitshift(x, -2)),  0x0f0f0f0f0f0f0f0f);
    x = bitand(bitxor(x, bitshift(x, -4)),  0x00ff00ff00ff00ff);
    x = bitand(bitxor(x, bitshift(x, -8)),  0x0000ffff0000ffff);
    x = bitand(bitxor(x, bitshift(x, -16)), 0x00000000ffffffff);
end
