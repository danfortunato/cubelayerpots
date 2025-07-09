function ibflag = mkibflag(nside)
%
%   INPUT:
%   nside        number of square patches in each linear dimension
%
%   OUTPUT:
%   ibflag(it,is,ifacet,ifaces) = -2  if patches are not neighbors
%   ibflag(it,is,ifacet,ifaces) = (1,..,9) if patches are neighbors
%                                           on same face in
%                                           standard ordering
%   ibflag(it,is,ifacet,ifaces) = 0 if patches share a side
%                                           on different faces
%   ibflag(it,is,ifacet,ifaces) = -1 if "it" patch is shifted down
%                                    by one on different face w.r.t
%                                    "is" patch
%   ibflag(it,is,ifacet,ifaces) =  1 if "it" patch is shifted up
%                                    by one on different face w.r.t
%                                    "is" patch
%   Note that
%    ibflag(it,is,ifacet,ifaces) = -ibflag(is,it,ifaces,ifacet)
%   for this situation: neighboring patches on different faces.


%     ORDERING CONVENTION
%
%             ____
%            /|F6 /|      F6 = top                  y
%           /_|_ / |      F5 = bottom           z| /
%        F1|  |_ |_|F2    F4 = (y=1 face)        |/___ x
%          | /   | /      F3 = (y=-1 face)
%          |/_F5_|/       F2 = (x=1 face)
%                         F1 = (x=-1 face)
%
%
%     Near neighbors on same face:
%
%              7 8 9            slow
%              4 5 6            /|\       on any give face
%              1 2 3  -> fast    |
%
%
%
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
ibflag = zeros(nsq_per_face, nsq_per_face, 6, 6);

for ifacet = 1:6
    for ifaces = 1:6
        for jsq = 1:nsq_per_face
            for jtq = 1:nsq_per_face
                ibflag(jtq,jsq,ifacet,ifaces) = -2;
            end
        end
        if (ifacet == ifaces)
            for jsq = 1:nsq_per_face
                nfasts = mod(jsq, nside);
                if (nfasts == 0)
                    nfasts = nside;
                end
                nslows = (jsq-nfasts)/nside + 1;

                nslowt = nslows;
                if ((nfasts-1) >= 1)
                    nfastt = nfasts-1;
                    jtq = (nslowt-1)*nside + nfastt;
                    ibflag(jtq,jsq,ifacet,ifaces) = 4;
                end

                nfastt = nfasts;
                jtq = (nslowt-1)*nside + nfastt;
                ibflag(jtq,jsq,ifacet,ifaces) = 5;

                if ((nfasts+1) <= nside)
                    nfastt = nfasts+1;
                    jtq = (nslowt-1)*nside + nfastt;
                    ibflag(jtq,jsq,ifacet,ifaces) = 6;
                end

                if ((nslows-1) >= 1)
                    nslowt = nslows-1;
                    if ((nfasts-1) >= 1)
                        nfastt = nfasts-1;
                        jtq = (nslowt-1)*nside + nfastt;
                        ibflag(jtq,jsq,ifacet,ifaces) = 1;
                    end

                    nfastt = nfasts;
                    jtq = (nslowt-1)*nside + nfastt;
                    ibflag(jtq,jsq,ifacet,ifaces) = 2;

                    if ((nfasts+1) <= nside)
                        nfastt = nfasts+1;
                        jtq = (nslowt-1)*nside + nfastt;
                        ibflag(jtq,jsq,ifacet,ifaces) = 3;
                    end
                end
                if ((nslows+1) <= nside)
                    nslowt = nslows+1;
                    if ((nfasts-1) >= 1)
                        nfastt = nfasts-1;
                        jtq = (nslowt-1)*nside + nfastt;
                        ibflag(jtq,jsq,ifacet,ifaces) = 7;
                    end

                    nfastt = nfasts;
                    jtq = (nslowt-1)*nside + nfastt;
                    ibflag(jtq,jsq,ifacet,ifaces) = 8;

                    if ((nfasts+1) <= nside)
                        nfastt = nfasts+1;
                        jtq = (nslowt-1)*nside + nfastt;
                        ibflag(jtq,jsq,ifacet,ifaces) = 9;
                    end
                end
            end
        end
    end
end

% 13 edge:
ifaces = 1;
ifacet = 3;
is1 = 1;
if3 = 1;
for if1 = 1:nside
    jsq = (is1-1)*nside + if1;
    for is3 = 1:nside
        if (is3 == if1)
            jtq = (is3-1)*nside + if3;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (is3 == if1+1)
            jtq = (is3-1)*nside + if3;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (is3 == if1-1)
            jtq = (is3-1)*nside + if3;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

% 14 edge:
ifaces = 1;
ifacet = 4;
is1 = nside;
is4 = 1;
for if1 = 1:nside
    jsq = (is1-1)*nside + if1;
    for if4 = 1:nside
        if (if4 == if1)
            jtq = (is4-1)*nside + if4;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (if4 == if1+1)
            jtq = (is4-1)*nside + if4;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (if4 == if1-1)
            jtq = (is4-1)*nside + if4;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

% 15 edge:
ifaces = 1;
ifacet = 5;
if1 = 1;
is5 = 1;
for is1 = 1:nside
    jsq = (is1-1)*nside + if1;
    for if5 = 1:nside
        if (if5 == is1)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (if5 == is1+1)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (if5 == is1-1)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

% 16 edge:
ifaces = 1;
ifacet = 6;
if1 = nside;
if6 = 1;
for is1 = 1:nside
    jsq = (is1-1)*nside + if1;
    for is6 = 1:nside
        if (is6 == is1)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (is6 == is1+1)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (is6 == is1-1)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end


% 23 edge:
ifaces = 2;
ifacet = 3;
if2 = 1;
if3 = nside;
for is2 = 1:nside
    jsq = (is2-1)*nside + if2;
    for is3 = 1:nside
        if (is3 == is2)
            jtq = (is3-1)*nside + if3;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (is3 == is2+1)
            jtq = (is3-1)*nside + if3;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (is3 == is2-1)
            jtq = (is3-1)*nside + if3;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

% 24 edge:
ifaces = 2;
ifacet = 4;
if2 = nside;
is4 = nside;
for is2 = 1:nside
    jsq = (is2-1)*nside + if2;
    for if4 = 1:nside
        if (if4 == is2)
            jtq = (is4-1)*nside + if4;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (if4 == is2+1)
            jtq = (is4-1)*nside + if4;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (if4 == is2-1)
            jtq = (is4-1)*nside + if4;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

% 25 edge:
ifaces = 2;
ifacet = 5;
is2 = 1;
is5 = nside;
for if2 = 1:nside
    jsq = (is2-1)*nside + if2;
    for if5 = 1:nside
        if (if5 == if2)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (if5 == if2+1)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (if5 == if2-1)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

% 26 edge:
ifaces = 2;
ifacet = 6;
is2 = nside;
if6 = nside;
for if2 = 1:nside
    jsq = (is2-1)*nside + if2;
    for is6 = 1:nside
        if (is6 == if2)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (is6 == if2+1)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (is6 == if2-1)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

% 36 edge:
ifaces = 3;
ifacet = 6;
is3 = nside;
is6 = 1;
for if3 = 1:nside
    jsq = (is3-1)*nside + if3;
    for if6 = 1:nside
        if (if6 == if3)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (if6 == if3+1)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (if6 == if3-1)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

% 35 edge:
ifaces = 3;
ifacet = 5;
is3 = 1;
if5 = 1;
for if3 = 1:nside
    jsq = (is3-1)*nside + if3;
    for is5 = 1:nside
        if (is5 == if3)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (is5 == if3+1)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (is5 == if3-1)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

% 46 edge:
ifaces = 4;
ifacet = 6;
if4 = nside;
is6 = nside;
for is4 = 1:nside
    jsq = (is4-1)*nside + if4;
    for if6 = 1:nside
        if (if6 == is4)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (if6 == is4+1)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (if6 == is4-1)
            jtq = (is6-1)*nside + if6;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

% 45 edge:
ifaces = 4;
ifacet = 5;
if4 = 1;
if5 = nside;
for is4 = 1:nside
    jsq = (is4-1)*nside + if4;
    for is5 = 1:nside
        if (is5 == is4)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) = 0;
            ibflag(jsq,jtq,ifaces,ifacet) = 0;
        end
        if (is5 == is4+1)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) =  1;
            ibflag(jsq,jtq,ifaces,ifacet) = -1;
        end
        if (is5 == is4-1)
            jtq = (is5-1)*nside + if5;
            ibflag(jtq,jsq,ifacet,ifaces) = -1;
            ibflag(jsq,jtq,ifaces,ifacet) =  1;
        end
    end
end

end
