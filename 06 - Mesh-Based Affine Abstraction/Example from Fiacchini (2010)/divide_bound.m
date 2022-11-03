function sub = divide_bound(bound)
    sub = cell(6,1);

    xp = bound(1,1);
    xq = bound(1,2);
    yp = bound(2,1);
    yq = bound(2,2);
    zp = bound(3,1);
    zq = bound(3,2);
    
    sub{1} = [xp (xp+xq)/2;yp (yp+yq)/2;zp (zp+yq)/2];
    sub{2} = [(xp+xq)/2 xq;yp (yp+yq)/2;zp (zp+yq)/2];
    sub{3} = [xp (xp+xq)/2;(yp+yq)/2 yq;zp (zp+yq)/2];
    sub{4} = [xp (xp+xq)/2;yp (yp+yq)/2;(zp+zq)/2 zq];
    sub{5} = [(xp+xq)/2 xq;(yp+yq)/2 yq;zp (zp+yq)/2];
    sub{6} = [(xp+xq)/2 xq;(yp+yq)/2 yq;(zp+zq)/2 zq];
    
end