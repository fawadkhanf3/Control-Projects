function sub = divide_bound(bound)
    sub = cell(4,1);

    xp = bound(1,1);
    xq = bound(1,2);
    yp = bound(2,1);
    yq = bound(2,2);
    
    sub{1} = [xp (xp+xq)/2;yp (yp+yq)/2];
    sub{2} = [(xp+xq)/2 xq;yp (yp+yq)/2];
    sub{3} = [xp (xp+xq)/2;(yp+yq)/2 yq];
    sub{4} = [(xp+xq)/2 xq;(yp+yq)/2 yq];
end