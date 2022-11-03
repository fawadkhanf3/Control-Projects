function sub = divide_bound(bound)

    sub = cell(2,1);

    xp = bound(1,1);
    xq = bound(1,2);
    
    sub{1} = [xp (xp+xq)/2];
    sub{2} = [(xp+xq)/2 xq];

end