a_abst = cell(size(final,1),2);
h_abst = cell(size(final,1),2);
bounds = cell(size(final,1),1);

for i = 1:size(final,1)
    
    Au = final{i,1};
    hu = real(final{i,2});
    
    Ab = final{i,3};
    hb = real(final{i,4});
    
    a_abst{i,1} = Au;
    h_abst{i,1} = hu;
    
    a_abst{i,2} = Ab;
    h_abst{i,2} = hb;
    
    subBound = final{i,6};
    xp = subBound(1,1);
    xq = subBound(1,2);
    
    bounds{i,1} = [xp,xq];
    
end

save coeffs_1D.mat a_abst h_abst bounds