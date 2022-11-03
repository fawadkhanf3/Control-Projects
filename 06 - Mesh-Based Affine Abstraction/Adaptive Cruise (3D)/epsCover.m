function hyperplanes = epsCover(bound,param)

des_error = param.des_error;

[error,planes] = abstract(bound,param); % Linear program for Affine Abstraction in Theorem 1

if(error>des_error)
    
    sub_bounds = divide_bound(bound);
    
    cell_1 = epsCover(sub_bounds{1},param);
    cell_2 = epsCover(sub_bounds{2},param);
    cell_3 = epsCover(sub_bounds{3},param);
    cell_4 = epsCover(sub_bounds{4},param);
    
    s_1 = size(cell_1,1);
    s_2 = size(cell_2,1);
    s_3 = size(cell_3,1);
    s_4 = size(cell_4,1);
    
    size_f = s_1+s_2+s_3+s_4;
    
    hyperplanes = cell(size_f,1);
    k = 1;
    
    for i = 1:s_1
        hyperplanes{k,1} = cell_1{i,1};
        k = k+1;
    end
    for i = 1:s_2
        hyperplanes{k,1} = cell_2{i,1};
        k = k+1;
    end
    for i = 1:s_3
        hyperplanes{k,1} = cell_3{i,1};
        k = k+1;
    end
    for i = 1:s_4
        hyperplanes{k,1} = cell_4{i,1};
        k = k+1;
    end
    
else
    hyperplanes = cell(1);
    hyperplanes{1,1} = planes;
end

end