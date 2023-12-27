gaussleg = struct;
% point number list
gaussleg.point_index = [1 2 3 4 5 6 7 8 9 10 ...
                   20 30 40 50 60 70 80 90 100];

% cell for saving the parameters
gaussleg.para = cell(1,size(gaussleg.point_index,2));

for i = 1: size(gaussleg.point_index,2)
    gaussleg.para{1,i}= struct;  
    [zerosave,weight,~] = gengausslegquadrule(gaussleg.point_index(i));  
    gaussleg.para{1,i}.zero   = zerosave;
    gaussleg.para{1,i}.weight = weight;
end

save('gaussleg_para.mat','gaussleg');


