function [scale] = search_effscale(invec)

scale = zeros(2,1);

[~,address] = find(~isnan(invec(2,:)));

if size(address,2) >=1
    scale(1) = min(invec(1,address));
    scale(2) = max(invec(1,address));
else
    scale = NaN(2,1);
end

end

