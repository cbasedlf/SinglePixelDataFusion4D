function [ y ] = timeDeInt( x, tempRes )
%specDeInt Temporal De-integration
%   Generate new temporal columns (channels) replicating available one but
%   weighted by the global variable wt
    global wt
    y = zeros([size(x,1), tempRes, size(x,3), size(x,4)]);
    for i=1:size(x,3)
        for j=1:size(x,4)
            y(:,:,i,j) = squeeze(x(:,:,i,j))*wt(:,:,i,j);
        end
    end
end