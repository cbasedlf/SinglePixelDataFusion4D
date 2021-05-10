function [ y ] = timeInt( x )
%timeInt Temporal Integration
%   Integrates the time dimension of the 4D tensor and stores the
%   weight of each temporal column (channel) in a global variable wt
    y = sum(x,2);
    global wt
    wt = zeros([1, size(x,2), size(x,3), size(x,4)]);
    for i=1:size(x,3)
        for j=1:size(x,4)
            temp = squeeze(x(:,:,i,j));
            wt(:,:,i,j) = sum(temp,1)/sum(temp(:));
        end
    end
end