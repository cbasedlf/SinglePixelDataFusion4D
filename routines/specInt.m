function [ y ] = specInt( x )
%specInt Spectral Integration
%   Integrates the spectral dimension of the 4D tensor and stores the
%   weight of each spectral row (channel) in a global variable ws
    y = sum(x,1);
    global ws
    ws = zeros([size(x,1), 1, size(x,3), size(x,4)]);
    for i=1:size(x,3)
        for j=1:size(x,4)
            temp = squeeze(x(:,:,i,j));
            ws(:,:,i,j) = sum(temp,2)/sum(temp(:));
        end
    end
end

