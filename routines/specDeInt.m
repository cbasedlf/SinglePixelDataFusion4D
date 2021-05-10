function [ y ] = specDeInt( x, specRes )
%specDeInt Spectral De-integration
%   Generate new spectral rows (channels) replicating available one but
%   weighted by the global variable ws
    global ws
    y = zeros([specRes, size(x,2), size(x,3), size(x,4)]);
    for i=1:size(x,3)
        for j=1:size(x,4)
            y(:,:,i,j) = ws(:,:,i,j)*squeeze(x(:,:,i,j));
        end
    end
end

