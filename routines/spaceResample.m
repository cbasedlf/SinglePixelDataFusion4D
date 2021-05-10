function [ hypcube ] = spaceResample( input, newSize )
%Given a 4D (lambda,t,x,y) tensor [input], spaceResample provides a
%spatially downsampled/upsampled 4D tensor with a size of [newSize]. Also
%normalizes the 4D hyperspectral tensor so norm(hypcube(:)) = 1;

%preallocating
hypcube = zeros(size(input,1),size(input,2),newSize,newSize);
%downsample process
for i=1:size(hypcube,1)
    for j=1:size(hypcube,2)
        temp = imresize(squeeze(input(i,j,:,:)),[newSize newSize],'nearest'); %resize
        temp2 = norm(tens2vec(input(i,j,:,:))); %old norm
        temp3 = temp2*temp/norm(temp(:)); %make new norm = old norm
        hypcube(i,j,:,:) = temp3;
    end
end

end