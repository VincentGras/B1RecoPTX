function [out, chunk] = map2vect(msk, data)

% Reconstruct 3d map from data
% Usage:
%   out = mat2vect(msk, data)
%   data : N*M
%   mask : 3D logical, with count(data==true) = N 
% See also: vect2map


m = numel(msk);
n = numel(data)/m;

if (isempty(msk))
    
   assert(ndims(data) <= 4, 'expecting 4-D input !'); 
    
   m = size(data,1)*size(data,2)*size(data,3); 
   n = numel(data)/m;
   out = reshape(data, m, n); 
   
   return;
    
end

mskv = msk(:)>0;
data = reshape(data, numel(mskv), numel(data)/numel(mskv)); 
out = zeros(sum(mskv(:)), size(data, 2), class(data));


for i = 1:n
        
    out(:,i) = data(mskv,i);
    
end

chunk = ones(size(data, 2)+1,1);

for i = 1:size(msk, 4)
    
    chunk(i+1) = chunk(i)+nnz(msk(:,:,:,i));
        
end
