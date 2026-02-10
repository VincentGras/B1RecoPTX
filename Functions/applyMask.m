function data = applyMask(Mask, data, fillval)

if (nargin < 3)
    fillval = NaN;
end



if (~islogical(Mask))
    data = pTXUtils.applyMask(Mask>0, data, fillval);
else
    sz = size(data);
    nvol = numel(data)/numel(Mask);
    assert(rem(nvol,1)==0, 'Mask size is not compatible with data size');
    data = reshape(data, numel(Mask), numel(data)/numel(Mask));
    MaskC = reshape(~Mask, numel(Mask), 1);
    
    
%     n = max(ndims(data), ndims(Mask)+1);
%     s = repmat({':'}, 1, n-1);
    
    
    
    for i = 1:nvol
        
%         s_i = [s, i];        
%         data_ = data(s_i{:});
%         data_(~Mask) = fillval;
%         data(s_i{:}) = data_;        

        data(MaskC, i) = fillval; 

    end
    
    data = reshape(data, sz);
end


