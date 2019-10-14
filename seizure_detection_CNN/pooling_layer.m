function f1x = pooling_layer(f1, layer)
% Average pooling layer
% Limited - average pooling only
%         - only 1D pooling

A = size(f1); K = floor(A(1)/layer.Stride(1))-floor(layer.PoolSize(1)/layer.Stride(1))+1; 
f1x = zeros(K, A(2));
for ii = 1:A(2) 
    for jj = 1:K
        f1x(jj,ii) = mean(f1((jj-1)*layer.Stride(1)+1:(jj-1)*layer.Stride(2)+layer.PoolSize(1),ii));
    end
end
