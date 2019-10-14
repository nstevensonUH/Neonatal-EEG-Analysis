function f1x = conv_layer(f1, layer)
% Convolutional layer
% f1 is the input matrix (initially 1xN, but MxN for many later stages - M is the number of filters)
% f1x is the output (MxN-K), K is the filter length-1
% Limited - 2D structures
%           stride length set to 1

A = size(layer.Weights); B = size(f1);
f1x = single(zeros(B(1)-(A(1)-1), length(layer.Weights(1,1,1,:))));
for ii = 1:A(4)
    val = zeros(1,B(1)-(A(1)-1));
if A(3)==1
    for jj = 1:B(1)-(A(1)-1)        
        val(jj) = sum(f1(jj:(A(1)-1)+jj,:).*layer.Weights(:,1,:,ii))-layer.Bias(ii);
    end 
else
    for jj = 1:B(1)-(A(1)-1)        
        val(jj) = sum(sum(f1(jj:(A(1)-1)+jj,:).*reshape(layer.Weights(:,1, :, ii), A(1), A(3))))-layer.Bias(ii);
    end 
end
    f1x(:,ii) = val';
end

