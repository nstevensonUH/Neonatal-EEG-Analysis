function f1x = normalization_layer(f1, layer)
% Normalization layer
%    simple zscore type normalization

A = layer.NumChannels; 
f1x = zeros(size(f1));
for ii = 1:A
    f1x(:,ii) = (f1(:,ii)-layer.TrainedMean(ii))./sqrt(layer.TrainedVariance(ii)+layer.Epsilon);
    f1x(:,ii) = f1x(:,ii).*layer.Scale(ii)+layer.Offset(ii);
end