% DEMO
% Seizure detection in the neonatal EEG using deep convolutional neural
% convnetworks
%
% General idea is that the original time series is sequentiall process by a
% number of processing layers (defined by the convnetwork in trained_cnn_all_sgdm_v1.mat)
%
% Currently very messy implementation; an attempt to de-Matlab the
% implementation
%
% NOTES
% -  the sub-functions have not been generalized the work for the
%    seizure detection implementation in trained_cnn_all_sgdm_v1, there are
%    no guarantees this will work for other network strucutres (minor tweaks may make it work)
% -  data are processed as single-precision floats as our NVIDIA graphics card
%    is in single precision (works just fine).
%
% Nathan Stevenson
% QIMR Berghofer
% 14/10/2019

addpath('C:\QIMRBerghofer_Stevenson\Document_and_Code\papers\seizure_detection_cnn\unpack_cnn') % or whatever you have

load trained_cnn_all_sgdm_v1 % This is the trained network in the standard Matlab data structure, it contains all the INFO required to implement the CNN on an icoming time series

% Data assumed to be 16s of EEG from a single channel sampled at 256Hz
% This generates a dummy piece of EEG using a random number generator (randn)
[b,a]=cheby2(6,80,[0.5 32]/128);  % I can send on the actual filter coefficients is this helps
dum = resample(filter(b,a, randn(1,8*512)*250), 1, 8);
tx = uint16(zeros(512,1,1,1));
dum(dum>200) = 200; dum(dum<-200)=-200;
dum = dum./400*(2^16-1)+2^15+0.5;
fn = single(uint16(dum))-32768; % all this mucking about here is due to a silly original assumption I had the Matlab likes images for its CNNs so I tried to re-create a 1D, 16-bit image

% the input to the network is variable fn
% this ridiculously long-winded implementation is because I already know
% the structure of the network, this could be greatly simplified by reading
% the network strucutre as part of the implementation
fn = conv_layer(fn', convnet.Layers(2,1));  % Convolution Layer
fn = normalization_layer(fn, convnet.Layers(3,1)); % Batch Normalization Layer
fn(fn<0)=0; % Rectified Linear Unit Layer
fn = conv_layer(fn, convnet.Layers(5,1));
fn = normalization_layer(fn, convnet.Layers(6,1));
fn(fn<0)=0;
fn = conv_layer(fn, convnet.Layers(8,1));
fn = normalization_layer(fn, convnet.Layers(9,1));
fn(fn<0)=0;
fn = pooling_layer(fn, convnet.Layers(11,1));
fn = conv_layer(fn, convnet.Layers(12,1));
fn = normalization_layer(fn, convnet.Layers(13,1));
fn(fn<0)=0;
fn = conv_layer(fn, convnet.Layers(15,1));
fn = normalization_layer(fn, convnet.Layers(16,1));
fn(fn<0)=0;
fn = conv_layer(fn, convnet.Layers(18,1));
fn = normalization_layer(fn, convnet.Layers(19,1));
fn(fn<0)=0;
fn = pooling_layer(fn, convnet.Layers(21,1));  % Average Pooling layer
fn = conv_layer(fn, convnet.Layers(22,1));
fn = normalization_layer(fn, convnet.Layers(23,1));
fn(fn<0)=0;
fn = conv_layer(fn, convnet.Layers(25,1));
fn = normalization_layer(fn, convnet.Layers(26,1));
fn(fn<0)=0;
fn = conv_layer(fn, convnet.Layers(28,1));
fn = normalization_layer(fn, convnet.Layers(29,1));
fn(fn<0)=0;
fn = pooling_layer(fn, convnet.Layers(31,1));
fn = conv_layer(fn, convnet.Layers(32,1));
fn = normalization_layer(fn, convnet.Layers(33,1));
fn(fn<0)=0;
fn = conv_layer(fn, convnet.Layers(35,1));
fn = normalization_layer(fn, convnet.Layers(36,1));
fn(fn<0)=0; 
fn = convnet.Layers(38,1).Weights*fn(:)+convnet.Layers(38,1).Bias; % fully connected layer - the output of this layer is what is usually analyzed
fn = exp(fn)./sum(exp(fn)); % softmax



