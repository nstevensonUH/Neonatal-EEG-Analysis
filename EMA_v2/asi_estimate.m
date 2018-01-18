function asi = asi_estimate(n1, n2, art, fs);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%n1(art==1)=0; n2(art==1) = 0;
n1 = filter([1 -0.95],1, n1);
n2 = filter([1 -0.95],1, n2);
epl = 2*fs; olap = 0.125*fs; M = length(n1);
block_no = floor(M/olap)-15; env1 = zeros(1,block_no); env2 = env1; aa = env1;
for z1 = 1:block_no;
    r1 = 1+(z1-1)*olap; r2 = r1+epl-1;
    X1 = abs(fft(n1(r1:r2).*hamming(epl)'))/epl;
    X2 = abs(fft(n2(r1:r2).*hamming(epl)'))/epl;
    env1(z1)=sum(X1(4:40));
    env2(z1)=sum(X2(4:40));
    aa(z1) = sum(art(r1:r2));
end

env1 = log(env1(aa==0)); env2 = log(env2(aa==0));
bin_no = floor(sqrt(sqrt(length(env1))));
range = [floor(min([env1 env2])) ceil(max([env1 env2]))];
dr = diff(range);
xe = linspace(range(1), range(2), bin_no+1); ye = xe;
xc = xe(1:end-1)+(xe(2)-xe(1))/2; yc = xc;

lags = -40:40; edtf = zeros(1,length(lags));
for kk = 1:length(lags);
    x = circshift(env2', lags(kk));
    histmat  = hist2(env1, x', xe, ye)+1;
    histmat = histmat(1:end-1, 1:end-1);
    histmat = histmat./sum(sum(histmat))/dr^2;
    histx = sum(histmat')*dr;
    histy = sum(histmat)*dr;
    edtf(kk) = sum(sum((xc'*yc)  .*  (histmat.^2) ./ (histx'*histy))) ;
end

en = (edtf-min(edtf));
asi = en(lags==0)/mean(en);

end



