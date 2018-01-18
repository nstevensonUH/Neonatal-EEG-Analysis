function va = estimate_vigilance_state_new(bursts, ch_left, ch_right, M, art)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fs2 = 64;
Nl = length(ch_left); Nr = length(ch_right);
block_no = floor(length(bursts(ch_left(1),:))/(M/5))-4;
SATper = zeros(Nl+Nr,block_no);% AMPper = SATper;

c1 = 1; 
for z1 = 1:Nl;
    ba = bursts(ch_left(z1),:);
    %nd = ndm(ch_left(z1), :);
     for z2 = 1:block_no;
            r1 = (z2-1)*M/5+1; r2 = r1+M-1;
           SATper(c1, z2) = sum(ba(r1:r2))/M*100;
           a(z2) = sum(art(r1:r2));
         %  AMPper(c1, z2) = quantile(abs(hilbert(nd(r1:r2))), 0.1);
      end
    c1 = c1+1;
end


for z1 = 1:Nr;
    ba = bursts(ch_right(z1),:);
    %nd = ndm(ch_right(z1), :);
     for z2 = 1:block_no;
            r1 = (z2-1)*M/5+1; r2 = r1+M-1;
           SATper(c1, z2) = sum(ba(r1:r2))/M*100;
           %AMPper(c1, z2) = quantile(abs(hilbert(nd(r1:r2))), 0.1);
     end
    c1=c1+1;
end

err = zeros(1, block_no);
for z2 = 1:block_no;
       r1 = (z2-1)*M/5+1; r2 = r1+M-1;
       if sum(art(r1:r2))>0.1*M; err(z2) = NaN; end
 end


x2 = detrend(median(SATper)');
x1 = x2; x1(isnan(err)==1)=NaN;
rng = quantile(x1, [0.05 0.95]);
th = linspace(rng(2), rng(1), 100);
va1 = zeros(1, length(x1)); 

if sum(isnan(err))<26
va1 = zeros(1, length(x1)); c1 = 1;
while length(find(va1(isnan(err)==0) == 0))>13;
    va1 = x1>th(c1);
    va1 = double(process_ba_1ch(va1', 1, 5, 5, 0));
    va1(isnan(err)==1)=NaN;
    c1 = c1+1;
end
end
va1 = x1>th(c1-2);
va1 = double(process_ba_1ch(va1', 1, 5, 5, 0));
va1(isnan(err)==1)=NaN;

   va = zeros(1,length(bursts(ch_left(1),:)));
     for z2 = 1:block_no;
            r1 = (z2-1)*M/5+1; r2 = r1+M-1;
           va(r1:r2) = va1(z2);
     end
     va(r2:end) = va1(end);


end

