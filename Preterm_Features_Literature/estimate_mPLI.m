function mPLIx = estimate_mPLI(data, fs, art)

% assume incoming data has already been HPFed at 0.5Hz
A = size(data);
[Num1, Den1] = butter(8, 8/fs);
[Num2, Den2] = butter(4, 1/fs, 'high');
for ii = 1:A(1)
    data(ii,:) = filter(Num1, Den1, filter(Num2, Den2, data(ii,:)));
end

epoch = 16*fs;
block_no = floor(A(2)/epoch);
mPLI = zeros(1, block_no);
for zz = 1:block_no
    r1 = (zz-1)*epoch+1; r2 = r1+epoch-1;
    pli = zeros(A(1), A(1));
    for ii = 1:A(1)
       for jj = 1:A(1)
           if ii~=jj
            phase1 = unwrap(angle(hilbert(data(ii,r1:r2))));
            phase2 = unwrap(angle(hilbert(data(jj,r1:r2))));
            pli(ii,jj) = abs(mean(sign(sin(phase1-phase2))));
           end
       end
    end
    mPLI(zz) = sum(sum(pli))/(length(pli(:))-A(1));
end

epoch = 16*fs;
block_no = floor(A(2)/epoch);
nart = zeros(1, block_no);
for zz = 1:block_no
    r1 = (zz-1)*epoch+1; r2 = r1+epoch-1;
    nart(zz) = sum(isnan(art(r1:r2)));
end
mPLIx = median(mPLI(nart==0));

