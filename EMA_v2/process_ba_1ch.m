function burst_anno = process_ba_1ch(burst_anno, fs2, tl1, tl2, flag);
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

lim1 = fs2*tl1; lim2 = fs2*tl2;

if flag == 1;
    
    r1 = find(diff([0 burst_anno 0]) == 1);
    r2 = find(diff([0 burst_anno 0]) == -1);
    tst2 = r2-r1;
    rf = find(tst2<lim1);
    for z2 = 1:length(rf)
        burst_anno(r1(rf(z2)):r2(rf(z2)))=0;
    end
    
    r1 = find(diff([0 burst_anno 0]) == 1);
    r2 = find(diff([0 burst_anno 0]) == -1);
    r1d = r1(2:end); r2d = r2(1:end-1);
    tst1 = r1d-r2d;
    rf = find(tst1<lim2);
    for z2 = 1:length(rf)
        burst_anno(r2d(rf(z2)):r1d(rf(z2)))=1;
    end
    
else

    r1 = find(diff([0 burst_anno 0]) == 1);
    r2 = find(diff([0 burst_anno 0]) == -1);
    %if r1(1) == 1; r1 = r1(2:end); r2 = r2(2:end); end
    %if r2(end) == length(burst_anno)+1; r1 = r1(1:end-1); r2 = r2(1:end-1); end
    r1d = r1(2:end); r2d = r2(1:end-1);
    tst1 = r1d-r2d;
    rf = find(tst1<lim2);
    for z2 = 1:length(rf)
        burst_anno(r2d(rf(z2)):r1d(rf(z2)))=1;
    end

    r1 = find(diff([0 burst_anno 0]) == 1);
    r2 = find(diff([0 burst_anno 0]) == -1);
    %if r1(1) == 1; r1 = r1(2:end); r2 = r2(2:end); end
    if r2(end) == length(burst_anno)+1; r1 = r1(1:end-1); r2 = r2(1:end-1); end
    tst2 = r2-r1;
    rf = find(tst2<lim1);
    for z2 = 1:length(rf)
        burst_anno(r1(rf(z2)):r2(rf(z2)))=0;
    end

end

