function [asi1, asi2] = estimate_ASI(data, fs, href, art)
%
% ASI1 global synchrony
% ASI2 hemispheric synchrony
%
A = size(data); B = size(href);
epoch = 2.5*60*fs;
block_no = floor(A(2)/epoch);
a = ones(A(1), A(1));
asi1x = NaN*ones(1, block_no); asi2x = asi1x;
for zz = 1:block_no
    r1 = (zz-1)*epoch+1; r2 = r1+epoch-1;
    if sum(isnan(art(r1:r2)))==0
    ASIval = zeros(A(1), A(1));
    for ii = 1:A(1)
       for jj = ii+1:A(1)
           if ii~=jj
            [ASIval(ii,jj),~, ~] = getASI(data([ii jj],r1:r2)',fs);   
           end
       end
    end
    asi1x(zz) = median(ASIval(triu(a,1)==1));
    
    ASIval = zeros(1, B(1));
    for ii = 1:B(1)
        [ASIval(ii), ~, ~] = getASI(data([href(ii,1), href(ii,2)],r1:r2)',fs);
    end
    asi2x(zz) = median(ASIval);
    end
    
end
asi1 = median(asi1x(isnan(asi1x)==0));
asi2 = median(asi2x(isnan(asi2x)==0));
