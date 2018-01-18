function art = find_artefacts_amp(ndm, bursts, ch_left, ch_right);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fs = 64;
Nl = length(ch_left); Nr = length(ch_right);
cc = [ch_left ch_right];
art = zeros(Nl+Nr, length(bursts{1, cc(1)}));
for z1 = 1:Nl;
    ba = bursts{1,ch_left(z1)};
    r1 = ceil(find(diff([0 ba 0])==1));
    r2 = floor((find(diff([0 ba 0])==-1)-1));
    c1 = 1; be = zeros(1,length(r1));
     for z2 = 1:length(r1);
         be(c1) = max(abs(ndm(ch_left(z1), r1(z2):r2(z2))));
         c1 = c1+1;
     end
     
     clear r1 r2

     %ref = find(be>1.5*diff(quantile(be, [0.25 0.75]))+quantile(be, 0.75));
     ref = find(be>500);
     r1 = find(diff([0 ba 0])==1);
     r2 = find(diff([0 ba 0])==-1)-1;
     for z3 = 1:length(ref);
         art(ch_left(z1), r1(ref(z3)):r2(ref(z3))) = 1;
     end
     
     dt = abs(diff(ndm(ch_left(z1),:)));
     rr = find(dt>200);
     for z4 = 1:length(rr);
         s1 = rr(z4)-fs; s2 = s1+2*fs-1;
         if s1<1; s1 = 1; end; if s2>length(ndm); s2 = length(ndm); end
        art(ch_left(z1), s1:s2) = 1;
     end
     
     clear be r1 r2 ref
     
  %   dum = abs(diff())
     
end



for z1 = 1:Nr;
    ba = bursts{1,ch_right(z1)};
    r1 = ceil(find(diff([0 ba 0])==1));
    r2 = floor((find(diff([0 ba 0])==-1)-1));
    c1 = 1; be = zeros(1,length(r1));
     for z2 = 1:length(r1);
         be(c1) = max(abs(ndm(ch_right(z1), r1(z2):r2(z2))));
         c1 = c1+1;
     end
     
     clear r1 r2

     %ref = find(be>1.5*diff(quantile(be, [0.25 0.75]))+quantile(be, 0.75));
     ref = find(be>500);
     r1 = find(diff([0 ba 0])==1);
     r2 = find(diff([0 ba 0])==-1)-1;
     for z3 = 1:length(ref);
         art(ch_right(z1), r1(ref(z3)):r2(ref(z3))) = 1;
     end
     
     dt = abs(diff(ndm(ch_right(z1),:)));
     rr = find(dt>200);
     for z4 = 1:length(rr);
         s1 = rr(z4)-fs; s2 = s1+2*fs-1;
         if s1<1; s1 = 1; end; if s2>length(ndm); s2 = length(ndm); end
        art(ch_right(z1), s1:s2) = 1;
     end
     
     clear be r1 r2 ref
     
end





end

