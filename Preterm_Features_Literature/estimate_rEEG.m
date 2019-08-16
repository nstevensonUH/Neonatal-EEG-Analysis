function reeg = estimate_rEEG(nd1, fs2);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

epl = fs2*2; N = length(nd1);
block_no = N/epl;
for ii = 1:block_no;
    q1 = (ii-1)*epl+1; q2 = q1+epl-1;
    reeg(ii) = max(nd1(q1:q2))-min(nd1(q1:q2));
end

end

