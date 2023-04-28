function [bcB] = bc_pairs(xs_b,QB)

nSamp = size(xs_b,1);
for s=1:nSamp

    e1 = xs_b(s,:);e2 = QB(s,:);
    f  =  f_braycurtis([e1;e2]');f = f(1,2);
    bcB(s) = f;
end


end