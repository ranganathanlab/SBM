function [W, N_eff] = CalcWeights(align,theta)
% calculate reweighting weights per sequence and effective number of seqs
    W = (1./(1+sum(squareform(pdist(align, 'hamming')<theta))));
    N_eff=sum(W);
end

