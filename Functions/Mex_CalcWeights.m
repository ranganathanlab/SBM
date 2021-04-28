function [W, N_eff] = Mex_CalcWeights(align,theta)
% calculate reweighting weights per sequence and effective number of seqs
    out=C_InverseWeights(int32(align),theta);
    W=1./out;
    N_eff=sum(W);
end

