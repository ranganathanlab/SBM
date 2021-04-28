function [align] = Mex_IndepMonteCarlo(N, w, q, delta_t)
% use a cpp compiled mex file to create an alignment of N seqs given a
% parameter vector w
% creation of MSA is done with N independent MC chains each with burn time of
% delta_t
    q = double(q);
    L = int32(((q^2 - 2 * q)+sqrt((2*q-q^2)^2 + 8 * size(w, 1) * q * q))/2/q/q);
    seed = randi(intmax-N);
    align = double(transpose(C_MSA_MonteCarlo(int32(seed), w, int32(N), int32(q), L, int32(delta_t))));

end

