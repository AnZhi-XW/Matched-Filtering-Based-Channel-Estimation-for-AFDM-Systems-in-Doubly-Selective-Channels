function [h_w]=hw_matrix_DAFT(M,l_tau,k_v,chan_coef,c1,c2)
%%% 计算DAFT域有效信道矩阵
%   M：              子载波数
%   l_tau：          归一化时延
%   k_v：            归一化多普勒
%   chan_coef：      信道增益
%   c1,c2：          DAFT参数

h_w=0;
idx_Y=0:M-1;
idx_X=(0:M-1);
for itap=1:length(l_tau)
    l=l_tau(itap);
    k=k_v(itap);

    vec_phase1=-M*c1.*(l.^2)+(idx_X+k).*l + M*c2*((idx_Y.').^2-idx_X.^2);
    alpha1=exp( -1i*(2*pi/M).*vec_phase1 );
    alpha21=( exp(-1i*2*pi*(idx_Y.'-idx_X+2*M*c1*l-k+eps))-1);
    alpha22=( exp(-(1i*2*pi/M)*(idx_Y.'-idx_X+2*M*c1*l-k+eps))-1);
    alpha2=alpha21./alpha22./M;
    alpha=alpha1.*alpha2;

    h_w=h_w+chan_coef(itap)*alpha;
end

end

