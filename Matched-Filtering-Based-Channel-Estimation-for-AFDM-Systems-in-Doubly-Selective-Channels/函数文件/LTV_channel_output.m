function [r,noise] = LTV_channel_output(M,L_CP,taps,l_tau,f,chan_coef,s,sigma2)
%% 双选信道 和 高斯白噪声
s_chan=0;

for itao = 1 : taps
    r_itao=[ s.*exp(1i*(2*pi/M)*f(itao).* (-L_CP:-L_CP+length(s)-1) ).' ; zeros(max(l_tau),1) ];
    s_chan = s_chan + chan_coef( itao ) .* circshift( r_itao, l_tau(itao) );
end

noise=sqrt(sigma2/2)*(randn(size(s_chan)) + 1i*randn(size(s_chan)));
r=s_chan+noise;
r=r(L_CP+1:L_CP+M);
end



