function [l_tau_est_MF,k_v_est_MF,chan_coef_est_MF] = MF_estimator( M,Y,Lt,Kv,idx_start,idx_YT,guard_length,mp,G_l,G_r,max_iter,n_v,tol,c1,c2,pilot_trans_power)
%%% MF信道估计

l_tau_est_MF=[];
k_v_est_MF=[];
chan_coef_est_MF=[];

avg_power=(norm(Y,"fro").^2)./M;

for i=1:max_iter
    Y_old=Y;

    YT=Y(mp-G_l:mp+G_r);

    tau_grid1=Lt:-1:0;
    fd_grid1=-Kv-guard_length:Kv+guard_length;
    [xq_fix,yq_fix]=meshgrid(tau_grid1,fd_grid1);
    vec_DAFT_tau1=reshape(xq_fix,[],1);
    vec_DAFT_v1=reshape(yq_fix,[],1);

    [h,pos]=sort(YT,"descend");
    % path_taps=length(find(abs(YT)>Ty));

    l_est=vec_DAFT_tau1(pos(1));
    k_est=vec_DAFT_v1(pos(1));
    idx=idx_start+pos(1);
    % plot(abs(Y))
    l_tau_est_MF=[l_tau_est_MF;l_est];


    Lamda_k=k_est+(-0.5:1/n_v:0.5).';
    vec_phase1=-M*c1.*(l_est.^2)+((mp-1)+Lamda_k).*l_est + M*c2*(((idx_YT-1)).^2-(mp-1)^2);
    alpha1=exp( -1i*(2*pi/M).*vec_phase1.' );
    vec_phase2=( ((idx_YT-1)-(mp-1)+2*M*c1*l_est-Lamda_k+eps)/M).';
    C = num2cell(vec_phase2.', 2);
    alpha21=arrayfun(@(x) sum( exp(-1i*2*pi*(0:M-1)'* cell2mat(x)) ).'/M, C, 'UniformOutput', 0)'; % 分数信道的相偏
    alpha2=cell2mat(alpha21);
    alpha=alpha1.*alpha2;

    tau_c_l=abs(alpha'*YT).^2;
    [~,idx_v]=max(tau_c_l);
    k_est=Lamda_k(idx_v);
    k_v_est_MF=[k_v_est_MF;k_est];


    vec_phase1=-M*c1.*(l_est.^2)+((mp-1)+k_est).*l_est+M*c2*(((0:M-1).').^2-(mp-1)^2);
    alpha1=exp( -1i*(2*pi/M).*vec_phase1 );
    alpha2=sum( exp(-1i*2*pi*(0:M-1)'*( ((0:M-1)-(mp-1)+2*M*c1*l_est-k_est+eps)/M) ) ).'/M;  % 分数多普勒的相偏

    % % CCP的相偏
    % alpha31=exp( -1i*2*pi*c1*(M^2-2*M*l_est) );
    % alpha3=ones(M,1);
    % alpha3(1:l_est+1)=alpha31;

    alpha=alpha1.*alpha2;
    chan_coef_est1=alpha'*Y/sqrt(pilot_trans_power);
    chan_coef_est_MF=[chan_coef_est_MF;chan_coef_est1];

    yi=chan_coef_est1*alpha*sqrt(pilot_trans_power);
    Y=Y-yi;

    error= norm( ((Y'*Y)/avg_power) - ((Y_old'*Y_old)/avg_power) );
    if error<=tol
        break
    end

end
end


