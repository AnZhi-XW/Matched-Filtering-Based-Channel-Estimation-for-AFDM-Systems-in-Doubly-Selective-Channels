%% 嵌入式单脉冲信道估计
%% 分数多普勒信道
%%% 多径
%%% 斐波那契搜索(FS)

clear
clc
close all
rng("shuffle")
addpath("函数文件\")

%%% 比较不同调制方案的BER性能
%%% 时变多径信道


%% parameters
M =256;             % Number of Chirp subcarriers                      
N =1;             % Number of symbols

M_mod = 4;   % size of constellation
M_bits = log2(M_mod);

Fc=4e9;             % 载波频率4GHZ
dlta_f=1e3;        % 载波间隔

Tf=1/dlta_f;
B = M*dlta_f;
Ts = 1/B;
c=3e8;             % 光速

tau_max=1.56e-05;%2.6e-06;   % 最大时延 8.2e-06;%
SPE=540;            % 速度 km/h
v_max=SPE/3.6;  %通信与雷达不同：“/2”
fd_max=v_max*Fc/c;         % 多普勒频移

% AFDM 参数
Lt=ceil(tau_max*B);   % 最大归一化时延
Kv=fd_max/dlta_f;   % 最大归一化多普勒
guard_length = 4;     % AFDM中每条路径的保护间隔长度
c1=(2*(Kv+guard_length)+1)/(2*M);
c2=1/(M^2*pi);   % c2需要是远小于1/(2*M)的有理数或任意无理数

mp=floor(M/2);    % 单脉冲导频位置
Q=(Lt+1)*(2*(Kv+guard_length)+1)-1;   % 保护间隔大小
DAFT_grid=ones(M,1);
DAFT_grid(mp-Q:mp+Q)=0;

N_zeros=2*Q; % 零保护符号数
N_syms_perfram = M-N_zeros-1;  % 每帧的数据符号数
N_bits_perfram = N_syms_perfram*M_bits;  % 每帧的数据bit数

% 信道估计区域长度
G_r=Kv+guard_length;   % 导频右边
G_l=Q-G_r; % 导频左边


L1 = diag(exp(-1j*2*pi*c1*((0:M-1).^2)));
L = diag(exp(-1j*2*pi*c2*((0:M-1).^2)));
F = dftmtx(M)/sqrt(M);
A = L*F*L1;

% 导频功率
SNRp_dB = 30;
SNRp = 10.^(SNRp_dB/10);
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
pilot_trans_power=SNRp*abs(eng_sqrt)^2;

if (2*(Kv + guard_length)*(Lt + 1)) + (Lt) > M
    fprintf("AFDM的正交性不满足!");
end

idx_start=mp-G_l-1;
idx_YT=(idx_start+1:idx_start+Q+1);
% idx_YT=(mp-G_l:mp+G_r);



%% 构造测量矩阵
idx_start=mp-G_l-1;
idx_YT=(idx_start+1:idx_start+Q+1);

dl=1;
dk=2;
tau_grid=Lt:-1/dl:0;
fd_grid=-Kv-guard_length:1/dk:Kv+guard_length;
[xq_fix,yq_fix]=meshgrid(tau_grid,fd_grid);
vec_DAFT_tau=reshape(xq_fix,[],1);
vec_DAFT_kv=reshape(yq_fix,[],1);

Phi=zeros(length(idx_YT),length(vec_DAFT_tau));
for ii=1:length(vec_DAFT_tau)
    l_tau=vec_DAFT_tau(ii);
    k_v=vec_DAFT_kv(ii);
    vec_phase1=-M*c1.*(l_tau.^2)+((mp-1)+k_v).*l_tau + M*c2*(((idx_YT-1)).^2-(mp-1)^2);
    alpha1=exp( -1i*(2*pi/M).*vec_phase1.' );

    alpha21=exp( -1i*2*pi*(0:M-1)'* ((idx_YT-1)-(mp-1)+2*M*c1*l_tau-k_v+eps)/M );
    alpha2=sum( alpha21 ).'/M;

    Phi(:,ii)=alpha1.*alpha2;
end


%%
test_num=1e2;%500;
test_SNR=0:5:30;
ber_SNR_Nfram=zeros(test_num,length(test_SNR));

for iesn0=1:length(test_SNR)
    %% Parameter
    SNR_dB = test_SNR(iesn0);%10;
    SNR=10.^(SNR_dB/10);
    sigma2 = (abs(eng_sqrt).^2)./SNR;
    Ty=3*sqrt(sigma2);                 % 信道估计阈值
    parfor ifram=1:test_num 
        %% 生成发送数据矩阵
        data_info_bit = randi([0,1],N_bits_perfram,1); % 随机生成数据比特以得到数据符号
        data = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
        data = qammod(data,M_mod,'gray');

        X=zeros(M,1);   % 将数据符号排列在DD网格中
        [data_pos,~]=find(DAFT_grid>0);
        X(data_pos)=data;
        X(mp)=sqrt(pilot_trans_power);

        %% 生成信道
        P=5;
        normDelays=randi([0,Lt],P,1);
        normDelays=sort(normDelays);
        max_f = Kv/M*B;
        Dopplers = (max_f*cos(2*pi*rand(1,P)));  % 多普勒考虑Jake's谱
        normDopplers = (Dopplers/B*M).';
        pow_prof = (1/P) * (ones(P,1));
        chan_coef = sqrt(pow_prof).*(sqrt(1/2) * (randn(P,1)+1i*randn(P,1)));

        %% AFDM调制
        s_mat=IDAFT(X,M,c1,c2);

        %% 过信道
        % 加CPP(Chirp-CP)
        L_CP=Lt;
        idx_CP=(-L_CP:-1).';
        CP_phi=repmat( exp( -1i*2*pi*c1*(M.^2+2*M*idx_CP) ) , 1,N );
        CP_temp=s_mat(end-L_CP+1:end,:);
        CCP=CP_temp.*CP_phi;
        s_mat_CP=[CCP;s_mat];
        s=s_mat_CP(:);

        % sigma2=0;
        [r,noise] = LTV_channel_output(M,L_CP,P,normDelays,normDopplers,chan_coef,s,sigma2);

        %% AFDM解调
        Y=DAFT(r,M,c1,c2);

        %% 传统LS方法
        tau_grid=Lt:-1:0;
        fd_grid=-Kv-guard_length:Kv+guard_length;
        [xq_fix,yq_fix]=meshgrid(tau_grid,fd_grid);
        vec_DAFT_tau1=reshape(xq_fix,[],1);
        vec_DAFT_v1=reshape(yq_fix,[],1);
    
        YT=Y(mp-G_l:mp+G_r);
        [h,pos]=sort(YT,"descend");
        % path_taps=length(find(abs(YT)>Ty));
        Q1=2;
        path_taps=P*(Q1+1);

        l_tau_est=vec_DAFT_tau1(pos(1:path_taps));
        k_v_est=vec_DAFT_v1(pos(1:path_taps));
        idx=idx_start+pos(1:path_taps);

        vec_phase1=-M*c1.*(l_tau_est.^2)+((mp-1)+k_v_est).*l_tau_est + M*c2*(((idx-1)).^2-(mp-1)^2);
        alpha1=exp( -1i*(2*pi/M).*vec_phase1 );
        alpha2=sum( exp(-1i*2*pi*(0:M-1)'*( ((idx-1)-(mp-1)+2*M*c1*l_tau_est-k_v_est+eps)/M).' ) ).'/M;  % 分数多普勒的相偏
        chan_phase=alpha1.*alpha2;
        chan_coef_est=h(1:path_taps)./sqrt(pilot_trans_power)./chan_phase;
        [H_est_LS]=hw_matrix_DAFT(M,l_tau_est, k_v_est, chan_coef_est,c1,c2);
        

        %% MF-GFS搜索(GFS搜索)
        max_iter=path_taps;
        tol=1e-3;
        tol_error=1e-4;
        n_order_F=15;     % 广义斐波那契数的阶数
        [l_tau_est_MF,k_v_est_MF,chan_coef_est_MF] = MF_GFS_estimator( M,Y,Lt,Kv,idx_start,idx_YT,guard_length,mp,G_l,G_r,max_iter,tol,tol_error,c1,c2,pilot_trans_power,n_order_F);
        [H_est_MF_GFS]=hw_matrix_DAFT(M,l_tau_est_MF, k_v_est_MF, chan_coef_est_MF,c1,c2);
        

        %% 计算NMSE
        [H_true]=hw_matrix_DAFT(M,normDelays, normDopplers, chan_coef,c1,c2);
        
        nmse_LS(ifram,iesn0)=norm(H_true-H_est_LS,'fro')^2/norm(H_true,'fro')^2;
        nmse_MF_GFS(ifram,iesn0)=norm(H_true-H_est_MF_GFS,'fro')^2/norm(H_true,'fro')^2;

        fprintf('AFDM-MF: %d   %d\n',iesn0,ifram);
    end
end
NMSE_LS=mean(nmse_LS);
NMSE_MF_GFS=mean(nmse_MF_GFS);

%% Compare
figure(1)
plot(test_SNR, 10*log10(NMSE_LS),'--g',"LineWidth",1);
hold on
plot(test_SNR, 10*log10(NMSE_MF_GFS),'--Or',"LineWidth",1);
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
grid on
legend("LS","MF-GFS")

rmpath("函数文件\")

