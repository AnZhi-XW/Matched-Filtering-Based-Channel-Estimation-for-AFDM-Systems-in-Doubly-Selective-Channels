function [k_est] = GFS( YT,M,n_order_F,a,b,p,q,k_est,idx_YT,mp,c1,c2,l_est,tol_error)
% GFS搜索

GF_set=[];
for ii=1:n_order_F
    if ii==1
        GF_set(ii)=a;
    elseif ii==2
        GF_set(ii)=b;
    else
        GF_set(ii)=p*GF_set(ii-1)+q*GF_set(ii-2);
    end
end
x_start = k_est-0.5;
x_finish = k_est+0.5;
d = x_finish - x_start;
L = (q*GF_set(n_order_F-2)/GF_set(n_order_F))*d;
x1 = x_start + L;
x2 = x_finish - L;
% x1 = (p*GF_set(3)/GF_set(4))*x_start + (q*GF_set(2)/GF_set(4))*x_finish;
% x2 = (q*GF_set(2)/GF_set(4))*x_start + (p*GF_set(3)/GF_set(4))*x_finish;

idx_Y=idx_YT-1;
idx_X=mp-1;
func = @(x) -ObjectiveFunction(M,c1,c2,l_est,x,idx_Y,idx_X,YT);
f1=func(x1);
f2=func(x2);

for k = 0:n_order_F-3
    if func(x1) < func(x2)
        x_finish = x2;
        d = x_finish - x_start;
        L = (q*GF_set(n_order_F-k-2)/GF_set(n_order_F-k))*d;
        x2 = x1;
        x1 = x_start + L;
        f2=f1;
        f1=func(x1);

    elseif func(x1) == func(x2)
        x_start = x1;
        x_finish = x2;

        d = x_finish - x_start;
        L = (q*GF_set(n_order_F-k-2)/GF_set(n_order_F-k))*d;
        x1 = x_start + L;
        x2 = x_finish - L;
        f1=func(x1);
        f2=func(x2);

    else %func(x1) > func(x2)
        x_start = x1;

        d = x_finish - x_start;
        L = (q*GF_set(n_order_F-k-2)/GF_set(n_order_F-k))*d;
        x1 = x2;
        x2 = x_finish - L;
        f1=f2;
        f2=func(x2);
    end
    % x1 = (p*GF_set(3+k)/GF_set(4+k))*x_start + (q*GF_set(2+k)/GF_set(4+k))*x_finish;
    % x2 = (q*GF_set(2+k)/GF_set(4+k))*x_start + (p*GF_set(3+k)/GF_set(4+k))*x_finish;

    error = x_finish - x_start;
    if error<=tol_error
        break
    end
end
k_est = (x_finish + x_start)/2;
end



function h_w = ObjectiveFunction(M,c1,c2,l_tau_est,k_v_est,idx_Y,idx_X,YT)

l=l_tau_est;
k=k_v_est;

vec_phase1=-M*c1.*(l.^2)+(idx_X+k).*l + M*c2*((idx_Y.').^2-idx_X.^2);
alpha1=exp( -1i*(2*pi/M).*vec_phase1 );
alpha21=( exp(-1i*2*pi*(idx_Y.'-idx_X+2*M*c1*l-k+eps))-1);
alpha22=( exp(-(1i*2*pi/M)*(idx_Y.'-idx_X+2*M*c1*l-k+eps))-1);
alpha2=alpha21./alpha22./M;
alpha=alpha1.*alpha2;

h_w=(abs(alpha'*YT)).^2;

end


