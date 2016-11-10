function [ pre_connect, phi_x_zhat, delta_zstar_zhat ] = CXSL_cal_MatchLift_object( w, vars_ADN_F, phi_x_z, loss, n_person_cum )

% 解包
vars = vars_ADN_F.vars;
F = vars_ADN_F.F;

%% map图中，边的个数（注意这个边是2个对象的边）此处为2个场景的连接
% 见章节2.1
n_seq = numel(n_person_cum)-1;
E = n_seq*(n_seq-1)/2; % (Si,Sj)∈ε，不重复计算
if ~exist('lambda_l1', 'var')
    lambda_l1 = sqrt(E)/2/n_seq; % 正则参数lambda_l1
end

disp('      开始求解损失增强型预测问题...')
obj = dot(w, phi_x_z) - lambda_l1*sum(vars(:)) + loss; % 正则为L1，即X_out中1的个数
options = sdpsettings('verbose', 0, 'solver','mosek'); % 商业求解器 mosek 比 sedumi 快很多
sol = solvesdp( F, -obj, options ); % checkset(F)

%% 输出得到的各个变量的值
if sol.problem == 0      
    phi_x_zhat = value(phi_x_z);
    delta_zstar_zhat = value(loss);
else
    sol.info
    yalmiperror(sol.problem)
end

X_out = value(vars);
X_out_int = round(X_out);
pre_connect = mat2bolckmat(X_out_int, n_person_cum );

