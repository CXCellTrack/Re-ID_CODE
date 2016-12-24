function [ pre_connect, ADN_value, phi_x_z_hat, delta_zstar_zhat ] = CXSL_cal_GMMCP_object( w, vars_ADN_F, Cd, phi_x_z, loss )

% 解包
vars = vars_ADN_F.vars;
ADN = vars_ADN_F.ADN;
F = vars_ADN_F.F;

% 注意在非speedup版本中，连接到dummy node的权值为Cd（将其设为0）
% 在speedup版本中，ADN不是边变量，而是节点变量，节点权值为Cd/2（不妨设为0.01）
% Cd = 0.01;
h_c = size(vars,1);

% 2、计算目标函数并进行求解（损失增强型预测）
OBJ = dot(w, phi_x_z) + Cd/2*sum(ADN) + loss; % 把ADN的值*权重也算上

% 求解BILP
disp('    开始求解损失增强型预测问题...')
options = sdpsettings('verbose', 0, 'solver', 'cplex');
sol = solvesdp( F, -OBJ, options );
% 输出得到的各个变量的值
if sol.problem == 0      
    phi_x_z_hat = value(phi_x_z);
    delta_zstar_zhat = value(loss);
else
    sol.info
    yalmiperror(sol.problem)
end

% 当前求出的预测匹配
pre_connect = cell(size(vars));
for i1=1:h_c-1
    for i2=i1+1:h_c
        pre_connect{i1,i2} = value(vars{i1,i2});
    end
end

ADN_value = value(ADN);



end

