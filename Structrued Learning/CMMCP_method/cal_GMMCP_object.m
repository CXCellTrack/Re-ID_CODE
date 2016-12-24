function [ pre_connect, ADN_value ] = cal_GMMCP_object( Kval, vars_ADN_F, Cd )

% 解包
vars = vars_ADN_F.vars;
ADN = vars_ADN_F.ADN;
F = vars_ADN_F.F;

% 注意在非speedup版本中，连接到dummy node的权值为Cd（将其设为0）
% 在speedup版本中，ADN不是边变量，而是节点变量，节点权值为Cd/2（不妨设为0.01）
% Cd = 0.01;
weight = Kval;
n_seq = size(vars,1);

% 5、计算目标函数并进行求解
OBJ = 0;
for i1=1:n_seq-1
    for i2=i1+1:n_seq
        OBJ = OBJ + sum(sum( vars{i1,i2}.*weight{i1,i2} ));
    end
end
OBJ = OBJ + Cd/2*sum(ADN); % 把ADN的值*权重也算上

% 求解BILP
disp('开始求解整数混合规划...')
options = sdpsettings('verbose',0,'solver','cplex');
sol = solvesdp( F, -OBJ, options )
if sol.problem == 0      
    % 当前求出的预测匹配
    pre_connect = cell(size(vars));
    for i1=1:n_seq-1
        for i2=i1+1:n_seq
            pre_connect{i1,i2} = value(vars{i1,i2});
        end
    end
else
    sol.info
    yalmiperror(sol.problem)
end

ADN_value = value(ADN);



end

