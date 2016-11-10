

% w_best

%% 2、特征合成
% 选择测试集特征
Port = 'P2L';
disp(['choose ', Port]);
feature_test = getfeild(feature, Port);
[ feature_in_1_cam_comb, ~ ] = feature_merge( feature_test, i_cam, numClusters, cluster_centers );
% 滤掉其中的nan，并获得人脸标签
[ feature_in_1_cam_comb, P_labels ] = filter_nan_and_get_label( feature_in_1_cam_comb );
n_p_all = cellfun(@(x) size(x,1), feature_in_1_cam_comb,'un',1)';

% 设定要使用的seq数目
n_seq = 3;
P_fe_to_use = feature_in_1_cam_comb(1:n_seq);
P_n_person_in_seq = n_p_all(1:n_seq); % 每个seq中人数
P_labels = P_labels(1:n_seq);

% 截取部分进行匹配
if 0
    fe_s = 5; 
    fe_e = 20;
    P_fe_to_use = cellfun(@(x) x(fe_s:fe_e,:), P_fe_to_use,'un',0);
    P_n_person_in_seq = ones(size(P_n_person_in_seq))*(fe_e-fe_s+1);
    P_labels = cellfun(@(x) x(fe_s:fe_e,:), P_labels,'un',0);
end

[ gt_connect, truth_num ] = get_gt_match_result( P_labels );


%%
vars_ADN_F = cal_GMMCP_constraints( P_n_person_in_seq );  
[ phi_x_z, phi_x_zstar, loss ] = CXSL_cal_phi_and_loss( vars_ADN_F, P_fe_to_use, gt_connect, truth_num );
% 注意在非speedup版本中，连接到dummy node的权值为Cd（将其设为0）
% 在speedup版本中，ADN不是边变量，而是节点变量，节点权值为Cd/2（不妨设为0.01）
% 解包
vars = vars_ADN_F.vars;
ADN = vars_ADN_F.ADN;
F = vars_ADN_F.F;

Cd = 0;
n_seq = numel(P_n_person_in_seq);

%% 2、计算目标函数并进行求解（损失增强型预测）
OBJ = dot(w_best, phi_x_z) + Cd/2*sum(ADN); % 把ADN的值*权重也算上

% 求解BILP
disp('开始求解预测问题...')
options = sdpsettings('verbose',0,'solver','gurobi');
sol = solvesdp( F, -OBJ, options )
% 输出得到的各个变量的值
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

%% 计算精度
ADN_value = value(ADN);
[ ACC, TP_cell ] = cal_precision( pre_connect, gt_connect, truth_num );
  
