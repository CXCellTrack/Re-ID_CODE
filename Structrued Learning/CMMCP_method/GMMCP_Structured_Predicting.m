

% w_best

%% 2�������ϳ�
% ѡ����Լ�����
Port = 'P2L';
disp(['choose ', Port]);
feature_test = getfeild(feature, Port);
[ feature_in_1_cam_comb, ~ ] = feature_merge( feature_test, i_cam, numClusters, cluster_centers );
% �˵����е�nan�������������ǩ
[ feature_in_1_cam_comb, P_labels ] = filter_nan_and_get_label( feature_in_1_cam_comb );
n_p_all = cellfun(@(x) size(x,1), feature_in_1_cam_comb,'un',1)';

% �趨Ҫʹ�õ�seq��Ŀ
n_seq = 3;
P_fe_to_use = feature_in_1_cam_comb(1:n_seq);
P_n_person_in_seq = n_p_all(1:n_seq); % ÿ��seq������
P_labels = P_labels(1:n_seq);

% ��ȡ���ֽ���ƥ��
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
% ע���ڷ�speedup�汾�У����ӵ�dummy node��ȨֵΪCd��������Ϊ0��
% ��speedup�汾�У�ADN���Ǳ߱��������ǽڵ�������ڵ�ȨֵΪCd/2��������Ϊ0.01��
% ���
vars = vars_ADN_F.vars;
ADN = vars_ADN_F.ADN;
F = vars_ADN_F.F;

Cd = 0;
n_seq = numel(P_n_person_in_seq);

%% 2������Ŀ�꺯����������⣨��ʧ��ǿ��Ԥ�⣩
OBJ = dot(w_best, phi_x_z) + Cd/2*sum(ADN); % ��ADN��ֵ*Ȩ��Ҳ����

% ���BILP
disp('��ʼ���Ԥ������...')
options = sdpsettings('verbose',0,'solver','gurobi');
sol = solvesdp( F, -OBJ, options )
% ����õ��ĸ���������ֵ
if sol.problem == 0      
    % ��ǰ�����Ԥ��ƥ��
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

%% ���㾫��
ADN_value = value(ADN);
[ ACC, TP_cell ] = cal_precision( pre_connect, gt_connect, truth_num );
  
