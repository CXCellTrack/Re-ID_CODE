clear, close all


% �Ƚ��ڲ�ѧϰ/pspsp��ලѧϰ�£�GMMCP/MATCHLIST/NAIVE����������
% ------------------------- Main file ------------------------- % 

%% 1��������ȡ
% ������ GMMCP_Structured_Learning.m �м���
load('feature.mat');

%% 2�������ϳ�

i_cam = 1;
numClusters = 1000;
load('cluster_centers-1000.mat');

% ѡ����Լ�����
feature_test = feature.P1L;
[ feature_in_1_cam_comb, ~ ] = feature_merge( feature_test, i_cam, numClusters, cluster_centers );
% �˵����е�nan�������������ǩ
[ feature_in_1_cam_comb, labels ] = filter_nan_and_get_label( feature_in_1_cam_comb );
n_p_all = cellfun(@(x) size(x,1), feature_in_1_cam_comb,'un',1)';

% �趨Ҫʹ�õ�seq��Ŀ
n_seq = 3;
fe_to_use = feature_in_1_cam_comb(1:n_seq);
n_person_in_seq = n_p_all(1:n_seq); % ÿ��seq������
n_person_cum = [0, cumsum(n_person_in_seq)]; % �������ۼӣ��������ֿ飩
labels = labels(1:n_seq);
[ gt_connect, truth_num ] = get_gt_match_result( labels );
% ��������������W ��Խ�������ƣ�
W  = cal_W( fe_to_use, n_person_in_seq );

% �Ƿ�������GMMCP����������ÿ��cliqueΪһ���ˣ�����1���PCPSP��Լ������Ӱ�죩
global each_person_a_clique;
each_person_a_clique = 0; % ��Ч���Ͽ�GMMCP��ÿ��cliqueֻ��һ�����ǲ��õģ�

%% 3��ʹ��PCPSP���������ڻ�Kval����Ϊ���������
if 1
    % ֱ��ʹ��ŷ�Ͼ�����ΪKval
    Kval = mat2bolckmat(W, n_person_cum);
else
    max_must_link = 0;
    max_cannot_link = 0; % ����������һ��mustlink��cannotlink������kvalȫΪ1
    [Kval, seq_mustlink, seq_cannotlink] = PCPSP( W, n_person_in_seq, labels, max_must_link, max_cannot_link );
end

if each_person_a_clique
    W = log(W./(1-W)); % ��W�任Ϊ�����и���������Ϊȫ����
%     W(W==-Inf) = 0; % W�Խ�����Ϊ0
    Kval = num2cell(W);
end

%% 4��ʹ��GMMCP�������й켣����

% Ϊ�˱������в�ͬ��Ŀlink�µĽ����������Լ�����������ֿ�
% [vars, ADN] = GMMCP(Kval, n_person_in_seq);
vars_ADN_F = cal_GMMCP_constraints( n_person_in_seq );  

%% 5�����GMMCP�����о��ȼ���
Cd = 0.00; % CdΪ��ADN�ĸ�ֵ��Խ��˵��Խ����ADN�ĳ��֣�each_person_a_clique=1ʱӦ��Cd����
[ pre_connect, ADN_value ] = cal_GMMCP_object( Kval, vars_ADN_F, Cd );
% ����ÿ2��cliqueƥ��ľ���
[ ACC, TP_cell ] = cal_precision( pre_connect, gt_connect, truth_num );

%% 6��������ⷨ������Ȩֵ22����ƥ�䣨�������㷨�ȣ�
pre_connect = naive_match( Kval );
% ��22ƥ���Ч����ȫ�ֵ�Ҫ��һЩ
[ ACC, TP_cell ] = cal_precision( pre_connect, gt_connect, truth_num );

%% 7��ʹ��matchlift��5��7����ȫ��ƥ��ķ�����Ӧ����1ʹ�ã�
% ���� vars_value Ϊ�������Xin��ÿ��������Ϊһ������2��������ƥ���ϵ��Ϊһ��map
% ����һ���ֿ����
% �Ծ����������W��Ϊ����
[ X_in_mat, X_out_val_int ] = MatchLift( Kval, n_person_in_seq );
X_out = mat2bolckmat( X_out_val_int, n_person_cum );
% ��Ч����
[ ACC, TP_cell ] = cal_precision( X_out, gt_connect, truth_num );
















