
% GMMCP�Ľṹ��ѧϰ���
clear; close all

diary('learning_log.txt')
diary on 

%% 1��������ȡ
dataaddr = 'E:\����ʶ�����ݼ�\ChokePoint\cropped_face\P1E_cropped_face\P1E\';
if exist('feature.mat', 'file')
    load('feature.mat');
else
    feature = get_feature_wrapper(dataaddr);
    save('feature.mat', 'feature');
end

%% 2�������ϳ�

% ��һ���˵���Ƶ���кϳ�һ������
i_cam = 1;
[ feature_in_1_cam_comb ] = feature_merge( feature, i_cam );
% �˵����е�nan�������������ǩ
[ feature_in_1_cam_comb, labels ] = filter_nan_and_get_label( feature_in_1_cam_comb );
n_p_all = cellfun(@(x) size(x,1), feature_in_1_cam_comb,'un',1)';

% �趨Ҫʹ�õ�seq��Ŀ
n_seq = 3;
n_person_in_seq = n_p_all(1:n_seq); % ÿ��seq������
n_person_cum = [0, cumsum(n_person_in_seq)]; % �������ۼӣ��������ֿ飩
labels = labels(1:n_seq);
% ��ȡgt�����ӷ�ʽ
[ gt_connect, truth_num ] = get_gt_match_result( labels );

% ��ѡ����seq����pca��ά
im_mat = cell2mat(feature_in_1_cam_comb(1:n_seq));
[X, mean_face, COEFF] = pca_eig( im_mat, 0.975 ); % imshow(reshape(mean_face,96,96))
feature_in_1_cam_comb_pca = mat2cell(X', n_person_in_seq, size(X',2));

fe_to_use = feature_in_1_cam_comb_pca(1:n_seq);

%% ʹ������ƥ��Ľ��������m
% ��������������W ��Խ�������ƣ�
W  = cal_W( fe_to_use, n_person_in_seq );
Kval = mat2bolckmat(W, n_person_cum);
pre_connect = naive_match( Kval );
% ��������ƥ��Ľ��������m
[ m_hat, X_in_mat ] = estimate_m_by_inputX( pre_connect, n_person_in_seq );

%% �������������Լ����phi
N = 1; % ��ʱֻ��һ��������4��������

vars_ADN_F = cal_MatchLift_constraints( X_in_mat, m_hat, n_person_in_seq );
[ phi_x_z, phi_x_zstar, loss ] = MatchLift_cal_phi_and_loss( vars_ADN_F, fe_to_use, gt_connect, truth_num, n_person_cum );

%% ��������ռ䣬�������
rng(0)
w = (rand(size(phi_x_z))-0.5)*10;
w = zeros(size(phi_x_z));

iter = 500;
W = cell(iter,1); % W ����ۺ�Ȩֵw
Wavg = cell(iter,1);
Wi = cell(iter,N); % Wi�������Ȩֵw
Wavg{1} = w;
W{1} = w; % ȫ��������W��Ҫ�趨��ֵw
for i=1:N
    Wi{1,i} = w; % Wi ���ÿ��ѭ�����ض��������º��Wi
end
L = zeros(iter,1); % L ����ۺ�L
Li = zeros(iter,N);
for i=1:N
    Li(1,i) = 0;
end

lambda = 1e0;
pre_connect = cell(iter,1); % ÿ��Ԥ����
gamma = zeros(iter,1); % ����gamma
time_consume = zeros(iter,1); % ��¼ÿ��ѭ�����õ�ʱ��
aver_loss = zeros(iter,1); % ��¼ÿһ����������ʧ������ֵ
linesearch = 0; % �Ƿ�ʹ������������
t = 0;


%% ��ʼѭ��
while t<iter
    t = t + 1;
    tic;
    disp('  ==========================');
    disp(['  ��ʼ�� ', num2str(t), ' ��ѭ��...']);
    % ѡ�е�ind��������Ϊѵ������
    ind = randi(N); 

    % �����ٶ�������dot(w,phi)���𣿣�����˲�ʹ��������ʽ������ʱ����<w,feature>�������ΪĿ�꺯��
    [ pre_connect{t}, phi_x_zhat, delta_zstar_zhat ] = CXSL_cal_MatchLift_object( Wavg{t}, vars_ADN_F, phi_x_z, loss, n_person_cum );
    disp('      ����Ȩ���� w...');
    
    U_x_zstar_zhat = phi_x_zstar - phi_x_zhat;
    Ws = U_x_zstar_zhat/(lambda*N);
    ls = delta_zstar_zhat/N;
    aver_loss(t) = delta_zstar_zhat;
    fprintf('      ��ǰ������ʧ������(z*,z^):\t%f\n', delta_zstar_zhat);
    
    % ����gap��gap��ֵ��������lambda����仯���޷�ȷ����������˻�����loss��gap�ȽϺ���
    if linesearch
        tmp = lambda*(Wi{t,ind}- Ws)'*W{t}- Li(t,ind)+ ls;
        gamma(t) = tmp/(lambda*norm(Wi{t,ind}- Ws)^2);
        gamma(t) = max([0, min([gamma(t),1])]);
    else
        gamma(t) = 2*N/(2*N + t-1);
    end

    % ���� wi��Li�������º��w������ W{t+1,ind}��
    % ����Wi��˵��������һ�β�һ��ѡ����i�����Wi{t,ind}����Ϊ�յģ��ᵼ�¸��µ�ʱ���ƺ������⣿
    % 2015.10.23 ������ʵû�б�ѡ�е�������ֱ�ӽ�Wi���������֣�
    % ֱ�ӽ�����������Wi������һ��
    for ii=1:N
        Wi{t+1,ii} = Wi{t,ii}; 
        Li(t+1,ii) = Li(t,ii);
    end
    % �ٸ��µ���������Wi
    Wi{t+1,ind} = (1- gamma(t))*Wi{t,ind}+ gamma(t)*Ws; 
    Li(t+1,ind) = (1- gamma(t))*Li(t,ind)+ gamma(t)*ls;

    % ���� w��L�������º��w������ W{t+1,N+1}��
    W{t+1} = W{t} + Wi{t+1,ind} - Wi{t,ind};
    L(t+1) = L(t) + Li(t+1,ind) - Li(t,ind); % ����L�ƺ�ûʲô�ã�
    
    % ����Wavg
    Wavg{t+1} = (t-1)/(t+1)*Wavg{t} + 2/(t+1)*W{t+1};
    
    % ��¼ʱ��
    time_consume(t) = toc;
    fprintf('      ʱ�仨��:\t%1.2f s\n', time_consume(t)); 
end

%% ��������
% ������ʧ����
plot(aver_loss, '-*');
t_best = find(aver_loss==min(aver_loss));
[ ACC, TP_cell ] = cal_precision( pre_connect{t_best}, gt_connect, truth_num );


save('aver_loss.mat', 'aver_loss')
diary off;
system('shutdown /s');


