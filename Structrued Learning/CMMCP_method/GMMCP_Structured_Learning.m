
% GMMCP�Ľṹ��ѧϰ���
clear; close all

%% 1��������ȡ
baseaddr = 'E:\����ʶ�����ݼ�\ChokePoint\cropped_face\';
savename = 'feature.mat';
if exist(savename, 'file')
    load(savename);
else
    feature = struct;
    is_front = struct;
    for ds={'P1E','P1L','P2E','P2L'}
        disp(['process daatset ', ds{1}]);
        dataaddr = [baseaddr, ds{1}, '_cropped_face\', ds{1}, '\'];
        feature = setfield( feature, ds{1}, get_feature_wrapper(dataaddr) );
        is_front = setfield( is_front, ds{1}, check_if_front_face(dataaddr) );
    end
    save(savename, 'feature', 'is_front');
end

% ȥ�����������
feature = remove_non_front_face(feature, is_front);

%% 2�������ϳ�

% �ϲ��ü��������µ�����
feature_train = {feature.P1E, feature.P1L};
f_train_final = make_f_train_final(feature_train);

savename = 'cluster_centers-1000.mat';
i_cam = 1;
numClusters = 1000;
if exist(savename, 'file')
    load(savename);
else
    [ ~, cluster_centers ] = feature_merge( f_train_final, i_cam, numClusters );
    save(savename, 'cluster_centers');
end

[ feature_in_1_cam_comb, ~ ] = feature_merge( f_train_final, i_cam, numClusters, cluster_centers );
% �˵����е�nan�������������ǩ
[ feature_in_1_cam_comb, labels ] = filter_nan_and_get_label( feature_in_1_cam_comb );
% ÿ��seq�е�����
n_p_all = cellfun(@(x) size(x,1), feature_in_1_cam_comb,'un',1)';
% �趨Ҫʹ�õ�seq��Ŀ
n_seq = 3;
n_person_in_seq = n_p_all(1:n_seq); % ÿ��seq������
n_person_cum = [0, cumsum(n_person_in_seq)]; % �������ۼӣ��������ֿ飩
labels = labels(1:n_seq);

% ----------------------------------------------------------------------- %
% ��ѡ����seq����pca��ά
if 0
    im_mat = cell2mat(feature_in_1_cam_comb(1:n_seq));
    [X, mean_face, COEFF] = pca_eig( im_mat, 0.975 ); % imshow(reshape(mean_face,96,96))
    feature_in_1_cam_comb_pca = mat2cell(X', n_person_in_seq, size(X',2));
    fe_to_use = feature_in_1_cam_comb_pca(1:n_seq);
else
    fe_to_use = feature_in_1_cam_comb(1:n_seq);
end
    
%% ��ȡgt�����ӷ�ʽ
% ����������Ŀ���涨�����зַ�ʽ
N = 5;
nn = cell(1,n_seq);
nn_cum = cell(1,n_seq);
for ii=1:n_seq
    % ��ÿ��n_person_in_seq��ΪN��
    nn{ii} = [];
    cut = n_person_in_seq(ii);
    for i_N=1:N-1
        nn{ii} = [nn{ii}, floor(cut/N)];
    end
    nn{ii} = [nn{ii}, cut-sum(nn{ii})];
    nn_cum{ii} = [0, cumsum(nn{ii})];
end

% �������������Ӧ�ı���
N_p_in_seq = cell(N,1);
N_fe_to_use = cell(N,1);
N_labels = cell(N,1);
N_gt_connect = cell(N,1);
N_truth_num = cell(N,1);
for ii=1:N
    for jj=1:n_seq
        N_p_in_seq{ii}(jj) = nn{jj}(ii);
        fe_s = nn_cum{jj}(ii)+1;
        fe_e = nn_cum{jj}(ii+1);
        N_fe_to_use{ii}{jj,1} = fe_to_use{jj}(fe_s:fe_e,:);
        N_labels{ii}{jj,1} = labels{jj}(fe_s:fe_e,:);
    end
    [ N_gt_connect{ii}, N_truth_num{ii} ] = get_gt_match_result( N_labels{ii} );
end

%% �������������Լ����phi

N_vars_ADN_F = cell(N,1);
N_phi_x_z = cell(N,1);
N_phi_x_zstar = cell(N,1);
N_loss = cell(N,1);

for ii=1:N
    fprintf('�������� %d ��Լ��...\n', ii);
    N_vars_ADN_F{ii} = cal_GMMCP_constraints( N_p_in_seq{ii} );  
    [ N_phi_x_z{ii}, N_phi_x_zstar{ii}, N_loss{ii} ] = ...
    CXSL_cal_phi_and_loss( N_vars_ADN_F{ii}, N_fe_to_use{ii}, N_gt_connect{ii}, N_truth_num{ii} );
end

%% ��������ռ䣬�������
rng(0)
iter = 200;
w = zeros(size(N_phi_x_z{1}));
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

lambda = 1e-3;
gamma = zeros(iter,1); % ����gamma
time_consume = zeros(iter,1); % ��¼ÿ��ѭ�����õ�ʱ��
sample_loss = zeros(iter,N); % ��¼ÿһ����������ʧ������ֵ
Cd = 0;
linesearch = 0; % �Ƿ�ʹ������������
t = 0;
random = 0;
ind = 0;

%% ��ʼѭ��
while t<iter
    
    %% ���Ŀ�꺯��
    t = t + 1;
    tic;
    disp('  ==========================');
    disp(['  ��ʼ�� ', num2str(t), ' ��ѭ��...']);
    % ѡ�е�ind��������Ϊѵ������
    if random
        ind = randi(N);
    else
        ind = ind + 1;
        ind(ind==(N+1)) = 1;
    end
        
    disp(['  ѡ������ ', num2str(ind)]);
    % �����ٶ�������dot(w,phi)���𣿣�����˲�ʹ��������ʽ������ʱ����<w,feature>�������ΪĿ�꺯��
    [ pre_connect, ADN_value, phi_x_zhat, delta_zstar_zhat ] = ...
        CXSL_cal_GMMCP_object( Wavg{t}, N_vars_ADN_F{ind}, Cd, N_phi_x_z{ind}, N_loss{ind} );
    
    %% ����w����
    disp('      ����Ȩ���� w...'); 
    U_x_zstar_zhat = N_phi_x_zstar{ind} - phi_x_zhat;
    Ws = U_x_zstar_zhat/(lambda*N);
    ls = delta_zstar_zhat/N;
    sample_loss(t, ind) = delta_zstar_zhat;
    fprintf('      ��ǰ������ʧ������(z*,z^):\t%f\n', delta_zstar_zhat);
    
    % ����gap��gap��ֵ��������lambda����仯���޷�ȷ����������˻�����loss��gap�ȽϺ���
    if linesearch
        tmp = lambda*(Wi{t,ind}- Ws)'*W{t}- Li(t,ind)+ ls;
        gamma(t) = tmp/(lambda*norm(Wi{t,ind}- Ws)^2);
        gamma(t) = max([0, min([gamma(t),1])]);
    else
        gamma(t) = 2*N/(2*N + t-1);
    end

    %% ���� wi��Li�������º��w������ W{t+1,ind}��
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
aver_loss = sum(sample_loss,2);
plot(aver_loss, '-*');
fprintf('min loss = %f\n', min(aver_loss));
t_best = find(aver_loss==min(aver_loss));
t_best
w_best = Wavg{t_best};






