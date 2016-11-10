function [ X_in_mat, X_out_val_int ] = MatchLift( X_in, n_person_in_seq, lambda_l1 )

%% ������cell X_in��Ϊ������ʽ
n_seq = size(X_in,1);

for ii=1:n_seq
    X_in{ii,ii} = eye(n_person_in_seq(ii));
    if ii==n_seq
        break
    end
    % ���°����ԳƲ���
    for jj=ii+1:n_seq
        X_in{jj,ii} = X_in{ii,jj}';
    end
end
X_in_mat = cell2mat(X_in); % X_in�ľ�����ʽ

%% �㷨1������m��ֵ�����������Ϊ����people�ĸ�����
lambda = sort(eig(X_in_mat), 'descend');
mi = n_person_in_seq; % mi=|Si|����ÿ�������е�Ԫ����Ŀ���˴�Ϊÿ�������е�Ŀ������
N = sum(mi); % �����¹�ʽ��2���·�
M = max(2, max(mi));
% �����㷨1���m_hat
m_hat = 0;
max_sub = 0;
for i=M:N-1
    tmp = lambda(i) - lambda(i+1);
    if tmp>max_sub
        max_sub = tmp;
        m_hat = i; % ע����argmax�������������m_hat����������ֵ������ֵ
    end
end

% mapͼ�У��ߵĸ�����ע���������2������ıߣ��˴�Ϊ2������������
% ���½�2.1
E = n_seq*(n_seq-1)/2; % (Si,Sj)�ʦţ����ظ�����
if ~exist('lambda_l1', 'var')
    lambda_l1 = sqrt(E)/2/n_seq; % �������lambda_l1
end

%% ����Ŀ�꺯��
disp('����Լ������...');
X_out = sdpvar(size(X_in_mat,1));
obj = trace(X_in_mat'*X_out) - lambda_l1*sum(X_out(:)); % ����ΪL1����X_out��1�ĸ���

n_person_cum = [0, cumsum(n_person_in_seq)];
% ��ѡԼ����X��ÿ��ÿ�� 0<=��<=1
% ע�����Լ��������Ϊ0-1����ʱ����ʡȥ���Լӿ�����ٶ�
% ��������ΪȨֵʱ����Ҫ���ϣ����������Լ�����ģ���
F0 = [];
for ii=1:n_seq-1
    for jj=ii+1:n_seq
        ii_s = n_person_cum(ii)+1;
        ii_e = n_person_cum(ii+1);
        jj_s = n_person_cum(jj)+1;
        jj_e = n_person_cum(jj+1);
        for ih=ii_s:ii_e % �ֿ��к�Լ��
            F0 = [ F0, sum(X_out(ih, jj_s:jj_e))<=1, sum(X_out(ih, jj_s:jj_e))>=0 ];
        end
        for iw=jj_s:jj_e % �ֿ��к�Լ��
            F0 = [ F0, sum(X_out(ii_s:ii_e, iw))<=1, sum(X_out(ii_s:ii_e, iw))>=0 ];
        end
    end
end

% Լ��1��X_out�ĶԽǷֿ�Ϊ��λ��
F1 = [];
for ii=1:size(X_out,1)
    F1 = [F1, X_out(ii,ii)==1];
end

% Լ��2��X>=0������X�Գƣ�ֻ��Ū1�뼴��
F2 = [];
for i=1:size(X_out,1)
    for j=i:size(X_out,1)
        F2 = [F2, X_out(i,j)>=0];
    end
end

% Լ��3��������
X_out_ex = [m_hat, ones(1,size(X_out,1));
            ones(size(X_out,1),1), X_out];
F3 = [ X_out_ex>=0 ];

%% �������
disp('��ʼ���...');
F = [ F0, F1, F2, F3 ];
options = sdpsettings('verbose', 0, 'solver','mosek'); % ��ҵ����� mosek �� sedumi ��ܶ�
sol = solvesdp( F, -obj, options ) % checkset(F)
X_out_val = value(X_out);

%% round���ԣ��㷨2
% [Vec, D] = eig(X_out_val);
% [Y, I] = sort(diag(D), 'descend');
% r = m_hat; % rΪ�����е���Ŀ�Ĺ��ƣ���ͬ��m����
% SIGMA = diag(Y(1:r));
% U = Vec(:,I(1:r));
% V = U*SIGMA^(1/2);
% 
% % ���и���
% n_row = size(X_out_val,1);
% eye_mat = eye(size(V,2));
% 
% % ʣ���vj��������2��
% j_remain = 1:n_row;
% i = 1;
% 
%     vi = V(i,:)'; % V�ĵ�i�У�has not been fixed.��
%     ei = eye_mat(:,i); % ��λ��
%     % ---- �����u����O ---- %
%     orth_vecs = null(vi'); % ��vi����������ÿ����һ����������
%     O = [orth_vecs'; vi'];
%     O(end,:) = O(i,:); O(i,:) = vi'; % ����vi����Ӧλ��
% 	% O*vi
%     % ---------------------- %
%     V = V*O';
%     % ��ʣ�µ�ÿһ��vj���в�����
%     dotvv = [];
%     for j=i+1:n_row
%         sj = find(j<=n_person_cum, 1,'first')-1;
%         vj = V(j,:)';
%         dotvv(sj,j) = dot(vj,vi);
%     end
%     index = [];
%     for nn=1:n_seq
%         [maxvalue, ind] = max(dotvv(nn,:));
%         if maxvalue>0.5
%             index(nn,1) = ind;
%         end
%     end
%     % indexΪ�գ���˵��û��<vi,v1>����0.5����δ���
%     V_remain = V(j_remain,:);
    
    
    
% ���������round���
X_out_val_int = round(X_out_val); 

%% ���㾫��
% n_bingo1 = trace(X_in_mat'*X_out_val_int);
% n_in1 = sum(sum(X_in_mat));
% fprintf('��ԭ����� %d ��1���������� %d ��1�� ������ %f\n', n_bingo1, n_in1, n_bingo1/n_in1);

end

