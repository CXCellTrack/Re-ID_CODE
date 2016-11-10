function [ m_hat, X_in_mat ] = estimate_m_by_inputX( X_in, n_person_in_seq )

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