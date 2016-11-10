function [ m_hat, X_in_mat ] = estimate_m_by_inputX( X_in, n_person_in_seq )

%% 将输入cell X_in变为矩阵形式
n_seq = size(X_in,1);

for ii=1:n_seq
    X_in{ii,ii} = eye(n_person_in_seq(ii));
    if ii==n_seq
        break
    end
    % 把下半区对称补上
    for jj=ii+1:n_seq
        X_in{jj,ii} = X_in{ii,jj}';
    end
end
X_in_mat = cell2mat(X_in); % X_in的矩阵形式

%% 算法1：估计m的值（在这里表现为估计people的个数）
lambda = sort(eig(X_in_mat), 'descend');
mi = n_person_in_seq; % mi=|Si|，即每个对象中的元素数目（此处为每个场景中的目标数）
N = sum(mi); % 见文章公式（2）下方
M = max(2, max(mi));
% 依照算法1求出m_hat
m_hat = 0;
max_sub = 0;
for i=M:N-1
    tmp = lambda(i) - lambda(i+1);
    if tmp>max_sub
        max_sub = tmp;
        m_hat = i; % 注意是argmax，因此最后求的是m_hat，而非特征值差的最大值
    end
end