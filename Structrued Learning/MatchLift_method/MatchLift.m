function [ X_in_mat, X_out_val_int ] = MatchLift( X_in, n_person_in_seq, lambda_l1 )

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

% map图中，边的个数（注意这个边是2个对象的边）此处为2个场景的连接
% 见章节2.1
E = n_seq*(n_seq-1)/2; % (Si,Sj)∈ε，不重复计算
if ~exist('lambda_l1', 'var')
    lambda_l1 = sqrt(E)/2/n_seq; % 正则参数lambda_l1
end

%% 构造目标函数
disp('构造约束条件...');
X_out = sdpvar(size(X_in_mat,1));
obj = trace(X_in_mat'*X_out) - lambda_l1*sum(X_out(:)); % 正则为L1，即X_out中1的个数

n_person_cum = [0, cumsum(n_person_in_seq)];
% 可选约束：X的每行每列 0<=和<=1
% 注意这个约束在输入为0-1矩阵时可以省去，以加快求解速度
% 但在输入为权值时，需要加上（这种是我自己试验的！）
F0 = [];
for ii=1:n_seq-1
    for jj=ii+1:n_seq
        ii_s = n_person_cum(ii)+1;
        ii_e = n_person_cum(ii+1);
        jj_s = n_person_cum(jj)+1;
        jj_e = n_person_cum(jj+1);
        for ih=ii_s:ii_e % 分块行和约束
            F0 = [ F0, sum(X_out(ih, jj_s:jj_e))<=1, sum(X_out(ih, jj_s:jj_e))>=0 ];
        end
        for iw=jj_s:jj_e % 分块列和约束
            F0 = [ F0, sum(X_out(ii_s:ii_e, iw))<=1, sum(X_out(ii_s:ii_e, iw))>=0 ];
        end
    end
end

% 约束1：X_out的对角分块为单位阵
F1 = [];
for ii=1:size(X_out,1)
    F1 = [F1, X_out(ii,ii)==1];
end

% 约束2：X>=0，由于X对称，只需弄1半即可
F2 = [];
for i=1:size(X_out,1)
    for j=i:size(X_out,1)
        F2 = [F2, X_out(i,j)>=0];
    end
end

% 约束3：半正定
X_out_ex = [m_hat, ones(1,size(X_out,1));
            ones(size(X_out,1),1), X_out];
F3 = [ X_out_ex>=0 ];

%% 进行求解
disp('开始求解...');
F = [ F0, F1, F2, F3 ];
options = sdpsettings('verbose', 0, 'solver','mosek'); % 商业求解器 mosek 比 sedumi 快很多
sol = solvesdp( F, -obj, options ) % checkset(F)
X_out_val = value(X_out);

%% round策略：算法2
% [Vec, D] = eig(X_out_val);
% [Y, I] = sort(diag(D), 'descend');
% r = m_hat; % r为对所有点数目的估计（等同于m？）
% SIGMA = diag(Y(1:r));
% U = Vec(:,I(1:r));
% V = U*SIGMA^(1/2);
% 
% % 逐行更新
% n_row = size(X_out_val,1);
% eye_mat = eye(size(V,2));
% 
% % 剩余的vj（见步骤2）
% j_remain = 1:n_row;
% i = 1;
% 
%     vi = V(i,:)'; % V的第i行（has not been fixed.）
%     ei = eye_mat(:,i); % 单位阵
%     % ---- 构造出u矩阵O ---- %
%     orth_vecs = null(vi'); % 求vi的正交基（每列是一个正交基）
%     O = [orth_vecs'; vi'];
%     O(end,:) = O(i,:); O(i,:) = vi'; % 交换vi到对应位置
% 	% O*vi
%     % ---------------------- %
%     V = V*O';
%     % 对剩下的每一行vj进行操作：
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
%     % index为空，，说明没有<vi,v1>超过0.5，如何处理？
%     V_remain = V(j_remain,:);
    
    
    
% 先用这个简单round替代
X_out_val_int = round(X_out_val); 

%% 计算精度
% n_bingo1 = trace(X_in_mat'*X_out_val_int);
% n_in1 = sum(sum(X_in_mat));
% fprintf('复原结果有 %d 个1，输入结果有 %d 个1， 命中率 %f\n', n_bingo1, n_in1, n_bingo1/n_in1);

end

