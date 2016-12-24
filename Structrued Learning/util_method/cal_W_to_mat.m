function [ W ] = cal_W_to_mat( fe_to_use, n_person_in_seq )


%（W是代表2个点之间的相似度，那么肯定值越大相似度越高，则W为距离的反比）公式24

n_all = sum(n_person_in_seq);
% 合成整个 feature_all 矩阵，并求出距离
feature_all = cell2mat(fe_to_use);
distance = squareform( pdist(feature_all, 'euclidean') ); % 求||xi-xj||

% 求r，r代表所有节点到其最近的第10个节点之间的平均距离
n_tenth = 10; % 这里因为点比较少，因此取4
tenth = zeros(n_all,1);
for i=1:n_all
    [Y,~]= sort(distance(i,:));
    tenth(i) = Y(n_tenth+1);
end
r = mean(tenth);

% S(0.1r, r, 5) 并 S(r, 10r, 5)
S_r5 = linspace(0.1*r,r,5);
S_10r5 = linspace(r,10*r,5);
sigma_range = [S_r5(1:4), S_10r5]; % 求出sigma的范围

sigma = sigma_range(5); % 先取一个作为sigma的值
W = zeros(n_all);
for i=1:n_all-1
    for j=i+1:n_all
        W(i,j) = exp(-(distance(i,j)^2/(2*sigma^2))); % W的公式可以自己设计，满足w越大，2者越相似即可
        W(j,i) = W(i,j); % 利用对称生成W
    end
end

end

