clear
close all

%% 制造数据
n1 = 5;
n2 = 10;
n3 = 40;
r1 = 10;
r2 = 30;
[data_n1, data_n2, data_n3] = make_data(n1, n2, n3, r1, r2, 1);

%% 1、求W矩阵（W是代表2个点之间的相似度，那么肯定值越大相似度越高，则W为距离的反比）公式24
n_all = n1 + n2 + n3;
data_all = [data_n1; data_n2; data_n3];
distance = squareform( pdist(data_all, 'euclidean') ); % 求||xi-xj||

% 求r，r代表所有节点到其最近的第10个节点之间的平均距离
n_tenth = 4; % 这里因为点比较少，因此取4
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
        W(i,j) = exp(-(distance(i,j)^2/(2*sigma^2)));
        W(j,i) = W(i,j); % 利用对称生成W
    end
end

%% 2、生成L矩阵（normalized graph Laplacian）
D = zeros(n_all);
for i=1:n_all
    D(i,i) = sum(W(i,:));
end
L = D - W;
L_ = D^(-0.5)*L*D^(-0.5); % L_ 公式（8） 
% L_是对称半正定的，特征值在0-2之间（可以用eig(L_)检查一下）

%% 求解半正定规划
label = [ones(n1,1); ones(n2,1)*2; ones(n3,1)*3]; % 没有用到label
% 对3种进行取样
n_sample = 5;
sample1 = randi(n1, n_sample,1);
sample2 = randi(n2, n_sample,1)+n1;
sample3 = randi(n3, n_sample,1)+n2+n1;
% 进行must-link和cannot-link连接
MustLink = [sample1, sample3]; % MustLink = []; % 没有mustlink和cannotlink时，解得K所有都为1
CannotLink = [sample1, sample2]; % CannotLink = [];
%
% 绘制连线
for i=1:size(MustLink,1)
    p1 = data_all(MustLink(i,1),:);
    p2 = data_all(MustLink(i,2),:);
    line([p1(1),p2(1)], [p1(2),p2(2)], 'linestyle', '-','color','r');
end
for i=1:size(CannotLink,1)
    p1 = data_all(CannotLink(i,1),:);
    p2 = data_all(CannotLink(i,2),:);
    line([p1(1),p2(1)], [p1(2),p2(2)], 'linestyle', '--', 'color','k');
end
hold off

% ----------------------------------- %
% 进行变量分配、求解
K = sdpvar(n_all);
% kij =<phi(xi),phi(xj)>F. Then the matrix K =[kij ]n×n
%  is symmetric and positive semidefinite, denoted
% by K>=0：K为对称半正定，因此可作为一个核函数
OBJ = sum(sum(L_.*K));

% 构造3条约束
F1 = [];
F2 = [];
F3 = [];

disp('构造约束条件...');
tic;
for i=1:n_all
%     E = zeros(n_all);
%     E(i,i) = 1;
%     F1 = [ F1, sum(sum(E.*K))==1 ]; % 自身为1
    F1 = [ F1, abs( K(i,i)-1 )<=1e-6 ]; % 自身为1
end
for ii=1:size(MustLink,1) % MustLink为1
%     E = zeros(n_all);
%     E(MustLink(ii,1), MustLink(ii,2)) = 1;
%     F2 = [ F2, sum(sum(E.*K))==1 ];
    F2 = [ F2, abs( K(MustLink(ii,1), MustLink(ii,2))-1 )<=1e-6 ];
end
for ii=1:size(CannotLink,1) % CannotLink为0
%     E = zeros(n_all);
%     E(CannotLink(ii,1), CannotLink(ii,2)) = 1;
%     F3 = [ F3, sum(sum(E.*K))==0 ];
    F3 = [ F3, abs( K(CannotLink(ii,1), CannotLink(ii,2)) )<=1e-6 ];
end

F = [ K>=0, F1, F2, F3 ]; 
toc;
% 开始求解
disp('开始求解...');
options = sdpsettings('verbose', 0, 'debug',1, 'solver','SEDUMI');
sol = solvesdp( F, OBJ, options ) % checkset(F)
Kval = value(K);













    
    
    
    