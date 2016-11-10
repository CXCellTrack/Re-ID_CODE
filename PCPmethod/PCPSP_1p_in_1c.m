function [weight, seq_mustlink, seq_cannotlink] = PCPSP_1p_in_1c( W, n_person_in_seq, labels, max_must_link, max_cannot_link )

% 来自论文 Pairwise Constraint Propagation by Semidefinite Programming
% feature 为 n*1 的形式
% n为场景的个数
% feature{i,1}中存放的是场景i下所有人各自的特征（融合结果）
% feature{i,1}{p}代表人物p的特征，若为空，说明该人在此处不存在
% 即p也包含了人物的标签信息

%% 1、求W矩阵（W是代表2个点之间的相似度，那么肯定值越大相似度越高，则W为距离的反比）公式24
n_seq = numel(n_person_in_seq); % 序列的个数

% 每个序列中的人数
n_all = sum(n_person_in_seq);

%% 2、生成L矩阵（normalized graph Laplacian）
D = zeros(n_all);
for i=1:n_all
    D(i,i) = sum(W(i,:));
end
L = D - W;
L_ = D^(-0.5)*L*D^(-0.5); % L_ 公式（8） 
% L_是对称半正定的，特征值在0-2之间（可以用eig(L_)检查一下）

%% 设定must-link和cannot-link约束
% 人数累加
n_person_cum = cumsum(n_person_in_seq);
n_person_cum = [0, n_person_cum]; % 补0方便计算第一个

% 规定每2组中must_link的个数（增多则更加精确）
% max_must_link = 4;
% max_cannot_link = 4;
rng(0)

MustLink = [];
CannotLink = [];

for ii=1:n_seq-1
    for jj=ii+1:n_seq
        n_must_link = 0;
        n_cannot_link = 0;
        % -------- 随机生成一些index给MustLink（保证label相同）
        index = randsample(n_person_in_seq(ii), max_must_link)';
        for ii_ind=index 
            final_i = n_person_cum(ii) + ii_ind;
            % 在ii中找到这些index对应的标签，看看jj中是否有相同标签的，将其设置为mustlink
            label_ii = labels{ii}(ii_ind);
            if any(labels{jj}==label_ii) && n_must_link<max_must_link
                jj_ind = find(labels{jj}==label_ii);
                % 给出MustLink中的最终索引
                final_j = n_person_cum(jj) + jj_ind;
                MustLink = [MustLink; [final_i, final_j]];
                n_must_link = n_must_link + 1;
            end
        end
        
        % -------- 随机生成一些index给CannotLink（保证label相同）
        index = randsample(n_person_in_seq(ii), max_cannot_link)';
        for ii_ind=index 
            label_ii = labels{ii}(ii_ind);
            final_i = n_person_cum(ii) + ii_ind;
            not_equal_ind = find(labels{jj}~=label_ii);
            jj_ind = not_equal_ind( randsample(numel(not_equal_ind),1) );
            % 给出MustLink中的最终索引
            final_j = n_person_cum(jj) + jj_ind;
            if n_cannot_link<max_cannot_link
                CannotLink = [CannotLink; [final_i, final_j]];
                n_cannot_link = n_cannot_link + 1;
            end
        end
    end
end

% 检查mustlink和cannotlink是否矛盾
% 一旦出现矛盾，则问题变为不可解问题
seq_mustlink = cell(n_seq);
for ii=1:size(MustLink,1)
    p1 = MustLink(ii,1);
    p2 = MustLink(ii,2);
    seq1 = find(p1-n_person_cum>0, 1, 'last');
    seq2 = find(p2-n_person_cum>0, 1, 'last');
    sp1 = p1-n_person_cum(seq1);
    sp2 = p2-n_person_cum(seq2);
    % tmp_mustlink为区分seq的mustlink方案，便于检查
    seq_mustlink{seq1,seq2} = [seq_mustlink{seq1,seq2}; sp1,sp2];
    if labels{seq1}(sp1)~=labels{seq2}(sp2)
        error('MustLink的第%d行错误：%d和%d不是mustlink关系！', ii,p1,p2);
    end
end
seq_cannotlink = cell(n_seq);
for ii=1:size(CannotLink,1)
    p1 = CannotLink(ii,1);
    p2 = CannotLink(ii,2);
    seq1 = find(p1-n_person_cum>0, 1, 'last');
    seq2 = find(p2-n_person_cum>0, 1, 'last');
    sp1 = p1-n_person_cum(seq1);
    sp2 = p2-n_person_cum(seq2);
    % tmp_mustlink为区分seq的mustlink方案，便于检查
    seq_cannotlink{seq1,seq2} = [seq_cannotlink{seq1,seq2}; sp1,sp2];
    if labels{seq1}(sp1)==labels{seq2}(sp2)
        error('CannotLink的第%d行错误：%d和%d不是cannotlink关系！', ii,p1,p2);
    end
end

%% 定义半正定规划的约束
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
gap = 1e-6; % 使用松弛的=1和=0，不然会出现 numerical problem
exact_euqal = 1;
for i=1:n_all
    if exact_euqal
        E = zeros(n_all);
        E(i,i) = 1;
        F1 = [ F1, sum(sum(E.*K))==1 ]; % 自身为1
    else
        F1 = [ F1, abs( K(i,i) -1 )<=gap]; % 自身为1
    end
end
for ii=1:size(MustLink,1) % MustLink为1
    if exact_euqal
        E = zeros(n_all);
        E(MustLink(ii,1), MustLink(ii,2)) = 1;
        F2 = [ F2, sum(sum(E.*K))==1 ];
    else
        F2 = [ F2, abs( K(MustLink(ii,1), MustLink(ii,2)) -1 )<=gap ];
    end
end
for ii=1:size(CannotLink,1) % CannotLink为0
    if exact_euqal
        E = zeros(n_all);
        E(CannotLink(ii,1), CannotLink(ii,2)) = 1;
        F3 = [ F3, sum(sum(E.*K))==0 ];
    else
        F3 = [ F3, abs( K(CannotLink(ii,1), CannotLink(ii,2)) )<=gap ];
    end
end
% 额外添加同一seq中各项必须cannotlink
for ii=1:n_seq
    for p1 = n_person_cum(ii)+1:n_person_cum(ii+1)-1
        for p2 = p1+1:n_person_cum(ii+1)
            F3 = [F3, K(p1,p2)==0];
        end
    end
end
         
toc;
F = [ K>=0, F1, F2, F3 ]; 

%% 求解半正定规划

% 开始求解
disp('开始求解...');
options = sdpsettings('verbose', 0, 'debug',1, 'solver','mosek'); % 商业求解器 mosek 比 sedumi 快很多
sol = solvesdp( F, OBJ, options ) % checkset(F)
Kval = value(K);


% 将kval转换后输出（按seq分开，便于观察）
weight = cell(n_seq);
for i1=1:n_seq-1
    for i2=i1+1:n_seq
        i1_s = n_person_cum(i1)+1;
        i1_e = n_person_cum(i1+1);
        i2_s = n_person_cum(i2)+1;
        i2_e = n_person_cum(i2+1);
        % 将Kval的值作为权值分配下来
        weight{i1,i2} = Kval(i1_s:i1_e, i2_s:i2_e);
    end
end










    
    
    
    