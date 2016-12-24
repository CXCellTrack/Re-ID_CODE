function vars_ADN_F = cal_GMMCP_constraints( n_person_in_seq )

% 实现下面这篇文章
% GMMCP Tracker: Globally Optimal Generalized Maximum Multi Clique Problem for Multiple Object Tracking

%% 1、加入聚合哑节点（ADN）（在后面直接intvar即可）
global each_person_a_clique;
if each_person_a_clique
    % 是否把每个人都当作一个clique（这样就无需考虑每个clique中对象不能是一个人这个条件）
    n_nodes = ones(1, sum(n_person_in_seq));
else
    n_nodes = n_person_in_seq; % 每个cluster中的节点数目
end

% 每个cluster中加入Aggregated Dummy Nodes (ADN)
% 其值为整数，这样减少了变量的数目
% ADN的加入在建立整数规划的时候再体现!

h_c = numel(n_nodes); % h_c 为cluster的个数

%% 2、计算约束条件
vars = cell(h_c);
for i1=1:h_c-1 % 分配变量，每个变量对应一个权值
    for i2=i1+1:h_c % vars指的是边的变量（0、1代表是否连接上），论文中还有节点变量Vij，似乎不参与运作
        vars{i1,i2} = binvar(n_nodes(i1), n_nodes(i2), 'full'); 
        vars{i2,i1} = vars{i1,i2}'; % 镜像复制
    end
end 
% 分配 ADN 整型变量
ADN = intvar(1,h_c,'full');

% 计算约束条件1（对于每个cluster i，进来的边的个数+ADN(i)的个数 = (h-1)*K），
% K为希望生成的clique的个数)
F1 = [];
for i=1:h_c
    sum_in = 0;
    for j=1:h_c
        sum_in = sum_in + sum(sum(vars{i,j}));
    end
    F1 = [F1, sum_in + ADN(i) == (h_c-1)*(10+max(n_nodes))]; % 估计K值为 10+max(n_nodes)
end
        
% 计算约束条件2（任意一个节点的有效边为h-1个，h为cluster个数）在ADN版本中，改为<=1
F2 = [];
for i1=1:h_c-1 % weight 为上三角阵
    for i2=i1+1:h_c
        for j1=1:n_nodes(i1)
            F2 = [F2, sum(vars{i1,i2}(j1,:))<=1]; % i1中的j1至多连接到i2中的某一个
        end
        for j2=1:n_nodes(i2)
            F2 = [F2, sum(vars{i1,i2}(:,j2))<=1]; % i2中的j2至多连接到i1中的某一个
        end
    end
end

% 计算约束条件3（点a与点b连，点b与点c连，则a要与c连）
F3 = [];
tic
for i1=1:h_c-1 % 非常非常非常慢 eij的总数目为 h*(h-1)/2*K^2 在这里为 6*14^2 = 1176
    for i2=i1+1:h_c % 每3个变量之间就有3个约束3，因此约束3条数为: Ch_2*(h-2)*K^3 = 6*2*14^3 = 32928
        for i3=1:h_c % i1和i2交换位置，其实是同一条约束（因此不妨规定i2>i1）
            % 三个cluter必须不相同
            if i1==i2 || i2==i3 || i3==i1
                continue;
            end
            fprintf('计算 {%d} + {%d} <= 1 + {%d} 的约束...\n', i1,i2,i3);
            for j1=1:n_nodes(i1)
                for j2=1:n_nodes(i2)
                    for j3=1:n_nodes(i3)
                    	F3 = [F3, vars{i1,i2}(j1,j2)+vars{i2,i3}(j2,j3)<=1+vars{i1,i3}(j1,j3)];
                    end
                end
            end
        end
    end
end
toc

F = [F1, F2, F3];

vars_ADN_F = struct();
vars_ADN_F.vars = vars;
vars_ADN_F.ADN = ADN;
vars_ADN_F.F = F;


end

