clear, close all
 
% 实现下面这篇文章的demo
% GMMCP Tracker: Globally Optimal Generalized Maximum Multi Clique Problem for Multiple Object Tracking

%% 1、生成数据
% 仿照图1，生成4段数据，以平面上的点来表示
rng(0)
n = 35; r = 10;
theta = linspace(0,360,n);
data = [r*cosd(theta)', r*sind(theta)']; % 在一个圆周上生成一些数据
noisy = (rand(n,2)*2-1)*2;
data = data + noisy;

figure
plot(data(:,1),data(:,2),'bo');hold on;
line([-15,15],[0, 0],'color','r','linewidth',2);  line([0,0],[-15,15],'color','r','linewidth',2);

% 将数据data 分入4个cluster中
cluster = cell(4,1);
for h=1:size(data,1)
    if data(h,1)>0 && data(h,2)>0 % 第1象限
        cluster{1} = [cluster{1}; data(h,:)];
    elseif data(h,1)<0 && data(h,2)>0 % 第2象限
        cluster{2} = [cluster{2}; data(h,:)];
    elseif data(h,1)<0 && data(h,2)<0 % 第3象限
        cluster{3} = [cluster{3}; data(h,:)];       
    elseif data(h,1)>0 && data(h,2)<0 % 第4象限
        cluster{4} = [cluster{4}; data(h,:)];
    end
end

%% 2、加入哑节点
n_nodes = zeros(size(cluster)); % 每个cluster中的节点数目
for i=1:numel(cluster)
    n_nodes(i) = size(cluster{i},1);
end
% 每个cluster需要加入的哑节点数目的上界 Nd_ub
Nd_ub = zeros(size(n_nodes));
for i=1:numel(Nd_ub)
    Nd_ub(i) = sum(n_nodes(1:i)) + sum(n_nodes(i:end));
end
% 按照作者的说法，取1/3大小的上界就够用了
% K = round(Nd_ub(1)/10) + n_nodes(1);
K = max(n_nodes);
% 将哑节点加入cluster中，为了便于计算权重，设其值为nan
for i=1:numel(cluster)
    cluster{i} = [cluster{i}; nan*zeros(K-n_nodes(i),2)];
end

%% 3、计算边的权值
weight = cell(4,4);
for i1=1:3
    for i2=i1+1:4
        weight{i1,i2} = zeros(K,K); 
    end
end
% weight{i1,i2}(j1,j2)为第i1个cluster中第j1个节点到第i2个cluster中第j2个节点权重

for i1=1:3 % weight 为上三角阵
    for i2=i1+1:4
        for j1=1:K
            for j2=1:K
                if isnan(cluster{i1}(j1,1)) || isnan(cluster{i2}(j2,1))
                    weight{i1,i2}(j1,j2) = 0; % 如果有一个为哑节点，则wegith为0
                else % 否则weigh为距离的倒数，越近weight越大
                    weight{i1,i2}(j1,j2) = 1/norm(cluster{i1}(j1,:)-cluster{i2}(j2,:));
                end
            end
        end
        weight{i2,i1} = weight{i1,i2}'; % 对称
    end
end

%% 4、计算约束条件
vars = cell(4,4);
for i1=1:3
    for i2=i1+1:4 % vars指的是边的变量（0、1代表是否连接上），论文中还有节点变量Vij，似乎不参与运作
        vars{i1,i2} = binvar(K, K, 'full'); 
        vars{i2,i1} = vars{i1,i2}'; % 镜像复制
    end
end % 分配变量，每个变量对应一个权值

% 计算约束条件1（每个cluster中节点的个数必须为k，这样才能最终形成k个clique，自动成立！）
% 计算约束条件2（任意一个节点的有效边为h-1个，h为cluster个数）
F2 = [];
for i1=1:3 % weight 为上三角阵
    for i2=i1+1:4
        for j1=1:K
            F2 = [F2, sum(vars{i1,i2}(j1,:))==1, sum(vars{i1,i2}(:,j1))==1]; % i1中的j1必须连接到i2中的某一个
        end
    end
end

% 计算约束条件3（点a与点b连，点b与点c连，则a要与c连）
F3 = [];
tic
for i1=1:4 % 非常非常非常慢 eij的总数目为 h*(h-1)/2*K^2 在这里为 6*14^2 = 1176
    for i2=i1+1:4 % 每3个变量之间就有3个约束3，因此约束3条数为: C_h_2*(h-2)*K^3 = 6*2*14^3 = 32928
        % 因为是 i1和i2的和，因此2个没有必须都从1到4（否则会使约束数目翻倍，而结果并不变）
        for i3=1:4 % i1和i2交换位置，其实是同一条约束（因此不妨规定i2>i1）
            % 三个cluter必须不相同
            if i1==i2 || i2==i3 || i3==i1
                continue;
            end
            fprintf('计算 %d — %d — %d 的约束...\n', i1,i2,i3);
            for j1=1:K
                for j2=1:K
                    for j3=1:K
                         % ============================================== %
%                         fprintf('计算 %d：%d — %d：%d — %d：%d 的约束\n', i1,j1,i2,j2,i3,j3);
%                         % 只接受后面大前面小的组合（因为此前例如var{3,1}为空，必须先转换为var{1,3}）
%                         i1a = i1; i1b = i1; i1c = i1; j1a = j1; j1b = j1; j1c = j1; 
%                         i2a = i2; i2b = i2; i2c = i2; j2a = j2; j2b = j2; j2c = j2; 
%                         i3a = i3; i3b = i3; i3c = i3; j3a = j3; j3b = j3; j3c = j3; 
%                         if i1>i2 
%                             i1a = i2; i2a = i1; j1a = j2; j2a = j1;
%                         end
%                         if i2>i3 
%                             i2b = i3; i3b = i2; j2b = j3; j3b = j2;
%                         end
%                         if i1>i3 
%                             i1c = i3; i3c = i1; j1c = j3; j3c = j1;
%                         end
%                         
%                         F3 = [F3, vars{i1a,i2a}(j1a,j2a)+vars{i2b,i3b}(j2b,j3b)<=1+vars{i1c,i3c}(j1c,j3c)];
                        % =============================================== %
                        % 2016.4.15 将var填充为满矩阵后，就可以直接计算了
                        F3 = [F3, vars{i1,i2}(j1,j2)+vars{i2,i3}(j2,j3)<=1+vars{i1,i3}(j1,j3)];
                    end
                end
            end
        end
    end
end
tic

%% 5、计算目标函数并进行求解
F = [F2, F3];
OBJ = 0;
for i1=1:3
    for i2=i1+1:4
        OBJ = OBJ + sum(sum( vars{i1,i2}.*weight{i1,i2} ));
    end
end

% 求解BILP
options = sdpsettings('verbose',0,'solver','gurobi');
sol = solvesdp( F, -OBJ, options )
for i1=1:3
    for i2=i1+1:4
        vars{i1,i2} = value(vars{i1,i2});   
    end
end

%% 6、画出连接图
cluster_dummy = cluster; % 将dummy转变为顶角坐标
for i=1:4
    switch i % 将nan投射到顶角上
        case 1
            xy = [10 10];
        case 2
            xy = [-10 10];
        case 3
            xy = [-10 -10];
        case 4
            xy = [10 -10];
    end
    plot(xy(1),xy(2), 'rs');
    n_nan = sum(sum(isnan(cluster{i}(:,1))));
    text(xy(1)+1, xy(2)+1, ['dummy nodes: ', num2str(n_nan)]);
    cluster_dummy{i}(isnan(cluster{i}(:,1)),1) = xy(1);
    cluster_dummy{i}(isnan(cluster{i}(:,2)),2) = xy(2);
end

color = colormap(hsv); % 选择颜色
color = color(randperm(size(color,1)),:);
color_cluster = cell(4,1);
color_cluster{1} = color(1:K,:);

for i1=1:3
    for i2=i1+1:4
        if abs(i1-i2)==2
            continue; % 不画斜线
        end
        for j1=1:K
            j2 = find(vars{i1,i2}(j1,:));
            if isnan(cluster{i1}(j1,1)) || isnan(cluster{i2}(j2,1))
                color_cluster{i1}(j1,:) = [0 0 0];
            end
            tc = color_cluster{i1}(j1,:);
            color_cluster{i2}(j2,:) = tc;
            line([cluster_dummy{i1}(j1,1), cluster_dummy{i2}(j2,1)],...
                [cluster_dummy{i1}(j1,2), cluster_dummy{i2}(j2,2)],'color',tc);
        end
    end
end
hold off;




