clear, %close all
 
% 实现下面这篇文章的demo
% GMMCP Tracker: Globally Optimal Generalized Maximum Multi Clique Problem for Multiple Object Tracking

%% 1、生成数据
% 仿照图1，生成4段数据，以平面上的点来表示
rng(0)
n = 40; r = 10;
theta = linspace(0,360,n);
data = [r*cosd(theta)', r*sind(theta)']; % 在一个圆周上生成一些数据
noisy = (rand(n,2)*2-1)*2;
data = data + noisy;

figure
plot(data(:,1),data(:,2),'bo');hold on;
line([-15,15],[0, 0],'color','r','linewidth',2);  line([0,0],[-15,15],'color','r','linewidth',2);

% 将数据data 分入4个cluster中
h_c = 3; % h_c 为cluster的个数
cluster = cell(h_c,1);
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

%% 2、加入聚合哑节点（ADN）
n_nodes = zeros(size(cluster)); % 每个cluster中的节点数目
for i=1:numel(cluster)
    n_nodes(i) = size(cluster{i},1);
end

% 每个cluster中加入Aggregated Dummy Nodes (ADN)
% 其值为整数，这样减少了变量的数目
% ADN的加入在建立整数规划的时候再体现!

%% 3、计算边的权值

weight = cell(h_c,h_c);
for i1=1:h_c-1
    for i2=i1+1:h_c
        weight{i1,i2} = zeros(size(cluster{i1},1), size(cluster{i2},1)); 
    end
end
% weight{i1,i2}(j1,j2)为第i1个cluster中第j1个节点到第i2个cluster中第j2个节点权重

for i1=1:h_c-1 % weight 为上三角阵
    for i2=i1+1:h_c
        for j1=1:size(cluster{i1},1)
            for j2=1:size(cluster{i2},1)
                if isnan(cluster{i1}(j1,1)) || isnan(cluster{i2}(j2,1))
                    weight{i1,i2}(j1,j2) = 0; % 如果有一个为哑节点，则wegith为0
                else % 否则weigh为距离的倒数，越近weight越大
                    weight{i1,i2}(j1,j2) = 1/norm(cluster{i1}(j1,:)-cluster{i2}(j2,:));
%                     weight{i2,i1}(j2,j1) = weight{i1,i2}(j1,j2); % 对称
                end
            end
        end
    end
end

% 注意在非speedup版本中，连接到dummy node的权值为Cd（将其设为0）
% 在speedup版本中，ADN不是边变量，而是节点变量，节点权值为Cd/2（不妨设为0.01）
Cd = 0.00;

%% 4、计算约束条件
vars = cell(h_c,h_c);
for i1=1:h_c-1 % 分配变量，每个变量对应一个权值
    for i2=i1+1:h_c % vars指的是边的变量（0、1代表是否连接上），论文中还有节点变量Vij，似乎不参与运作
        vars{i1,i2} = binvar(size(cluster{i1},1), size(cluster{i2},1), 'full'); 
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
        sum_in = sum_in + sum(sum(vars{j,i}));
    end
    F1 = [F1, sum_in + ADN(i) == (h_c-1)*max(n_nodes)];
end
        
% 计算约束条件2（任意一个节点的有效边为h-1个，h为cluster个数）在ADN版本中，改为<=1
F2 = [];
for i1=1:h_c-1 % weight 为上三角阵
    for i2=i1+1:h_c
        for j1=1:size(cluster{i1},1)
            F2 = [F2, sum(vars{i1,i2}(j1,:))<=1]; % i1中的j1至多连接到i2中的某一个
        end
        for j2=1:size(cluster{i2},1)
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
            fprintf('计算 %d — %d — %d 的约束...\n', i1,i2,i3);
            for j1=1:size(cluster{i1},1)
                for j2=1:size(cluster{i2},1)
                    for j3=1:size(cluster{i3},1)
                    	F3 = [F3, vars{i1,i2}(j1,j2)+vars{i2,i3}(j2,j3)<=1+vars{i1,i3}(j1,j3)];
                    end
                end
            end
        end
    end
end
toc

%% 5、计算目标函数并进行求解
F = [F1, F2, F3];
OBJ = 0;
for i1=1:h_c-1
    for i2=i1+1:h_c
        OBJ = OBJ + sum(sum( vars{i1,i2}.*weight{i1,i2} ));
    end
end
OBJ = OBJ + Cd/2*sum(ADN); % 把ADN的值*权重也算上


% 求解BILP
options = sdpsettings('verbose',0,'solver','cplex');
sol = solvesdp( F, -OBJ, options )
for i1=1:h_c-1
    for i2=i1+1:h_c
        vars{i1,i2} = value(vars{i1,i2});
        vars{i2,i1} = vars{i1,i2}';
    end
end

ADN = value(ADN);

%% 6、画出连接图
ADN_xy = [];
for i=1:4
    switch i % 将nan投射到顶角上
        case 1
            xy = [10 10];
            ADN_xy = [ADN_xy; xy];
        case 2
            xy = [-10 10];
            ADN_xy = [ADN_xy; xy];
        case 3
            xy = [-10 -10];
            ADN_xy = [ADN_xy; xy];
        case 4
            xy = [10 -10];
            ADN_xy = [ADN_xy; xy];
    end
    plot(xy(1),xy(2), 'rs');
    text(xy(1)+1, xy(2)+1, ['ADN: ', num2str(ADN(i))]);
end

color = colormap(hsv); % 选择颜色
color = color(randperm(size(color,1)),:);
color_cluster = cell(4,1);
color_cluster{1} = color(1:size(cluster{1},1),:);

for i1=1:4
    for i2=1:4
        if abs(i1-i2)==2 || i1==i2
            continue; % 不画斜线
        end
        for j1=1:size(cluster{i1},1)
            j2 = find(vars{i1,i2}(j1,:));
            if isempty(j2) % i1中的j1与ADN相连
                color_cluster{i1}(j1,:) = [0 0 0];
                line([cluster{i1}(j1,1), ADN_xy(i2,1)],...
                    [cluster{i1}(j1,2), ADN_xy(i2,2)],'color',[0,0,0]);
            else
                tc = color_cluster{i1}(j1,:);
                color_cluster{i2}(j2,:) = tc;
                line([cluster{i1}(j1,1), cluster{i2}(j2,1)],...
                    [cluster{i1}(j1,2), cluster{i2}(j2,2)],'color',tc);
            end
        end
    end
end
hold off;




