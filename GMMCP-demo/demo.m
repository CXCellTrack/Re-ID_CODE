clear, close all
 
% ʵ��������ƪ���µ�demo
% GMMCP Tracker: Globally Optimal Generalized Maximum Multi Clique Problem for Multiple Object Tracking

%% 1����������
% ����ͼ1������4�����ݣ���ƽ���ϵĵ�����ʾ
rng(0)
n = 35; r = 10;
theta = linspace(0,360,n);
data = [r*cosd(theta)', r*sind(theta)']; % ��һ��Բ��������һЩ����
noisy = (rand(n,2)*2-1)*2;
data = data + noisy;

figure
plot(data(:,1),data(:,2),'bo');hold on;
line([-15,15],[0, 0],'color','r','linewidth',2);  line([0,0],[-15,15],'color','r','linewidth',2);

% ������data ����4��cluster��
cluster = cell(4,1);
for h=1:size(data,1)
    if data(h,1)>0 && data(h,2)>0 % ��1����
        cluster{1} = [cluster{1}; data(h,:)];
    elseif data(h,1)<0 && data(h,2)>0 % ��2����
        cluster{2} = [cluster{2}; data(h,:)];
    elseif data(h,1)<0 && data(h,2)<0 % ��3����
        cluster{3} = [cluster{3}; data(h,:)];       
    elseif data(h,1)>0 && data(h,2)<0 % ��4����
        cluster{4} = [cluster{4}; data(h,:)];
    end
end

%% 2�������ƽڵ�
n_nodes = zeros(size(cluster)); % ÿ��cluster�еĽڵ���Ŀ
for i=1:numel(cluster)
    n_nodes(i) = size(cluster{i},1);
end
% ÿ��cluster��Ҫ������ƽڵ���Ŀ���Ͻ� Nd_ub
Nd_ub = zeros(size(n_nodes));
for i=1:numel(Nd_ub)
    Nd_ub(i) = sum(n_nodes(1:i)) + sum(n_nodes(i:end));
end
% �������ߵ�˵����ȡ1/3��С���Ͻ�͹�����
% K = round(Nd_ub(1)/10) + n_nodes(1);
K = max(n_nodes);
% ���ƽڵ����cluster�У�Ϊ�˱��ڼ���Ȩ�أ�����ֵΪnan
for i=1:numel(cluster)
    cluster{i} = [cluster{i}; nan*zeros(K-n_nodes(i),2)];
end

%% 3������ߵ�Ȩֵ
weight = cell(4,4);
for i1=1:3
    for i2=i1+1:4
        weight{i1,i2} = zeros(K,K); 
    end
end
% weight{i1,i2}(j1,j2)Ϊ��i1��cluster�е�j1���ڵ㵽��i2��cluster�е�j2���ڵ�Ȩ��

for i1=1:3 % weight Ϊ��������
    for i2=i1+1:4
        for j1=1:K
            for j2=1:K
                if isnan(cluster{i1}(j1,1)) || isnan(cluster{i2}(j2,1))
                    weight{i1,i2}(j1,j2) = 0; % �����һ��Ϊ�ƽڵ㣬��wegithΪ0
                else % ����weighΪ����ĵ�����Խ��weightԽ��
                    weight{i1,i2}(j1,j2) = 1/norm(cluster{i1}(j1,:)-cluster{i2}(j2,:));
                end
            end
        end
        weight{i2,i1} = weight{i1,i2}'; % �Գ�
    end
end

%% 4������Լ������
vars = cell(4,4);
for i1=1:3
    for i2=i1+1:4 % varsָ���Ǳߵı�����0��1�����Ƿ������ϣ��������л��нڵ����Vij���ƺ�����������
        vars{i1,i2} = binvar(K, K, 'full'); 
        vars{i2,i1} = vars{i1,i2}'; % ������
    end
end % ���������ÿ��������Ӧһ��Ȩֵ

% ����Լ������1��ÿ��cluster�нڵ�ĸ�������Ϊk���������������γ�k��clique���Զ���������
% ����Լ������2������һ���ڵ����Ч��Ϊh-1����hΪcluster������
F2 = [];
for i1=1:3 % weight Ϊ��������
    for i2=i1+1:4
        for j1=1:K
            F2 = [F2, sum(vars{i1,i2}(j1,:))==1, sum(vars{i1,i2}(:,j1))==1]; % i1�е�j1�������ӵ�i2�е�ĳһ��
        end
    end
end

% ����Լ������3����a���b������b���c������aҪ��c����
F3 = [];
tic
for i1=1:4 % �ǳ��ǳ��ǳ��� eij������ĿΪ h*(h-1)/2*K^2 ������Ϊ 6*14^2 = 1176
    for i2=i1+1:4 % ÿ3������֮�����3��Լ��3�����Լ��3����Ϊ: C_h_2*(h-2)*K^3 = 6*2*14^3 = 32928
        % ��Ϊ�� i1��i2�ĺͣ����2��û�б��붼��1��4�������ʹԼ����Ŀ����������������䣩
        for i3=1:4 % i1��i2����λ�ã���ʵ��ͬһ��Լ������˲����涨i2>i1��
            % ����cluter���벻��ͬ
            if i1==i2 || i2==i3 || i3==i1
                continue;
            end
            fprintf('���� %d �� %d �� %d ��Լ��...\n', i1,i2,i3);
            for j1=1:K
                for j2=1:K
                    for j3=1:K
                         % ============================================== %
%                         fprintf('���� %d��%d �� %d��%d �� %d��%d ��Լ��\n', i1,j1,i2,j2,i3,j3);
%                         % ֻ���ܺ����ǰ��С����ϣ���Ϊ��ǰ����var{3,1}Ϊ�գ�������ת��Ϊvar{1,3}��
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
                        % 2016.4.15 ��var���Ϊ������󣬾Ϳ���ֱ�Ӽ�����
                        F3 = [F3, vars{i1,i2}(j1,j2)+vars{i2,i3}(j2,j3)<=1+vars{i1,i3}(j1,j3)];
                    end
                end
            end
        end
    end
end
tic

%% 5������Ŀ�꺯�����������
F = [F2, F3];
OBJ = 0;
for i1=1:3
    for i2=i1+1:4
        OBJ = OBJ + sum(sum( vars{i1,i2}.*weight{i1,i2} ));
    end
end

% ���BILP
options = sdpsettings('verbose',0,'solver','gurobi');
sol = solvesdp( F, -OBJ, options )
for i1=1:3
    for i2=i1+1:4
        vars{i1,i2} = value(vars{i1,i2});   
    end
end

%% 6����������ͼ
cluster_dummy = cluster; % ��dummyת��Ϊ��������
for i=1:4
    switch i % ��nanͶ�䵽������
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

color = colormap(hsv); % ѡ����ɫ
color = color(randperm(size(color,1)),:);
color_cluster = cell(4,1);
color_cluster{1} = color(1:K,:);

for i1=1:3
    for i2=i1+1:4
        if abs(i1-i2)==2
            continue; % ����б��
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




