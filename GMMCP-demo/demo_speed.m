clear, %close all
 
% ʵ��������ƪ���µ�demo
% GMMCP Tracker: Globally Optimal Generalized Maximum Multi Clique Problem for Multiple Object Tracking

%% 1����������
% ����ͼ1������4�����ݣ���ƽ���ϵĵ�����ʾ
rng(0)
n = 40; r = 10;
theta = linspace(0,360,n);
data = [r*cosd(theta)', r*sind(theta)']; % ��һ��Բ��������һЩ����
noisy = (rand(n,2)*2-1)*2;
data = data + noisy;

figure
plot(data(:,1),data(:,2),'bo');hold on;
line([-15,15],[0, 0],'color','r','linewidth',2);  line([0,0],[-15,15],'color','r','linewidth',2);

% ������data ����4��cluster��
h_c = 3; % h_c Ϊcluster�ĸ���
cluster = cell(h_c,1);
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

%% 2������ۺ��ƽڵ㣨ADN��
n_nodes = zeros(size(cluster)); % ÿ��cluster�еĽڵ���Ŀ
for i=1:numel(cluster)
    n_nodes(i) = size(cluster{i},1);
end

% ÿ��cluster�м���Aggregated Dummy Nodes (ADN)
% ��ֵΪ���������������˱�������Ŀ
% ADN�ļ����ڽ��������滮��ʱ��������!

%% 3������ߵ�Ȩֵ

weight = cell(h_c,h_c);
for i1=1:h_c-1
    for i2=i1+1:h_c
        weight{i1,i2} = zeros(size(cluster{i1},1), size(cluster{i2},1)); 
    end
end
% weight{i1,i2}(j1,j2)Ϊ��i1��cluster�е�j1���ڵ㵽��i2��cluster�е�j2���ڵ�Ȩ��

for i1=1:h_c-1 % weight Ϊ��������
    for i2=i1+1:h_c
        for j1=1:size(cluster{i1},1)
            for j2=1:size(cluster{i2},1)
                if isnan(cluster{i1}(j1,1)) || isnan(cluster{i2}(j2,1))
                    weight{i1,i2}(j1,j2) = 0; % �����һ��Ϊ�ƽڵ㣬��wegithΪ0
                else % ����weighΪ����ĵ�����Խ��weightԽ��
                    weight{i1,i2}(j1,j2) = 1/norm(cluster{i1}(j1,:)-cluster{i2}(j2,:));
%                     weight{i2,i1}(j2,j1) = weight{i1,i2}(j1,j2); % �Գ�
                end
            end
        end
    end
end

% ע���ڷ�speedup�汾�У����ӵ�dummy node��ȨֵΪCd��������Ϊ0��
% ��speedup�汾�У�ADN���Ǳ߱��������ǽڵ�������ڵ�ȨֵΪCd/2��������Ϊ0.01��
Cd = 0.00;

%% 4������Լ������
vars = cell(h_c,h_c);
for i1=1:h_c-1 % ���������ÿ��������Ӧһ��Ȩֵ
    for i2=i1+1:h_c % varsָ���Ǳߵı�����0��1�����Ƿ������ϣ��������л��нڵ����Vij���ƺ�����������
        vars{i1,i2} = binvar(size(cluster{i1},1), size(cluster{i2},1), 'full'); 
        vars{i2,i1} = vars{i1,i2}'; % ������
    end
end 
% ���� ADN ���ͱ���
ADN = intvar(1,h_c,'full');

% ����Լ������1������ÿ��cluster i�������ıߵĸ���+ADN(i)�ĸ��� = (h-1)*K����
% KΪϣ�����ɵ�clique�ĸ���)
F1 = [];
for i=1:h_c
    sum_in = 0;
    for j=1:h_c
        sum_in = sum_in + sum(sum(vars{j,i}));
    end
    F1 = [F1, sum_in + ADN(i) == (h_c-1)*max(n_nodes)];
end
        
% ����Լ������2������һ���ڵ����Ч��Ϊh-1����hΪcluster��������ADN�汾�У���Ϊ<=1
F2 = [];
for i1=1:h_c-1 % weight Ϊ��������
    for i2=i1+1:h_c
        for j1=1:size(cluster{i1},1)
            F2 = [F2, sum(vars{i1,i2}(j1,:))<=1]; % i1�е�j1�������ӵ�i2�е�ĳһ��
        end
        for j2=1:size(cluster{i2},1)
            F2 = [F2, sum(vars{i1,i2}(:,j2))<=1]; % i2�е�j2�������ӵ�i1�е�ĳһ��
        end
    end
end

% ����Լ������3����a���b������b���c������aҪ��c����
F3 = [];
tic
for i1=1:h_c-1 % �ǳ��ǳ��ǳ��� eij������ĿΪ h*(h-1)/2*K^2 ������Ϊ 6*14^2 = 1176
    for i2=i1+1:h_c % ÿ3������֮�����3��Լ��3�����Լ��3����Ϊ: Ch_2*(h-2)*K^3 = 6*2*14^3 = 32928
        for i3=1:h_c % i1��i2����λ�ã���ʵ��ͬһ��Լ������˲����涨i2>i1��
            % ����cluter���벻��ͬ
            if i1==i2 || i2==i3 || i3==i1
                continue;
            end
            fprintf('���� %d �� %d �� %d ��Լ��...\n', i1,i2,i3);
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

%% 5������Ŀ�꺯�����������
F = [F1, F2, F3];
OBJ = 0;
for i1=1:h_c-1
    for i2=i1+1:h_c
        OBJ = OBJ + sum(sum( vars{i1,i2}.*weight{i1,i2} ));
    end
end
OBJ = OBJ + Cd/2*sum(ADN); % ��ADN��ֵ*Ȩ��Ҳ����


% ���BILP
options = sdpsettings('verbose',0,'solver','cplex');
sol = solvesdp( F, -OBJ, options )
for i1=1:h_c-1
    for i2=i1+1:h_c
        vars{i1,i2} = value(vars{i1,i2});
        vars{i2,i1} = vars{i1,i2}';
    end
end

ADN = value(ADN);

%% 6����������ͼ
ADN_xy = [];
for i=1:4
    switch i % ��nanͶ�䵽������
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

color = colormap(hsv); % ѡ����ɫ
color = color(randperm(size(color,1)),:);
color_cluster = cell(4,1);
color_cluster{1} = color(1:size(cluster{1},1),:);

for i1=1:4
    for i2=1:4
        if abs(i1-i2)==2 || i1==i2
            continue; % ����б��
        end
        for j1=1:size(cluster{i1},1)
            j2 = find(vars{i1,i2}(j1,:));
            if isempty(j2) % i1�е�j1��ADN����
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




