function vars_ADN_F = cal_GMMCP_constraints( n_person_in_seq )

% ʵ��������ƪ����
% GMMCP Tracker: Globally Optimal Generalized Maximum Multi Clique Problem for Multiple Object Tracking

%% 1������ۺ��ƽڵ㣨ADN�����ں���ֱ��intvar���ɣ�
global each_person_a_clique;
if each_person_a_clique
    % �Ƿ��ÿ���˶�����һ��clique�����������迼��ÿ��clique�ж�������һ�������������
    n_nodes = ones(1, sum(n_person_in_seq));
else
    n_nodes = n_person_in_seq; % ÿ��cluster�еĽڵ���Ŀ
end

% ÿ��cluster�м���Aggregated Dummy Nodes (ADN)
% ��ֵΪ���������������˱�������Ŀ
% ADN�ļ����ڽ��������滮��ʱ��������!

h_c = numel(n_nodes); % h_c Ϊcluster�ĸ���

%% 2������Լ������
vars = cell(h_c);
for i1=1:h_c-1 % ���������ÿ��������Ӧһ��Ȩֵ
    for i2=i1+1:h_c % varsָ���Ǳߵı�����0��1�����Ƿ������ϣ��������л��нڵ����Vij���ƺ�����������
        vars{i1,i2} = binvar(n_nodes(i1), n_nodes(i2), 'full'); 
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
        sum_in = sum_in + sum(sum(vars{i,j}));
    end
    F1 = [F1, sum_in + ADN(i) == (h_c-1)*(10+max(n_nodes))]; % ����KֵΪ 10+max(n_nodes)
end
        
% ����Լ������2������һ���ڵ����Ч��Ϊh-1����hΪcluster��������ADN�汾�У���Ϊ<=1
F2 = [];
for i1=1:h_c-1 % weight Ϊ��������
    for i2=i1+1:h_c
        for j1=1:n_nodes(i1)
            F2 = [F2, sum(vars{i1,i2}(j1,:))<=1]; % i1�е�j1�������ӵ�i2�е�ĳһ��
        end
        for j2=1:n_nodes(i2)
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
            fprintf('���� {%d} + {%d} <= 1 + {%d} ��Լ��...\n', i1,i2,i3);
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

