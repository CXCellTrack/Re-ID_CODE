function [ pre_connect, ADN_value ] = cal_GMMCP_object( Kval, vars_ADN_F, Cd )

% ���
vars = vars_ADN_F.vars;
ADN = vars_ADN_F.ADN;
F = vars_ADN_F.F;

% ע���ڷ�speedup�汾�У����ӵ�dummy node��ȨֵΪCd��������Ϊ0��
% ��speedup�汾�У�ADN���Ǳ߱��������ǽڵ�������ڵ�ȨֵΪCd/2��������Ϊ0.01��
% Cd = 0.01;
weight = Kval;
h_c = size(vars,1);

% 5������Ŀ�꺯�����������
OBJ = 0;
for i1=1:h_c-1
    for i2=i1+1:h_c
        OBJ = OBJ + sum(sum( vars{i1,i2}.*weight{i1,i2} ));
    end
end
OBJ = OBJ + Cd/2*sum(ADN); % ��ADN��ֵ*Ȩ��Ҳ����

% ���BILP
disp('��ʼ���������Ϲ滮...')
options = sdpsettings('verbose',0,'solver','gurobi');
sol = solvesdp( F, -OBJ, options )

pre_connect = cell(size(vars));
for i1=1:h_c-1
    for i2=i1+1:h_c
        pre_connect{i1,i2} = value(vars{i1,i2});
    end
end

ADN_value = value(ADN);



end

