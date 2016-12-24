function [ pre_connect, ADN_value ] = cal_GMMCP_object( Kval, vars_ADN_F, Cd )

% ���
vars = vars_ADN_F.vars;
ADN = vars_ADN_F.ADN;
F = vars_ADN_F.F;

% ע���ڷ�speedup�汾�У����ӵ�dummy node��ȨֵΪCd��������Ϊ0��
% ��speedup�汾�У�ADN���Ǳ߱��������ǽڵ�������ڵ�ȨֵΪCd/2��������Ϊ0.01��
% Cd = 0.01;
weight = Kval;
n_seq = size(vars,1);

% 5������Ŀ�꺯�����������
OBJ = 0;
for i1=1:n_seq-1
    for i2=i1+1:n_seq
        OBJ = OBJ + sum(sum( vars{i1,i2}.*weight{i1,i2} ));
    end
end
OBJ = OBJ + Cd/2*sum(ADN); % ��ADN��ֵ*Ȩ��Ҳ����

% ���BILP
disp('��ʼ���������Ϲ滮...')
options = sdpsettings('verbose',0,'solver','cplex');
sol = solvesdp( F, -OBJ, options )
if sol.problem == 0      
    % ��ǰ�����Ԥ��ƥ��
    pre_connect = cell(size(vars));
    for i1=1:n_seq-1
        for i2=i1+1:n_seq
            pre_connect{i1,i2} = value(vars{i1,i2});
        end
    end
else
    sol.info
    yalmiperror(sol.problem)
end

ADN_value = value(ADN);



end

