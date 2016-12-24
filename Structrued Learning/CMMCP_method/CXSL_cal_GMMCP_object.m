function [ pre_connect, ADN_value, phi_x_z_hat, delta_zstar_zhat ] = CXSL_cal_GMMCP_object( w, vars_ADN_F, Cd, phi_x_z, loss )

% ���
vars = vars_ADN_F.vars;
ADN = vars_ADN_F.ADN;
F = vars_ADN_F.F;

% ע���ڷ�speedup�汾�У����ӵ�dummy node��ȨֵΪCd��������Ϊ0��
% ��speedup�汾�У�ADN���Ǳ߱��������ǽڵ�������ڵ�ȨֵΪCd/2��������Ϊ0.01��
% Cd = 0.01;
h_c = size(vars,1);

% 2������Ŀ�꺯����������⣨��ʧ��ǿ��Ԥ�⣩
OBJ = dot(w, phi_x_z) + Cd/2*sum(ADN) + loss; % ��ADN��ֵ*Ȩ��Ҳ����

% ���BILP
disp('    ��ʼ�����ʧ��ǿ��Ԥ������...')
options = sdpsettings('verbose', 0, 'solver', 'cplex');
sol = solvesdp( F, -OBJ, options );
% ����õ��ĸ���������ֵ
if sol.problem == 0      
    phi_x_z_hat = value(phi_x_z);
    delta_zstar_zhat = value(loss);
else
    sol.info
    yalmiperror(sol.problem)
end

% ��ǰ�����Ԥ��ƥ��
pre_connect = cell(size(vars));
for i1=1:h_c-1
    for i2=i1+1:h_c
        pre_connect{i1,i2} = value(vars{i1,i2});
    end
end

ADN_value = value(ADN);



end

