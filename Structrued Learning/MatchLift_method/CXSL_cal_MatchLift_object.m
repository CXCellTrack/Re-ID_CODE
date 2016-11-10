function [ pre_connect, phi_x_zhat, delta_zstar_zhat ] = CXSL_cal_MatchLift_object( w, vars_ADN_F, phi_x_z, loss, n_person_cum )

% ���
vars = vars_ADN_F.vars;
F = vars_ADN_F.F;

%% mapͼ�У��ߵĸ�����ע���������2������ıߣ��˴�Ϊ2������������
% ���½�2.1
n_seq = numel(n_person_cum)-1;
E = n_seq*(n_seq-1)/2; % (Si,Sj)�ʦţ����ظ�����
if ~exist('lambda_l1', 'var')
    lambda_l1 = sqrt(E)/2/n_seq; % �������lambda_l1
end

disp('      ��ʼ�����ʧ��ǿ��Ԥ������...')
obj = dot(w, phi_x_z) - lambda_l1*sum(vars(:)) + loss; % ����ΪL1����X_out��1�ĸ���
options = sdpsettings('verbose', 0, 'solver','mosek'); % ��ҵ����� mosek �� sedumi ��ܶ�
sol = solvesdp( F, -obj, options ); % checkset(F)

%% ����õ��ĸ���������ֵ
if sol.problem == 0      
    phi_x_zhat = value(phi_x_z);
    delta_zstar_zhat = value(loss);
else
    sol.info
    yalmiperror(sol.problem)
end

X_out = value(vars);
X_out_int = round(X_out);
pre_connect = mat2bolckmat(X_out_int, n_person_cum );

