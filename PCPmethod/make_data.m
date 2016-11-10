function [data_n1, data_n2, data_n3] = make_data(n1, n2, n3, r1, r2, isplot)

% 制造一些数据

% 中间点
center = [0,0];
% n1 = 50;
rng(0);
data_n1 =  repmat(center,n1,1) + rand(n1,2)*2-1;

% 内环
% r1 = 10;
% n2 = 100;
theta = linspace(0,360,n2);
data_n2 = [r1*cosd(theta)', r1*sind(theta)'];
noisy = (rand(n2,2)*2-1)*1;
data_n2 = data_n2 + noisy;

% 外环
% r2 = 20;
% n3 = 200;
theta = linspace(0,360,n3);
data_n3 = [r2*cosd(theta)', r2*sind(theta)'];
noisy = (rand(n3,2)*2-1)*1;
data_n3 = data_n3 + noisy;

% 数据显示
if isplot
    plot(data_n1(:,1),data_n1(:,2),'b*');hold on;
    plot(data_n2(:,1),data_n2(:,2),'go');
    plot(data_n3(:,1),data_n3(:,2),'rx');
end



