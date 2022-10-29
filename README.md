# ElectromagneticWavesMATLAB
MATLAB电磁场与电磁波：极化波的合成与分解，有耗介质中电磁波的传播特性。（Just a final assignment for electromagnetic wave）
代码见下，注意是多个文件，每个代码之间有间隔，别一下全复制到一个文件里去跑了。

Matlab程序代码




理想介质中静态极化波与动态圆极化波的仿真实现：
clear all;
clc;
% 用两个线极化波合成线极化波、圆极化波和椭圆极化波，并实现
% ==========参数========== %
N = 10000; % 单个极化波由N个样点构成
f = 100000000; %电磁波频率为100MHz
k = f/3e8 * 2 * pi; %波数为k
t = 0; %时间
xscale = 2 * pi / k * 5; % X轴显示的长度，单位为m，可以显示五个周期
dx = xscale / N; % 每个样点间的距离
x = 0 : dx : xscale-dx;
figure;
%======生成线极化波1======%
E1 = zeros(N,3); %每一列存放一个方向的坐标值
A1 = 0.5; %幅值为0.5
pha1 = 0; %初始相位
E1(:,2) = A1 * cos(2 * pi * f * t + pha1 - k * x); %沿x轴正方向传播
E1(:,1) = x; %传播方向为x轴
subplot(2,2,1);
scatter3(E1(:,1),E1(:,2),E1(:,3),'.');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======生成线极化波2======%
E2 = zeros(N,3); %每一列存放一个方向的坐标值
A2 = 0.5; %幅值为0.5
pha2 = 0; %初始相位
E2(:,3) = A2 * cos(2 * pi * f * t + pha2 - k * x);
E2(:,1) = x; %传播方向为x轴
subplot(2,2,2);
scatter3(E2(:,1),E2(:,2),E2(:,3),'.');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======生成线极化波3======%
E3 = zeros(N,3); %每一列存放一个方向的坐标值
A3 = 0.5; %幅值为0.5
pha3 = pi/2; %初始相位
E3(:,3) = A3 * cos(2 * pi * f * t + pha3 - k * x);
E3(:,1) = x; %传播方向为x轴
subplot(2,2,3);
scatter3(E3(:,1),E3(:,2),E3(:,3),'.');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======生成线极化波4======%
E4 = zeros(N,3); %每一列存放一个方向的坐标值
A4 = 0.8; %幅值为0.8
pha4 = pi/2; %初始相位
E4(:,2) = A4 * cos(2 * pi * f * t + pha4 - k * x);
E4(:,1) = x; %传播方向为x轴
subplot(2,2,4);
scatter3(E4(:,1),E4(:,2),E4(:,3),'.');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======合成线极化波======%
Ec1 = zeros(N,3);
Ec1(:,1) = x;
Ec1(:,2) = E1(:,2);
Ec1(:,3) = E2(:,3);
figure;
scatter3(Ec1(:,1),Ec1(:,2),Ec1(:,3),'.','g');
hold on;
scatter3(E1(:,1),E1(:,2),E1(:,3),'.','b');
hold on;
scatter3(E2(:,1),E2(:,2),E2(:,3),'.','b');
%======合成圆极化波======%
Ec2 = zeros(N,3);
Ec2(:,1) = x;
Ec2(:,2) = E1(:,2);
Ec2(:,3) = E3(:,3);
figure;
scatter3(Ec2(:,1),Ec2(:,2),Ec2(:,3),'.','g');
hold on;
scatter3(E1(:,1),E1(:,2),E1(:,3),'.','b');
hold on;
scatter3(E3(:,1),E3(:,2),E3(:,3),'.','b');
%======合成椭圆极化波(幅值不同)======%
Ec3 = zeros(N,3);
Ec3(:,1) = x;
Ec3(:,2) = E4(:,2);
Ec3(:,3) = E2(:,3);
figure;
%
scatter3(Ec3(:,1),Ec3(:,2),Ec3(:,3),'.','MarkerFaceAlpha',.0,'MarkerEdgeAlpha',.
0);
scatter3(Ec3(:,1),Ec3(:,2),Ec3(:,3),'.','g');
hold on;
%
scatter3(E2(:,1),E2(:,2),E2(:,3),'.','b','MarkerFaceAlpha',.0,'MarkerEdgeAlpha',
.0);
scatter3(E2(:,1),E2(:,2),E2(:,3),'.','b');
hold on;
%
scatter3(E4(:,1),E4(:,2),E4(:,3),'.','b','MarkerFaceAlpha',.0,'MarkerEdgeAlpha',
.0);
scatter3(E4(:,1),E4(:,2),E4(:,3),'.','b');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======电磁波传播动态演示======%
E_dy = zeros(3*N,3); %实现三个波形同时运动
c = zeros(3*N,3);
c(1:N,2) = 1; %c用于后续绘图区分三种波形的颜色
c(N+1:2*N,1) = 1;
c(2*N+1:3*N,3) = 1;
figure;
for t = 0 : 1 / (50*f) : 3 / f %演示时长与帧数
E1(:,2) = A1 * cos(2 * pi * f * t + pha1 - k * x); %沿x轴正方向传播
E3(:,3) = A3 * cos(2 * pi * f * t + pha3 - k * x);
Ec2(:,2) = E1(:,2);
Ec2(:,3) = E3(:,3);
E_dy(1:N,:) = E1;
E_dy(N+1:N*2,:) = E3;
E_dy(N*2+1:3*N,:) = Ec2;
scatter3(E_dy(:,1),E_dy(:,2),E_dy(:,3),10,c);
pause(0.1);
end







良介质中静态极化波与动态圆极化波的仿真实现：
clear all;
clc;
% 用两个线极化波合成线极化波、圆极化波和椭圆极化波
% ==========参数========== %
N = 10000; % 单个极化波由N个样点构成
f = 100000000; %电磁波频率为100MHz
k = f/3e8 * 2 * pi; %相位常数为k
a = 0.005; %衰减因子为a
t = 0; %时间
xscale = 2 * pi / k * 30; % X轴显示的长度，单位为m，可以显示五个周期
dx = xscale / N; % 每个样点间的距离
x = 0 : dx : xscale-dx;
figure;
%======生成线极化波1======%
E1 = zeros(N,3); %每一列存放一个方向的坐标值
A1 = 0.5*exp(-a*x); %幅值
pha1 = 0; %初始相位
E1(:,2) = A1 .* cos(2 * pi * f * t + pha1 - k * x); %沿x轴正方向传播
E1(:,1) = x; %传播方向为x轴
subplot(2,2,1);
scatter3(E1(:,1),E1(:,2),E1(:,3),'.');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======生成线极化波2======%
E2 = zeros(N,3); %每一列存放一个方向的坐标值
A2 = 0.5*exp(-a*x); %幅值为0.5
pha2 = 0; %初始相位
E2(:,3) = A2 .* cos(2 * pi * f * t + pha2 - k * x);
E2(:,1) = x; %传播方向为x轴
subplot(2,2,2);
scatter3(E2(:,1),E2(:,2),E2(:,3),'.');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======生成线极化波3======%
E3 = zeros(N,3); %每一列存放一个方向的坐标值
A3 = 0.5*exp(-a*x); %幅值为0.5
pha3 = pi/2; %初始相位
E3(:,3) = A3 .* cos(2 * pi * f * t + pha3 - k * x);
E3(:,1) = x; %传播方向为x轴
subplot(2,2,3);
scatter3(E3(:,1),E3(:,2),E3(:,3),'.');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======生成线极化波4======%
E4 = zeros(N,3); %每一列存放一个方向的坐标值
A4 = 0.8*exp(-a*x); %幅值为0.8
pha4 = pi/2; %初始相位
E4(:,2) = A4 .* cos(2 * pi * f * t + pha4 - k * x);
E4(:,1) = x; %传播方向为x轴
subplot(2,2,4);
scatter3(E4(:,1),E4(:,2),E4(:,3),'.');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======合成线极化波======%
Ec1 = zeros(N,3);
Ec1(:,1) = x;
Ec1(:,2) = E1(:,2);
Ec1(:,3) = E2(:,3);
figure;
scatter3(Ec1(:,1),Ec1(:,2),Ec1(:,3),'.','g');
hold on;
scatter3(E1(:,1),E1(:,2),E1(:,3),'.','b');
hold on;
scatter3(E2(:,1),E2(:,2),E2(:,3),'.','b');
%======合成圆极化波======%
Ec2 = zeros(N,3);
Ec2(:,1) = x;
Ec2(:,2) = E1(:,2);
Ec2(:,3) = E3(:,3);
figure;
scatter3(Ec2(:,1),Ec2(:,2),Ec2(:,3),'.','g');
hold on;
scatter3(E1(:,1),E1(:,2),E1(:,3),'.','b');
hold on;
scatter3(E3(:,1),E3(:,2),E3(:,3),'.','b');
%======合成椭圆极化波(幅值不同)======%
Ec3 = zeros(N,3);
Ec3(:,1) = x;
Ec3(:,2) = E4(:,2);
Ec3(:,3) = E2(:,3);
figure;
%
scatter3(Ec3(:,1),Ec3(:,2),Ec3(:,3),'.','MarkerFaceAlpha',.0,'MarkerEdgeAlpha',.
0);
scatter3(Ec3(:,1),Ec3(:,2),Ec3(:,3),'.','g');
hold on;
%
scatter3(E2(:,1),E2(:,2),E2(:,3),'.','b','MarkerFaceAlpha',.0,'MarkerEdgeAlpha',
.0);
scatter3(E2(:,1),E2(:,2),E2(:,3),'.','b');
hold on;
%
scatter3(E4(:,1),E4(:,2),E4(:,3),'.','b','MarkerFaceAlpha',.0,'MarkerEdgeAlpha',
.0);
scatter3(E4(:,1),E4(:,2),E4(:,3),'.','b');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======电磁波传播动态演示======%
E_dy = zeros(3*N,3); %实现三个波形同时运动
c = zeros(3*N,3);
c(1:N,2) = 1; %c用于后续绘图区分三种波形的颜色
c(N+1:2*N,1) = 1;
c(2*N+1:3*N,3) = 1;
figure;
for t = 0 : 1 / (50*f) : 3 / f %演示时长与帧数
E1(:,2) = A1 .* cos(2 * pi * f * t + pha1 - k * x); %沿x轴正方向传播
E3(:,3) = A3 .* cos(2 * pi * f * t + pha3 - k * x);
Ec2(:,2) = E1(:,2);
Ec2(:,3) = E3(:,3);
E_dy(1:N,:) = E1;
E_dy(N+1:N*2,:) = E3;
E_dy(N*2+1:3*N,:) = Ec2;
scatter3(E_dy(:,1),E_dy(:,2),E_dy(:,3),1,c);
pause(0.1);
end










良导体中电磁波的传播
clear all;
clc;
% ==========参数========== %
N = 10000; % 单个极化波由N个样点构成
f1 = 100000000; %电磁波1频率为100MHz
f2 = 150000000; %电磁波2频率为150MHz
S = 1.2; %取电导率为1.2
mu = 4 * pi * 1e-7; %磁导率
k1 = sqrt(pi * f1 * mu * S); %电磁波1相位常数
k2 = sqrt(pi * f2 * mu * S); %电磁波2相位常数
a1 = k1; %电磁波1衰减因子
a2 = k2; %电磁波2衰减因子
t = 0; %时间
xscale = 2 * pi / (k1-20) * 2; % X轴显示的长度，单位为m，可以显示2个周期
dx = xscale / N; % 每个样点间的距离
x = 0 : dx : xscale-dx;
figure;
%======良导体中电磁波1======%
E1 = zeros(N,3); %每一列存放一个方向的坐标值
A1 = 0.5*exp(-a1*x); %幅值
pha1 = 0; %初始相位
E1(:,3) = A1 .* cos(2 * pi * f1 * t + pha1 - k1 * x); %沿x轴正方向传播
E1(:,1) = x; %传播方向为x轴
E1(:,2) = -0.2;
% subplot(2,2,1);
scatter3(E1(:,1),E1(:,2),E1(:,3),'.');
hold on;
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======良导体中电磁波2======%
E2 = zeros(N,3); %每一列存放一个方向的坐标值
A2 = 0.5*exp(-a2*x); %幅值为
pha2 = 0; %初始相位
E2(:,3) = A2 .* cos(2 * pi * f2 * t + pha2 - k2 * x);
E2(:,1) = x; %传播方向为x轴
E2(:,2) = 0.2;
% subplot(2,2,2);
scatter3(E2(:,1),E2(:,2),E2(:,3),'.');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);







半导体中电磁波的传播
clear all;
clc;
% ==========参数========== %
N = 100000; % 单个极化波由N个样点构成
f1 = 600000000; %电磁波1频率为600MHz
f2 = 900000000; %电磁波2频率为900MHz
S = 0.0013 %取电导率为0.0013
si = 1 / 36 / pi * 1e-9; %介电常数
mu = 4 * pi * 1e-7; %磁导率
K1 = 2 * pi * f1 * sqrt(mu * si *(1 - j * S / f1 / si));
K2 = 2 * pi * f2 * sqrt(mu * si *(1 - j * S / f2 / si));
k1 = real(K1); %电磁波1相位常数
k2 = real(K2); %电磁波2相位常数
a1 = -imag(K1); %电磁波1衰减因子
a2 = -imag(K2); %电磁波2衰减因子
t = 0; %时间
xscale = 2 * pi / (k1-11) * 6; % X轴显示的长度，单位为m，可以显示2个周期
dx = xscale / N; % 每个样点间的距离
x = 0 : dx : xscale-dx;
figure;
%======半导体中电磁波1======%
E1 = zeros(N,3); %每一列存放一个方向的坐标值
A1 = 0.5*exp(-a1*x); %幅值
pha1 = 0; %初始相位
E1(:,3) = A1 .* cos(2 * pi * f1 * t + pha1 - k1 * x); %沿x轴正方向传播
E1(:,1) = x; %传播方向为x轴
E1(:,2) = -0.2;
% subplot(2,2,1);
scatter3(E1(:,1),E1(:,2),E1(:,3),'.');
hold on;
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);
%======半导体中电磁波2======%
E2 = zeros(N,3); %每一列存放一个方向的坐标值
A2 = 0.5*exp(-a2*x); %幅值为
pha2 = 0; %初始相位
E2(:,3) = A2 .* cos(2 * pi * f2 * t + pha2 - k2 * x);
E2(:,1) = x; %传播方向为x轴
E2(:,2) = 0.2;
% subplot(2,2,2);
scatter3(E2(:,1),E2(:,2),E2(:,3),'.');
axis([-inf,+inf,-0.8,0.8,-0.8,0.8]);





半导体代码2：
clear;
clc;
grid on;
x=[0:0.1:2000];
% x=[0:1:100];
zero = 0*ones(size(x));
E = ones(size(x));
H = ones(size(x));
t=0;
w=600*pi;
u0=4*pi*1e-7*1000;%真空中的磁导率
e0=1e-9/(36*pi)*100;%真空中的介电质常数
a=w*e0*sqrt(u0/e0);
b=w*sqrt(e0*u0);
% E = exp(-0.5*x).*cos(20*pi*t-0.5*x);
% H = (1/120/pi)*exp(-0.5*x).*cos(20*pi*t-0.5*x-pi/4);
E = exp(-a*x).*cos(w*t-b*x);
H = (1/120/pi)*exp(-a*x).*cos(w*t-b*x-pi/4);
quiver3(x,zero,zero,zero,zero,E,'R');
hold on;
quiver3(x,zero,zero,zero,H,zero,'B');
w1=1000*pi;
u0=4*pi*1e-7*1000;%真空中的磁导率
e0=1e-9/(36*pi)*100;%真空中的介电质常数
a1=w1*e0*sqrt(u0/e0);
b1=w1*sqrt(e0*u0);
% E = exp(-0.5*x).*cos(20*pi*t-0.5*x);
% H = (1/120/pi)*exp(-0.5*x).*cos(20*pi*t-0.5*x-pi/4);
E1 = exp(-a1*x).*cos(w1*t-b1*x);
H1 = (1/120/pi)*exp(-a1*x).*cos(w*t-b1*x-pi/4);
quiver3(x,zero+100,zero,zero,zero,E1,'g');
hold on;
quiver3(x,zero+100,zero,zero,H1,zero,'r');
ti = title('电磁波在半导体中传播','color','k','fontsize',12);
ylabel('磁场变化方向','fontsize',12);
zlabel('电场变化方向','fontsize',12);
% zlim([-0.05,0.1]);
% ylim([-0.05,0.1]);
hold off;
