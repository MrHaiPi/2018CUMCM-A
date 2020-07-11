%古典隐式差分求解二阶多层热传导方程
clc
clear all
close all

dt=1;%时间间隔
dx=10^(-3)*0.1;%距离间隔
L=10^(-3)*[0.6 6 3.6 5];
p=[300 862 74.2 1.18];%密度
c=[1377 2100 1726 1005];%比热容
k=[0.082 0.37 0.045 0.028];%热传导率
r=(dt/dx^2)*k./(c.*p);%稳定性系数

T_environment=75;%环境温度
T_skin_obj=xlsread('CUMCM-2018-Problem-A-Chinese-Appendix.xlsx');
T_skin=T_skin_obj(:,2);%体表温度
T_skin=fliplr(T_skin');
T_skin=T_skin';
T_initial=37;%初始温度

T=zeros(5400/dt+1,sum(L)/dx+1);%初始化温度分布
T(end,:)=T_initial;
T(:,1)=T_environment;
T(:,end)=T_skin;

%T(:,end)=37;

x_L=L/dx+1;%界面位置
x_L(1)=round(x_L(1));
x_L(2)=x_L(1)+x_L(2);
x_L(3)=x_L(3)+x_L(2);
x_L(4)=x_L(3)+x_L(4);

A=zeros(size(T,2)-2);
v_1=zeros(1,size(A,2));
v_2=zeros(1,size(A,2));
v_3=zeros(1,size(A,2));

for i=1:size(A,1) 
    
    if i<x_L(1)-1
    v_1(i)=-r(1);
    v_2(i)=1+2*r(1);
    v_3(i)=-r(1);
    end
    
    if x_L(1)+1<i&&i<x_L(2)-1
    v_1(i)=-r(2);
    v_2(i)=1+2*r(2);
    v_3(i)=-r(2);
    end
    
    if x_L(2)+1<i&&i<x_L(3)-1
    v_1(i)=-r(3);
    v_2(i)=1+2*r(3);
    v_3(i)=-r(3);
    end
    
    if x_L(3)+1<i&&i<x_L(4)-1
    v_1(i)=-r(4);
    v_2(i)=1+2*r(4);
    v_3(i)=-r(4);
    end
    
        %交界1处理
    if x_L(1)-1<=i&&i<=x_L(1)+1
        if i==x_L(1)-1
        v_2(i)=1+2*r(1)-r(1)*(k(1)/(k(1)+k(2)));
        v_3(i)=-k(1)/(k(1)+k(2));
        v_1(i)=0;
        A(x_L(1)-1,x_L(1)+1)=-r(1)*k(2)/(k(1)+k(2));
        end
        
        if i==x_L(1)
        v_1(i)=-k(2)/(k(1)+k(2));
        v_2(i)=1;
        v_3(i)=0;
        end
        
        if i==x_L(1)+1  
        v_1(i)=-r(2);
        v_2(i)=-k(2)/(k(1)+k(2))*r(2)+1+2*r(2);
        v_3(i)=-r(2);
        
        A(x_L(1)+1,x_L(1)-1)=-r(2)*k(1)/(k(1)+k(2));
        end
    end
    
     %交界2处理
    if x_L(2)-1<=i&&i<=x_L(2)+1
        if i==x_L(2)-1
        v_1(i)=0;
        v_2(i)=1+2*r(2)-r(2)*(k(2)/(k(2)+k(3)));
        v_3(i)=-k(2)/(k(2)+k(3));
        
        A(x_L(2)-1,x_L(2)+1)=-r(2)*k(3)/(k(2)+k(3));
        end
        
        if i==x_L(2)
        v_1(i)=-k(3)/(k(2)+k(3));
        v_2(i)=1;
        v_3(i)=0;
        end
        
        if i==x_L(2)+1  
        v_1(i)=-r(3);
        v_2(i)=-k(3)/(k(2)+k(3))*r(3)+1+2*r(3);
        v_3(i)=-r(3);
        
        A(x_L(2)+1,x_L(2)-1)=-r(3)*k(2)/(k(2)+k(3));
        end
    end
    
      %交界3处理
     if x_L(3)-1<=i&&i<=x_L(3)+1
        if i==x_L(3)-1
        v_1(i)=0;
        v_2(i)=1+2*r(3)-r(3)*(k(3)/(k(3)+k(4)));
        v_3(i)=-k(3)/(k(3)+k(4));
        
        A(x_L(3)-1,x_L(3)+1)=-r(3)*k(4)/(k(3)+k(4));
        end
        
        if i==x_L(3)
        v_1(i)=-k(4)/(k(3)+k(4));
        v_2(i)=1;
        v_3(i)=0;
        end
        
        if i==x_L(3)+1  
        v_1(i)=-r(4);
        v_2(i)=-k(4)/(k(3)+k(4))*r(4)+1+2*r(4);
        v_3(i)=-r(4);
        
        A(x_L(3)+1,x_L(3)-1)=-r(4)*k(3)/(k(3)+k(4));
        end
     end
     
     A(i,i)=v_2(i);
     if i<size(A,1)
     A(i,i+1)=v_1(i);
     A(i+1,i)=v_3(i);
     end
     
end

%A(end,end)=A(end,end)-r(4);

for k=1:5400
b=T(end-(k-1),2:end-1);
b=b';
b(1)=b(1)+r(1)*T_environment;
b(end)=b(end)+r(4)*T_skin(end-k);
%b(end)=b(end)+r(4)*37;
b(x_L(1:3))=0;

T_obj=A\b;

T(end-k,2:end-1)=T_obj';
fprintf('The %d layers have been solved\n',k+1)
end

T=flipud(T);
% filname='problem_1.xls';
% xlswrite(filname,T,1);

[X, Y]=meshgrid(1:size(T,2),1:size(T,1));
surf(X,Y,T)
shading interp

figure(1);
title('高温防护服—空气层热传导模型'); 
x1=xlabel(['距离x/', num2str(dx),'m']);       
x2=ylabel('时间t/s');        
x3=zlabel('温度T/℃');        
set(x1,'Rotation',30);    
set(x2,'Rotation',-30); 



X=T(:,x_L(1));
Y=T(:,x_L(2));
Z=T(:,x_L(3));

x=1:size(T,1);
figure(2);
plot(x,X);
hold on
plot(x,Y);
hold on
plot(x,Z);








