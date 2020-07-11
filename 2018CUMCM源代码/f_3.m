function ObjVal = f_3(Chrom)

for iter=1:size(Chrom,1)
   
dt=1;%时间间隔
dn=[30,100,180,250,250];%各层分割点数
L=10^(-3)*[0.6 Chrom(iter,1) 3.6 Chrom(iter,2) 250];
dx=L./dn;%距离间隔
p=[300 862 74.2 1.18 1000];%密度
c=[1377 2100 1726 1005 2600];%比热容
k=[0.082 0.37 0.045 0.028 0.190];%热传导率
r=(dt./dx.^2).*k./(c.*p);%稳定性系数

T_environment=80;%环境温度
% T_skin_obj=xlsread('C:\Users\16001\Desktop\数学建模\2018CUMCM\CUMCM-2018-Problem-A-Chinese-Appendix.xlsx');
% T_skin=T_skin_obj(:,2);%体表温度
% T_skin=fliplr(T_skin');
% T_skin=T_skin';
T_initial=37;%初始温度
T_skin=37;

T=zeros(1800/dt+1,sum(dn)+4);%初始化温度分布
T(end,:)=T_initial;
T(:,1)=T_environment;
T(:,end)=T_skin;

x_L=dn+1;%界面位置
x_L(2)=x_L(1)+x_L(2);
x_L(3)=x_L(3)+x_L(2);
x_L(4)=x_L(3)+x_L(4);
x_L(5)=x_L(5)+x_L(4);

A=zeros(size(T,2)-2);
v_1=zeros(1,size(A,2));
v_2=zeros(1,size(A,2));
v_3=zeros(1,size(A,2));

k=k./dx;%相对热传导率
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
    
    if x_L(4)+1<i&&i<x_L(5)-1
    v_1(i)=-r(5);
    v_2(i)=1+2*r(5);
    v_3(i)=-r(5);
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
     
        %交界4处理
     if x_L(4)-1<=i&&i<=x_L(4)+1
        if i==x_L(4)-1
        v_1(i)=0;
        v_2(i)=1+2*r(4)-r(4)*(k(4)/(k(4)+k(5)));
        v_3(i)=-k(4)/(k(4)+k(5));
        
        A(x_L(4)-1,x_L(4)+1)=-r(4)*k(5)/(k(5)+k(4));
        end
        
        if i==x_L(4)
        v_1(i)=-k(5)/(k(5)+k(4));
        v_2(i)=1;
        v_3(i)=0;
        end
        
        if i==x_L(4)+1  
        v_1(i)=-r(5);
        v_2(i)=-k(5)/(k(4)+k(5))*r(5)+1+2*r(5);
        v_3(i)=-r(5);
        
        A(x_L(4)+1,x_L(4)-1)=-r(5)*k(4)/(k(5)+k(4));
        end
     end
     
     A(i,i)=v_2(i);
     if i<size(A,1)
     A(i,i+1)=v_1(i);
     A(i+1,i)=v_3(i);
     end  
end


% A(end-1,end-1)=1+2*r(4)-k(4)/(k(4)+k(5))*r(4);
% A(end,end-1)=-k(4)/(k(4)+k(5));
% A(end-1,end)=0;
%A(end,end)=A(end,end)-r(5);

for m=1:1800
b=T(end-(m-1),2:end-1);
b=b';
b(1)=b(1)+r(1)*T_environment;
b(end)=b(end)+r(5)*T_skin;

% b(end-1)=b(end-1)-37*k(5)/(k(4)+k(5));
% b(end)=k(5)/(k(4)+k(5))*37;
b(x_L(1:4))=0;

T_obj=A\b;

if T_obj(x_L(4))>47
    T_obj(x_L(4))=T_environment;
    break;
end

if 1<=m&&m<=1500
   if T_obj(x_L(4))>44
      T_obj(x_L(4))=T_environment;
      break;
   end  
end

T(end-m,2:end-1)=T_obj';
%fprintf('The %d layers have been solved\n',m+1)
end

ObjVal(iter)=T_obj(x_L(4));

end
