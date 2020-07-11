clear all;close all;clc

%/*ABC算法参数设置*/
NP=10;%采蜜蜂的个数或观察蜂的个数与蜜源的数量相等
maxCycle=100; %算法最大迭代次数
objfun='f_2'; %目标函数
ub=[25*10^-3]; %参数的取值上限[x1,x2,···]
lb=[0.6*10^-3];%参数的取值下限[x1,x2,···]
model=0;%模式设置：0代表实数规划，1代表整数规划
u=1;%领域搜索收敛因子，取值越小，收敛越慢，取值为(0,1]
du=200;%收敛梯度
%/*ABC算法参数设置*/

FoodNumber=NP/2; %食物数量
D=size(ub,2); %参数个数
limit=NP*D;
runtime=1;%算法运行次数
GlobalMins=zeros(1,runtime);

for r=1:runtime   
Range=ub-lb;
Lower=lb;
Foods = rand(FoodNumber,D);%rand生成FoodNumber*D大小的矩阵，元素为0：1
for i=1:D
    if model==1
    Foods(:,i) = round(Foods(:,i) .* Range(i) + Lower(i));
    else
    Foods(:,i) = Foods(:,i) .* Range(i) + Lower(i);
    end
end

ObjVal=feval(objfun,Foods);%feval用来计算指定函数在某点的函数值
Fitness=calculateFitness(ObjVal);

trial=zeros(1,FoodNumber);%初始化每个食物的trial为0

%/*The best food source is memorized*/
BestInd=find(ObjVal==min(ObjVal));%find函数找到ObjVal数组中的最小值的位置
BestInd=BestInd(end);%取ObjVal数组中最后的最小值的位置
GlobalMin=ObjVal(BestInd);%记录函数的最小值
GlobalParams=Foods(BestInd,:);%记录函数最小值对应的参数（x,y）

iter=1;
while (iter <= maxCycle)
    deta=u^(iter/du);%引入领域搜索收敛因子
%雇佣蜂先采一批蜜源回去供观察蜂观察
%%%%%%%%% EMPLOYED BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:(FoodNumber)
        
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;%在[1:D]随机生成一个整数
        
        %/*A randomly chosen solution is used in producing a mutant 
        %solution of the solution i*/
        neighbour=fix(rand*FoodNumber)+1;%即在所有食物中随机挑选食物
        %/*Randomly selected solution must be different from the solution i        
            while(neighbour==i)
                neighbour=fix(rand*FoodNumber)+1;
            end
        
       sol=Foods(i,:);
       
      if model==0
       sol(Param2Change)=Foods(i,Param2Change)+...
           (Foods(i,Param2Change)-Foods(neighbour,Param2Change))...
           *(rand-0.5)*2*deta;%实数规划
       else
       sol(Param2Change)=Foods(i,Param2Change)+...
           round((Foods(i,Param2Change)-Foods(neighbour,Param2Change))...
           *(rand-0.5)*2*deta);%整数规划 
       end
       %  /*if generated parameter value is out of boundaries, it is 
       %shifted onto the boundaries*/
        ind=find(sol<lb);
        sol(ind)=lb(ind);
        ind=find(sol>ub);
        sol(ind)=ub(ind);
        
        %evaluate new solution
        ObjValSol=feval(objfun,sol);
        FitnessSol=calculateFitness(ObjValSol);
        
       % /*a greedy selection is applied between the current solution i
       %and its mutant 
       if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than 
           %the current solution i, replace the solution with the mutant 
           %and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; 
       end
    end

%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%
prob=Fitness./sum(Fitness);
%%%%%%%%%%%%%%%%%%%%%%%% ONLOOKER BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))
        t=t+1;
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;
        
        %/*A randomly chosen solution is used in producing a mutant 
        %solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i       
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end
        
       sol=Foods(i,:);
       
       if model==0
       sol(Param2Change)=Foods(i,Param2Change)+...
           (Foods(i,Param2Change)-Foods(neighbour,Param2Change))...
           *(rand-0.5)*2*deta;%实数规划
       else
       sol(Param2Change)=Foods(i,Param2Change)+...
           round((Foods(i,Param2Change)-Foods(neighbour,Param2Change))...
           *(rand-0.5)*2*deta);%整数规划 
       end
       
        ind=find(sol<lb);
        sol(ind)=lb(ind);
        ind=find(sol>ub);
        sol(ind)=ub(ind);
        
        %evaluate new solution
        ObjValSol=feval(objfun,sol);
        FitnessSol=calculateFitness(ObjValSol);
        
       % /*a greedy selection is applied between the current solution i and
       %its mutant*/
       if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than 
       %the current solution i, replace the solution with the mutant and 
       %reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; 
       end
    end
    
    i=i+1;
     if (i==(FoodNumber)+1) %确保每个蜜源至少被观察1次
         i=1;
     end 
end

%/*The best food source is memorized*/
         ind=find(ObjVal==min(ObjVal));
         ind=ind(end);
         if (ObjVal(ind)<GlobalMin)
         GlobalMin=ObjVal(ind);
         GlobalParams=Foods(ind,:);%储存最小函数值对应参数
         end
%%%%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    Bas(ind)=0;
    
    if model==0
    sol=(ub-lb).*rand(1,D)+lb;%实数规划
    else
    sol=round((ub-lb).*rand(1,D))+lb;%整数规划
    end
    
    ObjValSol=feval(objfun,sol);
    FitnessSol=calculateFitness(ObjValSol);
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
end

fprintf('代数=%d GlobalMin=%g\n',iter,GlobalMin);
iter=iter+1;
X(iter)=iter;
Y(iter)=GlobalMin;
end% End of ABC
GlobalMins(r)=GlobalMin;
end %end of runs

  figure(1);
  plot(X(2:end),Y(2:end));
  xlabel('迭代次数');
  ylabel('函数取值');
GlobalParams %输出参数
save all

