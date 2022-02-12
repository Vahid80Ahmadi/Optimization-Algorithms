clc;
close all;
clear all;
%%
%Problem Parameters
global NFE;
NFE=0;

model=CreateModel();
cost=@(x) Cost(x,model);

nvar=model.Ns;
varmax=1;
varmin=0;
varsize=([1 nvar]);
%%
%GA Parameters

maxit=200;
npop=200;

pc=0.8;
nc=2*round((pc*npop)/2);

pm=0.3;
nm=round(pm*npop);

beta=5.5;

TournamentSize=3;

mu=0.05;
%%
%Initiaalization
chromosome.position=[];
chromosome.cost=[];

c=repmat(chromosome,npop,1);

for i=1:npop
   
   c(i).position=CreateRandomSolution(model);
   c(i).cost=cost(c(i).position);
end

Costs=[c.cost];
[Costs,SortOrder]=sort(Costs);
c=c(SortOrder);

bestsol=c(1);
worstcost=c(end).cost;
nfe=zeros(1,maxit);
bestcost=zeros(1,maxit);
%%
%GA Main Loop
for it=1:maxit
   %Probability 
    P=exp(-beta*(Costs/worstcost))./sum(exp(-beta*Costs/worstcost));
    
   %Making Chilren & Crossover
   popc=repmat(chromosome,nc/2,2);
   for k=1:nc/2
      
      i1=RWS(P); 
      i2=RWS(P);
      
      c1=c(i1);
      c2=c(i2);
       
      [popc(k,1).position,popc(k,2).position]=Crossover(c1.position,c2.position);
        
      popc(k,1).cost=cost(popc(k,1).position);
      popc(k,2).cost=cost(popc(k,2).position);
         
   end
   popc=popc(:);
   
   %Mutation
   popm=repmat(chromosome,nm,1);
   for j=1:nm
      
       j1=randi([1 npop]);
               
       popm(j).position=Mutation(c(j1).position,mu);
       popm(j).cost=cost(popm(j).position);
   end
   
   c=[c
      popc
      popm];
  
   %Sort
   Costs=[c.cost];
   [Costs,SortOrder]=sort(Costs);
   c=c(SortOrder);
   
   
    %Truncate
    c=c(1:npop);
    Costs=Costs(1:npop);
    
    %Stor Best Solution
    bestsol=c(1);
    
    %Stor Best & Worst Cost
    bestcost(it)=bestsol.cost;
    worstcost=max(c(end).cost,worstcost);
    nfe(it)=NFE;
   
    
    disp(['Iteratin ' num2str(it) ' Best Cost ' num2str(bestcost(it)) ' NFE ' num2str(nfe(it))]);

end
%%
%Results
figure;
semilogy(nfe,bestcost,'LineWidth',2);
xlabel('NFE');
ylabel('Best Cost');
grid on;

%%
%Functions
function model=CreateModel()
    Nc=40;
    Ns=20;
    
    %Customers
    xc=[82 91 12 92 63 9 28 55 96 97 15 98 96 49 80 14 42 92 80 96 66 3 85 94 68 76 75 39 66 17 71 3 27 4 9 83 70 32 95 3];
    yc=[44 38 77 80 18 49 45 65 71 76 27 68 66 16 12 50 96 34 59 22 75 25 51 70 89 96 55 14 15 26 84 25 82 24 93 35 19 25 62 47];
    di=[21 43 31 30 47 18 39 39 22 31 8 7 29 40 47 10 31 26 5 20 12 41 19 29 12 32 17 35 36 39 25 8 15 47 12 42 29 50 8 25];
    
    %Servers
    c=[11243 8647 11029 8745 10192 10921 10411 9060 10430 8635 8127 11813 9090 10691 8081 8459 10938 11797 11954 9444];
    xs=[78 38 24 40 9 13 94 95 57 5 23 35 82 1 4 16 64 73 64 45];
    ys=[54 29 74 18 68 18 36 62 78 8 92 77 48 43 44 30 50 51 81 79];
    
    D=zeros(Nc,Ns);
    
    for i=1:Nc
        for j=1:Ns
            %D(i,j)=sqrt((xc(i)-xs(j))^2+(yc(i)-ys(j))^2);
            D(i,j)=norm([xc(i)-xs(j) yc(i)-ys(j)],1);
        end
    end
    
    model.Nc=Nc;
    model.Ns=Ns;
    model.xc=xc;
    model.yc=yc;
    model.di=di;
    model.c=c;
    model.xs=xs;
    model.ys=ys;
    model.D=D;
    
end

function f=CreateRandomSolution(model)
    
    Ns=model.Ns;
    
    f=randi([0 1],1,Ns);
    
end

function z=Cost(f,model)
    
    global NFE;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;
    
    Nc=model.Nc;
    D=model.D;
    di=model.di;
    c=model.c;
    %for i=1:Nc
    %    D(i:end)=D(i:end).*f;
    %end
    if all(f==0)
        z=inf;
        return;
    end
    
    Dmin=zeros(1,Nc);
    for i=1:Nc
        Dmin(i)=min(D(i,f==1));
    end
    
    z=sum(di.*Dmin)+0.001*sum(f.*c);
    
end

function [y1,y2]=Crossover(c1,c2)
    
    Selection=randi([1 3]);
    
    switch Selection
        case 1
            [y1,y2]=SinglePointCrossover(c1,c2);
        case 2
            [y1,y2]=DoublePointCrossover(c1,c2);
        case 3
            [y1,y2]=UniformCrossover(c1,c2);
    end
end

function i=RWS(p)
    
    select=rand;
    c=cumsum(p);
    i=find(select<c,1,'first');
    
end

function i=TournamentSelection(c,m)
    
    n=numel(c);
    S=randsample(n,m);
    
    SPop=c(S);
    
    sc=[SPop.cost];
    
    [~,j]=min(sc);
    i=S(j);
end

function [y1,y2]=SinglePointCrossover(c1,c2)
    
    n=numel(c1);
    cutpoint=randi([1,n-1]);
    y1=[c1(1:cutpoint) c2(cutpoint+1:end)];
    y2=[c2(1:cutpoint) c1(cutpoint+1:end)];
    
end

function [y1,y2]=DoublePointCrossover(c1,c2)
    
    n=numel(c1);
    cutpoint=randsample(n-1,2);
    cc1=min(cutpoint);
    cc2=max(cutpoint);
    
    y1=[c1(1:cc1) c2(cc1+1:cc2) c1(cc2+1:end)];
    y2=[c2(1:cc1) c1(cc1+1:cc2) c2(cc2+1:end)];
    
end

function [y1,y2]=UniformCrossover(c1,c2)
    
    n=numel(c1);
    
    alfa=randi([0 1],1,n);
    
   y1=alfa.*c1+(1-alfa).*c2;
   y2=(1-alfa).*c1+alfa.*c2;

end

function y=Mutation(x,mu)
    
    pim=mu;
    n=numel(x);
    
    SelectedN=ceil(n*pim);
    
    List=randsample(n,SelectedN);
   
    x(List)=1-x(List);
    y=x;
end
