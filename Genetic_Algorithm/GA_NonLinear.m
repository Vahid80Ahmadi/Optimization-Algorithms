clc;
clear all;
close all;
%%
global NFE;
NFE=0;
model=CreateModel();
cost=@(x) Cost(x,model);

nvar=4;

varmax=2;
varmin=0;
varsize=[1 nvar];
%%
%GA Parameters

maxit=100;
npop=100;

pc=0.8;
nc=2*round((pc*npop)/2);

pm=0.3;
nm=round(pm*npop);

beta=5.5;

TournamentSize=3;

mu=0.05;
%%
%Initialization
chromosome.position=[];
chromosome.cost=[];
chromosome.sol=[];
c=repmat(chromosome,npop,1);

for i=1:npop
    c(i).position=unifrnd(varmin,varmax,1,nvar);
    [c(i).cost,c(i).sol]=cost(c(i).position);
end

Costs=[c.cost];
[Costs,SortOrder]=sort(Costs);
c=c(SortOrder);

bestsol=c(1);

bestcost=zeros(maxit,1);
worstcost=c(end).cost;
nfe=zeros(maxit,1);
%%
%GA Main Loop
for it=1:maxit
    %Calculate Selction Probabities
    P=exp(-beta*(Costs/worstcost));
    P=P/sum(P);
    
    %Crossover
    popc=repmat(chromosome,nc/2,2);
    for k=1:nc/2
         
        i1=RWS(P);
        i2=RWS(P);
        %i1=TournamentSelection(c,TournamentSize);
        %i2=TournamentSelection(c,TournamentSize);
        
        
        %Copy
        c1=c(i1);
        c2=c(i2);
        
        [popc(k,1).position,popc(k,2).position]=Crossover(c1.position,c2.position,varmin,varmax);
        
        [popc(k,1).cost,popc(k,1).sol]=cost(popc(k,1).position);
        [popc(k,2).cost,popc(k,2).sol]=cost(popc(k,2).position);
       
    end
    popc=popc(:);
    
    %Mutation
    popm=repmat(chromosome,nm,1);
    for k=1:nm
        
        i=randi([1 npop]);
        
        popm(k).position=Mutation(c(i).position,mu,varmin,varmax);
        [popm(k).cost,popm(k).sol]=cost(popm(k).position);
    end
    
    %Merge
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
    
    figure(1);
    PlotSolution(model,bestsol.sol);
    
    disp(['Iteratin ' num2str(it) ' Best Cost ' num2str(bestcost(it)) ' NFE ' num2str(nfe(it))]);
end

%%
%Functions
function model=CreateModel()
    
    model.T=20;
    model.dt=0.001;
    model.x02=1;
    model.x01=2;
    model.a=1.2;
    model.b=0.6;
    model.c=0.8;
    model.d=0.3;
    
    out=SimulateModel(model);
    
    model.t=out.t;
    model.x1=out.x1;
    model.x2=out.x2;
    
end

function out=SimulateModel(params)

    a=params.a;
    b=params.b;
    c=params.c;
    d=params.d;
    
    x01=params.x01;
    x02=params.x02;
    
    T=params.T;
    dt=params.dt;
    
    t=0:dt:T;
    
    K=numel(t);
    
    x1=zeros(size(t));
    x2=zeros(size(t));
    
    x1(1)=x01;
    x2(1)=x02;
    
    for k=1:K-1
        x1(k+1)=x1(k)+(a*x1(k)-b*x1(k)*x2(k))*dt;
        x2(k+1)=x2(k)+(-c*x2(k)+d*x1(k)*x2(k))*dt;
    end
    
    out.t=t;
    out.x1=x1;
    out.x2=x2;
    
end

function [z,sol]=Cost(x,model)
    
    global NFE;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;
    
    p.T=model.T;
    p.dt=model.dt;
    p.x01=model.x01;
    p.x02=model.x02;
    p.a=x(1);
    p.b=x(2);
    p.c=x(3);
    p.d=x(4);
    
    out=SimulateModel(p);
    
    x1=out.x1;
    x2=out.x2;
    
    e1=model.x1-x1;
    e2=model.x2-x2;
    
    z=mean(e1.^2)+mean(e2.^2);
    
    sol.e1=e1;
    sol.e2=e2;
    sol.a=x(1);
    sol.b=x(2);
    sol.c=x(3);
    sol.d=x(4);
    sol.x1=x1;
    sol.x2=x2;
    sol.t=out.t;
end


function [y1,y2]=Crossover(c1,c2,varmin,varmax)
    
    gama=0.1;
    alfa=unifrnd(-gama,1+gama,size(c1));
   
    y1=alfa.*c1+(1-alfa).*c2;
    y2=(1-alfa).*c1+alfa.*c2;
    
    y1=min(max(y1,varmin),varmax);
    y2=min(max(y2,varmin),varmax);
end

function y=Mutation(x,mu,varmin,varmax)
    
    
    n=numel(x);
    
    SelectedN=ceil(n*mu);
    
    List=randsample(n,SelectedN);
    sigma=0.1*(varmax-varmin);
    y=x;
    y(List)=x(List)+sigma*randn(size(List));
    y=min(max(y,varmin),varmax);
    
end

function i=RWS(P) 
    
    r=rand;
    c=cumsum(P);
    i=min(find(r<c,1,'first'));
    
end

function i=TournamentSelection(c,m)
    
    n=numel(c);
    S=randsample(n,m);
    
    SPop=c(S);
    
    sc=[SPop.cost];
    
    [~,j]=min(sc);
    i=S(j);
end

function PlotSolution(model,c)
    
    %x1
    subplot(2,2,1);
    plot(model.t,model.x1,'r','LineWidth',2);
    hold on;
    plot(c.t,c.x1,'b:','LineWidth',2);
    legend('Main X1','Tried X1');
    xlabel('Time');
    ylabel('X1');
    hold off;
    
    %x2
    subplot(2,2,3);
    plot(model.t,model.x2,'r','LineWidth',2);
    hold on;
    plot(c.t,c.x2,'b:','LineWidth',2);
    legend('Main X2','Tried X2');
    xlabel('Time');
    ylabel('X2');
    hold off;
    
    %x1,x2
    subplot(2,2,[2,4]);
    plot(model.x1,model.x2,'r','LineWidth',2);
    hold on;
    plot(c.x1,c.x2,'b:','LineWidth',2);
    legend('Main X1,X2','Tried X1,X2');
    xlabel('x1');
    ylabel('X2');
    hold off;
    
    
end