clc;
clear all;
close all;
%%%Problem Parameters
global NFE;
NFE=0;

model=CreateModel();
cost=@(x) Cost(x,model);

nvar=model.m;
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

c=repmat(chromosome,npop,1);

for i=1:npop
    c(i).position=CreateRandomSolution(model);
    c(i).cost=cost(c(i).position);
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
        
        
        c1=c(i1);
        c2=c(i2);
        
        [popc(k,1).position,popc(k,2).position]=SinglePointCrossover(c1.position,c2.position);
        
        popc(k,1).cost=cost(popc(k,1).position);
        popc(k,2).cost=cost(popc(k,2).position);
       
    end
    popc=popc(:);
    
    %Mutation
    popm=repmat(chromosome,nm,1);
    for k=1:nm
        
        i=randi([1 npop]);
        
        popm(k).position=Mutation(c(i).position);
        popm(k).cost=cost(popm(k).position);
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
   
    
    disp(['Iteratin ' num2str(it) ' Best Cost ' num2str(bestcost(it)) ' NFE ' num2str(nfe(it))]);
end
%%
%Results
figure;
plot(nfe,bestcost,'LineWidth',2);
xlabel('NFE');
ylabel('Best Cost');
grid on;


%%
%Functions 
function model=CreateModel()
   
    n=20;
    m=30;
    
    xp=[58 54 87 26 32 12 94 65 48 64 55 65 54 72 52 100 22 10 11 6 40 45 36 77 63 77 94 98 19 14];
    yp=[70 9 53 53 86 48 39 67 74 52 35 15 59 26 4 76 24 44 69 36 74 39 69 71 44 1 33 42 27 19];

    
    w=[1  5  6  7  4  7  3  4  4  4  2  1  6  5  2  7  7  6  1  6
       5  9  3  6  6  7  5  3  6  8  6  7  6  4  2  1  7  4  7  4
       6  3  3  4  6  4  3  5  6  4  2  5  5  9  6  6  2  8  5  4
       7  6  4  4  2  3  6  5  7  3  1  6  9  4  1  3  5  9  4  7
       4  6  6  2  7  3  7  8  5  5  8  4  4  3  7  5  3  4  6  5
       7  7  4  3  3  9  9  7  4  6  2  4  5  3  9  5  3  5  4  6
       3  5  3  6  7  9  4  3  3  5  7  4  6  6  5  4  4  5  3  4
       4  3  5  5  8  7  3  6  5  7  9  5  8  6  4  3  6  2  6  6
       4  6  6  7  5  4  3  5  1  8  4  7  5  5  7  5  7  5  5  4
       4  8  4  3  5  6  5  7  8  6  2  4  9  5  3  2  9  6  8  5
       2  6  2  1  8  2  7  9  4  2  7  4  7  7  9  5  5  9  1  3
       1  7  5  6  4  4  4  5  7  4  4  4  4  7  5  4  7  1  7  5
       6  6  5  9  4  5  6  8  5  9  7  4  7  4  8  6  4  3  6  5
       5  4  9  4  3  3  6  6  5  5  7  7  4  4  4  5  8  3  4  5
       2  2  6  1  7  9  5  4  7  3  9  5  8  4  6  4  3  4  1  4
       7  1  6  3  5  5  4  3  5  2  5  4  6  5  4  9  6  6  4  5
       7  7  2  5  3  3  4  6  7  9  5  7  4  8  3  6  3  5  7  5
       6  4  8  9  4  5  5  2  5  6  9  1  3  3  4  6  5  8  6  5
       1  7  5  4  6  4  3  6  5  8  1  7  6  4  1  4  7  6  4  5
       6  4  4  7  5  6  4  6  4  5  3  5  5  5  4  5  5  5  5  6];
  
    dpq=zeros(m,m);
    for p=1:m
       for q=p:m
        dpq(p,q)=sqrt((xp(p)-xp(q))^2+(yp(p)-yp(q))^2);
        dpq(q,p)=dpq(p,q);
       end
    end
    
    model.w=w;
    model.dpq=dpq;
    model.n=n;
    model.m=m;
    model.xp=xp;
    model.yp=yp;
    
end   

function sol1=CreateRandomSolution(model)

    m=model.m;
    
    sol1=randperm(m);
    
end

function p=ParseSolution(s,model)
    
    n=model.n;
    
    p=s(1:n);
    
end

function z=Cost(s,model)
        
    global NFE;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;
    
    p=ParseSolution(s,model);
    
    n=model.n;
    z=0;
    w=model.w;
    d=model.dpq;
    for i=1:n
        for j=i:n
            z=z+w(i,j)*d(p(i),p(j));
        end
    end
    
    
end

function [y1,y2]=SinglePointCrossover(x1,x2)
      
    n=numel(x1);
    cutpoint=randi([1,n-1]);
    
    x11=x1(1:cutpoint);
    x12=x1(cutpoint+1:end);
    
    x21=x2(1:cutpoint);
    x22=x2(cutpoint+1:end);
    
    r1=intersect(x11,x22);
    r2=intersect(x12,x21);
    
    x11(ismember(x11,r1))=r2;
    x21(ismember(x21,r2))=r1;
    
    y1=[x11 x22];
    y2=[x21 x12];
    
      
end

function t=Mutation(s)
    
    %tour=ParseSolution(s,model);
    
    M=randi([1,3]);
    
    switch M
        case 1
            t=mswap(s);
        case 2
            t=mrev(s);
        case 3
            t=minser(s);
    end
    
end

function t=mswap(tour)
    
    n=numel(tour);
    
    i=randsample(n,2);
    i1=i(1);
    i2=i(2);
    
    t=tour;
    t([i1 i2])=tour([i2 i1]);
    
end

function t=mrev(tour)
    
    n=numel(tour);
    
    i=randsample(n,2);
    i1=min(i);
    i2=max(i);
    
    t=tour;
    t(i1:i2)=tour(i2:-1:i1);

end

function t=minser(tour)
    
    n=numel(tour);
    
    i=randsample(n,2);
    i1=i(1);
    i2=i(2);
    
    if i1>i2       
        t=[tour(1:i2-1) tour(i2+1:i1) tour(i2) tour(i1+1:end)];
    else
        t=[tour(1:i1) tour(i2) tour(i1+1:i2-1) tour(i2+1:end)];
    end

end

function i=RWS(P) 
    
    r=rand;
    c=cumsum(P);
    i=min(find(r<c,1,'first'));
    
end


