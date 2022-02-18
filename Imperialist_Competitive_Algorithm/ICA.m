clc;
clear;
close all;
%% Problem Parameters

cost=@(x) Sphere(x);

nvar=5;
varsize=[1 nvar];

varmin=-10;
varmax= 10;

%% ICA Parameters

maxit=500;

npop=50;
nEmp=round(0.4*npop);   % number of imperialists/empires
nCol=npop-nEmp;         % total number of colonies

beta=2;      % assimilation coefficient

pRev=0.2;    % revolution probability
mu=0.05;     % revolution rate

zeta=0.1;    % colonies cost coefficient
alpha=1;     % selection pressure   
pRevolution=0.5;

%% Share Settings

global ProblemSettings;
ProblemSettings.cost=cost;
ProblemSettings.nvar=nvar;
ProblemSettings.varsize=varsize;
ProblemSettings.varmin=varmin;
ProblemSettings.varmax=varmax;

global ICAsettings;
ICAsettings.maxit=maxit;
ICAsettings.npop=npop;
ICAsettings.nEmp=nEmp;
ICAsettings.nCol=nCol;
ICAsettings.beta=beta;
ICAsettings.pRev=pRev;
ICAsettings.mu=mu;
ICAsettings.zeta=zeta;
ICAsettings.alpha=alpha;
ICAsettings.pRevolution=pRevolution;

%% Setup ICA

emp=CreateEmpires();

bestcost=zeros(maxit,1);

%% Main Loop

for it=1:maxit
    
    % Assimilation
    emp=AssimilationColonies(emp);
    
    % Revolution
    emp=Revolution(emp);
    
    % Intra-Empire Competition
    emp=IntraEmpireCompetition(emp);
    
    % Update Total Cost of Empire
    emp=UpdateTotalCost(emp);
    
    % Inter-Empire Competition
    emp=InterEmpireCompetition(emp);
        
    % Update Best Solution Ever Found
    imp=emp.Imp;
    [~,i]=min([imp.cost]);
    
    bestsol=emp(i).Imp;
    
    % Update Best Cost
    bestcost(it)=bestsol.cost;
    
    % Display
    disp(['Iteration: ' num2str(it) ' Cost: ' num2str(bestcost(it))]);
    
end

%% Functions

function z=Sphere(x)
    z=sum(x.^2);
end

function pop=Sorting(pop)
    
    Costs=[pop.cost];
    [~,so]=sort(Costs);
    pop=pop(so);
    
end

function i=Roulett(p)
    
    r=rand;
    c=cumsum(p);
    
    i=find(r<=c,1,'first');
end

function emp=UpdateTotalCost(emp)

    global ICAsettings;
    zeta=ICAsettings.zeta;
    nEmp=numel(emp);
    
    for i=1:nEmp
        
        if emp(i).nCol>0
            emp(i).TotalCost=emp(i).Imp.cost+zeta*mean([emp(i).Col.cost]);
        else
            emp(i).TotalCost=emp(i).Imp.cost;
        end
        
    end
    
end

function emp=CreateEmpires()
    
    global ICAsettings;
    global ProblemSettings;
    
    cost=ProblemSettings.cost;
    varsize=ProblemSettings.varsize;
    varmin=ProblemSettings.varmin;
    varmax=ProblemSettings.varmax;
    
    npop=ICAsettings.npop;
    nEmp=ICAsettings.nEmp;
    nCol=ICAsettings.nCol;
    alpha=ICAsettings.alpha;
    
    country.position=[];
    country.cost=[];
    
    pop=repmat(country,npop,1);
    
    for i=1:npop
        
        pop(i).position=unifrnd(varmin,varmax,varsize);
        pop(i).cost=cost(pop(i).position);
             
    end
    
    pop=Sorting(pop);
    
    % Dividing Countries
    imp=pop(1:nEmp);
    col=pop(nEmp+1:end);
    
    Empire.Imp=[];
    Empire.Col=repmat(country,0,1);
    Empire.nCol=0;
    Empire.TotalCost=[];
    
    emp=repmat(Empire,nEmp,1);
    
    % Assign Imperialists
    for k=1:nEmp    
        emp(k).Imp=imp(k);
    end
    
    % Assign Colonies
    imp_costs=[imp.cost];
    p=exp(-alpha*imp_costs/max(imp_costs));
    p=p/sum(p);
    for j=1:nCol
        
        k=Roulett(p);
        emp(k).Col=[emp(k).Col
                    col(j)];
                
        emp(k).nCol=numel(emp(k).Col);
        
    end
    
    % Update Total Cost
    emp=UpdateTotalCost(emp);
    
    
    
end

function emp=AssimilationColonies(emp)
    
    global ICAsettings;
    global ProblemSettings;
    
    cost=ProblemSettings.cost;
    varsize=ProblemSettings.varsize;
    varmin=ProblemSettings.varmin;
    varmax=ProblemSettings.varmax;
    
    beta=ICAsettings.beta;
    nEmp=numel(emp);
    
    % Assimilation
    for i=1:nEmp  
       for j=1:emp(i).nCol
           
           emp(i).Col(j).position=emp(i).Col(j).position+...
               beta*rand(varsize).*(emp(i).Imp.position-emp(i).Col(j).position);
           
           emp(i).Col(j).position=min(max(emp(i).Col(j).position,varmin),varmax);
           
           emp(i).Col(j).cost=cost(emp(i).Col(j).position);
           
       end 
    end
 
end

function emp=Revolution(emp)

    global ICAsettings;
    global ProblemSettings;
    
    cost=ProblemSettings.cost;
    nvar=ProblemSettings.nvar;
    varsize=ProblemSettings.varsize;
    varmin=ProblemSettings.varmin;
    varmax=ProblemSettings.varmax;
    
    pRevolution=ICAsettings.pRevolution;
    mu=ICAsettings.mu;
    nEmp=numel(emp);
    
    sigma=0.03*(varmax-varmin);
    Nmu=ceil(mu*nvar);
    
    for i=1:nEmp  
       
        jj=randsample(nvar,Nmu)';
        NewPos=emp(i).Imp.position+sigma*randn(varsize);
        
        NewImp=emp(i).Imp;
        NewImp.position(jj)=NewPos(jj);
        NewImp.cost=cost(NewImp.position);
        
        if NewImp.cost<emp(i).Imp.cost
            
            emp(i).Imp.position=NewImp.position;

            emp(i).Imp.position=min(max(emp(i).Imp.position,varmin),varmax);

            emp(i).Imp.cost=cost(emp(i).Imp.position);

        end
        
        for j=1:emp(i).nCol        
            if rand<=pRevolution
       
               jj=randsample(nvar,Nmu)';
               NewPos=emp(i).Col(j).position+sigma*randn(varsize);
               
               emp(i).Col(j).position(jj)=NewPos(jj);
               
               emp(i).Col(j).position=min(max(emp(i).Col(j).position,varmin),varmax);

               emp(i).Col(j).cost=cost(emp(i).Col(j).position);

            end
        end
    end

end

function emp=IntraEmpireCompetition(emp)

    nEmp=numel(emp);

    for k=1:nEmp
        
        for i=1:emp(k).nCol
            
            if emp(k).Col(i).cost<emp(k).Imp.cost
                
                imp=emp(k).Imp;
                col=emp(k).Col(i);
                
                emp(k).Imp=col;
                emp(k).Col(i)=imp;
            end
            
        end
        
    end
 
end

function emp=InterEmpireCompetition(emp)
   
    if numel(emp)== 1
        return;
    end
    
    global ICAsettings;
    alpha=ICAsettings.alpha;
    
    TotalCost=[emp.TotalCost];
    
    [~,i]=max(TotalCost);
    Weakest=emp(i);
    
    [~,j]=max([Weakest.Col.cost]);
    WeakestCol=Weakest.Col(j);
    
    p=exp(-alpha*TotalCost/max(TotalCost));
    p(j)=0;
    p=p/sum(p);
    
    % Winner
    k=Roulett(p);
    
    emp(k).Col=[emp(k).Col
                WeakestCol];
    emp(k).nCol=numel(emp(k).Col);
    
    emp(i).Col(j)=[];
    emp(i).nCol=numel(emp(i).Col);
    
    if emp(i).nCol==0
        
        p=exp(-alpha*TotalCost/max(TotalCost));
        p(j)=0;
        p=p/sum(p);
        
        % Winner
        k=Roulett(p);
    
        emp(k).Col=[emp(k).Col
                    emp(i).Imp];
        emp(k).nCol=numel(emp(k).Col);
        
        emp(i)=[];
        
    end
   
    
end