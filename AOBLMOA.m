function [bestCost,cg_curve,bestPos,position_history,fitness_history,Trajectories]=AOBLMOA(nPop,nPopf,nc,MaxIt,LowerBound,UpperBound,ProblemSize,fobj)
g=0.9;                      % 惯性权重 Inertia Weight
gmax = 0.9;                 % 惯性权重最大值 Maximum Inertia Weight
gmin = 0.4;                 % 惯性权重最小值 Minimum Inertia Weight
gdamp=0.9;                  % 惯性权重阻尼比 Inertia Weight Damping Ratio
a1=1.0;                     % 个人学习因子 Personal Learning Coefficient
a2=1.5; a3=1.5;             % 全局学习因子 Global Learning Coefficient
beta=2;                     % 视距因子 Distance sight Coefficient
dance=5;                    % 交配舞蹈行为 Nuptial Dance
fl=1;                       % 随机飞行系数 Random flight
dance_damp=0.8;             % 舞蹈阻尼比 Damping Ratio
fl_damp=0.99;               % 配合参数 Mating Parameters
dim = ProblemSize;
ProblemSize = [1 ProblemSize];
alpha=0.1;
delta=0.1;
position_history=zeros(nPop,MaxIt,dim);
Trajectories=zeros(nPop,MaxIt);
fitness_history=zeros(1,MaxIt);
% 速度限制 Velocity Limits
VelMax=0.1*(UpperBound-LowerBound); VelMin=-VelMax;
%% 种群初始化 Initialization
empty_mayfly.Position=[];
empty_mayfly.Cost=[];
empty_mayfly.Velocity=[];
empty_mayfly.Best.Position=[];
empty_mayfly.Best.Cost=[];
Mayfly=repmat(empty_mayfly,nPop,1);   % Males
Mayflyf=repmat(empty_mayfly,nPopf,1); % Females
GlobalBest.Cost=inf;
funccount=0;
for i=1:nPop
    % 初始化父性群体 Initialize Position of Males
    Mayfly(i).Position=unifrnd(LowerBound,UpperBound,ProblemSize);
    % 初始化速度 Initialize Velocity
    Mayfly(i).Velocity=zeros(ProblemSize);
    % 估价 Evaluation
    Mayfly(i).Cost=fobj(Mayfly(i).Position);
    % 更新个人最优位置 Update Personal Best
    Mayfly(i).Best.Position=Mayfly(i).Position;
    Mayfly(i).Best.Cost=Mayfly(i).Cost;
    funccount=funccount+1;
    % 更新全局最优位置 Update Global Best
    if Mayfly(i).Best.Cost<GlobalBest.Cost
        GlobalBest=Mayfly(i).Best;
    end
end
for i=1:nPopf
    % 初始化母性群体位置 Initialize Position of Females
    Mayflyf(i).Position=unifrnd(LowerBound,UpperBound,ProblemSize);
    Mayflyf(i).Velocity=zeros(ProblemSize);
    Mayflyf(i).Cost=fobj(Mayflyf(i).Position);
    funccount=funccount+1;
end
BestSolution=zeros(MaxIt,1);
%% Mayfly Main Loop
for it=1:MaxIt
    G2=2*rand()-1; % Eq.26
    G1=2*(1-(it/MaxIt));  % Eq.25
    to = 1:dim;
    u = .0265;
    r0 = 10;
    r = r0 +u*to;
    omega = .005;
    phi0 = 3*pi/2;% Eq.21
    phi = -omega*to+phi0;% Eq.20
    x = r .* sin(phi);  % Eq.17
    y = r .* cos(phi); % Eq.18
    QF=it^((2*rand()-1)/(1-MaxIt)^2); % Eq.24
    for i = 1:nPopf
        xmf(i,:) = Mayflyf(i).Position;
    end
    for i = 1:nPop
        xm(i,:) = Mayfly(i).Position;
    end
    
    for i=1:nPopf
        % 更新母性位置 Update Females
        e=unifrnd(-1,+1,ProblemSize);
        rmf=(Mayfly(i).Position-Mayflyf(i).Position);
        if Mayflyf(i).Cost>Mayfly(i).Cost
            Mayflyf(i).Velocity = g*Mayflyf(i).Velocity ...
                +a3*exp(-beta.*rmf.^2).*(Mayfly(i).Position-Mayflyf(i).Position);
            Mayflyf(i).Velocity = max(Mayflyf(i).Velocity,VelMin);
            Mayflyf(i).Velocity = min(Mayflyf(i).Velocity,VelMax);
            Mayflyf(i).Position = Mayflyf(i).Position + Mayflyf(i).Velocity;
        else
            if it<=(2/3)*MaxIt
                Mayflyf(i).Position = GlobalBest.Position * (1- it/MaxIt) + (mean(Mayflyf(i).Position) - GlobalBest.Position * rand()); 
            else
                Mayflyf(i).Position = (GlobalBest.Position - mean(Mayflyf(i).Position))*alpha - rand + ((UpperBound-LowerBound)*rand+LowerBound)*delta;
            end
        end

        % 位置限制 Position Limits
        Mayflyf(i).Position = max(Mayflyf(i).Position,LowerBound);
        Mayflyf(i).Position = min(Mayflyf(i).Position,UpperBound);
        % 估价 Evaluation
        Mayflyf(i).Cost = fobj(Mayflyf(i).Position);
        funccount=funccount+1;
    end
    for i=1:nPop
        % 更新父性群体 Update Males
        rpbest=(Mayfly(i).Best.Position-Mayfly(i).Position);
        rgbest=(GlobalBest.Position-Mayfly(i).Position);
        e=unifrnd(-1,+1,ProblemSize);
        % 更新速度 Update Velocity
        if Mayfly(i).Cost>GlobalBest.Cost
            Mayfly(i).Velocity = g*Mayfly(i).Velocity ...
                +a1*exp(-beta.*rpbest.^2).*(Mayfly(i).Best.Position-Mayfly(i).Position) ...
                +a2*exp(-beta.*rgbest.^2).*(GlobalBest.Position-Mayfly(i).Position);
        % 应用速度限制 Apply Velocity Limits
            Mayfly(i).Velocity = max(Mayfly(i).Velocity,VelMin);
            Mayfly(i).Velocity = min(Mayfly(i).Velocity,VelMax);
        % 更新位置 Update Position
            Mayfly(i).Position = Mayfly(i).Position + Mayfly(i).Velocity;
        else
            if it<=(2/3)*MaxIt
                Mayfly(i).Position = GlobalBest.Position .* Levy(dim) + Mayfly((floor(nPop*rand()+1))).Position + (y-x)*rand;
            else
                Mayfly(i).Position = QF*GlobalBest.Position - (G2*Mayfly(i).Position*rand()-G1.*Levy(dim)+rand*G2);
            end
        end
        % 位置限制 Position Limits
        Mayfly(i).Position = max(Mayfly(i).Position,LowerBound);
        Mayfly(i).Position = min(Mayfly(i).Position,UpperBound);
        % 估价 Evaluation
        Mayfly(i).Cost = fobj(Mayfly(i).Position);
        funccount=funccount+1;
        % 更新个体最优 Update Personal Best
        if Mayfly(i).Cost<Mayfly(i).Best.Cost
            Mayfly(i).Best.Position=Mayfly(i).Position;
            Mayfly(i).Best.Cost=Mayfly(i).Cost;
            % 更新全局最优 Update Global Best
            if Mayfly(i).Best.Cost<GlobalBest.Cost
                GlobalBest=Mayfly(i).Best;
            end
        end
    end
    [~, SortMayflies]=sort([Mayfly.Cost]);
    Mayfly=Mayfly(SortMayflies);
    [~, SortMayflies]=sort([Mayflyf.Cost]);
    Mayflyf=Mayflyf(SortMayflies);
    % MATE
    MayflyOffspring=repmat(empty_mayfly,nc/2,2);
    for k=1:nc/2
        % 选择父母 Select Parents
        i1=k;
        i2=k;
        p1=Mayfly(i1);
        p2=Mayflyf(i2);
        % 应用交叉 Apply Crossover
        [MayflyOffspring(k,1).Position, MayflyOffspring(k,2).Position]=Crossover(p1.Position,p2.Position,LowerBound,UpperBound);
        % 评估后代 Evaluate Offsprings
        MayflyOffspring(k,1).Cost=fobj(MayflyOffspring(k,1).Position);
        antix =  rand().*((UpperBound + LowerBound)-MayflyOffspring(k,1).Position);
        % 位置限制 Position Limits
        antix=max(antix,LowerBound);
        antix=min(antix,UpperBound);
        if fobj(antix)<MayflyOffspring(k,1).Cost
            MayflyOffspring(k,1).Cost = fobj(antix);
            MayflyOffspring(k,1).Position = antix;
            if MayflyOffspring(k,1).Cost<GlobalBest.Cost
                GlobalBest=MayflyOffspring(k,1);
            end
        end
        MayflyOffspring(k,2).Cost=fobj(MayflyOffspring(k,2).Position);
        antix =  rand().*((UpperBound + LowerBound)-MayflyOffspring(k,2).Position);
        % 位置限制 Position Limits
        antix=max(antix,LowerBound);
        antix=min(antix,UpperBound);
        if fobj(antix)<MayflyOffspring(k,2).Cost
            MayflyOffspring(k,2).Cost = fobj(antix);
            MayflyOffspring(k,2).Position = antix;
            if MayflyOffspring(k,2).Cost<GlobalBest.Cost
                GlobalBest=MayflyOffspring(k,2);
            end
        end
        MayflyOffspring(k,1).Best.Position = MayflyOffspring(k,1).Position;
        MayflyOffspring(k,1).Best.Cost = MayflyOffspring(k,1).Cost;
        MayflyOffspring(k,1).Velocity= zeros(ProblemSize);
        MayflyOffspring(k,2).Best.Position = MayflyOffspring(k,2).Position;
        MayflyOffspring(k,2).Best.Cost = MayflyOffspring(k,2).Cost;
        MayflyOffspring(k,2).Velocity= zeros(ProblemSize);
    end
    MayflyOffspring=MayflyOffspring(:);
    % 创建合并群体 Create Merged Population
    split=round((size(MayflyOffspring,1))/2);
    newmayflies=MayflyOffspring(1:split);
    Mayfly=[Mayfly
        newmayflies];
    newmayflies=MayflyOffspring(split+1:size(MayflyOffspring,1));
    Mayflyf=[Mayflyf
        newmayflies]; 
    [~, SortMayflies]=sort([Mayfly.Cost]);
    Mayfly=Mayfly(SortMayflies);
    Mayfly=Mayfly(1:nPop); % 保持父性最优 Keep best males
    [~, SortMayflies]=sort([Mayflyf.Cost]);
    Mayflyf=Mayflyf(SortMayflies);
    Mayflyf=Mayflyf(1:nPopf); % 保持母性最优 Keep best females
    BestSolution(it)=GlobalBest.Cost;
    g=gmax - ((gmax-gmin)/MaxIt)*it;
    dance = dance*dance_damp;
    fl = fl*fl_damp;
    for i = 1:nPop
        Positions(i,:) = Mayfly(i).Position;
    end
    for i = 1:size(Positions,1)
        position_history(i,it,:) = Positions(i,:);
        Trajectories(:,it) = Positions(:,1);
        fitness_history(i,it) = fobj(Positions(i,:)); 
    end
    cg_curve(it) = GlobalBest.Cost;
    bestCost = GlobalBest.Cost;
    bestPos = GlobalBest.Position;
end
function [off1, off2]=Crossover(x1,x2,LowerBound,UpperBound)
    L=unifrnd(0,1,size(x1));
    off1=L.*x1+(1-L).*x2;
    off2=L.*x2+(1-L).*x1;
    % 位置限制 Position Limits
    off1=max(off1,LowerBound); off1=min(off1,UpperBound);
    off2=max(off2,LowerBound); off2=min(off2,UpperBound);
end
%%
function y=Mutate(x,mu,LowerBound,UpperBound)
    nVar=numel(x);
    nmu=ceil(mu*nVar);
    j=randsample(nVar,nmu);
    sigma(1:nVar)=0.1*(UpperBound-LowerBound);
    y=x;
    y(j)=x(j)+sigma(j).*(randn(size(j))');
    y=max(y,LowerBound); y=min(y,UpperBound);
end
function o=Levy(d)
    bet=1.5;
    sigma=(gamma(1+bet)*sin(pi*bet/2)/(gamma((1+bet)/2)*bet*2^((bet-1)/2)))^(1/bet);
    w=randn(1,d)*sigma;v=randn(1,d);step=w./abs(v).^(1/bet);
    o=step;
end

end
