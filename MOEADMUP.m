function MOEADMUP(Global)
% <algorithm> <H-N>
% MOEA/D-MUP
% operator      --- DE

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    delta=0.9;
    nr=20;
    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    T  = 20;

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);

    %% Generate random population
    Population = Global.Initialization();
    Z          = min(Population.objs,[],1);
    W1=W';
    [value,~]=min(W1);
    L=W1-repmat(value,Global.M,1);
    zz=L';
    %% Optimization
    while Global.NotTermination(Population)
        if Global.M==2
            mi1=10000;
            mi2=10000;
            for i=1:Global.N
                x=Population(i).obj-Z;
                v1=x(2);
                if v1<mi1
                   mi1=v1;
                   index1=i;
                end
                v2=x(1);
                if v2<mi2
                   mi2=v2; 
                   index2=i;
                end   
            end
            p1=Population(index1).obj;
            p2=Population(index2).obj;
            pp=[p1;p2];
            Z=min(pp);
            Znd=max(pp);
        else
            mi1=10000;
            mi2=10000;
            mi3=10000;
            for i=1:Global.N
                x=Population(i).obj;
                v1=sqrt(x(2)*x(2)+x(3)*x(3));
                if v1<mi1
                   mi1=v1;
                   index1=i;
                end
                v2=sqrt(x(1)*x(1)+x(3)*x(3));
                if v2<mi2
                   mi2=v2; 
                   index2=i;
                end   
                v3=sqrt(x(1)*x(1)+x(2)*x(2));
                if v3<mi3
                   mi3=v3; 
                   index3=i;
                end  
            end
            p1=Population(index1).obj;
            p2=Population(index2).obj;
            p3=Population(index3).obj;
            pp=[p1;p2;p3];
            Z=min(pp);
            Znd=max(pp);
        end

        for i = 1 : Global.N
            % Choose the parents
            s=rand;
            if s < delta
                P = B(i,randperm(size(B,2)));
            else
                P = randperm(Global.N);
            end

            % Generate an offspring
            Offspring = Global.Variation(Population([i P(1:2)]),1,@DE);
            Z = min(Z,Offspring.obj);
            
            if s<delta
                S=repmat(Z,T,1)+zz(P,:).*(Znd-Z);
                g_old = max((Population(P).objs-S)./W(P,:),[],2);                 
                g_new = max((repmat(Offspring.obj,T,1)-S)./W(P,:),[],2);
            else
                S=repmat(Z,Global.N,1)+zz(P,:).*(Znd-Z);
                g_old = max((Population(P).objs-S)./W(P,:),[],2);                 
                g_new = max((repmat(Offspring.obj,Global.N,1)-S)./W(P,:),[],2);
            end
            Population(P(find(g_old>=g_new,nr))) = Offspring;
        end
    end
end