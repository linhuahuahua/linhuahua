function varargout = F1(Operation,Global,input)
% <problem> <F1_16>
% Scalable Test Problems for Evolutionary Multi-Objective Optimization
% operator --- DE   

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 31;
            Global.lower    = [0 ones(1,Global.D-1)*(-1)];
            Global.upper    = ones(1,Global.D);
            Global.operator = @DE; 
            
            PopDec    = unifrnd(repmat(Global.lower,input,1),repmat(Global.upper,input,1));
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,n]  = size(PopDec);
            m  = Global.M;
            A=1;
            B=4;
            C=2;
            D=4;
            E=3;
            F=1;
      
            x1=PopDec(:,1:m-1);
            x2=PopDec;
            
            y(:,1)=1-x1(:,1);
            for i=2:m-1
                y(:,i)=(1-x1(:,i)).*prod(x1(:,1:i-1),2);
            end
            y(:,m)=prod(x1(:,1:m-1),2);
            
            h=y.^F;

            L=distance_ratio(y);
            b_B=position_dependent_scale(m,L,B);
            t=x2-0.9*b_B.*cos(E*pi*L+(n+2)*(1:n)*pi/2/n);
            tt=zeros(N,m);
            for i=1:m
               J=m+i-1:m:floor((n+1-i)/m)*m-1+i; 
               tt(:,i)=sum(abs(t(:,J)).^C,2);
            end
            
            b_D=position_dependent_scale(m,L,D);
            g=(A*b_D+1)/size(J,2).*tt;

            PopObj = g+h;
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            input=500;
            f = UniformPoint(input,Global.M);
            varargout = {f};
    end
end

function Output=distance_ratio(y)
    [n,m]=size(y);
    R=1/power(m*(m-1),1/2);
    N=ones(m,m)*(1/(m*(m-1)));
    N(logical(eye(m)))=-1/m;
    L=zeros(n,1);
    for i=1:n
        L(i)=max(N*y(i,:)')/R/R;
    end
    
    Output=L;
end

function Output=position_dependent_scale(m,L,a)
    b=power(sin(pi/2*power(L,m-1)),a);
    Output=b;
end
