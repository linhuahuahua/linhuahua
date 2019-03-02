function varargout = F14(Operation,Global,input)
% <problem> <F1_16>
% operator --- EAreal    
    switch Operation
        case 'init'
            Global.M        = 3;
            Global.D        = 32;
            Global.lower    = [zeros(1,Global.M-1) ones(1,Global.D-2)*(-1)];
            Global.upper    = ones(1,Global.D);
            Global.operator = @DE;
            
            PopDec    = unifrnd(repmat(Global.lower,input,1),repmat(Global.upper,input,1));
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,n]  = size(PopDec);
            m  = Global.M;
            A=4;
            B=1;
            C=2;
            D=1;
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
            input=1000;
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
