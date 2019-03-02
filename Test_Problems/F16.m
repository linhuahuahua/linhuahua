function varargout = F16(Operation,Global,input)
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
            A=2;
            B=2;
            C=2;
            D=2;
            E=3;
            F=2;
            
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
            %PF采样
            input=44;
            t1=0:1/(input-1):1;
            s=0;
            for i=1:length(t1)
                t2=0:1/(input-1):1-t1(i);
                t3=1-t1(i)-t2;
                f(s+1:s+length(t2),1)=repmat((t1(i)^2),length(t2),1)';
                f(s+1:s+length(t2),2)=(t2.^2)';
                f(s+1:s+length(t2),3)=(t3.^2)';
                s=s+length(t2);
            end
            
            %在x-y面上等间距选取
%             input=75;
%             x=0:1/(input-1):1;
%             a=[];
%             b=[];
%             c=[];
%             for i=1:input
%                y=0: 1/(input-1):(1-x(i)^0.5)^2;
%                z=(1-x(i)^0.5-y.^0.5).^2;
%                a=[a;repmat(x(i),length(y),1)];
%                b=[b;y'];
%                c=[c;z'];
%             end
%             f(:,1)=a;
%             f(:,2)=b;
%             f(:,3)=c;
            
            %近似PF的一种方式
%             input=126;
%             x=0:1/(input/3-1):1;
%             a=[];
%             b=[];
%             c=[];
%             for i=1:input/3
%                y=0: 1/(input/3-1):(1-x(i)^0.5)^2;
%                z=(1-x(i)^0.5-y.^0.5).^2;
%                a=[a;repmat(x(i),length(y),1)];
%                b=[b;y'];
%                c=[c;z'];
%             end
%             s=length(a);
%             f(1:s,1)=a;
%             f(1:s,2)=b;
%             f(1:s,3)=c;
%             
%             y=0:1/(input/3-1):1;
%             a=[];
%             b=[];
%             c=[];
%             for i=1:input/3
%                z=0: 1/(input/3-1):(1-y(i)^0.5)^2;
%                x=(1-y(i)^0.5-z.^0.5).^2;
%                b=[b;repmat(y(i),length(z),1)];
%                c=[c;z'];
%                a=[a;x'];
%             end
%             ss=length(a);
%             f(s+1:s+ss,1)=a;
%             f(s+1:s+ss,2)=b;
%             f(s+1:s+ss,3)=c;
%             
%             z=0:1/(input/3-1):1;
%             a=[];
%             b=[];
%             c=[];
%             for i=1:input/3
%                x=0: 1/(input/3-1):(1-z(i)^0.5)^2;
%                y=(1-z(i)^0.5-x.^0.5).^2;
%                c=[c;repmat(z(i),length(x),1)];
%                a=[a;x'];
%                b=[b;y'];
%             end
%             sss=length(a);
%             f(s+ss+1:s+ss+sss,1)=a;
%             f(s+ss+1:s+ss+sss,2)=b;
%             f(s+ss+1:s+ss+sss,3)=c;
            
            varargout = {f};
            scatter3(f(:,1),f(:,2),f(:,3),'.');
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

