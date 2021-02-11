function [S]=calMatrixS(q,dim,structure,beta)
switch structure.geometry
    case 'plane'
        a=beta(2);
        B=eye(dim)*sqrt(2)/a;
        B(1,1)=1/a;
        C=eye(dim)*sqrt(2/a);
        C(1,1)=sqrt(1/a);
        for i=1:dim
            for j=1:dim
                if (i==j)
                    S(i,j)=getSs(i-1,q,a);
                else
                    S(i,j)=0;
                end
               
            end
        end
        S=B*S;
        
        return;
    case 'sphere'
        beta(2)=beta(2)/2;
        a=beta(2);
        V_total=(a.^3)*pi*4/3;       
        [roots,alpha,index_row]=mati.MULTIPRO.calRoots(dim,beta);
        S=zeros(dim,dim);
       
        for k=1:dim
            bb=2*pi*q*a;
            aa=roots(k);
            nn=index_row(k)-1;
            S(k,k)=4*pi*(sqrt(-1).^nn)*(V_total.^(-0.5))*alpha(k)*(a.^3)/...
                (aa.^2-bb.^2)*(aa*mati.MULTIPRO.besselj_convert(nn,bb)*mati.MULTIPRO.besselj_convert(nn+1,aa)-bb*mati.MULTIPRO.besselj_convert(nn,aa)*mati.MULTIPRO.besselj_convert(nn+1,bb));
        end       
       
        return;
end

end

function Ss=getSs(k,q,a)
if k==0
    Ss=sin(q*a*pi)*exp(sqrt(-1)*q*a*pi)/pi/q;
elseif mod(k,2)==0 % k~=0,even
    Ss=2*a*exp(sqrt(-1)*q*a*pi)*(2*q*a*pi)*sin(pi*q*a)/(2*pi*q*a-k*pi)/(2*pi*q*a+k*pi);
else % k~=0, odd
    Ss=sqrt(-1)*2*a*exp(sqrt(-1)*q*a*pi)*(2*q*a*pi)*cos(pi*q*a)/(2*pi*q*a-k*pi)/(2*pi*q*a+k*pi); 
end
end






