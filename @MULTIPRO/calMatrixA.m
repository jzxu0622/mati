function [A]=calMatrixA(q,dim,structure,beta)
range=25;  % set up the intergal range for the geometry of sphere
switch structure.geometry
    case 'plane'
        a=beta(2);
        B=eye(dim)*sqrt(2)/a;
        B(1,1)=1/a;
        C=eye(dim)*sqrt(2/a);
        C(1,1)=sqrt(1/a);
        for i=1:dim
            for j=1:dim              
                A(i,j)=getSs(abs(i-j),q,a)/2+getSs(i+j-2,q,a)/2;
            end
        end
        
        A=ctranspose(C)*A*C;
     
        return;
    case 'sphere'    
        beta(2)=beta(2)/2;
        a=beta(2);
        [roots,alpha,index_row]=mati.MULTIPRO.calRoots(dim,beta);   
        A=zeros(dim,dim); 
     
        for k=1:dim            
            for m=1:dim
                if k>m
                    A(k,m)=A(m,k);
                    continue;
                end
                i_n=index_row(k)-1;
                i_k=index_row(m)-1;                
                A(k,m)=2*pi*alpha(k)*alpha(m)*(roots(k).^i_n)*(roots(m).^i_k)*(pi.^1.5)*(a.^3);
                t1=0;
                for l=abs(i_n-i_k):(i_n+i_k)  
                    if mod((l+i_n+i_k),2)==1
                        continue;
                    end
                    t1=t1+cal_l(q,a,l,i_n,i_k,roots,k,m,range);                    
                end
                A(k,m)=A(k,m)*t1;
              
            end            
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

function l=cal_l(q,a,l,n,k,roots,u,v,range)
c=(n+k+l)/2;
if l==n+k
    al=gamma(n+k+1).^2*gamma(2*n+1)*gamma(2*k+1)/gamma(2*n+2*k+1)/gamma(n+1)/gamma(n+1)/gamma(k+1)/gamma(k+1);
else
    al=(2*l+1)*gamma(c+1)*gamma(c-k+0.5)*gamma(c-n+0.5)*gamma(c-l+0.5)/...
    gamma(c+1.5)/gamma(c-n+1)/gamma(c-k+1)/gamma(c-l+1)/2/pi;
end
l=sqrt(-1).^l*((0.5).^(2*c+2))*((2*pi*q*a).^l)*al*cal_p(q,a,l,n,k,roots,u,v,range);

end

function p=cal_p(q,a,l,n,k,roots,u,v,range)
matrix_dim=0:range;
[i_p,i_m,i_s]=meshgrid(matrix_dim,matrix_dim,matrix_dim);
m_p=(-0.25*(roots(u).^2)).^i_p./gamma(i_p+1)./gamma(n+i_p+1.5).*...
        (-0.25*(roots(v).^2)).^i_m./gamma(i_m+1)./gamma(k+i_m+1.5).*...
        (-0.25*((2*pi*q*a).^2)).^i_s./gamma(i_s+1)./gamma(l+i_s+1.5)./...
        (n+k+l+2.*i_p+2.*i_m+2.*i_s+3);

p=sum(m_p(:));
end




