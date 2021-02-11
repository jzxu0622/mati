function [roots,alpha,index_row]=calRoots(dim,beta)
if exist(['roots_' num2str(dim) '.mat'],'file')
    load(['roots_' num2str(dim) '.mat'],'roots','index_row');
else
    tol = 0.001 ;
    xmax = 100 ;
    lmax = 20 ;
    x = 0*tol:tol:xmax ; 
    roots=zeros(lmax,lmax);
    for n=0:(lmax-1)
        y=besselj_convert_d(n,x);
        x=x+tol;
        y_plus=besselj_convert_d(n,x);
        x_zero = x(y.*y_plus<0) ; 
        if length(x_zero) > lmax
            roots(n+1,1:lmax) = x_zero(1:lmax) ; 
        else
            roots(n+1,1:length(x_zero))=x_zero;
        end
        index_row(n+1,1:lmax)=n+1;        
    end
    roots=roots(:);
    index_row=index_row(:);    
    [roots,I]=sort(roots);
    roots=[0; roots];
    index_row=index_row(I);
    index_row=[1; index_row];
    save(['roots_' num2str(dim) '.mat'],'roots','index_row');
end
roots=roots(1:dim);
index_row=index_row(1:dim);
for i=1:dim
    alpha(i)=1./sqrt((2.*pi./(2.*(index_row(i)-1)+1)).*(beta(2).^3).*(mati.MULTIPRO.besselj_convert(index_row(i)-1,roots(i)).^2-mati.MULTIPRO.besselj_convert(index_row(i)-2,roots(i)).*mati.MULTIPRO.besselj_convert(index_row(i),roots(i))));
end
alpha(1)=(4*pi*(beta(2).^3)/3).^(-0.5);
alpha=alpha';
end

function output=besselj_convert_d(n,x)
output=n./(2*n+1).*mati.MULTIPRO.besselj_convert(n-1,x)-(n+1)./(2*n+1).*mati.MULTIPRO.besselj_convert(n+1,x);
end