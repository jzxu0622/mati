function y=besselj_convert(n,x)
if length(n)>1

for i=1:length(n)
    y(i)=sqrt(pi./x(i)/2).*besselj(n(i)+0.5,x(i));
    if n(i)==0 && x(i)==0
        y(i)=1;
    elseif n(i)>0 && x(i)==0
        y(i)=0;
    end
end
else
    y=sqrt(pi./x/2).*besselj(n+0.5,x);
    i_zeros=find(x==0);
    if n==0
        y(i_zeros)=1;
    elseif n>0
        y(i_zeros)=0;
    end
end
end