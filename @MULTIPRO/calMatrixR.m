function [R]=calMatrixR(t,dim,structure,beta)
switch structure.geometry
    case 'plane'
        alpham=0:(dim-1);
        R=eye(dim,dim).*exp(-alpham.*alpham*pi*pi*beta(4)*t/beta(2)/beta(2));
        return;
    case 'sphere'           
        beta(2)=beta(2)/2;
        [roots,alpha,index_row]=mati.MULTIPRO.calRoots(dim,beta);       
        R=eye(dim,dim).*exp(-roots.*roots*beta(4)*t/beta(2)/beta(2));
        return;
end

end




