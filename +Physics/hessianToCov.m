% reference: "Easy computation of the Bayes factor to fully quantify Occam’s razor in least‑squares ftting and to guide actions"
function [cov,ncov,minflag]=hessianToCov(hessian,fval,num_signal,num_par)
if ~isempty(hessian)|| sum(nonzeros(isnan(hessian))) || sum(nonzeros(isinf(hessian)))
    if det(-hessian)>0
        cov=inv(hessian);
%         [U,D,V]=svd(hessian);
%         inv_D=inv(D);
%         inv_D(inv_D>1e3)=0;
%         cov2=(V*inv_D*U');
    elseif abs(det(-hessian))<1e-2
        cov=pinv(hessian);
    else
        [L , DMC , P , D ] = more_sorensen (hessian); %P*(A + E)*P ' = L*DMC*L' 
        inv_D=inv(DMC);
        cov=(L*inv_D*L');        
    end  
    
%    cov=cov*fval/(num_signal-num_par);
    minflag=1;
    for ii=1:num_par
        for jj=1:num_par
            if cov(ii,ii)<0||cov(jj,jj)<0
                minflag=0;
                ncov(ii,jj)=0;
            else
                ncov(ii,jj)=cov(ii,jj)/sqrt(cov(ii,ii)*cov(jj,jj));
            end
        end
    end
else
    cov=eye(num_par);
    ncov=eye(num_par);
    minflag=0;
end
end

function [L , DMC , P , D ] = more_sorensen (A , delta )
% more_sorensen More and Sorensen modified Cholesky algorithm based on LDL ' factorization .
% [L D,P,D0] = more_sorensen (A, delta ) computes the modified Cholesky factorization P*(A + E)*P ' = L*D*L',
% where P is a permutation matrix , L is unit lower triangular , and D is block diagonal and positive
% definite with 1 -by -1 and 2 -by -2 diagonal blocks . Thus A+E is symmetric positive definite , but E is
% not explicitly computed . Also returned is a block diagonal D0 such that P*A*P ' = L*D0*L '. If A is
% sufficiently positive definite then E = 0 and D = D0.The algorithm sets the smallest eigenvalue of D
% to the tolerance delta , which defaults to eps.
% The LDL ' factorization is computed using a symmetric form of rook pivoting proposed by Ashcraft , Grimes
% and Lewis .
% This code is a modification of an existing code of Cheng and Higham , altered to use the modified
% Cholesky algorithm of More and Sorensen rather than that of Cheng and Higham ; the original code is
% available here : https :// github .com/ higham / modified -cholesky .
% Reference :
% J. J. More and D. C. Sorensen . On the use of directions of negative curvature in a modified Newton
% method . Mathematical Programming , 16(1) :1 -20 , 1979.
% Authors : Bobby Cheng and Nick Higham , 1996; revised 2015.
% Modified by Thomas McSweeney , 2017.
if ~ ishermitian ( A ) , error ( 'Must supply symmetric matrix . ') ,end
if nargin < 2, delta = eps ; end
n = max ( size ( A ) ) ;
[L ,D , p] = ldl (A ,'vector') ;
DMC = eye (n ) ;
% ( More and Sorensen ) modified Cholesky perturbations .
k = 1;
while k <= n
    if k == n || D (k , k +1) == 0 % 1-by -1 block
        if abs ( D (k ,k ) ) <= delta
            DMC (k , k ) = delta ;
        else
            DMC (k , k ) = abs ( D (k , k )) ;
        end
        k = k +1;
    else % 2 -by -2 block
        E = D ( k :k +1 , k : k +1) ;
        [U , T ] = eig ( E ) ;
        T = abs (T ) ;
        for ii = 1:2
            if T ( ii , ii ) <= delta
                T ( ii , ii ) = delta ;
            end
        end
        temp = U *T *U';
        DMC ( k : k +1 , k : k +1) = ( temp + temp')/2; % Ensure symmetric .
        k = k + 2;
    end
end
if nargout >= 3 , P = eye ( n ) ; P = P (p ,:) ;end    
end