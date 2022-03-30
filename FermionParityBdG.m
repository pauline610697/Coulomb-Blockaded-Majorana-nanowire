function [pf] = FermionParityBdG(A)
%modified from Wimmer to give only sign
%pfaffian_hessenberg: Compute the Pfaffian of a real skew-symmetric matrix
%
% pf = pfaffian_hessenberg(A) computes the Pfaffian of the skew-symmetric matrix A
% using the Hessenberg decomposition

    assert(ndims(A)==2, 'argument must be a matrix')
    assert(size(A,1)==size(A,2), 'argument is not skew-symmetric')
    assert(isreal(A), 'argument must be a real matrix')
    
    N=size(A,1);
    
    if( mod(N,2) == 1 )
        pf = 0.0;
        return;
    end
    
    % calculate number of spatial dimensions
    n=N/4;
    
    sy=[0 -1i;1i 0];% Paulimatrix sigmay
    syty=kron(sy,sy); % particle-hole matrix
    A=A*kron(eye(n),syty); % depends on how the PHM is implemented
    %make sure input is skew-symmetric
    assert(norm(A+A.')<1e-14*size(A,1), 'argument does not seem skew-symmetric')
    
    % antisymmetrize A
    A=0.5*(A-transpose(A));
    
    [Q,H] = hess(A);
    
    pf = 1.0;
    for i=1:2:N-1
        pf = pf * sign(H(i,i+1));
    end
    pf = pf * sign(det(Q));
    
end
