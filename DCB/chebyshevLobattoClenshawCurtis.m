function [x,D,w] = chebyshevLobattoClenshawCurtis(N)

    if N < 2
        error('N must be at least 2.');
    end

    n = N-1;
    theta = pi*(0:n)'/n;
    xDesc = cos(theta);

    c = [2;ones(n-1,1);2].*(-1).^(0:n)';
    X = repmat(xDesc,1,N);
    dX = X-X.';

    Ddesc = (c*(1./c).')./(dX+eye(N));
    Ddesc = Ddesc-diag(sum(Ddesc,2));

    wDesc = zeros(N,1);
    ii = 2:n;
    v = ones(n-1,1);

    if mod(n,2)==0
        wDesc(1) = 1/(n^2-1);
        wDesc(N) = wDesc(1);
        for k = 1:n/2-1
            v = v-2*cos(2*k*theta(ii))/(4*k^2-1);
        end
        v = v-cos(n*theta(ii))/(n^2-1);
    else
        wDesc(1) = 1/n^2;
        wDesc(N) = wDesc(1);
        for k = 1:(n-1)/2
            v = v-2*cos(2*k*theta(ii))/(4*k^2-1);
        end
    end

    wDesc(ii) = 2*v/n;

    perm = N:-1:1;
    x = xDesc(perm);
    D = Ddesc(perm,perm);
    w = wDesc(perm);
end
