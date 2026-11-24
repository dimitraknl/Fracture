function [x,Dx,wx,domainMap] = multidomainChebyshevLobatto(breaks,nodesPerDomain)
% Build a fixed multidomain Chebyshev-Lobatto grid with shared interface
% nodes. Differentiation is local to each subdomain. At a shared interface,
% the left and right derivative rows are averaged for the smooth healthy
% fields. Clenshaw-Curtis weights from neighboring domains are summed.

    nDom = numel(nodesPerDomain);
    if numel(breaks) ~= nDom+1
        error('breaks must contain one more entry than nodesPerDomain.');
    end

    x = [];
    domainMap = cell(nDom,1);
    localData = cell(nDom,1);

    for e = 1:nDom
        Ne = nodesPerDomain(e);
        [xi,Di,wi] = chebyshevLobattoClenshawCurtis(Ne);
        xa = breaks(e); xb = breaks(e+1);
        xe = 0.5*(xb-xa)*(xi+1)+xa;
        De = (2/(xb-xa))*Di;
        we = 0.5*(xb-xa)*wi;

        if e == 1
            idx = (1:Ne).';
            x = xe;
        else
            idx = [(numel(x)); (numel(x)+(1:Ne-1)).'];
            x = [x; xe(2:end)]; %#ok<AGROW>
        end
        domainMap{e} = idx;
        localData{e} = struct('D',De,'w',we);
    end

    Nx = numel(x);
    Dx = sparse(Nx,Nx);
    wx = zeros(Nx,1);
    rowCount = zeros(Nx,1);

    for e = 1:nDom
        idx = domainMap{e};
        De = localData{e}.D;
        we = localData{e}.w;
        Dx(idx,idx) = Dx(idx,idx) + De;
        wx(idx) = wx(idx) + we;
        rowCount(idx) = rowCount(idx) + 1;
    end

    shared = find(rowCount > 1);
    for r = shared.'
        Dx(r,:) = Dx(r,:)/rowCount(r);
    end

    if any(diff(x) <= 0)
        error('The assembled multidomain grid is not strictly increasing.');
    end
end
