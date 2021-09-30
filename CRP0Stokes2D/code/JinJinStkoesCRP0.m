function [soln,eqn, info] = JinJinStkoesCRP0(node,elem,bdFlag,pde)

% node 点
% elem 单元
% bdFlag 边界
% pde 偏微分方程

%% 数据结构
% elem2dof 在边edge上的点的全局坐标
[elem2dof, edge] = dofedge(elem);
NE = size(edge, 1); % 边的个数
NT = size(elem, 1); % 单元个数
Nu = NE; % 速度u的个数，自然等于NE，CRP0中u取的是每个边的中点
Np = NT; % 压力p的个数，自然等于单元的个数
N = size(node, 1); % 顶点个数

t = cputime; % 记录CPU时间

%% 计算局部基函数的几何量和梯度
% 通过重心坐标计算
[Dlambda, area] = gradbasis(node, elem);

%% 组装刚度矩阵，对于Laplace算子
A = sparse(Nu,Nu);
for i = 1:3
    for j = i:3
        % 局部节点到全局节点
        ii = double(elem2dof(:,i));
        jj = double(elem2dof(:,j));
        % 局部刚度矩阵
        Aij = 4*dot(Dlambda(:,:,i),Dlambda(:,:,j),2).*area;
        if (j==i)
            A = A + sparse(ii,jj,Aij,Nu,Nu);
        else
            A = A + sparse([ii,jj],[jj,ii],[Aij; Aij],Nu,Nu);
        end
    end
end
clear Aij
A = blkdiag(A,A);

%% 对于散度算子，组装刚度矩阵
d1 = -2.*Dlambda(:, :, 1) .* [area, area];
d2 = -2.*Dlambda(:, :, 2) .* [area, area];
d3 = -2.*Dlambda(:, :, 3) .* [area, area];
Dx = sparse(repmat((1:Np)', 3, 1), double(elem2dof(:)), ...
    [d1(:, 1); d2(:, 1); d3(:, 1)], Np, Nu);
Dy = sparse(repmat((1:Np)', 3, 1), double(elem2dof(:)), ...
    [d1(:, 2); d2(:, 2); d3(:, 2)], Np, Nu);
B = -[Dx Dy];
clear d1 d2 d3 B1 B2

%% 生成右端项
% 处理f
f1 = zeros(Nu, 1);
f2 = zeros(Nu, 1);
if ~isfield(pde, 'f') ||(isreal(pde.f) && (pde.f == 0)) % isfield函数是系统函数
    pde.f = [];
end
if ~isempty(pde.f) %如果右端项f非空
    mid1 = (node(elem(:,2),:) + node(elem(:,3),:))/2; %A2A3边的中点
    mid2 = (node(elem(:,3),:) + node(elem(:,1),:))/2;
    mid3 = (node(elem(:,1),:) + node(elem(:,2),:))/2;
    ft1 = repmat(area, 1, 2).* pde.f(mid1)/3; % 要除以3
    ft2 = repmat(area, 1, 2).* pde.f(mid2)/3;
    ft3 = repmat(area, 1, 2).* pde.f(mid3)/3;
    f1 = accumarray(elem2dof(:), [ft1(:,1); ft2(:,1); ft3(:,1)], [Nu 1]);
    f2 = accumarray(elem2dof(:), [ft1(:,2); ft2(:,2); ft3(:,2)], [Nu 1]);
end

% 处理h,预留随机接口
h = zeros(Np, 1);
if ~isfield(pde, 'h') ||(isreal(pde.h) && (pde.h == 0))
    pde.h = [];
end
if ~isempty(pde.h) %如果右端项h非空
    ht1 = repmat(area, 1, 2).* pde.h(mid1)/3; % 要除以3
    ht2 = repmat(area, 1, 2).* pde.h(mid1)/3;
    ht3 = repmat(area, 1, 2).* pde.h(mid1)/3;
    h = accumarray(repmat((1:Np)',3,1), [ht1(:,1); ht2(:,1); ht3(:,1)],...
        [Np 1]);
end

[AD, BD, f, g, u, p, ufreeDof, pDof] = getbdStokesCR;

%% Record assembeling time
assembleTime = cputime - t;

%% 求解线性方程组
if isempty(ufreeDof), return; end
t = cputime;
bigA = [AD, BD'; ...
        BD, sparse(Np,Np)];
bigF = [f; g];
bigu = [u; p];

% 对于只有Dirichlet边界，pDof = (1:Np-1)';???
bigFreeDof = [ufreeDof; 2*Nu+pDof]; 
bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
u = bigu(1:2*Nu);
p = bigu(2*Nu+1:end);
residual = norm(bigF - bigA*bigu);
info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);        

%% Post-process ???
if length(pDof) ~= Np % p is unique up to a constant
    % impose the condition int(p)=0
    c = sum(p.*area)/sum(area);
    p = p - c;
end

%% Output
soln = struct('u',u,'p',p);
eqn = struct('A',AD,'B',BD,'Lap',A,'f',f,'g',g,...
             'edge',edge,'ufreeDof',ufreeDof,'pDof',pDof);
info.assembleTime = assembleTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdStokesCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,BD,f,g,u,p,ufreeDof,pDof] = getbdStokesCR
        %% Boundary condition of Stokes equation: CR elements
        
        %% Initial set up
        f = [f1; f2];    % set in Neumann boundary condition
        g = zeros(Np,1); 
        u = zeros(2*Nu,1);
        p = zeros(Np,1);
        ufreeDof = (1:Nu)';
        pDof = (1:Np)';
        
        if ~exist('bdFlag','var'), bdFlag = []; end
        if ~isfield(pde,'g_D'), pde.g_D = []; end
        if ~isfield(pde,'g_N'), pde.g_N = []; end
        if ~isfield(pde,'g_R'), pde.g_R = []; end
        
        %% Part 1: Find Dirichlet dof and modify the matrix
        % Find Dirichlet boundary dof: fixedDof and pDof
        isFixedDof = false(Nu,1);
        if ~isempty(bdFlag)       % case: bdFlag is not empty
            isDirichlet(elem2dof(bdFlag(:)==1)) = true;
            isFixedDof(isDirichlet) = true;% dof on D-edges
            fixedDof = find(isFixedDof); % Dirichlet边界点
            ufreeDof = find(~isFixedDof); % 非Dirichlet边界点
        end
        % 如果bdFlag为空，只有Dirichlet边界点
        if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
            s = accumarray(elem2dof(:), 1, [Nu 1]);
            isFixedDof = (s==1);
            fixedDof = find(isFixedDof);
            ufreeDof = find(~isFixedDof);
        end
        if isempty(fixedDof) % 无Dirichlet边界点,即pure Neumann boundary condition
            % pde.g_N could be empty which is homogenous Neumann boundary condition
            fixedDof = 1;
            ufreeDof = 2:Nu;    % eliminate the kernel by enforcing u(1) = 0;
        end
        
        % Modify the matrix
        % Build Dirichlet boundary condition into the matrix AD by enforcing
        % AD(fixedDof,fixedDof)=I, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0.
        % BD(:,fixedDof) = 0 and thus BD'(fixedDof,:) = 0.
        bdidx = zeros(2*Nu,1);
        bdidx(fixedDof) = 1;
        bdidx(Nu+fixedDof) = 1;
        
        % 系统函数 spdiags
        %A = SPDIAGS(B,d,m,n) creates an m-by-n sparse matrix from the
        %    columns of B and places them along the diagonals specified by d.
       
        Tbd = spdiags(bdidx,0,2*Nu,2*Nu); % AD(fixedDof,fixedDof)=I
        % 1 - bdidx = 0, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0
        T = spdiags(1-bdidx,0,2*Nu,2*Nu); 
        AD = T*A*T + Tbd;
        BD = B*T; % BD(:,fixedDof) = 0
        
        %% Part 2: Find boundary edges and modify the right hand side f and g
        % Find boundary edges: Neumann and Robin
        Neumann = []; Robin = []; %#ok<*NASGU>
        if ~isempty(bdFlag)
            isNeumann(elem2dof((bdFlag(:)==2)|(bdFlag(:) == 3))) = true;
            isRobin(elem2dof(bdFlag(:)==3)) = true;
            Neumannidx = find(isNeumann);
            Neumann   = edge(isNeumann,:);
            Robin     = edge(isRobin,:);
        end
        if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
            % no bdFlag, only pde.g_N or pde.g_R is given in the input
            [tempvar,Neumann] = findboundary(elem);
            if ~isempty(pde.g_R)
                Robin = Neumann;
            end
        end
        
        % Neumann boundary condition
        if ~isempty(pde.g_N) && ~isempty(Neumann) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
            [lambda,w] = quadpts1(3);
            nQuad = size(lambda,1);
            % edge bases 
            bdphi = 2*(lambda(:,1)+lambda(:,2))-1;
            % length of edge
            ve = node(Neumann(:,1),:) - node(Neumann(:,2),:);
            edgeLength = sqrt(sum(ve.^2,2));
            % update RHS
            for pp = 1:nQuad
                pxy = lambda(pp,1)*node(Neumann(:,1),:)+lambda(pp,2)*node(Neumann(:,2),:);
                gp = pde.g_N(pxy);
                f1(Neumannidx) = f1(Neumannidx) + w(pp)*edgeLength.*gp(:,1).*bdphi(pp); % interior bubble
                f2(Neumannidx) = f2(Neumannidx) + w(pp)*edgeLength.*gp(:,2).*bdphi(pp); % interior bubble
            end
        end
        f = [f1; f2];
        % The case non-empty Neumann but g_N=[] corresponds to the zero flux
        % boundary condition on Neumann edges and no modification is needed.
        
        % Dirichlet boundary conditions
        % Dirichlet边界条件，处理右端项
        % u = g_D = 0 on (\Omega)
        if ~isempty(fixedDof) && ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && (pde.g_D == 0))
            u1 = zeros(Nu,1);
            u2 = zeros(Nu,1);
            bdEdgeMid = (node(edge(fixedDof,1),:)+node(edge(fixedDof,2),:))/2;
            uD = pde.g_D(bdEdgeMid);         % bd values at middle points of edges
            u1(fixedDof) = uD(:,1);
            u2(fixedDof) = uD(:,2);
            u = [u1; u2]; % Dirichlet bd condition is built into u
            
            % 处理边界条件g_D(u,t)中u涉及到非边界点的情况
            f = f - A*u;  % bring affect of nonhomgenous Dirichlet bd condition to
            
            % 处理右端项h
            g = h - B*u;
            % 无右端项h的情况
            % g = g - B*u;  % the right hand side
            
            % mean(g):所有元素的均值,系统函数
            g = g - mean(g); % impose the compatible condition
            
            f(fixedDof) = u1(fixedDof);
            f(fixedDof+Nu) = u2(fixedDof);
            u = [u1; u2]; % Dirichlet bd condition is built into u
        end
        % The case non-empty Dirichlet but g_D=[] corresponds to the zero Dirichlet
        % boundary condition and no modification is needed.
        
        % modfiy pressure dof for pure Dirichlet
        if isempty(Neumann)
            pDof = (1:Np-1)';
        end
        
        ufreeDof = [ufreeDof; Nu+ufreeDof];                
    end % end of function getbdStokesCR
end

