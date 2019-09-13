function [acc, ts, Q1, active1, Q2, active2, S_core] = Moana(A1, B1, H, alpha, maxiter, tol, num_levels, par, top_K, K, gnd)
    %% Multilevel network alignment algorithm
    % Description:
    %   This algorithm is to align networks not only at the node level, but
    %   also at other coarse levels to discover the correspondences between
    %   clusters at different resolutions.
    %
    % Inputs: 
    %   - A1, B1: the adjacency matrices of input networks with
    %   size n1-by-n1 and b2-by-n2 respectively
    %   - H: the n2-by-n1 prior similarity matrix
    %   - alpha: regularization parameter
    %   - maxiter: maximum number of iterations
    %   - tol: threshold to terminate the algorithm
    %   - num_levels: the number of levels to coarsen the networks
    %   - par: number of threads for parallelization
    %   - top_K: whether to use top_K approx as Line 13 in Algorithm 1
    %   - K: the K value for top-K preservation on S
    %   - gnd: node-level alignment ground-truth for evaluations
    %
    % Outputs:
    %   - acc: node-level alignment accuracy
    %   - ts: total time of the algorithm
    %   - Q1: interpolation matrices at different levels for graph 1
    %   - active1: active nodes at different levels for graph 1
    %   - Q2: interpolation matrices at different levels for graph 2
    %   - active2: active nodes at different levels for graph 2
    %   - S_core: alignment matrix obtained at the coarsest level
    %
    % Reference:
    %   Zhang, Si, Hanghang Tong, Ross Maciejewski, and Tina Eliassi-Rad.
    %   "Multilevel Network Alignment." In The World Wide Web Conference,
    %   pp. 2344-2354. ACM, 2019.
    
    
    if isempty(gcp('nocreate'))
        parpool(par);
    end
%     addpath(genpath('graclus1.2'));
    deg1 = sum(A1, 1); deg2 = sum(B1, 1);
    deg1 = deg1.^(-0.5); deg2 = deg2.^(-0.5);
    deg1(deg1 == Inf) = 0; deg2(deg2 == Inf) = 0;
    W1 = bsxfun(@times, deg1, bsxfun(@times, deg1', A1));
    W2 = bsxfun(@times, deg2, bsxfun(@times, deg2', B1));
    inputs = cell(2, 1); inputs{1} = W1; inputs{2} = W2; pmmf_error = zeros(2, 1);
    Qs = cell(2, 1); Hs = cell(2, 1); actives = cell(2, 1); true_levels = zeros(2, 1);
    MMFStart = tic;

    for i = 1: 2
    	[Qs{i}, Hs{i}, actives{i}, Aout{i}, pmmf_error(i), true_levels(i)] = pMMF(inputs{i}, num_levels, par);
    end
    
    tMMF = toc(MMFStart);
    Q1 = Qs{1}; H1 = Hs{1}; active1 = actives{1}; A1out = Aout{1};
    Q2 = Qs{2}; H2 = Hs{2}; active2 = actives{2}; B1out = Aout{2};
    levels = min(true_levels);
    Q1 = Q1(1: levels); Q2 = Q2(1: levels); 
    active1 = active1(1: levels); active2 = active2(1: levels);
    A1out = A1out(1: levels); B1out = B1out(1: levels);
    % in case that the off-core diagonal elements are beyond [-1,1] and make
    % the algorithm unable to converge.
    d1 = diag(H1); d2 = diag(H2);
    n1 = size(H1, 1); n2 = size(H2, 1);
    i1 = setdiff(find(d1>1), active1{levels}, 'stable'); i2 = setdiff(find(d2>1), active2{levels}, 'stable');
    j1 = setdiff(find(d1<-1), active1{levels}, 'stable'); j2 = setdiff(find(d2<-1), active2{levels}, 'stable');
    diff1 = [1-d1(i1); -1-d1(j1)]; diff2 = [1-d2(i2); -1-d2(j2)];
    D1 = sparse([i1;j1], [i1;j1], diff1, n1, n1); 
    D2 = sparse([i2;j2], [i2;j2], diff2, n2, n2);
    
    A = H1 + D1; B = H2 + D2;
    
    L = H; 
    innerStart = tic;
    for j = 1: levels  
        L = Q2{j}*L*Q1{j}';     
    end    
    tH = toc(innerStart);
    
    [S, tFINAL] = FINAL_core(A, B, active1{levels}, active2{levels}, L, alpha, maxiter, tol, par);
    if top_K
        [B1, I1] = maxk(S, K, 2, 'ComparisonMethod', 'abs'); I1_rows = repmat([1: n2]', [1, K]);
        [B2, I2] = maxk(S, K, 1, 'ComparisonMethod', 'abs'); I2_cols = repmat(1: n1, [K, 1]);
        idx1 = [I1_rows(:); I2(:)]; idx2 = [I1(:); I2_cols(:)]; vals = [B1(:); B2(:)];  
        Sdash = sparse(idx1, idx2, vals, n2, n1);
    else
        Sdash = S;
    end
    S_core = S;
    innerStart = tic;
    for j = levels: -1: 1
        Sdash = Q2{j}'*Sdash*Q1{j};
    end
    tS = toc(innerStart);
    % node-level alignment refinement
    postStart = tic;
    Sdash1 = alpha*W2*full(Sdash)*W1' + (1-alpha)*H;
    tPost = toc(postStart);
    ts = sum([tMMF, tH, tFINAL, tS, tPost]);
    
    %% evaluations

    M2 = greedy_match(Sdash1);
    [row, col] = find(M2);
    acc = length(intersect([col row], gnd, 'rows'))/length(gnd);
    
