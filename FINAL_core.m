function [S, tf] = FINAL_core(A, B, active1, active2, H, alpha, maxiter, tol, par)
    %% Coarsest-level alignment
    % Description:
    %   This algorithm is to find the alignment at the coarsest level based
    %   on the existing algorithm named FINAL.
    %
    % Inputs: 
    %   - A: the core-diagonal matrix of adjacency matrix of graph 1
    %   - B: the core-diagonal matrix of adjacency matrix of graph 2
    %   - active1: active nodes at the coarsest levels for graph 1
    %   - active2: active nodes at the coarsest levels for graph 2
    %   - H: rotated matrix of the prior alignment matrix at node level
    %   - alpha: regularization parameter
    %   - maxiter: maximum number of iterations
    %   - tol: threshold to terminate the algorithm
    %   - par: number of threads for parallelization
    %
    % Outputs:
    %   - S: alignment matrix at the coarsest level
    %   - tf: running time for alignment
    %
    % Reference:
    %   Zhang, Si, and Hanghang Tong. "Final: Fast attributed network alignment."
    %   KDD, pp. 1345-1354. ACM, 2016.


if isempty(gcp('nocreate'))
    parpool(par);
end


n1 = size(A, 1); n2 = size(B, 1);
id1 = setdiff(1: n1, active1, 'stable')';
A_core = A(active1, active1); A_off = A(id1, id1);
id2 = setdiff(1: n2, active2, 'stable')';
B_core = B(active2, active2); B_off = B(id2, id2);

H_all = cell(4, 1); input1 = cell(4, 1); input2 = cell(4, 1);
H_all{1} = H(active2, active1); H_all{2} = H(active2, id1);
H_all{3} = H(id2, active1); H_all{4} = sparse(H(id2, id1));
input1{1} = A_core; input1{2} = A_off; input1{3} = A_core; input1{4} = A_off;
input2{1} = B_core; input2{2} = B_core; input2{3} = B_off; input2{4} = B_off;
rows = cell(4, 1); cols = cell(4, 1); vals = cell(4, 1);
row_idx = cell(4, 1); col_idx = cell(4, 1);
row_idx{1} = active2; col_idx{1} = active1; row_idx{2} = active2; 
col_idx{2} = id1; row_idx{3} = id2; col_idx{3} = active1;
row_idx{4} = id2; col_idx{4} = id1;

t0 = tic;
fprintf('started distributed final_core.\n');

parfor (i = 1: 4, par)
    if i == 4
        dA = diag(input1{i});
        dB = diag(input2{i});
        % find nonzero indices from H4
        h = H_all{i}(:); [m, ~] = size(H_all{i});
        [id_h4, ~, v] = find(h);
        idx_A = floor((id_h4-1) / m) + 1;
        idx_B = rem(id_h4 - 1, m) + 1;
        val = dA(idx_A) .* dB(idx_B);
        val = (1-alpha*val).^(-1); val(val == Inf) = 0;
        val = (1 - alpha) * val .* v;
        s = sparse(id_h4, 1, val, length(dB)*length(dA), 1);
        [r, c, v] = find(reshape(s, [length(dB), length(dA)]));
        rows{i} = row_idx{i}(r); cols{i} = col_idx{i}(c); vals{i} = v;
    else
        fprintf('iteration begins for core %d.\n', i);
        S1 = rand(length(input2{i}), length(input1{i}));
        for j = 1: maxiter
            prev = S1;
            S1 = alpha*input2{i}*S1*input1{i} + (1-alpha)*H_all{i};
            delta = norm(S1 - prev, 'fro');
%             fprintf('iteration %d, delta = %f for core %d.\n', j, delta, i);

            if delta < tol, break; end
        end

        [r, c, v] = find(S1);
        rows{i} = row_idx{i}(r); cols{i} = col_idx{i}(c); vals{i} = v;
    end

end
fprintf('finished distributed final_core.\n');
tf = toc(t0);
%% compose outputs

row = [rows{1}(:); rows{2}(:); rows{3}(:); rows{4}(:)];
col = [cols{1}(:); cols{2}(:); cols{3}(:); cols{4}(:)];
val = [vals{1}(:); vals{2}(:); vals{3}(:); vals{4}(:)];

S = sparse(row, col, val, n2, n1);


