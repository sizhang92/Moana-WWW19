function [Q, H, active, Aout, error, L] = pMMF(A, L, par)
    %% Parallel multiresolution matrix factorization
    % Description:
    %   This algorithm is to factorize an input matrix such that the
    %   representations at different resolution of the matrix can be
    %   obtained.
    %
    % Inputs: 
    %   - A: the input matrix to be factorized
    %   - L: levels to be factorized (i.e., # of resolutions)
    %   - par: number of threads for parallelization
    %
    % Outputs:
    %   - Q: orthogonal rotation matrices at different levels 
    %   - H: core-diagonal matrix at the coarsest level
    %   - active: active nodes at different levels 
    %   - Aout: rotated matrices at different levels
    %   - error: difference at the coarsest level
    %   - L: number of levels
    %
    % Reference:
    %   Teneva, Nedelina, Pramod Kaushik Mudrakarta, and Risi Kondor. 
    %   "Multiresolution matrix compression." AISTATS, pp. 1441-1449. 2016.
    %
    %   Kondor, Risi, Nedelina Teneva, and Vikas Garg. 
    %   "Multiresolution matrix factorization." ICML, pp. 1620-1628. 2014.



if isempty(gcp('nocreate'))
    parpool(par);
end

rng(123);
eta = 0.5;
n = size(A, 1);
act = [1: n]';
active = cell(L, 1); Q = cell(L, 1); Aout = cell(L, 1);

for l = 1: L
    fprintf('level %d.\n', l);
    if l == 1
        m = max(2, round(sqrt(length(act))));
        B = A; B(B~=0) = 1; 
        G = graph(B);
        partition = graclus(B, m);
        deg = centrality(G, 'pagerank');
    else
        act = active{l-1};
        B = Aout{l-1}(act, act);
        B(B > 0) = 1;
        B(B < 0) = 0.001;       
        m = max(2, round(sqrt(length(act))));
        partition = graclus(B, m);
    end
    fprintf('finished clustering.\n');
    blocks = cell(m, 1); index = cell(m, 1);
    active1 = cell(m, 1); degrees = cell(m, 1);
    row_Qs = cell(m, 1); col_Qs = cell(m, 1); val_Qs = cell(m, 1);
    for i = 1: m
        index{i} = act(partition == i);
        degrees{i} = deg(partition == i);
        if l == 1
            blocks{i} = A(:, index{i});
        else
            blocks{i} = Aout{l-1}(:, index{i});
        end
    end
    
%     for i = 1: m
    parfor (i = 1: m, par)
        Qs = eye(length(index{i})); 
        if length(index{i}) <= 1, continue; end
        activeset = 1: length(index{i}); % index based on blocks
        w = degrees{i};
        
        gram = full(blocks{i}' * blocks{i});
        s = 1; pos = 1; s_total = floor(eta*length(activeset));
        while s <= s_total 
            if length(activeset) <= 1, break; end
            blks = blocks{i};
            % using approximate algorithm
            [~, w_idx] = sort(w, 'descend');
            if pos > length(w_idx), break; end

            idx = activeset(w_idx(pos));                           
            temp = full(sum(blocks{i}.^2, 1));
            score = abs(gram(idx, activeset))./sqrt(temp(activeset));
            score(activeset == idx) = -Inf;
            [val, id] = max(score);
            if val > 0
                id(activeset(id) == idx) = [];
                if isempty(id)
                    pos = pos + 1;
                    continue; 
                end                
                if length(id) > 1
                    [~, min_id] = min(w(id));
                    id = id(min_id(1));
                end
                id = activeset(id);
                theta = 0: 0.01*pi: 2*pi;
                i1 = index{i}(idx); i2 = index{i}(id);
                a1 = blks(i1, idx); a2 = blks(i1, id);
                a3 = blks(i2, idx); a4 = blks(i2, id);
                cost1 = ((a1-a4)*cos(theta).*sin(theta)+a3*cos(theta).^2-a2*sin(theta).^2).^2;
                cost2 = sin(theta')*blocks{i}(i1, activeset)+cos(theta')*blocks{i}(i2, activeset);
                cost2 = sum(cost2.^2, 2)';
                cost = cost1 + cost2;
                [~, min_theta] = min(cost); min_theta = theta(min_theta(1));
            else
                if l >= ceil(L/2)
                    theta = 0: 0.01*pi: 2*pi;
                    a1 = blks(index{i}(idx), idx) * ones(length(activeset), 1); 
                    a2 = blks(index{i}(idx), activeset)';
                    i1 = index{i}(activeset); i2 = activeset';
                    a4 = blks(sub2ind(size(blks), i1, i2));
                    cost1 = (0.5*(a1-a4)*sin(2*theta) + a2*cos(2*theta)).^2;
                    b1 = gram(idx, idx) * ones(length(activeset), 1); b2 = gram(idx, activeset)';
                    b4 = diag(gram); b4 = b4(activeset);
                    cost2 = b1*sin(theta).^2 + b4*cos(theta).^2 + 2*b2*sin(theta).*cos(theta);
                    cost3 = cost1 + cost2; 
                    cost3(a2 == 0, :) = Inf;
                    cost3(w_idx(pos), :) = Inf;
                    [~, min_idx] = min(cost3, [], 'all', 'linear');

                    [r, c] = ind2sub(size(cost3), min_idx);
                    if length(c) > 1
                        [~, min_r] = min(w(r));
                        r = r(min_r(1)); c = c(min_r(1));
                    end
                
                    min_theta = theta(c);
                    id = activeset(r);
                    i1 = index{i}(idx); i2 = index{i}(id);
                else
                    pos = pos + 1;
                    s = s + 1;
                    continue; 
                end
            end               
                  
            R = [cos(min_theta), -sin(min_theta); sin(min_theta), cos(min_theta)];
            
            blocks{i}(:, [idx, id]) = blocks{i}(:, [idx, id]) * R';
            temp = speye(n); temp([i1, i2], [i1, i2]) = R;
            blocks{i} = temp * blocks{i};
            Qs([idx, id], :) = R * Qs([idx, id], :);

            gram(:, [idx, id]) = gram(:, [idx, id]) * R';
            gram([idx, id], :) = R * gram([idx, id], :);
            
            a = norm(blocks{i}(:, idx))^2 - blocks{i}(i1, idx)^2;
            b = norm(blocks{i}(:, id))^2 - blocks{i}(i2, id)^2;
            if a < b
                w(activeset == idx) = [];
                activeset(activeset == idx) = [];
            else
                w(activeset == id) = [];
                activeset(activeset == id) = []; 
                pos = pos + 1;
            end
            
            s = s + 1;
        end
        [row1, col1, val1] = find(Qs);
        row_Qs{i} = index{i}(row1); col_Qs{i} = index{i}(col1);
        val_Qs{i} = val1; 

        active1{i} = index{i}(activeset);
    end
        
    % combine active index
    active{l} = unique(cell2mat(active1));
    row_Q = cell2mat(row_Qs); col_Q = cell2mat(col_Qs); val_Q = cell2mat(val_Qs);
    Q{l} = sparse(row_Q, col_Q, val_Q, n, n);
    dQ = diag(Q{l}); new_diag = ones(n, 1); new_diag(dQ ~= 0) = 0;
    Q{l} = Q{l} + diag(sparse(new_diag));
    if l == 1
        Aout{l} = Q{l}*A*Q{l}';
    else
        Aout{l} = Q{l}*Aout{l-1}*Q{l}';
    end
    % check if isolated nodes exist
    B = Aout{l}(active{l}, active{l}); B(B>0) = 1; B(B<0) = 0.001;
    d1 = sum(B, 1);
    iso_idx = find(d1 == 0);
    B(iso_idx, :) = []; B(:, iso_idx) = [];
    active{l}(iso_idx) = [];
    
    Aout{l} = max(Aout{l}, Aout{l}');
    
    if isempty(active{l})
        L = l - 1;
        break;
    end
end
n1 = length(active{L}); 
edgelist = zeros(n1*n1+n-n1, 3);

edgelist(1: n1^2, 1) = repmat(active{L}(:), [n1, 1]);
temp = repmat(1: n1, [n1, 1]);
edgelist(1: n1^2, 2) = active{L}(temp(:));
temp = Aout{L}(active{L}, active{L});
edgelist(1: n1^2, 3) = temp(:);

rest = setdiff(1: n, active{L}, 'stable');
edgelist(n1^2+1 : end, 1) = rest;
edgelist(n1^2+1 : end, 2) = rest;
edgelist(n1^2+1 : end, 3) = Aout{L}(sub2ind([n, n], rest, rest));

H = sparse(edgelist(:, 1), edgelist(:, 2), edgelist(:, 3), n, n);
error = norm(Aout{L} - H, 'fro');


end


            
            
            