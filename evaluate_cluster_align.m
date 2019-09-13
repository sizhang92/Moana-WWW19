function [acc_cluster, hits] = evaluate_cluster_align(S_core, Q1, Q2, active1, active2, ratio, cluster_corrs)
    %% Evaluate cluster alignment at different levels
    % Description:
    %   This algorithm is to evaluate cluster levels alignment.
    %
    % Inputs: 
    %   - S_core: alignment matrix obtained at the coarsest level
    %   - Q1: interpolation matrices at different levels for graph 1
    %   - Q2: interpolation matrices at different levels for graph 2
    %   - active1: active nodes at different levels for graph 1
    %   - active2: active nodes at different levels for graph 2
    %   - ratio: percentage to be considered as a hit (e.g., 0.15 in paper)
    %   - cluster_corrs: cluster correspondence at different levels
    %
    % Outputs:
    %   - acc_cluster: cluster alignment accuracy at different levels
    %   - hits: number of hits at different levels
    %
    
    
levels = length(Q1);
Sdash = S_core; hits = zeros(levels, 1);

for l = levels: -1: 1
    act1 = active1{l}; act2 = active2{l};
    if l < levels
        Sdash = Q2{l+1}'*Sdash*Q1{l+1};
        S_cluster = Sdash(act2, act1);
    else
        S_cluster = Sdash(act2, act1);
    end

    [~, I] = maxk(abs(S_cluster), round(ratio*length(act1)), 2);
    cluster_corr = cluster_corrs{l};
    for i = 1: length(act2)
        supernode2 = act2(i);
        if ismember(supernode2, cluster_corr(:,2))
            supernode1 = find(act1==cluster_corr(cluster_corr(:, 2)==supernode2, 1));
            if ismember(supernode1, I(i, :))
                hits(l) = hits(l) + 1;
            end
        end
    end
end

acc_cluster = zeros(levels, 1);
for i = 1: levels
    acc_cluster(i) = hits(i)/length(cluster_corrs{i});
end
