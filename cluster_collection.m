function [clusters1, clusters2, cluster_corrs] = cluster_collection(Q1, active1, Q2, active2, gnd)
    %% Collect cluster assignment at different levels 
    % Description:
    %   This algorithm is to collect cluster assignments at different
    %   levels of two graphs and form the cluster correspondence described
    %   in the paper.
    %
    % Inputs: 
    %   - Q1: interpolation matrices at different levels for graph 1
    %   - active1: active nodes at different levels for graph 1
    %   - Q2: interpolation matrices at different levels for graph 2
    %   - active2: active nodes at different levels for graph 2
    %   - gnd: node-level alignment ground-truth for evaluations
    %
    % Outputs:
    %   - clusters1: the cluster assignment at different levels of graph 1
    %   - clusters2: the cluster assignment at different levels of graph 2
    %   - cluster_corrs: cluster correspondence at different levels
    %

levels = length(Q1);
clusters1 = cell(levels, 1); clusters2 = cell(levels, 1);
locate1 = cell(levels, 1); locate2 = cell(levels, 1);
%% collect nodes for cluster at each level
for l = 1: levels
    supernodes1 = cell(length(active1{l}), 1);
    locate1{l} = zeros(length(active1{l}), 1);
    fprintf('collect graph1 clusters at level %d\n', l);
    for i = 1: length(active1{l})
        start_node = active1{l}(i);
        locate1{l}(i) = start_node;
        candidate = [start_node];
        cluster = [start_node];
        while ~isempty(candidate)
            removed_cand = candidate(1);
            new_nodes = find(Q1{l}(removed_cand, :));
            new_nodes = new_nodes(~ismember(new_nodes, cluster));
            candidate = candidate(2: end);
            candidate = [candidate, new_nodes];
            cluster = [cluster, new_nodes];           
        end
        if l == 1
            cluster_combined = cluster;
        else
            cluster_combined = [];
            for j = 1: length(cluster)
                cluster_combined = [cluster_combined, clusters1{l-1}{locate1{l-1} == cluster(j)}];
            end
        end
        supernodes1{i} = cluster_combined;
    end
    clusters1{l} = supernodes1;
    fprintf('collect graph2 clusters at level %d\n', l);    
    supernodes2 = cell(length(active2{l}), 1);
    locate2{l} = zeros(length(active1{l}), 1);
    for i = 1: length(active2{l})
        start_node = active2{l}(i);
        locate2{l}(i) = start_node;
        candidate = [start_node];
        cluster = [start_node];
        while ~isempty(candidate)
            removed_cand = candidate(1);
            new_nodes = find(Q2{l}(removed_cand, :));
            new_nodes = new_nodes(~ismember(new_nodes, cluster));
            candidate = candidate(2: end);
            candidate = [candidate, new_nodes];
            cluster = [cluster, new_nodes];           
        end
        if l == 1
            cluster_combined = cluster;
        else
            cluster_combined = [];
            for j = 1: length(cluster)
                cluster_combined = [cluster_combined, clusters2{l-1}{locate2{l-1} == cluster(j)}];
            end
        end
        supernodes2{i} = cluster_combined;
    end
    clusters2{l} = supernodes2;
end


%% generate cluster correspondence by overlapping nodes
cluster_corrs = cell(levels, 1);
for l = 1: levels
    cluster1 = clusters1{l}; cluster2 = clusters2{l};
    c1 = length(cluster1); c2 = length(cluster2);
    overlap_count = zeros(c1, c2);
    fprintf('collect cluster correspondence ground-truth at level %d\n', l);
    cluster_id2 = zeros(1, size(Q2{1}, 1));
    for i = 1: c2
        cluster_id2(cluster2{i}) = i;
    end
    for i = 1: c1
        nodes1 = cluster1{i};
        nodes1_mapped = gnd(nodes1, 2);
        cluster_id_mapped = cluster_id2(nodes1_mapped);
        val = unique(cluster_id_mapped,'stable');
        overlap_node_count = sum(bsxfun(@eq,val,cluster_id_mapped(:)));
        overlap_node_count(val == 0) = []; val(val == 0) = []; 
        overlap_count(i, val) = overlap_node_count;
    end
    M = greedy_match(overlap_count);
    [row, col] = find(M);
    cluster_corr = zeros(length(row), 2);
    for p = 1: length(row)
        cluster_corr(p, 1) = cluster1{row(p)}(1);
        cluster_corr(p, 2) = cluster2{col(p)}(1);
    end
    cluster_corrs{l} = cluster_corr;
end




