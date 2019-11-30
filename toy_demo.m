%   Note even the graclus I compiled worked on my machine, but it is not
%   guaranteed to work on all Windows users. I highly recommend users to run
%   the code on Linux. Windows users may have to compile graclus on their
%   own. The website of graclus is:
%   http://www.cs.utexas.edu/users/dml/Software/graclus.html

% For Linux users, use the Linux version of Graclus
addpath(genpath('graclus1.2(linux)'));
% For Windows users, use the Windows version of Graclus
% addpath(genpath('graclus1.2(windows)'));

%% CAGrQc
% run Moana
load('CAGrQc.mat');
[acc_moana, t_moana, Q1, active1, Q2, active2, S_core] = Moana(A1, B1, L, 0.5, 100, 1e-4, 5, 4, 0, 0, gnd);
% Here, acc_moana is the node-level alignment accuracy.
% For cluster-level alignment, one can use the following code
[clusters1, clusters2, cluster_corrs] = cluster_collection(Q1, active1, Q2, active2, gnd);
ratio = 0.15;
[acc_cluster, hits] = evaluate_cluster_align(S_core, Q1, Q2, active1, active2, ratio, cluster_corrs);

% run FINAL_P
[S, t_final] = FINAL_P(A1, B1, L, 0.5, 100, 1e-4);
M = greedy_match(S);
[row, col] = find(M);
acc_final = length(intersect([col row], gnd, 'rows'))/length(gnd);