function [S, t] = FINAL_P(A, B, H, alpha, max_iter, tol)

d1 = sum(A, 1); d2 = sum(B, 1);
d1 = d1.^(-0.5); d2 = d2.^(-0.5);
d1(d1 == Inf) = 0; d2(d2 == Inf) = 0;
W1 = bsxfun(@times, d1, bsxfun(@times, d1', A));
W2 = bsxfun(@times, d2, bsxfun(@times, d2', B));

S = H;
fprintf('iteration begins.\n');
t0 = tic;
for i=1:max_iter
    tic;
    prev = S;
    S = alpha*W2*S*W1 + (1-alpha)*H;
   
    delta = norm(S - prev, 'fro');
    fprintf('iteration %d, delta = %f, running time = %f\n', i, delta, toc);
    
    if delta < tol, break; end
end
t = toc(t0);

