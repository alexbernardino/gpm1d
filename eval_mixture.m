function loglik = eval_mixture(pts, w1, w2, m1, m2, v1, v2)

pdf_tmp = w1/sqrt(2*pi*v1)*exp(-0.5*(pts-m1).^2/v1) + ...
                  w2(i)/sqrt(2*pi*v2)*exp(-0.5*(pts-m2).^2/v2);
loglik = sum(log(pdf_tmp))
%alternative evaluation using matlab toolbox
mu = [m1, m2]';
var = reshape([v1 v2], 1, 1, 2); 
p = [w1 w2]; 
gm = gmdistribution(mu,var,p);
pdf_vals = pdf(gm, pts); 
loglik = sum(log(pdf_tmp))

