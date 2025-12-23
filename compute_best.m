
function [w1, w2, m1, m2, v1, v2] = compute_best( tab, pts )

loglik = [];
[Nw, Nm, Nv] = size(tab.v1);
for i = 1:Nw 
    for j = 1:Nm
        sqr1 = (pts-tab.m1(i,j)).^2;
        sqr2 = (pts-tab.m1(i,j)).^2;
        for k = 1:Nv
            pdf_tmp = tab.w1(i)/sqrt(2*pi*tab.v1(i,j,k))* ...
                  exp(-0.5*sqr1/tab.v1(i,j,k)) + ...
                  tab.w2(i)/sqrt(2*pi*tab.v2(i,j,k))* ...
                  exp(-0.5*sqr2/tab.v2(i,j,k));
                  loglik(i,j,k) = sum(log(pdf_tmp+realmin)); %add a small number to avoid inf ?
        end
    end
end
   
% choose hypothesis with max log likelihood
[maxVal, linearIdx] = max(loglik(:));
[i, j, k] = ind2sub(size(loglik), linearIdx);
w1 = tab.w1(i);
w2 = tab.w2(i);
m1 = tab.m1(i,j);
m2 = tab.m2(i,j);
v1 = tab.v1(i,j,k);
v2 = tab.v2(i,j,k);

return;

Mu = [m1, m2]';
Sigma(:,:,1)=v1;
Sigma(:,:,2)=v2;
P = [w1, w2];
S = struct('mu',Mu,'Sigma',Sigma,'ComponentProportion',P);
options = statset('MaxIter',50, 'TolFun',1e-3 );
fine_tune = fitgmdist(pts, 2, 'RegularizationValue', 0.0001, 'ProbabilityTolerance', 1e-6, 'Options',options, 'Start', S);
%fine_tune = fitgmdist(pts, 2, 'RegularizationValue', 0.0001, 'ProbabilityTolerance', 1e-6, 'Options',options);
w1 = fine_tune.ComponentProportion(1);
w2 = fine_tune.ComponentProportion(2);
m1 = fine_tune.mu(1);
m2 = fine_tune.mu(2);
v1 = fine_tune.Sigma(1);
v2 = fine_tune.Sigma(2);

