function exec_exps(numRepetitions, numPoints, Nw, Nm, Nv, saveBoxPlot2File, outputFileName, outputFolder)

tab = create_table(Nw, Nm, Nv, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%% Simulations

%allow user interaction
ui = 0;

% Testing arrays
ours_timer_arr = zeros(numRepetitions, 1, 'double');
em_mtc_arr = zeros(numRepetitions, 1, 'double');
gmmDist_ours_arr = zeros(numRepetitions, 1, 'double');
gmmDist_em_arr = zeros(numRepetitions, 1, 'double');
gmmLoglik_ours_arr = zeros(numRepetitions, 1, 'double');
gmmLoglik_em_arr = zeros(numRepetitions, 1, 'double');


% quantile values for statistics - do not change this (hardcoded)
q = [0.025 0.25 0.5 0.75 0.975]; %lowtail, lowquart, median, highquart, hightail
rng(1,'twister'); %reset seed and random number generator for repeatable results

for i = 1:numRepetitions    
    disp('Iteration: ');
    disp( i );
    % create a 1D gmm with 2 components
    mu = [20*rand-10,20*rand-10]'; % mean uniform between -10 and 10
    std = [rand*10+0.0001, rand*10+0.0001]; % std uniform between 0.0001 and 10
    sigma = reshape(std.^2,1,1,2);
    p = rand; %mixing proportions
    gmm = gmdistribution(mu,sigma,[p,1-p]);
    
    % generate numPoints sample points from the distribution
    pts = random(gmm, numPoints);
    
    % compute theoretical mean and covariance
    theo_mean = p*mu(1)+(1-p)*mu(2);
    theo_var = p*(std(1)^2+mu(1)^2)+(1-p)*(std(2)^2+mu(2)^2)-theo_mean^2;

    % compute sample mean and covariance
    samp_mean = mean(pts);
    samp_var = var(pts);
    samp_std = sqrt(samp_var);
    % normalize to zero mean unit norm
    pts_norm = (pts - samp_mean)/samp_std;
    mu_norm = (mu - samp_mean)/sqrt(samp_var);
    std_norm = std/samp_std;
    sigma_norm = reshape(std_norm.^2, 1, 1, 2);
    gmm_norm = gmdistribution(mu_norm,sigma_norm,[p,1-p]);
 
    % Standard EM fit to compare
    options = statset('MaxIter',1500, 'TolFun',1e-3 );
    tic;
    %gmm_fit = fitgmdist(pts2, 2);
    %gmm_fit = fitgmdist(pts_norm, 2, 'RegularizationValue', 0.001);
    gmm_fit = fitgmdist(pts_norm, 2, 'RegularizationValue', 0.0001, 'ProbabilityTolerance', 1e-6, 'Options',options);
    em_timer = toc;

    % Our algorithm
    tic;
    % cycle over hypotheses and evaluate
    loglik = [];
    [w1, w2, m1, m2, v1, v2] = compute_best(tab,pts_norm);

    % denormalize
    %std1_est = sqrt(v1_best)*samp_std;
    %std2_est = sqrt(v2_best)*samp_std;
    %mean1_est = mean1_best*samp_std+samp_mean;
    %mean2_est = mean2_best*samp_std+samp_mean;
    %w1_est = w1_best;
    %w2_est = w2_best;
    ours_timer = toc;       
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test time performance
    ours_timer_arr(i) = ours_timer;
    em_mtc_arr(i) = em_timer;  

    % Test precision performance
    mu_best = [m1 m2]';
    v_best = [v1 v2];
    sigma_best = reshape(v_best, 1, 1, 2);
    p_best = [w1 w2];
    gmm_best = gmdistribution(mu_best, sigma_best, p_best);
    
    gmmLoglik_em_arr(i) = sum(log(pdf(gmm_fit, pts_norm) + 0.000001));
    gmmLoglik_ours_arr(i) = sum(log(pdf(gmm_best, pts_norm) + 0.000001 ));

    gmmDist_ours_arr(i) = dist_kld(gmm_norm, gmm_best);
    gmmDist_em_arr(i) = dist_kld(gmm_norm, gmm_fit);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    if ui,
        ours_mtc = ours_timer/i;
        em_mtc = em_timer/i;
        X = ['Ours MTC: ',num2str(ours_mtc) , '    EM MTC: ', num2str(em_mtc)];
        disp(X);  
        
        if gmmDist_em_arr(i) > 1,
            figure(3);
            mu_best = [m1,m2]';
            std_best = [v1, v2];
            sigma_best = reshape(std_best,1,1,2);
            p_best = w1_best;
            gmm_best = gmdistribution(mu_best,sigma_best,[p_best,1-p_best]);
            maximo = max(pts);
            minimo = min(pts);
            range = minimo:(maximo-minimo)/1000:maximo;
            plot(range',pdf(gmm_norm,range'),'r-',range',pdf(gmm_best,range'),'g',range',pdf(gmm_fit,range'),'b');
            legend('True','Ours','EM');
            keyboard;
        end
    end
end


% Testing results

% Means
ours_timer_mean = sum(ours_timer_arr)/numRepetitions;
em_mtc_mean = sum(em_mtc_arr)/numRepetitions;
gmmDist_ours_mean = sum(gmmDist_ours_arr)/numRepetitions;
gmmDist_em_mean = sum(gmmDist_em_arr)/numRepetitions;
gmmLoglik_ours_mean = sum(gmmLoglik_ours_arr)/numRepetitions;
gmmLoglik_em_mean = sum(gmmLoglik_em_arr)/numRepetitions;

% Stdev
ours_timer_std_arr = zeros(numRepetitions, 1, 'double');
em_mtc_std_arr = zeros(numRepetitions, 1, 'double');
gmmDist_ours_std_arr = zeros(numRepetitions, 1, 'double');
gmmDist_em_std_arr = zeros(numRepetitions, 1, 'double');
gmmLoglik_ours_std_arr = zeros(numRepetitions, 1, 'double');
gmmLoglik_em_std_arr = zeros(numRepetitions, 1, 'double');

for i = 1:numRepetitions
    ours_timer_std_arr(i) = ours_timer_mean - ours_timer_arr(i);
    em_mtc_std_arr(i) = em_mtc_mean - em_mtc_arr(i);
    gmmDist_ours_std_arr(i) = gmmDist_ours_mean - gmmDist_ours_arr(i);
    gmmDist_em_std_arr(i) = gmmDist_em_mean - gmmDist_em_arr(i);
    gmmLoglik_ours_std_arr(i) = gmmLoglik_ours_mean - gmmLoglik_ours_arr(i);
    gmmLoglik_em_std_arr(i) = gmmLoglik_em_mean - gmmLoglik_em_arr(i);
end
ours_timer_std = sqrt(1/(numRepetitions-1)*sum(ours_timer_std_arr.^2));
em_mtc_std = sqrt(1/(numRepetitions-1)*sum(em_mtc_std_arr.^2));
gmmDist_ours_std = sqrt(1/(numRepetitions-1)*sum(gmmDist_ours_std_arr.^2));
gmmDist_em_std = sqrt(1/(numRepetitions-1)*sum(gmmDist_em_std_arr.^2));
gmmLoglik_ours_std = sqrt(1/(numRepetitions-1)*sum(gmmLoglik_ours_std_arr.^2));
gmmLoglik_em_std = sqrt(1/(numRepetitions-1)*sum(gmmLoglik_em_std_arr.^2));

% Max time
ours_timer_max = max(ours_timer_arr);
em_mtc_max = max(em_mtc_arr);

% Max dist
gmmDist_ours_max_dist = max(gmmDist_ours_arr);
gmmDist_em_max_dist = max(gmmDist_em_arr);

% Min loglik
gmmLoglik_ours_min = min(gmmLoglik_ours_arr);
gmmLoglik_em_min = min(gmmLoglik_em_arr);

%quartiles
qt_dist_ours = quantile(gmmDist_ours_arr,q);
qt_dist_em = quantile(gmmDist_em_arr,q);
qt_time_ours = quantile(ours_timer_arr,q);
qt_time_em = quantile(em_mtc_arr, q);
qt_llik_ours = quantile(gmmLoglik_ours_arr,q);
qt_llik_em = quantile(gmmLoglik_em_arr,q);

% Save output to screen
dump_results(1, numRepetitions, numPoints, Nw, Nm, Nv, ...
              ours_timer_mean, ours_timer_std, ours_timer_max, ...
              em_mtc_mean, em_mtc_std, em_mtc_max, ...
              gmmDist_ours_mean, gmmDist_ours_std, gmmDist_ours_max_dist, ...
              gmmDist_em_mean, gmmDist_em_std, gmmDist_em_max_dist, ...
              gmmLoglik_ours_mean, gmmLoglik_ours_std, gmmLoglik_ours_min, ...
              gmmLoglik_em_mean, gmmLoglik_em_std, gmmLoglik_em_min, ...
              qt_dist_ours, qt_dist_em, qt_time_ours, qt_time_em, ...
              qt_llik_ours, qt_llik_em);

% Save ouput to file
output_file_name = [outputFolder '/' outputFileName '.txt'];
fid = fopen(output_file_name, 'w');
dump_results(fid, numRepetitions, numPoints, Nw, Nm, Nv, ...
                ours_timer_mean, ours_timer_std, ours_timer_max, ...
                em_mtc_mean, em_mtc_std, em_mtc_max, ...
                gmmDist_ours_mean, gmmDist_ours_std, gmmDist_ours_max_dist, ...
                gmmDist_em_mean, gmmDist_em_std, gmmDist_em_max_dist, ...
                gmmLoglik_ours_mean, gmmLoglik_ours_std, gmmLoglik_ours_min, ...
                gmmLoglik_em_mean, gmmLoglik_em_std, gmmLoglik_em_min, ...
                qt_dist_ours, qt_dist_em, qt_time_ours, qt_time_em, ...
                qt_llik_ours, qt_llik_em);
fclose(fid);

% Draw boxplot kl divergence
fig = figure(1);
boxplot([gmmDist_ours_arr, gmmDist_em_arr]);
xlabel('Ours (1)                                                          EM (2)');
ylabel('KL-Divergence');
titleParameters = sprintf('Nw = %d; Nm = %d; Nv = %d', Nw, Nm, Nv);
title(['KL-Divergence D_{KL}(true || estim)   ' titleParameters])

% Save boxplot to file
if (saveBoxPlot2File == 1)
    resolution = '-r300';
    format = '-depsc'; %EPS color; or format = '-depsc2' for EPS Level 2 color
    print(fig, resolution, format, '-painters', [outputFolder '/' outputFileName '_kl.eps']);
    format = '-djpeg'; %jpeg-24 bit
    print(fig, resolution, format, '-painters', [outputFolder '/' outputFileName '_kl.jpg']);
end

% Draw boxplot loglikelihood
fig = figure(2);
boxplot([gmmLoglik_ours_arr, gmmLoglik_em_arr]);
xlabel('Ours (1)                                                          EM (2)');
ylabel('Log-Likelihood');
titleParameters = sprintf('Nw = %d; Nm = %d; Nv = %d', Nw, Nm, Nv);
title(['Data Loglikelihood   ' titleParameters])

% Save boxplot to file
if (saveBoxPlot2File == 1)
    resolution = '-r300';
    format = '-depsc'; %EPS color; or format = '-depsc2' for EPS Level 2 color
    print(fig, resolution, format, '-painters', [outputFolder '/' outputFileName '_ll.eps']);
    format = '-djpeg'; %jpeg-24 bit
    print(fig, resolution, format, '-painters', [outputFolder '/' outputFileName '_ll.jpg']);
end




