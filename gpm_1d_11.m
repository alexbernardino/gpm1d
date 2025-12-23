function gpm_1d_11(numRepetitions, numPoints, mean_bins, std_bins, saveBoxPlot2File, outputFileName, ouputFolder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%% Generate the hypotheses table
 %m1 hypotheses are positive. According to a numSigma rule, we define them
 %in the interval [-numSigma,0[
    numSigma = 3;
 
    m1_array = [];
    m2_array = [];
    w1_array = [];
    std1_array = [];
    std2_array = [];
    index_array = 0;
    
    m1 = -numSigma:numSigma/mean_bins:-numSigma/mean_bins;
    m2 = -m1;
    w = [];
    var2max = [];
    for j=1:length(m1);
        for k = 1:length(m2)
            if (m1(j)*m2(k) <= -1) % impossible - will make variance negative
                w(j,k) = 0;
                var2max(j,k) = 0;
            else
                w(j,k) = m2(k)/(m2(k)-m1(j));
                var2max(j,k)=(1-w(j,k)*m1(j)^2)/(1-w(j,k)) - m2(k)^2 ;
                std2max(j,k) = sqrt(var2max(j,k));
                %split the range of sigma2 into 10 bins
                range = std2max(j,k);
                std2 = [range/std_bins/2:range/std_bins:range-range/std_bins/2];
                var1 = ((m2(k)-m1(j))*(1+m1(j)*m2(k))+m1(j)*std2.^2)/m2(k);
                std1 = sqrt(var1);
                for l=1:length(std2),
                    index_array = index_array+1;
                    m1_array(index_array) = m1(j);
                    m2_array(index_array) = m2(k);
                    w1_array(index_array) = w(j,k);
                    std1_array(index_array) = std1(l);
                    std2_array(index_array) = std2(l);
                end;
                    
            end;
        end;
    end;
    
    num_hyp = length(m1_array)
    gmm_hyp = cell(1,num_hyp);
    for m = 1:num_hyp;
        mu_hyp = [m1_array(m), m2_array(m)]';
        var_hyp = reshape([std1_array(m)^2 std2_array(m)^2], 1, 1, 2); 
        p_hyp = w1_array(m); 
        gmm_hyp{m} = gmdistribution(mu_hyp,var_hyp,[p_hyp,1-p_hyp]);
    end;
    

    
    
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

for i = 1:numRepetitions,    
    % create a 1D gmm with 2 components
    %mu = [rand,rand]'; % mean uniform between 0 and 1
    %std = [rand*0.3, rand*0.3]; % std uniform between 0 and 0.3

    mu = [20*rand-10,20*rand-10]'; % mean uniform between -10 and 10
    std = [rand*10+0.0001, rand*10+0.0001]; % std uniform between 0.0001 and 10
    
    sigma = reshape(std.^2,1,1,2);

    p = rand; %mixing proportions
    gmm = gmdistribution(mu,sigma,[p,1-p]);
    
    % generate numPoints sample points from the distribution
    pts = random(gmm, numPoints);
    
    %compute theoretical mean and covariance
    theo_mean = p*mu(1)+(1-p)*mu(2);
    theo_var = p*(std(1)^2+mu(1)^2)+(1-p)*(std(2)^2+mu(2)^2)-theo_mean^2;
 
  
    % Our algorithm
    tic;
    % compute sample mean and covariance
    samp_mean = mean(pts);
    samp_var = var(pts);
    samp_std = sqrt(samp_var);
    % normalize to zero mean unit norm
    pts2 = (pts - samp_mean)/samp_std;
    mu_norm = (mu - samp_mean)/sqrt(samp_var);
    std_norm = std/samp_std;
    sigma_norm = reshape(std_norm.^2, 1, 1, 2);
    gmm_norm = gmdistribution(mu_norm,sigma_norm,[p,1-p]);
    
    % cycle over hypotheses and evaluate
    loglik = [];
    length_m1 = length(m1_array);
    for m = 1:num_hyp 
        %pdf_tmp = pdf(gmm_hyp{m}, pts2); % ci sono un sacco di ripetizioni
        pdf_tmp = w1_array(m)/std1_array(m)* ...
                  exp(-0.5*(pts2-m1_array(m)).^2/std1_array(m)^2) + ...
                  (1-w1_array(m))/std2_array(m)* ...
                  exp(-0.5*(pts2-m2_array(m)).^2/std2_array(m)^2);
        loglik(m) = sum(log(pdf_tmp));
    end;
    % choose hypothesis with max log likelihood
    [val_hyp, index_hyp] = max(loglik);
    mean1_best = m1_array(index_hyp);
    mean2_best = m2_array(index_hyp);
    std1_best = std1_array(index_hyp);
    std2_best = std2_array(index_hyp);
    w1_best = w1_array(index_hyp);
    % denormalize 
    std1_est = std1_best; %*samp_std;
    std2_est = std2_best; %*samp_std;
    mean1_est = mean1_best;%*samp_std+samp_mean;
    mean2_est = mean2_best;%*samp_std+samp_mean;
    w1_est = w1_best;
    ours_timer = toc;       
    
      % Standard EM fit to compare
    tic;
    %maxIterationNumber =  100;%6000;
    %minIncrementPerc_EM1D = 0.004;
    %usePrecalculatedMeansVars = 0;
    %[gmm_our_em_output, lkHood_em, ~] = executeWholeEM(pts', 2, usePrecalculatedMeansVars, p, mu, std, maxIterationNumber, minIncrementPerc_EM1D);
    options = statset('MaxIter',1500, 'TolFun',1e-3 );
    gmm_output = fitgmdist(pts2, 2,  'RegularizationValue', 0.0001, 'ProbabilityTolerance', 1e-6, 'Options',options);
    em_timer = toc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test time performance
    ours_timer_arr(i) = ours_timer;
    em_mtc_arr(i) = em_timer;  

    % Test precision performance
    mu_best = [mean1_est, mean2_est]';
    std_best = [std1_est, std2_est];
    sigma_best = reshape(std_best.^2, 1, 1, 2);
    p_best = w1_est;
    gmm_best = gmdistribution(mu_best, sigma_best, [p_best, 1-p_best]);
    
    gmmLoglik_em_arr(i) = sum(log(pdf(gmm_output, pts2) + 0.000001));
    gmmLoglik_ours_arr(i) = sum(log(pdf(gmm_best, pts2) + 0.000001));
    
    % convert the gmm distribution format for convenience
    gmm_input.numComponents = gmm.NumComponents;
    gmm_input.priorProbs(1) = gmm.ComponentProportion(1);
    gmm_input.priorProbs(2) = gmm.ComponentProportion(2);
    gmm_input.singleComponent(1).mean = gmm.mu(1);
    gmm_input.singleComponent(2).mean = gmm.mu(2);
    gmm_input.singleComponent(1).rxx = gmm.Sigma(1);
    gmm_input.singleComponent(2).rxx = gmm.Sigma(2);

    % gmm_em_output.numComponents = gmm_output.NumComponents;
    % gmm_em_output.priorProbs(1) = gmm_output.ComponentProportion(1);
    % gmm_em_output.priorProbs(2) = gmm_output.ComponentProportion(2);
    % gmm_em_output.singleComponent(1).mean = gmm_output.mu(1);
    % gmm_em_output.singleComponent(2).mean = gmm_output.mu(2);
    % gmm_em_output.singleComponent(1).rxx = gmm_output.Sigma(1);
    % gmm_em_output.singleComponent(2).rxx = gmm_output.Sigma(2);

    gmm_best_maxLik.numComponents = gmm_best.NumComponents;
    gmm_best_maxLik.priorProbs(1) = gmm_best.ComponentProportion(1);
    gmm_best_maxLik.priorProbs(2) = gmm_best.ComponentProportion(2);
    gmm_best_maxLik.singleComponent(1).mean = gmm_best.mu(1);
    gmm_best_maxLik.singleComponent(2).mean = gmm_best.mu(2);
    gmm_best_maxLik.singleComponent(1).rxx = gmm_best.Sigma(1);
    gmm_best_maxLik.singleComponent(2).rxx = gmm_best.Sigma(2);

    %gmmDist_ours_arr(i) = normalizedL2DistanceGMMs(gmm_input, gmm_best_maxLik);
    %gmmDist_em_arr(i) = normalizedL2DistanceGMMs(gmm_input, gmm_our_em_output);% gmm_em_output);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gmmDist_ours_arr(i) = dist_kld(gmm_norm, gmm_best);
    gmmDist_em_arr(i) = dist_kld(gmm_norm, gmm_output);
    
    
    if ui,
        ours_mtc = ours_timer/i;
        em_mtc = em_timer/i;
        X = ['Ours MTC: ',num2str(ours_mtc) , '    EM MTC: ', num2str(em_mtc)];
        disp(X);  
        
        if gmmDist_em_arr(i) > 1,
            figure(3);
            mu_best = [mean1_est,mean2_est]';
            std_best = [std1_est, std2_est];
            sigma_best = reshape(std_best.^2,1,1,2);
            p_best = w1_est;
            gmm_best = gmdistribution(mu_best,sigma_best,[p_best,1-p_best]);
            maximo = max(pts);
            minimo = min(pts);
            range = minimo:(maximo-minimo)/1000:maximo;
            plot(range',pdf(gmm,range'),'r-',range',pdf(gmm_best,range'),'g',range',pdf(gmm_fit,range'),'b');
            legend('True','Ours','EM');
            keyboard;
        end;
    end;
end;


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
gmmLoglik_em_min = max(gmmLoglik_em_arr);

%quartiles
qt_dist_ours = quantile(gmmDist_ours_arr,q);
qt_dist_em = quantile(gmmDist_em_arr,q);
qt_time_ours = quantile(ours_timer_arr,q);
qt_time_em = quantile(em_mtc_arr, q);
qt_llik_ours = quantile(gmmLoglik_ours_arr,q);
qt_llik_em = quantile(gmmLoglik_em_arr,q);



% Save output to screen
saveOutput(1, numRepetitions, numPoints, mean_bins, std_bins, ...
              ours_timer_mean, ours_timer_std, ours_timer_max, ...
              em_mtc_mean, em_mtc_std, em_mtc_max, ...
              gmmDist_ours_mean, gmmDist_ours_std, gmmDist_ours_max_dist, ...
              gmmDist_em_mean, gmmDist_em_std, gmmDist_em_max_dist, ...
              gmmLoglik_ours_mean, gmmLoglik_ours_std, gmmLoglik_ours_min, ...
              gmmLoglik_em_mean, gmmLoglik_em_std, gmmLoglik_em_min, ...
              qt_dist_ours, qt_dist_em, qt_time_ours, qt_time_em, ...
              qt_llik_ours, qt_llik_em);

% Save ouput to file
output_file_name = [ouputFolder '/' outputFileName '.txt'];
fid = fopen(output_file_name, 'w');
saveOutput(fid, numRepetitions, numPoints, mean_bins, std_bins, ...
                ours_timer_mean, ours_timer_std, ours_timer_max, ...
                em_mtc_mean, em_mtc_std, em_mtc_max, ...
                gmmDist_ours_mean, gmmDist_ours_std, gmmDist_ours_max_dist, ...
                gmmDist_em_mean, gmmDist_em_std, gmmDist_em_max_dist, ...
                gmmLoglik_ours_mean, gmmLoglik_ours_std, gmmLoglik_ours_min, ...
                gmmLoglik_em_mean, gmmLoglik_em_std, gmmLoglik_em_min, ...
                qt_dist_ours, qt_dist_em, qt_time_ours, qt_time_em, ...
                qt_llik_ours, qt_llik_em);
fclose(fid);

% Draw boxplot
fig = figure(1);
boxplot([gmmDist_ours_arr, gmmDist_em_arr]);
xlabel('Ours (1)                                                          EM (2)');
ylabel('normalizedL2DistanceGMMs');
titleParameters = sprintf('mean-bins = %d; std-bins = %d', mean_bins, std_bins);
title(['Normalized L2 distance to the ground thruth - ' titleParameters])

% Save boxplot to file
if (saveBoxPlot2File == 1)
    resolution = '-r300';
    format = '-depsc'; %EPS color; or format = '-depsc2' for EPS Level 2 color
    print(fig, resolution, format, '-painters', [ouputFolder '/' outputFileName '.eps']);
    format = '-djpeg'; %jpeg-24 bit
    print(fig, resolution, format, '-painters', [ouputFolder '/' outputFileName '.jpg']);
end

return;



function saveOutput(fileID, numRepetitions, numPoints, mean_bins, std_bins, ...
                    ours_timer_mean, ours_timer_std, ours_timer_max, ...
                    em_mtc_mean, em_mtc_std, em_mtc_max, ...
                    gmmDist_ours_mean, gmmDist_ours_std, gmmDist_ours_max_dist, ...
                    gmmDist_em_mean, gmmDist_em_std, gmmDist_em_max_dist, ...
                    gmmLoglik_ours_mean, gmmLoglik_ours_std, gmmLoglik_ours_min, ...
                    gmmLoglik_em_mean, gmmLoglik_em_std, gmmLoglik_em_min, ...
                    qt_dist_ours, qt_dist_em, qt_time_ours, qt_time_em, ...
                    qt_llik_ours, qt_llik_em)

fprintf(fileID, '\n');
fprintf(fileID, '-------------------------------------------');
fprintf(fileID, '\n');
fprintf(fileID, '\n');
 
fprintf(fileID, '--> Trials over %d runs:\n', numRepetitions);
fprintf(fileID, '- Using numPoints = %d\n', numPoints);
fprintf(fileID, '- Considering mean_bins = %d\n', mean_bins);
fprintf(fileID, '- Considering std_bins =  %d\n', std_bins);
fprintf(fileID, '\n');
fprintf(fileID, 'ours_timer = \t\t%f\n', ours_timer_mean);
fprintf(fileID, 'ours_timer_std = \t%f\n', ours_timer_std);
fprintf(fileID, 'ours_timer_max = \t%f\n', ours_timer_max);
fprintf(fileID, 'ours_timer_lowtail = \t%f\n', qt_time_ours(1));
fprintf(fileID, 'ours_timer_lowquart = \t%f\n', qt_time_ours(2));
fprintf(fileID, 'ours_timer_median = \t%f\n', qt_time_ours(3));
fprintf(fileID, 'ours_timer_highquart = \t%f\n', qt_time_ours(4));
fprintf(fileID, 'ours_timer_hightail = \t%f\n', qt_time_ours(5));
fprintf(fileID, '-\n');
fprintf(fileID, 'em_mtc = \t\t%f\n', em_mtc_mean);
fprintf(fileID, 'em_mtc_std = \t\t%f\n', em_mtc_std);
fprintf(fileID, 'em_mtc_max = \t\t%f\n', em_mtc_max);
fprintf(fileID, 'em_timer_lowtail = \t%f\n', qt_time_em(1));
fprintf(fileID, 'em_timer_lowquart = \t%f\n', qt_time_em(2));
fprintf(fileID, 'em_timer_median = \t%f\n', qt_time_em(3));
fprintf(fileID, 'em_timer_highquart = \t%f\n', qt_time_em(4));
fprintf(fileID, 'em_timer_hightail = \t%f\n', qt_time_em(5));
fprintf(fileID, '\n');
fprintf(fileID, '\n');
fprintf(fileID, 'gmmDist_ours = \t\t%f\n', gmmDist_ours_mean);
fprintf(fileID, 'gmmDist_ours_std = \t%f\n', gmmDist_ours_std);
fprintf(fileID, 'gmmDist_ours_max_dist = %f\n', gmmDist_ours_max_dist);
fprintf(fileID, 'ours_dist_lowtail = \t%f\n', qt_dist_ours(1));
fprintf(fileID, 'ours_dist_lowquart = \t%f\n', qt_dist_ours(2));
fprintf(fileID, 'ours_dist_median = \t%f\n', qt_dist_ours(3));
fprintf(fileID, 'ours_dist_highquart = \t%f\n', qt_dist_ours(4));
fprintf(fileID, 'ours_dist_hightail = \t%f\n', qt_dist_ours(5));
fprintf(fileID, '-\n');
fprintf(fileID, 'gmmDist_em = \t\t%f\n', gmmDist_em_mean);
fprintf(fileID, 'gmmDist_em_std = \t%f\n', gmmDist_em_std);
fprintf(fileID, 'gmmDist_em_max_dist = \t%f\n', gmmDist_em_max_dist);
fprintf(fileID, 'em_dist_lowtail = \t%f\n', qt_dist_em(1));
fprintf(fileID, 'em_dist_lowquart = \t%f\n', qt_dist_em(2));
fprintf(fileID, 'em_dist_median = \t%f\n', qt_dist_em(3));
fprintf(fileID, 'em_dist_highquart = \t%f\n', qt_dist_em(4));
fprintf(fileID, 'em_dist_hightail = \t%f\n', qt_dist_em(5));
fprintf(fileID, '\n');
fprintf(fileID, '\n');
fprintf(fileID, 'gmmLoglik_ours = \t\t%f\n', gmmLoglik_ours_mean);
fprintf(fileID, 'gmmLoglik_ours_std = \t%f\n', gmmLoglik_ours_std);
fprintf(fileID, 'gmmLoglik_ours_min = %f\n', gmmLoglik_ours_min);
fprintf(fileID, 'ours_llik_lowtail = \t%f\n', qt_llik_ours(1));
fprintf(fileID, 'ours_llik_lowquart = \t%f\n', qt_llik_ours(2));
fprintf(fileID, 'ours_llik_median = \t%f\n', qt_llik_ours(3));
fprintf(fileID, 'ours_llik_highquart = \t%f\n', qt_llik_ours(4));
fprintf(fileID, 'ours_llik_hightail = \t%f\n', qt_llik_ours(5));
fprintf(fileID, '-\n');
fprintf(fileID, 'gmmLoglik_em = \t\t%f\n', gmmLoglik_em_mean);
fprintf(fileID, 'gmmLoglik_em_std = \t%f\n', gmmLoglik_em_std);
fprintf(fileID, 'gmmLoglik_em_min = \t%f\n', gmmLoglik_em_min);
fprintf(fileID, 'em_llik_lowtail = \t%f\n', qt_llik_em(1));
fprintf(fileID, 'em_llik_lowquart = \t%f\n', qt_llik_em(2));
fprintf(fileID, 'em_llik_median = \t%f\n', qt_llik_em(3));
fprintf(fileID, 'em_llik_highquart = \t%f\n', qt_llik_em(4));
fprintf(fileID, 'em_llik_hightail = \t%f\n', qt_llik_em(5));
fprintf(fileID, '-------------------------------------------');
fprintf(fileID, '\n'); 
fprintf(fileID, '\n');


return;







