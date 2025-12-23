function dump_results(fileID, numRepetitions, numPoints, Nw, Nm, Nv, ...
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
fprintf(fileID, '- Considering Nw = %d\n', Nw);
fprintf(fileID, '- Considering Nm =  %d\n', Nm);
fprintf(fileID, '- Considering Nv =  %d\n', Nv);
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
