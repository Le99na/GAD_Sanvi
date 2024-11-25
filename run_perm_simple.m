% 
% dat_eng = readtable("theta_indiv_eng_betaCoeffs.csv")
% data_condition1 = table2array(dat_eng);
% 
% dat_dis = readtable("theta_indiv_diseng_betaCoeffs.csv")
% data_condition2 = table2array(dat_dis);

cor = load('rPFC_theta-gamma_pdist_correct.mat')
icor = load('rPFC_theta-gamma_pdist_icorrect.mat')

clear all; clc

[t_values2, p_values, sig_clusters] = cluster_permutation_test_Correlation2(mat_diff',behav_dat(:,4)',10000, 0.05)


addpath('C:\Users\magcam\Documents\code_matlab\permutation_analysis\permutation_analysis_MC\2D')



 a={sig_clusters.channels}
 channels = [a{1},a{2}]

rpfc = mean(resultMatrix([10,11,19,20],:),1)

x=1-pup_eng_dist
y=1-rpfc
% Scatter plot
figure;
scatter(x, y, 'filled', 'MarkerFaceAlpha', 0.6); % Semi-transparent points
hold on;

% Fit a linear regression model
mdl = fitlm(x, y);  
slope = mdl.Coefficients.Estimate(2);
intercept = mdl.Coefficients.Estimate(1);

% Regression line
x_line = linspace(min(x), max(x), 100);
y_line = slope * x_line + intercept;
plot(x_line, y_line, 'r', 'LineWidth', 2);

% Confidence intervals for the regression line
[ypred, yCI] = predict(mdl, x_line'); % yCI gives lower and upper bounds

% Shaded area for standard deviation around the line
fill([x_line, fliplr(x_line)], [yCI(:, 1)', fliplr(yCI(:, 2)')], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Spearman correlation
[r_s, p_value] = corr(x', y', 'Type', 'Spearman');

% Display Spearman correlation in title
xlabel('Subjective engagement');
ylabel('Theta power');
title(sprintf('Regression Line\nSpearman r = %.2f, p = %.3f', r_s, p_value));

% Legend
legend({'Data points', 'Regression line', '95% Confidence Interval'}, 'Location', 'best');
grid on;
hold off;