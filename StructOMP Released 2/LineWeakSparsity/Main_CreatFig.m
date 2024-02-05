clc;clear; close all;
load CompareResults;
load CompareResultsStd;

K=32; 
m_num=round(K*[2 2.5 3 3.5 4 4.5 5]);
ms=5; ts=12;

figure; hold on;
errorbar(m_num/K, CompareResults(1,:), CompareResultsStd(1,:), 'bv-', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(2,:), CompareResultsStd(2,:), 'gd-', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(3,:), CompareResultsStd(3,:), 'ch-', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(4,:), CompareResultsStd(4,:), 'mh-.', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(5,:), CompareResultsStd(5,:), 'kh--', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(6,:), CompareResultsStd(6,:), 'yh-', 'linewidth', 3,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(7,:), CompareResultsStd(7,:), 'rp-', 'linewidth', 2,  'MarkerSize',ms);
ylabel('Recovery Error'); xlabel('Sample Size Ratio (n / k)');box on;
legend('OMP', 'Lasso', 'GroupLasso, gs=2','GroupLasso, gs=4','GroupLasso, gs=8','GroupLasso, gs=16','StructOMP');
axis([m_num(1)/K-0.5 m_num(end)/K+0.5 0 0.4])
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts); 



figure; hold on;
errorbar(m_num/K, CompareResults(8,:), CompareResultsStd(8,:), 'bv-', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(9,:), CompareResultsStd(9,:), 'gd-', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(10,:), CompareResultsStd(10,:), 'ch-', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(11,:), CompareResultsStd(11,:), 'mh-.', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(12,:), CompareResultsStd(12,:), 'kh--', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(13,:), CompareResultsStd(13,:), 'yh-', 'linewidth', 3,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(14,:), CompareResultsStd(14,:), 'rp-', 'linewidth', 2,  'MarkerSize',ms);
ylabel('CPU Time (Second)');xlabel('Sample Size Ratio (n / k)');box on;
legend('OMP','Lasso', 'StructOMP');
axis([m_num(1)/K-0.5 m_num(end)/K+0.5 0 0.8])
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts); 

