clc;clear; close all;
load CompareResults;
load CompareResultsStd;
K=1133;
m_num=[1.25 1.5 1.75 2 2.25 2.5 2.75]*K; 
m_num=round(m_num);
ms=5; ts=12;


figure; hold on;
errorbar(m_num, CompareResults(1,:), CompareResultsStd(1,:), 'bv-.', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num, CompareResults(2,:), CompareResultsStd(2,:), 'gd--', 'linewidth', 2, 'MarkerSize',ms); 
errorbar(m_num, CompareResults(3,:), CompareResultsStd(3,:), 'rp-', 'linewidth', 2, 'MarkerSize',ms); 
ylabel('Recovery Error'); xlabel('Sample Size');box on;
legend('Group OMP', 'Group Lasso','StructOMP');
% axis([m_num(1)-round(0.25*K) m_num(end)+round(0.25*K) -0 0.3])
% axis([m_num(1)-0.25 m_num(end)+0.25 0 0.4])
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts); 

figure; hold on;
errorbar(m_num, CompareResults(4,:), CompareResultsStd(4,:), 'bv-.', 'linewidth', 2, 'MarkerSize',ms);
errorbar(m_num, CompareResults(5,:), CompareResultsStd(5,:), 'gd--', 'linewidth', 2, 'MarkerSize',ms);
errorbar(m_num, CompareResults(6,:), CompareResultsStd(6,:), 'rp-', 'linewidth', 2, 'MarkerSize',ms);
ylabel('CPU Time (Second)'); xlabel('Sample Size');box on;
legend('Group OMP','Group Lasso', 'StructOMP','Location','NorthWest');
% axis([m_num(1)-0.25 m_num(end)+0.25 0 80])
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts); 
