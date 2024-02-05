clc;clear; close all;
load CompareResults13;
load CompareResultsStd13;
CompareResults=CompareResults13;
CompareResultsStd=CompareResultsStd13;
load CompareResults2;
load CompareResultsStd2;
CompareResults(2,:)=CompareResults2(2,:);
CompareResults(5,:)=CompareResults2(5,:);
CompareResultsStd(2,:)=CompareResultsStd2(2,:);
CompareResultsStd(5,:)=CompareResultsStd2(5,:);



% CompareResults(:,1)=[];
% CompareResultsStd(:,1)=[];


K=297; 
m_num=[1.5 2 2.5 3 3.5 4];
ms=5; ts=12;


figure; hold on;
errorbar(m_num, CompareResults(1,:), CompareResultsStd(1,:), 'bv-.', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num, CompareResults(2,:), CompareResultsStd(2,:), 'gd--', 'linewidth', 2, 'MarkerSize',ms); 
errorbar(m_num, CompareResults(3,:), CompareResultsStd(3,:), 'rp-', 'linewidth', 2, 'MarkerSize',ms); 
ylabel('Recovery Error'); xlabel('Sample Size Ratio (n / k)');box on;
legend('Group OMP', 'Group Lasso','StructOMP');
% axis([m_num(1)-round(0.25*K) m_num(end)+round(0.25*K) -0 0.3])
axis([m_num(1)-0.25 m_num(end)+0.25 0 1.6])
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
ylabel('CPU Time (Second)'); xlabel('Sample Size Ratio (n / k)');box on;
legend('Group OMP','Group Lasso', 'StructOMP','Location','NorthWest');
axis([m_num(1)-0.25 m_num(end)+0.25 0 320])
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts); 
