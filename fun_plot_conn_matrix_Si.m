function fun_plot_conn_matrix_Si(connectivity_matrix,S_tfd)
% Drawing of connectivity matrix
% Inputs:
%   S_tfd - structure, structure is same as produced with EEGLab
%   connectivity_matrix - matrix of connection between observed electrodes
% Outputs:
%   figure with connectivity matrix representation
%


% ?verko, Z.; Sajovic, J.; Dreven?ek, G.; Vlahini?c, S.; Rogelj, P. Generation of Oscillatory Synthetic Signal Simulating Brain Network
% Dynamics. In Proceedings of the 2021 44th International Convention on Information, Communication and Electronic Technology
% (MIPRO), MEET¡ªMicroelectronics, Electronics and Electronic Technology, Opatija, Croatia, 27 September¨C1 October 2021.
% ---------------------------------------------------------------------- 
% Copyright (2021): Zoran ?verko
%-----------------------------------------------------------------------



h=imagesc(connectivity_matrix);
% set(gca,'clim',[0 1],'xtick',1:1:EEG.nbchan,'xticklabel',{EEG.chanlocs(1:1:end).labels},'ytick',1:1:EEG.nbchan,'yticklabel',{EEG.chanlocs(1:1:end).labels},'fontsize', 6);
set(gca,'clim',[0 0.5]);

ax = ancestor(h, 'axes');
yrule = ax.YAxis;
% Change properties of the axes
ax.YTick = 1:1:S_tfd.nbchan;
% ax.YTick.label = down;
ax.YTickLabel = {32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};
% set(gca,'ytick',1:1:32);
% set(gca,'yminortick','on');
% ylim([1 32]);
% Change properties of the ruler
yrule.FontSize = 8;
% % Change properties of the label
% yL.FontSize = 8;

xrule = ax.XAxis;
% Change properties of the axes
ax.XTick = 1:1:S_tfd.nbchan;
ax.XTickLabel = {1:1:S_tfd.nbchan};
% Change properties of the ruler
xrule.FontSize = 8;
% % Change properties of the label
% yL.FontSize = 8;
xtickangle(90)
axis square
colorbar
end