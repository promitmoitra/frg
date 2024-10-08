cd '/home/decision_lab/work/github/frg/data/spec_data/';
clear;

state_conditions = {'pre','pre_short','pre_long','post_short','post_long'};
spec_vars = {'times','freqs','ersp','exponents','offsets','bp','ap_bp'};
freq_bands = {'delta','theta','alpha','lbeta','hbeta','lgamma','hgamma'};

%%
subid = 43000;
chan_label = 'FZ';
spec_data = load(sprintf('spec_%d_%s',subid,chan_label));
spec_data = spec_data.spec_data;

%% Fit gev dist: Focusing on alpha for now:
clc;
fprintf([repmat('=',1,80),'\n',char(string(subid)),'\n',repmat('=',1,80),'\n'])
pre_alpha_pdf = fitdist(spec_data.pre.alpha_bp','gev')
pre_alpha_ap_pdf = fitdist(spec_data.pre.alpha_ap_bp','gev')
% pre_short_alpha_pd = fitdist(spec_data.pre_short.alpha_bp','gev')
% pre_long_alpha_pd = fitdist(spec_data.pre_long.alpha_bp','gev')
% post_short_alpha_pd = fitdist(spec_data.post_short.alpha_bp','gev')
% post_long_alpha_pd = fitdist(spec_data.post_long.alpha_bp','gev')
%%
% x=spec_data.pre.times;y=spec_data.pre.freqs;
% fig=figure('Position',[380 380 1280 720]);
% surf(x/1000,log10(y),log10(spec_data.pre.ersp),'FaceColor','interp','FaceAlpha',0.5,'EdgeColor','flat')
% xl = xlabel('$Time\ [s]$','Interpreter','latex','FontSize',14);
% set(xl,'Units','Normalized');set(xl,'Position',[0.75,0.03,0]);set(xl,'Rotation',30);
% yl = ylabel('$Frequency\ [log_{10} (Hz)]$','Interpreter','latex','FontSize',14);
% set(yl,'Units','Normalized');set(yl,'Position',[0.25,0.05,0]);set(yl,'Rotation',-9);
% zl = zlabel('$Power\ [log_{10}({\mu}V^2)]$','Interpreter','latex','FontSize',14);
% colorbar; axis tight;
% view(-240,40);

% figure; ax=gca; c=ax.ColorOrder;
% ax.XScale='log'; ax.YScale='log';
% hold on;
% 
% histogram(ax,spec_data.pre.alpha_bp,'BinMethod','fd','FaceColor',c(1,:),...
%           'Normalization','pdf','DisplayName','pre game');
% x = sort(spec_data.pre.alpha_bp);
% y = pdf(pre_alpha_pdf,x);
% area(ax,x,y,'FaceColor',c(1,:),'FaceAlpha',0.5,'HandleVisibility','on','DisplayName','pre game')
% 
% histogram(ax,spec_data.pre.alpha_ap_bp,'BinMethod','fd','FaceColor',c(2,:),...
%           'Normalization','pdf','DisplayName','pre game');
% x = sort(spec_data.pre.alpha_ap_bp);
% y = pdf(pre_alpha_ap_pdf,x);
% area(ax,x,y,'FaceColor',c(2,:),'FaceAlpha',0.5,'HandleVisibility','on','DisplayName','pre game ap')
%
% % histogram(ax,spec_data.pre_short.alpha_bp,'BinMethod','fd','FaceColor',c(2,:),...
% %           'Normalization','pdf','DisplayName','pre stress short');
% x = sort(spec_data.pre_short.alpha_bp);
% y = pdf(pre_short_alpha_pd,x);
% area(ax,x,y,'FaceColor',c(2,:),'FaceAlpha',0.5,'HandleVisibility','on','DisplayName','pre stress short')
% 
% % % histogram(ax,spec_data.pre_long.alpha_bp,'BinMethod','fd','FaceColor',c(2,:),...
% % %           'Normalization','pdf','DisplayName','pre stress long');
% x = sort(spec_data.pre_long.alpha_bp);
% y = pdf(pre_long_alpha_pd,x);
% area(ax,x,y,'FaceColor',c(3,:),'FaceAlpha',0.5,'HandleVisibility','on','DisplayName','pre stress long')
% % 
% % % histogram(ax,spec_data.post_short.alpha_bp,'BinMethod','fd','FaceColor',c(3,:),...
% % %           'Normalization','pdf','DisplayName','post stress short');
% x = sort(spec_data.post_short.alpha_bp);
% y = pdf(post_short_alpha_pd,x);
% area(ax,x,y,'FaceColor',c(4,:),'FaceAlpha',0.5,'HandleVisibility','on','DisplayName','post stress short')
% 
% % % histogram(ax,spec_data.post_long.alpha_bp,'BinMethod','fd','FaceColor',c(3,:),...
% % %           'Normalization','pdf','DisplayName','post stress long');
% x = sort(spec_data.post_long.alpha_bp);
% y = pdf(post_long_alpha_pd,x);
% area(ax,x,y,'FaceColor',c(5,:),'FaceAlpha',0.5,'HandleVisibility','on','DisplayName','post stress long')
%
% legend;
% hold off;

% figure; ax=gca; c=ax.ColorOrder; hold on;
% 
% plot(ax,spec_data.pre.times,spec_data.pre.alpha_bp,'Color',c(1,:),'LineWidth',1,'HandleVisibility','on','DisplayName','pre')
% plot(ax,spec_data.pre.times,spec_data.pre.alpha_ap_bp,'Color',c(1,:),'LineWidth',1,'LineStyle','--','HandleVisibility','on','DisplayName','pre ap')
% 
% plot(ax,spec_data.pre_short.times,spec_data.pre_short.alpha_bp,'Color',c(2,:),'LineWidth',1,'HandleVisibility','on','DisplayName','pre short')
% plot(ax,spec_data.pre_short.times,spec_data.pre_short.alpha_ap_bp,'Color',c(2,:),'LineWidth',1,'LineStyle','--','HandleVisibility','on','DisplayName','pre short ap')
% 
% plot(ax,spec_data.pre_long.times,spec_data.pre_long.alpha_bp,'Color',c(3,:),'LineWidth',1,'HandleVisibility','on','DisplayName','pre long')
% plot(ax,spec_data.pre_long.times,spec_data.pre_long.alpha_ap_bp,'Color',c(3,:),'LineWidth',1,'LineStyle','--','HandleVisibility','on','DisplayName','pre long ap')
% 
% plot(ax,spec_data.post_short.times,spec_data.post_short.alpha_bp,'Color',c(4,:),'LineWidth',1,'HandleVisibility','on','DisplayName','post short')
% plot(ax,spec_data.post_short.times,spec_data.post_short.alpha_ap_bp,'Color',c(4,:),'LineWidth',1,'LineStyle','--','HandleVisibility','on','DisplayName','post short ap')
% 
% plot(ax,spec_data.post_long.times,spec_data.post_long.alpha_bp,'Color',c(5,:),'LineWidth',1,'HandleVisibility','on','DisplayName','post long')
% plot(ax,spec_data.post_long.times,spec_data.post_long.alpha_ap_bp,'Color',c(5,:),'LineWidth',1,'LineStyle','--','HandleVisibility','on','DisplayName','post long ap')
% 
% legend;
% hold off;
