%% MAKE_FIGURES  — draws Fig. 1 to Fig. 7 and exports them at 600 dpi
%
%  Run run_simulation first. Figures are exported as PNG at 600 dpi, which is
%  what a journal expects. Do NOT screen-capture the figure windows: a capture
%  is 72-96 dpi and prints blurred.
%
%  If you prefer to paste directly into Word, use the figure window menu
%  Edit > Copy Figure after setting Edit > Copy Options to a white background
%  and metafile format. That pastes as vector art and stays sharp when scaled.

clear; clc; close all;
load('results.mat');
set(0,'DefaultFigureColorMap',gray);
set(groot,'defaultFigureColor','w', ...
    'defaultAxesColor','w', ...
    'defaultAxesXColor','k','defaultAxesYColor','k', ...
    'defaultTextColor','k','defaultAxesGridColor',[.5 .5 .5]);
T = cfg.T;  PL = cfg.P_list;  NMAX = cfg.NMAX;
GY = gray(numel(PL)+2);  GY = GY(1:numel(PL),:);
EXPORT = @(h,name) exportgraphics(h, name, 'Resolution', 600);
set(0,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');

%% ---------- Fig. 1 : Cho & Cho (2020) replot ----------
pl   = [1 2 3 4 6 8 12 24];
sp   = [24 12 8 6 4 3 2 1];
rmin = [0.99 1.98 4.71 0.54 4.64 6.27 4.32 0.47];
rmea = [15.94 16.58 16.21 17.68 16.99 17.18 17.30 16.98];
rmax = [1049.26 347.54 138.17 84.44 64.47 29.55 36.44 78.38];

h = figure('Color','w','Units','centimeters','Position',[2 2 14 8.6]);
hold on; box on; grid on;
plot(pl, rmax, 'ko-',  'LineWidth',cfg.lw, 'MarkerFaceColor','k','MarkerSize',5);
plot(pl, rmea, 'ks--', 'LineWidth',cfg.lw, 'MarkerFaceColor','w','MarkerSize',5);
plot(pl, rmin, 'k^:',  'LineWidth',1.0,    'MarkerFaceColor',[.6 .6 .6],'MarkerSize',4);
yline(30,'k-.','30 min','LineWidth',0.8,'FontSize',cfg.fs-1,'LabelHorizontalAlignment','left');
set(gca,'XScale','log','YScale','log','FontSize',cfg.fs,'XTick',pl, ...
        'XTickLabel',compose('%d\\times%d',pl',sp'),'GridAlpha',0.15);
xtickangle(45);
xlabel('Orbital planes \times satellites per plane  (24 satellites)','FontSize',cfg.fs);
ylabel('Revisit time (min)','FontSize',cfg.fs);
legend({'Maximum revisit time','Mean revisit time','Minimum revisit time'}, ...
       'Box','off','FontSize',cfg.fs-1,'Location','northeast');
hold off;  EXPORT(h,'fig1.png');

%% ---------- Fig. 2 : access geometry ----------
%  Ground-track spacing for a gamma = 15/1 repeat orbit is 360/15 = 24.0 deg
%  by definition. Compute it from the NODAL period and the NODAL day, not
%  from the Keplerian period and the sidereal day; the latter omits J2 nodal
%  regression and gives 23.67 deg, which is inconsistent with the repeat
%  condition the baseline orbit satisfies.
walk = info.track_walk_deg;
half = info.halfwin_deg;
fprintf('Fig.2  ground-track spacing %.2f deg, access window +-%.2f deg (ratio %.2f)\n', ...
        walk, half, walk/(2*half));

h = figure('Color','w','Units','centimeters','Position',[2 2 15.5 7.2]);
hold on;
lbl = containers.Map({-3,-1,1,3}, ...
      {'rev {\itn}-1','rev {\itn}','rev {\itn}+1','rev {\itn}+2'});
for m = [-5 -3 -1 1 3 5]
    x = m*walk/2;
    plot([x x],[0.46 0.94],'-','Color',[.58 .58 .58],'LineWidth',1.7);
    if isKey(lbl,m), text(x,0.99,lbl(m),'HorizontalAlignment','center','FontSize',cfg.fs); end
end
patch([-half half half -half],[0.40 0.40 0.79 0.79],[.13 .13 .13],'EdgeColor','none');
plot(0,0.70,'wp','MarkerSize',11,'MarkerFaceColor','w');

% dimension lines, both in axis coordinates
plot([-half half],[0.34 0.34],'k-','LineWidth',1.2);
plot(-half,0.34,'k<','MarkerSize',4,'MarkerFaceColor','k');
plot( half,0.34,'k>','MarkerSize',4,'MarkerFaceColor','k');
text(half+3, 0.34, sprintf('access window  %.2f\\circ', 2*half), ...
     'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',cfg.fs);

plot([-walk/2 walk/2],[0.17 0.17],'k-','LineWidth',1.2);
plot(-walk/2,0.17,'k<','MarkerSize',4,'MarkerFaceColor','k');
plot( walk/2,0.17,'k>','MarkerSize',4,'MarkerFaceColor','k');
text(walk/2+3, 0.17, sprintf('track spacing  %.1f\\circ', walk), ...
     'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',cfg.fs);

xlim([-2.9*walk 2.9*walk]); ylim([0.05 1.12]);
set(gca,'YTick',[],'FontSize',cfg.fs,'Box','off', ...
        'XTick',(-5:1:5)*walk/2, 'XTickLabel',compose('%.0f',(-5:1:5)'*walk/2));
xlabel('Relative longitude (deg)','FontSize',cfg.fs);
hold off;  EXPORT(h,'fig2.png');

%% ---------- Fig. 3 : single-satellite clustering ----------
h = figure('Color','w','Units','centimeters','Position',[2 2 15 9]);
subplot(2,1,1); hold on; box on;
area(t1, double(u1), 'FaceColor',[.2 .2 .2], 'EdgeColor','none');
for k = 1:5, xline(k*24,':','Color',[.7 .7 .7],'LineWidth',0.7); end
xlim([0 6*24]); ylim([0 1.3]);
set(gca,'YTick',[0 1],'YTickLabel',{'no','yes'},'FontSize',cfg.fs);
xlabel('Time (h)','FontSize',cfg.fs); ylabel('Access','FontSize',cfg.fs);
title('(a) Access timeline, single satellite','FontWeight','normal','FontSize',cfg.fs);
hold off;

subplot(2,1,2); hold on; box on; grid on;
histogram(iv1, 0:0.5:25, 'FaceColor',[.55 .55 .55], 'EdgeColor','k','LineWidth',0.4);
set(gca,'YScale','log','FontSize',cfg.fs,'GridAlpha',0.15);
xlabel('Interval between consecutive accesses (h)','FontSize',cfg.fs);
ylabel('Count','FontSize',cfg.fs);
title('(b) Bimodal interval distribution','FontWeight','normal','FontSize',cfg.fs);
hold off;  EXPORT(h,'fig3.png');

%% ---------- Fig. 4 : daily repeatability ----------
show = [find(PL==5) find(PL==10) find(PL==20)];
mk = {'o','s','^'};
h = figure('Color','w','Units','centimeters','Position',[2 2 14 7.6]);
hold on; box on; grid on;
lbl = cell(1,numel(show));
for i = 1:numel(show)
    ii = show(i);  dv = daily(ii,:);
    plot(1:cfg.days, dv, ['-' mk{i}], 'Color',GY(2*i-1,:), 'LineWidth',cfg.lw, ...
         'MarkerFaceColor',GY(2*i-1,:), 'MarkerSize',5);
    lbl{i} = sprintf('P=%d, s=%d  (spread %.1f%%)', peace(ii).P, peace(ii).s, ...
                     100*(max(dv)-min(dv))/min(dv));
end
set(gca,'YScale','log','FontSize',cfg.fs,'GridAlpha',0.15);
xlabel('Simulation day','FontSize',cfg.fs);
ylabel('Daily maximum observation gap (min)','FontSize',cfg.fs);
legend(lbl,'Box','off','FontSize',cfg.fs-1,'Location','best');
hold off;  EXPORT(h,'fig4.png');

%% ---------- Fig. 5 : degradation ----------
mks = 'osd^vp';
h = figure('Color','w','Units','centimeters','Position',[2 2 14 9]);
hold on; box on; grid on;
lbl = cell(1,numel(PL));
for ii = 1:numel(PL)
    plot(0:NMAX, G(ii,:,1), ['-' mks(ii)], 'Color',GY(ii,:), 'LineWidth',cfg.lw, ...
         'MarkerFaceColor',GY(ii,:), 'MarkerSize',4.5);
    lbl{ii} = sprintf('P=%d, s=%d', peace(ii).P, peace(ii).s);
end
yline(30,'k--','LineWidth',1.0);
text(0.1, 26, 'W = 30 min (stated target)','FontSize',cfg.fs-1);
set(gca,'YScale','log','FontSize',cfg.fs,'GridAlpha',0.15);
xlabel('Satellites lost, {\itn}','FontSize',cfg.fs);
ylabel('{\itG}_{max}({\itn})  (min)','FontSize',cfg.fs);
legend(lbl,'Box','off','FontSize',cfg.fs-1,'NumColumns',2,'Location','southeast');
hold off;  EXPORT(h,'fig5.png');

%% ---------- Fig. 6 : Pareto front ----------
jw = find(cfg.W_list == cfg.W_main);
x = arrayfun(@(p) max(p.mrts), peace);
y = Ncrit(:, jw)';
dom = false(1,numel(PL));
for i = 1:numel(PL)
    for j = 1:numel(PL)
        if i~=j && x(j)<=x(i) && y(j)>=y(i) && (x(j)<x(i) || y(j)>y(i))
            dom(i) = true; break
        end
    end
end
h = figure('Color','w','Units','centimeters','Position',[2 2 13 8.6]);
hold on; box on; grid on;
for i = 1:numel(PL)
    if dom(i), mfc = 'w'; ms = 8; else, mfc = 'k'; ms = 10; end
    plot(x(i), y(i), 'o','MarkerSize',ms,'MarkerFaceColor',mfc,'MarkerEdgeColor','k','LineWidth',1.2);
    text(x(i), y(i), sprintf('  P=%d\\times%d', peace(i).P, peace(i).s),'FontSize',cfg.fs-1);
end
px = x(~dom); py = y(~dom);
[px, o] = sort(px); py = py(o);
plot(px, py, 'k-','LineWidth',1.0);
set(gca,'XScale','log','FontSize',cfg.fs,'GridAlpha',0.15);
xlabel('Peacetime maximum observation gap (min)','FontSize',cfg.fs);
ylabel(sprintf('{\\itN}_{crit}  ({\\itW} = %d min)', cfg.W_main),'FontSize',cfg.fs);
hold off;  EXPORT(h,'fig6.png');

%% ---------- Fig. 7 : re-phasing ----------
RE = 6378.137;  a = RE + cfg.alt_km;
dV_budget = 30;                 % m/s  <-- CITE A SOURCE FOR THIS
Tp = logspace(log10(1/48), log10(60), 400);
h = figure('Color','w','Units','centimeters','Position',[2 2 14 8.6]);
hold on; box on; grid on;
S_list = sort(unique(T./PL));
lbl = cell(1,numel(S_list));
for i = 1:numel(S_list)
    S = S_list(i);  dth = deg2rad(360/S);
    plot(Tp, 2*a*dth./(3*Tp*86400)*1000, '-', 'Color',GY(min(i,size(GY,1)),:), 'LineWidth',cfg.lw);
    lbl{i} = sprintf('s=%d  (\\Delta\\theta=%.0f\\circ)', S, 360/S);
end
yline(dV_budget,'k--',sprintf('\\DeltaV budget %g m/s',dV_budget), ...
      'LineWidth',1.0,'FontSize',cfg.fs-1);
xline(64/24,'k:','TacRS 64 h','LineWidth',0.9,'FontSize',cfg.fs-2,'LabelOrientation','horizontal');
xline(121/24,'k:','TacRS 121 h','LineWidth',0.9,'FontSize',cfg.fs-2,'LabelOrientation','horizontal');
set(gca,'XScale','log','YScale','log','FontSize',cfg.fs,'GridAlpha',0.15);
xlabel('Recovery (drift) time (days)','FontSize',cfg.fs);
ylabel('Round-trip \DeltaV (m/s)','FontSize',cfg.fs);
legend(lbl,'Box','off','FontSize',cfg.fs-2,'NumColumns',2,'Location','northeast');
hold off;  EXPORT(h,'fig7.png');

fprintf('exported fig1.png ... fig7.png at 600 dpi\n');
