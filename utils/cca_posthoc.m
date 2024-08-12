function [pval,wald,pval_adj] = cca_posthoc(lme,periods,area2test,ar,fig2do)

emm = emmeans(lme,{'pair' 'period'});
% h=emmip(emm,'pair | period')

%- renaming area to show interaction
area2plot = area2test(~ismember(1:length(area2test),ar));
for i = 1 : length(area2plot)
    area2plot{i} = [area2test{ar}(2:end) '-' area2plot{i}(2:end)];
end
if fig2do
    figure;
end

clear pval wald pval_adj
for pe = 1 : length(periods)
    pval{pe}=NaN(length(area2test),length(area2test));
    wald{pe}=NaN(length(area2test),length(area2test));

    for ar1 = 1 : length(area2test)
        for ar2 = 1 : length(area2test)
            if ar1 > ar2 & ar1 ~= ar & ar2 ~= ar
                name1 = [area2test{ar} '-' area2test{ar1}];
                name1b = [area2test{ar1} '-' area2test{ar}];
                name2 = [area2test{ar} '-' area2test{ar2}];
                name2b = [area2test{ar2} '-' area2test{ar}];

                L_post = (strcmp(emm.table.period,periods{pe}) & ( strcmp(emm.table.pair,name1) |  strcmp(emm.table.pair,name1b))   )' ...
                    - (strcmp(emm.table.period,periods{pe}) & ( strcmp(emm.table.pair,name2) |  strcmp(emm.table.pair,name2b)) )';
                H0_post = contrasts_wald(lme, emm, L_post);
                pval{pe}(ar1,ar2) = H0_post.pVal;
                wald{pe}(ar1,ar2) = H0_post.Wald;
            end
        end
    end

    %- remove the area considered!
    pval{pe} = pval{pe}(~ismember(1:length(area2test),ar),~ismember(1:length(area2test),ar));
    wald{pe} = wald{pe}(~ismember(1:length(area2test),ar),~ismember(1:length(area2test),ar));

    %- FDR correction
    [h, crit_p, pval_adj{pe}]=fdr_bh(pval{pe},.05);
    
    pval_adj{pe}(pval_adj{pe}>1)=1;
    maxWald(pe)=max(max(wald{pe}));
end

max2plot = ceil(max(maxWald));

for pe = 1 : length(periods)

    % plot
    if fig2do
        colorWald = cbrewer('seq', 'Purples', 100);
        subplot(1,length(periods),pe)
        imagesc(wald{pe}',[0 max2plot]);colormap(colorWald(1:70,:));colorbar;
        for ar1 = 1 : length(area2plot)
            for ar2 = 1 : length(area2plot)
                if ~isnan(wald{pe}(ar1,ar2))
                    n=1;p=0;
                    while p == 0
                        n = n + 1;
                        p = round(pval_adj{pe}(ar1,ar2),n);
                    end
                    p = round(pval_adj{pe}(ar1,ar2),n+1);
                    if pval_adj{pe}(ar1,ar2)<.05
                        text(ar1,ar2,num2str(p),'HorizontalAlignment','center','FontWeight','bold')
                    else
                        text(ar1,ar2,num2str(p),'HorizontalAlignment','center')
                    end
                end
            end
        end
        set(gca,'Xtick',1:length(area2plot),'XtickLabel',area2plot,'Ytick',1:length(area2plot),'YtickLabel',area2plot,'FontSize',18)
        ylim([0.5 length(area2plot)-.5])
        xlim([1.5 length(area2plot)+.5])
    end

end


