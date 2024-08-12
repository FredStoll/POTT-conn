function [pval,wald,pval_adj] = cca_posthoc_periods(lme,periods,area2test,ar,fig2do)

emm = emmeans(lme,{'pair' 'period'});
% h=emmip(emm,'pair | period')

%- renaming area to show interaction
area2plot = area2test(~ismember(1:length(area2test),ar));
for i = 1 : length(area2plot)
    area2plot{i} = [area2test{ar}(2:end) '-' area2plot{i}(2:end)];
end

% find all the correlation with the area of interest
for i = 1 : length(emm.table.pair)
    dash = find(emm.table.pair{i}=='-');
    ar1{i} = emm.table.pair{i}(1:dash-1);
    ar2{i} = emm.table.pair{i}(dash+1:end);
end
takeme = ismember(ar1,area2test(ar)) | ismember(ar2,area2test(ar));

clear pval wald pval_adj
pval = NaN(length(periods),length(periods));
wald = NaN(length(periods),length(periods));
for pe1 = 1 : length(periods)
    for pe2 = 1 : length(periods)
        if pe1>pe2

            L_post = (strcmp(emm.table.period,periods{pe1})' & takeme) - (strcmp(emm.table.period,periods{pe2})' & takeme);

            H0_post = contrasts_wald(lme, emm, L_post);
            pval(pe1,pe2) = H0_post.pVal;
            wald(pe1,pe2) = H0_post.Wald;
        end
    end
end


%- FDR correction
[h, crit_p, pval_adj]=fdr_bh(pval,.05);

pval_adj(pval_adj>1)=1;


% plot
if fig2do
    figure;
    colorWald = cbrewer('seq', 'Purples', 100);
    imagesc(wald',[0 ceil(max(wald(:)))]);colormap(colorWald(1:70,:));colorbar;
    for pe1 = 1 : length(periods)
        for pe2 = 1 : length(periods)
            if ~isnan(wald(pe1,pe2))
                n=1;p=0;
                while p == 0
                    n = n + 1;
                    p = round(pval_adj(pe1,pe2),n);
                end
                p = round(pval_adj(pe1,pe2),n+1);
                if pval_adj(pe1,pe2)<.05
                    text(pe1,pe2,num2str(p),'HorizontalAlignment','center','FontWeight','bold')
                else
                    text(pe1,pe2,num2str(p),'HorizontalAlignment','center')
                end
            end
        end
    end
    set(gca,'Xtick',1:length(periods),'XtickLabel',periods,'Ytick',1:length(periods),'YtickLabel',periods,'FontSize',18)
    ylim([0.5 length(periods)-.5])
    xlim([1.5 length(periods)+.5])

    for pe1 = 1 : length(periods)
        avg_conn(pe1,:) = varfun(@mean, emm.table(takeme & strcmp(emm.table.period,periods{pe1})',{'Estimated_Marginal_Mean' 'SE'}))
    end
    figure;errorbar(avg_conn.mean_Estimated_Marginal_Mean,avg_conn.mean_SE);
    xlim([0 length(periods)+1])
    set(gca,'Xtick',1:length(periods),'XtickLabel',periods,'FontSize',18)
end




