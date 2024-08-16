%% Canonical Correlation Analysis on POTT dataset - Posthoc analyses
%-
%- Reproduce the analyses and figures from:
%- Stoll & Rudebeck (2024) Decision-making shapes dynamic inter-areal communication within macaque ventral frontal cortex
%- https://doi.org/10.1101/2024.07.05.602229
%-
%- Require data available at : 10.5281/zenodo.13306842
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Last updated: 2024.08

clear

%% choose your adventure!

if isunix
    path2go = '/home/fred/POTTconn/data/'; %- path where SPKpool files are!
elseif ispc
    path2go = 'R:\POTTconn\data\'; %- path where SPKpool files are!
end

ana2run = 'rewTr'; % 'main' / 'lesion' / 'alltime' / 'rewTr'
area2run = 'a12o'; %- only used for ana2run='lesion' or 'rewTr' / can be a12o/a11ml/LAI
period2run = 'Rew'; %- only used for ana2run='lesion' or 'rewTr' / can be 'Stim' 'Rew'
crosstemp = false; %- only used for ana2run='lesion' / can be true (FIGURE S5) or false (FIGURES 3/S4)
bothresid = false; %- when using residuals on both area (Figure S4A): true, otherwise false

%- colors for areas
colorsArea = cbrewer('qual', 'Paired', 12);
colorsArea = colorsArea([1:7 12 8:11],:);
nb_sub=[6 5]; %- nb of subplots

%- colors for periods
colors = cbrewer('qual','Set1',9);
colors = colors([9 3 1],:);
colors_light = cbrewer('qual','Pastel1',9);
colors_light = colors_light([9 3 1],:);

%- define time periods (always put one called BL == used as reference for some analyses)
if strcmp(ana2run,'alltime')
    periods =      { 'Fix'  'BL'       'Stim'  'Resp' 'Rew'};
    periods_bin =  { 1       2           2      4      6};
    periods_time = {[100 700] [-700 -100] [100 700] [100 700] [100 700]};
else
    periods =      { 'BL'       'Stim'   'Rew'}; %      { 'Fix'  'BL'       'Stim'  'Resp' 'Rew'};
    periods_bin =  { 2           2        6};    %      { 1       2           2      4      6};
    periods_time = {[-700 -100] [100 700] [100 700]}; % {[100 700] [-700 -100] [100 700] [100 700] [100 700]};
end

if strcmp(ana2run,'lesion') && strcmp(period2run,'Stim')
    periods = periods([1 2]);periods_bin = periods_bin([1 2]);periods_time = periods_time([1 2]);
elseif (strcmp(ana2run,'lesion') && strcmp(period2run,'Rew')) || strcmp(ana2run,'rewTr')
    periods = periods([1 3]);periods_bin = periods_bin([1 3]);periods_time = periods_time([1 3]);
end

%% load the correct data!
if strcmp(ana2run,'main')
    filename = [path2go 'FINAL_CCA_INS_8areas_perm_25sdf_25ms_raw.mat'];
    load(filename)

elseif strcmp(ana2run,'alltime')
    filename = [path2go 'FINAL_CCA_INS_8areas_perm_25sdf_25ms_raw_alltime.mat'];
    load(filename)

elseif strcmp(ana2run,'rewTr')
    LESION{1} = load([path2go 'FINAL_CCA_INS_8areas_perm_25sdf_25ms_raw_rewarded.mat']);
    LESION{2} = load([path2go 'FINAL_CCA_INS_8areas_perm_25sdf_25ms_raw_unrewarded.mat']);
    area2test = LESION{1}.area2test;
    param = LESION{1}.param;

    predic = {'rew' 'norew'};
    predic_stats = {'rewarded' 'not'};
    colors = cbrewer('qual','Dark2',8);
    colors_light = cbrewer('qual','Set2',8);
    colors = [colors(2,:) ; colors(6,:) ];
    colors_light = [colors_light(2,:) ; colors_light(6,:)];

elseif strcmp(ana2run,'lesion') && ~bothresid %- lesion analyses
    colors = cbrewer('qual','Set1',9);
    colors_light = cbrewer('qual','Pastel1',9);
    if strcmp(period2run,'Stim')
        predic = {'none' 'allstim' 'chosenproba' 'chosenjuice' 'chosenside' };
        predic_stats = {'none' 'all' 'proba' 'juice' 'side' }; %- same but fix problems with emmeans and 'interactions'!
        colors = [colors(3,:) ;  0 0 0 ; colors(7,:) ; colors(4,:) ; colors(8,:)] ;
        colors_light = [colors_light(3,:)  ; .6 .6 .6 ; colors_light(7,:) ; colors_light(4,:) ; colors_light(8,:)];
    elseif strcmp(period2run,'Rew')
        predic = {'none' 'allrew' 'rew' 'chosenjuicerew' 'chosenproba' 'chosenside' };
        predic_stats = {'none' 'all' 'rew' 'chosenjuice' 'proba' 'side' };
        colors = [colors(1,:) ; 0 0 0; colors(5,:) ; colors(4,:) ; colors(7,:) ;  colors(8,:)] ;
        colors_light = [colors_light(1,:)  ;  .6 .6 .6 ; colors_light(5,:) ; colors_light(4,:) ; colors_light(7,:) ; colors_light(8,:)];
    end

    for p = 1 : length(predic)
        if strcmp(predic{p},'none')
            LESION{p} = load([path2go 'FINAL_CCA_INS_8areas_perm_25sdf_25ms_raw.mat']);
        else
            filename = [path2go 'FINAL_CCA_INS_8areas_25sdf_25ms_resid' area2run(2:end) '_' predic{p} '.mat'];
            LESION{p} = load(filename);
        end
    end
    area2test = LESION{1}.area2test;
    param = LESION{1}.param;

elseif strcmp(ana2run,'lesion') && bothresid %-  lesion analyses for Figure S4A (control)

    colors = cbrewer('qual','Set1',9);
    colors_light = cbrewer('qual','Pastel1',9);
    if strcmp(period2run,'Stim')
        predic = {'none' 'allstim'};
        colors = [colors(3,:) ;  0 0 0 ] ;
        colors_light = [colors_light(3,:)  ; .6 .6 .6 ];
    elseif strcmp(period2run,'Rew')
        predic = {'none' 'allrew'};
        colors = [colors(1,:) ; 0 0 0; ] ;
        colors_light = [colors_light(1,:)  ;  .6 .6 .6 ];
    end

    for p = 1 : length(predic)
        if strcmp(predic{p},'none')
            LESION{p} = load([path2go 'FINAL_CCA_INS_8areas_perm_25sdf_25ms_raw.mat']);
        else
            filename = [path2go 'FINAL_CCA_INS_8areas_25sdf_25ms_resid' area2run(2:end) '_both_' predic{p} '.mat'];
            LESION{p} = load(filename);
        end
    end
    area2test = LESION{1}.area2test;
    param = LESION{1}.param;
end

%% MAIN ANALYSES (FIGURES 2/5/6/S2/S3)
if strcmp(ana2run,'main') || strcmp(ana2run,'alltime')

    %% Average CCA for each area pairs
    minVal = [];
    maxVal = [];
    CCA_sta=struct();
    for p = 1 : length(periods)

        %- bin 2 takes depending on event alignement and time windows
        bin2take = param.bin_name == periods_bin{p} & ( param.time >= periods_time{p}(1)  & param.time <= periods_time{p}(2) );

        %- extract average cca
        cca_map=[];
        cca_sta=[];
        cca_sem=[];
        for ar = 1 : length(area2test)
            for ar2 = 1 : length(area2test)
                if ar>ar2 && size(res_CCA.cvr{ar,ar2},1)>=minReplicates
                    temp = [];
                    temp_perm = [];
                    for j = 1 : size(res_CCA.cvr{ar,ar2},1)
                        temp(j,:) = diag(squeeze(res_CCA.cvr{ar,ar2}(j,:,:)));
                        temp_perm(j,:) = diag(squeeze(res_CCA.cvr_perm{ar,ar2}(j,:,:)));
                    end
                    temp = nanmean(temp(:,bin2take),2);
                    temp_perm = nanmean(temp_perm(:,bin2take),2);

                    [pval,h,stats] = signrank(temp,temp_perm,'tail','right','method','approximate');

                    % temp = nanmean(squeeze(res_CCA.cvr_perm{ar,ar2}(:,1,:)));
                    cca_map(ar,ar2)=nanmean(temp);
                    cca_sem(ar,ar2)=nanstd(temp)/sqrt(length(temp));
                    cca_sta(ar,ar2)=pval;
                else
                    cca_map(ar,ar2)=NaN;
                    cca_sem(ar,ar2)=NaN;
                    cca_sta(ar,ar2)=NaN;
                end
            end
        end

        CCA_avg.(periods{p})=cca_map; %- create a variable with the map
        CCA_sem.(periods{p})=cca_sem; %- create a variable with the map
        CCA_sta.(periods{p})=cca_sta; %- create a variable with the map

        maxVal= ceil(max(max([maxVal ; cca_map(:)]))*100)/100; %- for plotting purposes...
        minVal= floor(min(min([minVal ; cca_map(:)]))*100)/100;
    end

    %% FIGURE 2A - Brain plot
    figure;
    maximaxi = .16 ; % ceil(max(max([max_crosstemp{:}]))*100)/100;
    for p = 1 : length(periods)
        [h, crit_p, adj_p]=fdr_bh(CCA_sta.(periods{p}),.05); %- FDR correction
        temp = CCA_avg.(periods{p});
        temp(~h)=NaN; %- remove connectivity that doesn't pass threshold
        temp(isnan(temp))=0;
        temp = temp + temp';
        gr = graph(temp,area2test,'omitselfloops');
        fig = subplot(1, length(periods),p);nobrain_map(gr,area2test,[0 maximaxi],'Blues',fig)
        title([periods{p}])
    end

    %% FIGURE S3 - Graph representation (irrespective of anatomical coordinates)
    figure;
    maximaxi = .16 ; % ceil(max(max([max_crosstemp{:}]))*100)/100;
    for p = 1 : length(periods)
        temp = CCA_avg.(periods{p});
        temp(isnan(temp))=0;
        temp = temp + temp';
        gr = graph(temp,area2test,'omitselfloops');
        fig = subplot(1, length(periods),p);

        LWidths = 5*gr.Edges.Weight/max(gr.Edges.Weight);
        LWidths(LWidths<=0)=0.01;
        plot(gr,'Layout','force','WeightEffect','inverse','LineWidth',LWidths)
        title([periods{p}])
    end

    %% Derive significance from permutations
    if isfield(res_CCA,'cvr_perm') && strcmp(ana2run,'main')

        p_thr = 0.05;
        n_thr = 10;

        % one sided pvalues // step is a bit long (a minute or so..)
        res_CCA.pval = cell(length(area2test),length(area2test));
        res_CCA.psig = cell(length(area2test),length(area2test));
        for ar = 1 : length(area2test)
            for ar2 = 1 : length(area2test)
                if ar>ar2 && size(res_CCA.cvr{ar,ar2},1)>=minReplicates

                    perf = res_CCA.cvr{ar,ar2};
                    perm = res_CCA.cvr_perm{ar,ar2};

                    pval=[];
                    for p = 1 : size(perf,1)
                        perm_sub = squeeze(perm(p,:,:));
                        perf_sub = squeeze(perf(p,:,:));
                        perm_sub_nanfree = perm_sub(:);
                        perm_sub_nanfree(isnan(perm_sub_nanfree))=[];
                        allperm = repmat(perm_sub_nanfree,1,size(perf_sub,1));

                        for t = 1 : size(perf_sub,1)
                            pval(p,t,:) = (size(allperm,1) - nansum(perf_sub(t,:) > allperm))  / size(allperm,1) ;
                        end

                        %- cluster correction
                        p_value_sub = squeeze(pval(p,:,:));
                        CC = bwconncomp(p_value_sub < p_thr,4);
                        cMapPrimary = zeros(size(p_value_sub));
                        for i=1:CC.NumObjects
                            if length(CC.PixelIdxList{i}) >= n_thr
                                cMapPrimary(CC.PixelIdxList{i}) = true;
                            end
                        end
                        psig(p,:,:) = cMapPrimary;
                    end

                    res_CCA.pval{ar,ar2} = pval;
                    res_CCA.psig{ar,ar2} = psig;
                end
            end
        end

        mask = NaN(size(param.mask));
        mask(param.mask)= 1;
        mask = mask(1:size(perf_sub),1:size(perf_sub));

        %% FIGURE 5A - cross-temporal average rhos
        %- for subset of plots:
        showme = {'12l - 12o' ; '13l - 12o' ; '11ml - 12l' ; 'AI - 12o' ; '13l - 12l' ; '13m - 12l'}; nb_sub=[3 2];

        %- for all plots
        % showme = {}; nb_sub=[6 5];

        figure;x=0;
        for ar = 1 : length(area2test)
            for ar2 = 1 : length(area2test)
                pair2plot = [area2test{ar}(2:end) ' - ' area2test{ar2}(2:end)];
                if  (~isempty(showme) && ~isempty(find(ismember(showme,pair2plot))))  | (isempty(showme) && ar>ar2 && size(res_CCA.cvr{ar,ar2},1)>=minReplicates)

                    x = x +1 ;
                    subplot(nb_sub(1),nb_sub(2),x);
                    data2plot = squeeze(nanmean(res_CCA.cvr{ar,ar2},1));
                    imagesc(data2plot,'AlphaData',~isnan(data2plot));hold on;axis xy square
                    set(gca,"CLim",[0 .15],...
                        'XTick',1:length(param.time)/8:length(param.time),'XTickLabel',param.time(1:length(param.time)/8:length(param.time))-13,...
                        'YTick',1:length(param.time)/8:length(param.time),'YTickLabel',param.time(1:length(param.time)/8:length(param.time))-13)
                    contour(1:size(res_CCA.psig{ar,ar2},3),1:size(res_CCA.psig{ar,ar2},3),squeeze(nanmean(res_CCA.psig{ar,ar2},1)),[.25 .25],'LineWidth',1,'Color','y')
                    line([0 length(param.time)],[0 length(param.time)],'Color','w')
                    xline(find(param.time==-12)+.5,'Color','w');yline(find(param.time==-12)+.5,'Color','w');

                    %- put the rectangle around the periods of interest
                    for pp = 1 : length(periods)
                        bin2take = param.bin_name == periods_bin{pp} & ( param.time >= periods_time{pp}(1)  & param.time <= periods_time{pp}(2) );
                        rectangle('Position',[find(bin2take==true,1,'first') find(bin2take==true,1,'first') ...
                            sum(bin2take) sum(bin2take)],'EdgeColor',colors(pp,:),'LineWidth',1.5);

                    end
                    title([area2test{ar}(2:end) ' - ' area2test{ar2}(2:end)])
                end
            end
        end

        %% FIGURE 5B - right panels:  average correlation for the different
        figure;
        nb_sub=[6 5];
        bins_lag = -10 : 1 : 10;
        time_lag = bins_lag*(param.time(2)-param.time(1));
        x=0;
        max_crosstemp{1}=NaN(length(area2test),length(area2test));
        max_crosstemp{2}=NaN(length(area2test),length(area2test));
        max_crosstemp{3}=NaN(length(area2test),length(area2test));
        lag_crosstemp = max_crosstemp;
        for ar = 1 : length(area2test)
            for ar2 = 1 : length(area2test)
                if  ar>ar2 && size(res_CCA.cvr{ar,ar2},1)>=minReplicates
                    perf = squeeze(nanmean(res_CCA.cvr{ar,ar2},1));
                    x=x+1;
                    subplot(nb_sub(1),nb_sub(2),x);
                    for pp = 1 : length(periods)
                        bin2take = param.bin_name == periods_bin{pp} & ( param.time >= periods_time{pp}(1)  & param.time <= periods_time{pp}(2) );
                        perf1 = perf(bin2take, bin2take);
                        y=0;xx=[];xx_sem=[];
                        for d = bins_lag
                            y=y+1;
                            xx(y) = nanmean(diag(perf1,d));
                            xx_sem(y) = nanstd(diag(perf1,d))/sqrt(length(diag(perf1,d)));
                        end
                        plot(time_lag,xx,'Color',colors(pp,:));hold on;
                        plot(time_lag,xx+xx_sem,'Color',colors(pp,:));hold on;
                        plot(time_lag,xx-xx_sem,'Color',colors(pp,:));hold on;

                        max_crosstemp{pp}(ar,ar2)=max(xx); % at best time lag
                        [~,ii] = max(xx);
                        lag_crosstemp{pp}(ar,ar2)=time_lag(ii);
                    end

                    title([area2test{ar}(2:end) ' -> ' area2test{ar2}(2:end) '  |  ' area2test{ar2}(2:end) ' -> ' area2test{ar}(2:end)])
                    hold on
                    %line([0 0],[.025 .2],'Color','k')
                    xline(0)
                    ylim([0 .2])
                end
            end
        end

        %% FIGURE 6B - average lag per area and period
        x=0;
        lag_crosstemp_all{1}=NaN(100,length(area2test),length(area2test));
        lag_crosstemp_all{2}=NaN(100,length(area2test),length(area2test));
        lag_crosstemp_all{3}=NaN(100,length(area2test),length(area2test));
        max_crosstemp_all = lag_crosstemp_all;
        for ar = 1 : length(area2test)
            for ar2 = 1 : length(area2test)
                if ar>ar2 && size(res_CCA.cvr{ar,ar2},1)>=minReplicates
                    perf = res_CCA.cvr{ar,ar2};
                    x=x+1;
                    for r = 1 : size(perf,1)

                        for pp = 1 : length(periods)
                            bin2take = param.bin_name == periods_bin{pp} & ( param.time >= periods_time{pp}(1)  & param.time <= periods_time{pp}(2) );
                            perf1 = squeeze(perf(r,bin2take, bin2take));
                            y=0;xx=[];xx_sem=[];
                            for d = bins_lag
                                y=y+1;
                                xx(y) = nanmean(diag(perf1,d));
                                xx_sem(y) = nanstd(diag(perf1,d))/sqrt(length(diag(perf1,d)));
                            end

                            [~,ii] = max(xx);
                            max_crosstemp_all{pp}(r,ar,ar2)=max(xx);
                            lag_crosstemp_all{pp}(r,ar,ar2)=time_lag(ii);
                        end
                    end
                end
            end
        end

        figure;
        xgrid=linspace(-200,200,100);
        for p = 1 : length(periods)
            temp = lag_crosstemp_all{p};
            all_lag=[];all_val=[];
            subplot(1,3,p)
            for ar = 1 : length(area2test)
                temp2 = [squeeze(temp(:,ar,:)) ,-squeeze(temp(:,:,ar))];
                all_lag(ar,:) = temp2(:); %- I have to inverse cos area1-area2 vs area2-area1.. but weird that all points don't have a friend, which they should
                [f,ep]=ksdensity(all_lag(ar,:),xgrid); % remove the outputs to see a 3D plot of the distribution
                % plot(ep,f,'Color',colorsArea(ar,:));hold on
            end

            %- plot median lag and centrality
            temp_lag = prctile(all_lag',[50 40 60]);
            for ar = 1 : length(area2test)
                plot(temp_lag(1,ar),ar,'.','MarkerSize',40,'Color',colorsArea(ar,:));hold on
                line([temp_lag(2,ar) temp_lag(3,ar)],[ar ar],'Color',colorsArea(ar,:));hold on
            end
            xlim([-125 125]);ylim([0 length(area2test)+1]);axis ij
        end

        %% FIGURE 5B - extract FF/FB ratio
        FF_FB_ratio = cell(length(area2test),length(area2test));
        FF_FB_ratio_perm = cell(length(area2test),length(area2test));
        figure; x = 0;
        sig_diff = [];
        for ar = 1 : length(area2test)
            for ar2 = 1 : length(area2test)
                if ar>ar2 && size(res_CCA.cvr{ar,ar2},1)>=minReplicates
                    x=x+1;
                    subplot(nb_sub(1),nb_sub(2),x);
                    for p = 1 : size(res_CCA.cvr{ar,ar2},1)
                        perf = squeeze(res_CCA.cvr{ar,ar2}(p,:,:));
                        perm = squeeze(res_CCA.cvr_perm{ar,ar2}(p,:,:));

                        for pp = 1 : length(periods)
                            bin2take = param.bin_name == periods_bin{pp} & ( param.time >= periods_time{pp}(1)  & param.time <= periods_time{pp}(2) );

                            perf1 = perf(bin2take , bin2take);
                            perm1 = perm(bin2take , bin2take);
                            y=0;xx=[];xx_perm = [];
                            for d = bins_lag
                                y=y+1;
                                xx(y) = nanmean(diag(perf1,d));
                                xx_perm(y) = nanmean(diag(perm1,d));
                            end
                            minmax_xx = (xx-min(xx))/(max(xx)-min(xx)); %- min max normalization to account for potential negative correlations...
                            minmax_xx_perm = (xx_perm-min(xx_perm))/(max(xx_perm)-min(xx_perm)); %- min max normalization to account for potential negative correlations...
                            FF_FB_ratio{ar,ar2}(p,pp)=(sum(minmax_xx(time_lag>0)) - sum(minmax_xx(time_lag<0))) / sum(minmax_xx(time_lag~=0));
                            FF_FB_ratio_perm{ar,ar2}(p,pp)=(sum(minmax_xx_perm(time_lag>0)) - sum(minmax_xx_perm(time_lag<0))) / sum(minmax_xx_perm(time_lag~=0));
                        end
                    end
                    for pp = 1 : length(periods)
                        wdth = .5;
                        boxplot_ind(FF_FB_ratio{ar,ar2}(:,pp),pp,wdth,[colors_light(pp,:) ; colors(pp,:)])
                        plot(nanmean(FF_FB_ratio_perm{ar,ar2}(:,pp)),pp+.45,'x','Color',colors(pp,:)) ; hold on
                        [pval,~,sta]=signrank(FF_FB_ratio{ar,ar2}(:,pp),FF_FB_ratio_perm{ar,ar2}(:,pp),'method','approximate');
                        if pval<0.001 ;                  text(.9,pp,'***','FontSize',20)
                        elseif pval>=0.001 & pval<0.01 ; text(.9,pp,'**','FontSize',20)
                        elseif pval>=0.01 & pval<0.05 ;  text(.9,pp,'*','FontSize',20)
                        elseif pval>=0.05 & pval<0.1 ;   text(.9,pp,'.','FontSize',20)
                        end
                        if pval<0.1
                            sig_diff = [sig_diff ; ar ar2 pp pval sta.zval];
                        end
                    end
                    xlim([-1 1]);ylim([0 length(periods)+1])
                    set(gca,"yTick",1:length(periods),'YTickLabel',periods)
                    xline(0);box on;axis ij
                    title([area2test{ar}(2:end) ' -> ' area2test{ar2}(2:end) '  |  ' area2test{ar2}(2:end) ' -> ' area2test{ar}(2:end)])
                end
            end
        end

        %% FIGURE 6A - Brain plot based on max correlation at any lag (not like Figure 2 which is on average at lag 0)
        figure;
        maximaxi = ceil(max(max([max_crosstemp{:}]))*100)/100;
        for p = 1 : length(periods)
            temp = max_crosstemp{p};
            temp(isnan(temp))=0;
            temp = temp + temp';
            gr = graph(temp,area2test,'omitselfloops');
            fig = subplot(1, length(periods),p);nobrain_map(gr,area2test,[0 maximaxi],'Blues',fig,sig_diff(sig_diff(:,3)==p,[1 2 4 5]))
            title([periods{p}])
        end

    end

    %% FIGURE 2C - Time course of correlations
    cca_all = cell(length(area2test),length(area2test));
    for ar = 1 : length(area2test)
        for ar2 = 1 : length(area2test)
            cca_all{ar,ar2}=[res_CCA.cvr{ar,ar2} ; res_CCA.cvr{ar2,ar}];
        end
    end

    if strcmp(ana2run,'alltime')
        keeptime = {[-500 750] [-500 1000] [0 1000] [-500 500] [0 1000]};
    else
        keeptime = {[-500 1000] [-500 1000]} ;
    end

    figure;
    x=0;
    nb_bins = unique(param.bin_name);
    for ar = 1 : length(area2test)
        x = x + 1 ;

        subplot(2,length(area2test)/2,x)
        for ar2 = 1 : length(area2test)
            if size(cca_all{ar,ar2},1)>=minReplicates
                temp = diag(squeeze(nanmean(cca_all{ar,ar2})));

                p = find(strcmp(periods,'BL'));
                bin2take = param.bin_name == periods_bin{p} & ( param.time >= periods_time{p}(1)  & param.time <= periods_time{p}(2) );
                temp_norm = (temp - mean(temp(bin2take)))/std(temp(bin2take));

                y_smo_combined = [];x_sub_combined=[];
                for ti = 1 : length(nb_bins)
                    bin2show = param.time >= keeptime{ti}(1)  & param.time <= keeptime{ti}(2);
                    y_smoothed1 = smooth(temp_norm(param.bin_name==nb_bins(ti) & bin2show),5,'moving');
                    x_sub1 = param.time(param.bin_name==nb_bins(ti) & bin2show);
                    y_smoothed2 = smooth(temp_norm(param.bin_name==nb_bins(2) & bin2show),5,'moving');
                    x_sub2 = param.time(param.bin_name==nb_bins(2) & bin2show);
    
                    y_smo_combined = [y_smo_combined y_smoothed1' NaN(1,5)];
                    x_sub_combined = [x_sub_combined x_sub1 NaN(1,5)];

                end
                plot(y_smo_combined,'Color',colorsArea(ar2,:),'LineWidth',2);hold on;

            end
        end
        set(gca,"XTick",find(x_sub_combined==13),'XTickLabel',periods(2:end),'Color','none') %- wrong labels for alltime..
        title(area2test{ar}(2:end))
        ylim([-5 20])
        timebar = find(cumsum(diff(x_sub_combined))==500); %- 0.5s length
        line([length(x_sub_combined)-timebar length(x_sub_combined)-1],[-4 -4],'LineWidth',2,'Color','k')
        box off
    end

    %% FIGURE 2B - connectivity fingerprints (with a terrible implementation of fill, but the only one that truly work as intended!)

    %- assign order of areas based on loose anat!
    if length(area2test)==8
        neworder = [4 2 1 5 7 6 8 3];
    else
        neworder = 1 : length(area2test);
    end
    area2test_new = area2test(neworder);
    for ar = 1 : length(area2test_new)
        area2test_new{ar} = area2test_new{ar}(2:end);
    end

    figure;
    for p = 1 : length(periods)
        temp = CCA_avg.(periods{p});
        temp(isnan(temp))=0;
        temp = temp + temp';
        temp(temp==0)=NaN;

        temp2 = CCA_sem.(periods{p});
        temp2(isnan(temp2))=0;
        temp2 = temp2 + temp2';
        temp2(temp2==0)=NaN;

        temp = temp(neworder,neworder);
        temp2 = temp2(neworder,neworder);

        for ar = 1 : length(area2test_new)
            temp2plot = temp(ar,~isnan(temp(ar,:)));
            temp2plot_sem = temp2(ar,~isnan(temp(ar,:)));
            area2plot = area2test_new(~isnan(temp(ar,:)));
            polar_ax = 0:360/(length(area2test_new)):360 ;

            plot_min = [temp2plot temp2plot(1)]-[temp2plot_sem temp2plot_sem(1)];
            plot_max = [temp2plot temp2plot(1)]+[temp2plot_sem temp2plot_sem(1)];
            plot_all =[];
            for i = 1 : length(plot_min)
                plot_all(i,:) = linspace(plot_min(i),plot_max(i),50);
            end
            totake = ~isnan(temp(ar,:));
            polar_ax_sub = [deg2rad(polar_ax(totake)) deg2rad(polar_ax(find(totake==true,1,'first')))];

            subplot(2,length(area2test)/2,ar,'align');
            for i = 1 : size(plot_all,2)
                polarplot(polar_ax_sub, plot_all(:,i)','-','Color',[colors_light(p,:) 0.2]);hold on
            end
            polarplot(polar_ax_sub,    [temp2plot temp2plot(1)],'.-','Color',colors(p,:),'MarkerSize',20,'LineWidth',1.5);
            rlim([0 .2])
            % set(gca,'ThetaTick',polar_ax,'ThetaTickLabel',area2test_new) %- show all labels
            set(gca,'ThetaTick',rad2deg(polar_ax_sub(1:end-1)),'ThetaTickLabel',area2plot,'FontSize',16)
            hold on
            title(area2test_new{ar})
        end
    end
    set(gcf,'Color',[1 1 1])

    %% FIGURE S2 - individual monkey plots

    mks = {'M' 'X'};
    p=2;
    for m = 1 : length(mks)

        %- bin 2 takes depending on event alignement and time windows
        bin2take = param.bin_name == periods_bin{p} & ( param.time >= periods_time{p}(1)  & param.time <= periods_time{p}(2) );

        %- extract average cca
        cca_map=[];
        cca_sem=[];

        for ar = 1 : length(area2test)
            for ar2 = 1 : length(area2test)
                if ar>ar2 && sum(ismember(res_CCA.mk{ar,ar2},mks{m}))>=minReplicates
                    takeme = ismember(res_CCA.mk{ar,ar2},mks{m});
                    takeme = find(takeme==1);
                    temp = [];
                    for j = 1 : size(res_CCA.cvr{ar,ar2}(takeme,:,:),1)
                        temp(j,:) = diag(squeeze(res_CCA.cvr{ar,ar2}(takeme(j),:,:)));
                    end
                    temp = nanmean(temp(:,bin2take),2);
                    % temp = nanmean(squeeze(res_CCA.cvr_perm{ar,ar2}(:,1,:)));
                    cca_map(ar,ar2)=nanmean(temp);
                    cca_sem(ar,ar2)=nanstd(temp)/sqrt(length(temp));
                else
                    cca_map(ar,ar2)=NaN;
                    cca_sem(ar,ar2)=NaN;
                end

            end
        end
        CCA_avg.(mks{m})=cca_map; %- create a variable with the map
        CCA_sem.(mks{m})=cca_sem; %- create a variable with the map
    end

    colors_mk = [46 111 176 ; 245 130 68]/255;
    colors_mk_light = [164 194 209 ; 253 207 158]/255;

    figure;
    for m = 1 : length(mks)
        temp = CCA_avg.(mks{m});
        temp(isnan(temp))=0;
        temp = temp + temp';
        temp(temp==0)=NaN;

        temp2 = CCA_sem.(mks{m});
        temp2(isnan(temp2))=0;
        temp2 = temp2 + temp2';
        temp2(temp2==0)=NaN;

        temp = temp(neworder,neworder);
        temp2 = temp2(neworder,neworder);

        for ar = 1 : length(area2test_new)
            temp2plot = temp(ar,~isnan(temp(ar,:)));
            temp2plot_sem = temp2(ar,~isnan(temp(ar,:)));

            if ~isempty(temp2plot)
                area2plot = area2test_new(~isnan(temp(ar,:)));
                polar_ax = 0:360/(length(area2test_new)):360 ;

                plot_min = [temp2plot temp2plot(1)]-[temp2plot_sem temp2plot_sem(1)];
                plot_max = [temp2plot temp2plot(1)]+[temp2plot_sem temp2plot_sem(1)];
                plot_all =[];
                for i = 1 : length(plot_min)
                    plot_all(i,:) = linspace(plot_min(i),plot_max(i),50);
                end
                totake = ~isnan(temp(ar,:));
                polar_ax_sub = [deg2rad(polar_ax(totake)) deg2rad(polar_ax(find(totake==true,1,'first')))];

                subplot(2,length(area2test)/2,ar,'align');
                for i = 1 : size(plot_all,2)
                    polarplot(polar_ax_sub, plot_all(:,i)','-','Color',[colors_mk_light(m,:) 0.2]);hold on
                end
                polarplot(polar_ax_sub,    [temp2plot temp2plot(1)],'.-','Color',colors_mk(m,:),'MarkerSize',20,'LineWidth',1.5);
                rlim([0 .3])
                % set(gca,'ThetaTick',polar_ax,'ThetaTickLabel',area2test_new) %- show all labels
                set(gca,'ThetaTick',rad2deg(polar_ax_sub(1:end-1)),'ThetaTickLabel',area2plot,'FontSize',16)
                hold on
                title(area2test_new{ar})
            end
        end
    end
    set(gcf,'Color',[1 1 1])

    %% STAT - comparing between pairs of areas
    data_table = table();
    mks = {'M' 'X'};
    for p = 1 : 3
        for m = 1 : length(mks)

            %- bin 2 takes depending on event alignement and time windows
            bin2take = param.bin_name == periods_bin{p} & ( param.time >= periods_time{p}(1)  & param.time <= periods_time{p}(2) );

            %- extract average cca
            cca_map=[];cca_sem=[];
            for ar = 1 : length(area2test)
                for ar2 = 1 : length(area2test)
                    if ar>ar2 && sum(ismember(res_CCA.mk{ar,ar2},mks{m}))>=minReplicates
                        takeme = ismember(res_CCA.mk{ar,ar2},mks{m});
                        takeme = find(takeme==1);
                        temp = [];
                        for j = 1 : size(res_CCA.cvr{ar,ar2}(takeme,:,:),1)
                            temp(j,:) = diag(squeeze(res_CCA.cvr{ar,ar2}(takeme(j),:,:)));
                        end
                        temp = nanmean(temp(:,bin2take),2);

                        data_table = [data_table ; table(temp,repmat({[area2test{ar} '-' area2test{ar2}]},size(temp)),repmat(mks(m),size(temp)),repmat(periods(p),size(temp)),'VariableNames',{'corr','pair','mk','period'}    )];
                    end
                end
            end
        end
    end

    models_form = {'corr ~ 1 + pair*period  + (1|mk)'};
    lme = fitglme(data_table,models_form{1});

    [pval,wald,pval_adj] = cca_posthoc(lme,periods,area2test,3,true);
    [pval,wald,pval_adj] = cca_posthoc_periods(lme,periods,area2test,3,true);

%% LESION ANALYSES (FIGURES 3/S4)
elseif (strcmp(ana2run,'lesion') && ~crosstemp) || strcmp(ana2run,'rewTr') 
   
    %% Average CCA for each area pairs
    minVal = []; maxVal = [];
    for par = 1 : length(predic)
        for p = 1 : length(periods)

            %- bin 2 takes depending on event alignement and time windows
            bin2take = param.bin_name == periods_bin{p} & ( param.time >= periods_time{p}(1)  & param.time <= periods_time{p}(2) );

            %- extract average cca
            cca_map=[];cca_sem=[];
            for ar = 1 : length(area2test)
                for ar2 = 1 : length(area2test)
                    if ar>ar2 && size(LESION{par}.res_CCA.cvr{ar,ar2},1)>=LESION{par}.minReplicates
                        temp = [];
                        for j = 1 : size(LESION{par}.res_CCA.cvr{ar,ar2},1)
                            temp(j,:) = diag(squeeze(LESION{par}.res_CCA.cvr{ar,ar2}(j,:,:)));
                        end
                        temp = nanmean(temp(:,bin2take),2);

                        cca_map(ar,ar2)=nanmean(temp);
                        cca_sem(ar,ar2)=nanstd(temp)/sqrt(length(temp));
                    else
                        cca_map(ar,ar2)=NaN;
                        cca_sem(ar,ar2)=NaN;
                    end
                end
            end

            CCA_avg{par}.(periods{p})=cca_map; %- create a variable with the map
            CCA_sem{par}.(periods{p})=cca_sem; %- create a variable with the map

            maxVal= ceil(max(max([maxVal ; cca_map(:)]))*100)/100; %- for plotting purposes...
            minVal= floor(min(min([minVal ; cca_map(:)]))*100)/100;
        end

    end

    %- Time course of correlations
    for par = 1 : length(predic)
        cca_all{par} = cell(length(area2test),length(area2test));
        sess_all{par} = cell(length(area2test),length(area2test));
        for ar = 1 : length(area2test)
            for ar2 = 1 : length(area2test)
                cca_all{par}{ar,ar2}=[LESION{par}.res_CCA.cvr{ar,ar2} ; LESION{par}.res_CCA.cvr{ar2,ar}];
                sess_all{par}{ar,ar2}=[LESION{par}.res_CCA.sess{ar,ar2} ; LESION{par}.res_CCA.sess{ar2,ar}];
            end
        end
    end
    
    %% FIGURES 3 and S4 / Connectivity over time and lesion types
    figure;
    subplot(1,5,1:3)
    nb_bins = unique(param.bin_name);
    ar = find(ismember(area2test,area2run));

    for ar2 = 1 : length(area2test)
        if size(cca_all{2}{ar,ar2},1)>=LESION{1}.minReplicates %- check the lesion file cos only subselection of area
            temp = diag(squeeze(nanmean(cca_all{1}{ar,ar2})));
            p = find(strcmp(periods,'BL'));
            bin2take = param.bin_name == periods_bin{p} & ( param.time >= periods_time{p}(1)  & param.time <= periods_time{p}(2) );
            temp_norm = (temp - mean(temp(bin2take)))/std(temp(bin2take));
            bin2show = param.time >= -500  & param.time <= 1000;

            if strcmp(period2run,'Stim') ; xx = 1;
            else ; xx = 2;
            end

            y_smoothed1 = smooth(temp_norm(param.bin_name==nb_bins(xx) & bin2show),5,'moving');
            x_sub1 = param.time(param.bin_name==nb_bins(xx) & bin2show);

            y_smo_combined = y_smoothed1';
            x_sub_combined = x_sub1;
            for par = 2 : length(predic)
                temp_lesion = diag(squeeze(nanmean(cca_all{par}{ar,ar2})));
                temp_lesion_norm = (temp_lesion - mean(temp_lesion(bin2take)))/std(temp_lesion(bin2take));

                y_smoothed2 = smooth(temp_lesion_norm(param.bin_name==nb_bins(xx) & bin2show),5,'moving');
                x_sub2 = param.time(param.bin_name==nb_bins(xx) & bin2show);

                y_smo_combined = [y_smo_combined NaN(1,5) y_smoothed2'];
                x_sub_combined = [x_sub_combined NaN(1,5) x_sub2];
            end

            plot(y_smo_combined,'Color',colorsArea(ar2,:),'LineWidth',2);hold on;
        end
    end
    set(gca,"XTick",find(x_sub_combined==13),'XTickLabel',periods(2:end),'Color','none')
    title(area2test{ar}(2:end))
    ylim([-5 20])
    timebar = find(cumsum(diff(x_sub_combined))==500); %- 0.5s length
    line([length(x_sub_combined)-timebar length(x_sub_combined)-1],[-4 -4],'LineWidth',2,'Color','k')
    box off

    %% Connectivity fingerprints
    if length(area2test)==8
        neworder = [4 2 1 5 7 6 8 3];
    else
        neworder = 1 : length(area2test);
    end
    area2test_new = area2test(neworder);
    for ar = 1 : length(area2test_new)
        area2test_new{ar} = area2test_new{ar}(2:end);
    end
    ar = find(ismember(area2test_new,area2run(2:end)));

    subplot(1,5,4)
    p=2; %- ignore BL
    for m = 1 : length(predic)
        temp = CCA_avg{m}.(periods{p});
        temp(isnan(temp))=0;
        temp = temp + temp';
        temp(temp==0)=NaN;

        temp2 = CCA_sem{m}.(periods{p});
        temp2(isnan(temp2))=0;
        temp2 = temp2 + temp2';
        temp2(temp2==0)=NaN;

        temp = temp(neworder,neworder);
        temp2 = temp2(neworder,neworder);

        temp2plot = temp(ar,~isnan(temp(ar,:)));
        temp2plot_sem = temp2(ar,~isnan(temp(ar,:)));
        area2plot = area2test_new(~isnan(temp(ar,:)));
        polar_ax = 0:360/(length(area2test_new)):360 ;

        plot_min = [temp2plot temp2plot(1)]-[temp2plot_sem temp2plot_sem(1)];
        plot_max = [temp2plot temp2plot(1)]+[temp2plot_sem temp2plot_sem(1)];
        plot_all =[];
        for i = 1 : length(plot_min)
            plot_all(i,:) = linspace(plot_min(i),plot_max(i),50);
        end
        totake = ~isnan(temp(ar,:));
        polar_ax_sub = [deg2rad(polar_ax(totake)) deg2rad(polar_ax(find(totake==true,1,'first')))];

        for i = 1 : size(plot_all,2)
            polarplot(polar_ax_sub, plot_all(:,i)','-','Color',[colors_light(m,:) 0.2]);hold on
        end
        polarplot(polar_ax_sub,    [temp2plot temp2plot(1)],'.-','Color',colors(m,:),'MarkerSize',20,'LineWidth',1.5);
        rlim([0 .15])
        set(gca,'ThetaTick',rad2deg(polar_ax_sub(1:end-1)),'ThetaTickLabel',area2plot,'FontSize',16)
        hold on
    end
    set(gcf,'Color',[1 1 1])
    title(area2test_new{ar})

    %% stats
    % DO NOT WRITE 'chosenjuicerew' or emmeans freaking thinks its related to rew so dummy coding fails!!!!!!!!!!

    if ~bothresid & ~strcmp(ana2run,'rewTr')
        x=0;
        nb_bins = unique(param.bin_name);
        ar = find(ismember(area2test,area2run));
        bin2avg = param.bin_name == periods_bin{2} & ( param.time >= periods_time{2}(1)  & param.time <= periods_time{2}(2) );
        data_tab = table();
        for ar2 = 1 : length(area2test)
            if size(cca_all{1}{ar,ar2},1)>=LESION{1}.minReplicates

                for par = 1 : length(predic)
                    temp = cca_all{par}{ar,ar2};
                    temp_lesion=[];
                    for ii = 1 : size(temp,1)
                        temp_lesion(ii,:) = diag(squeeze(temp(ii,:,:)))';
                    end

                    conn_lesion_norm = (temp_lesion - repmat(mean(temp_lesion(:,bin2take),2),1,size(temp_lesion,2))    ) ./ repmat(   std(temp_lesion(:,bin2take)')',1,size(temp_lesion,2));
                    conn_lesion_norm_avg = nanmean(conn_lesion_norm(:,bin2avg),2);
                    temp_lesion_sess = sess_all{par}{ar,ar2};
                    temp_lesion_mk = cat(1,temp_lesion_sess{1,:});
                    temp_lesion_mk = temp_lesion_mk(:,1);

                    data_tab = [data_tab ; table(conn_lesion_norm_avg,...
                        repmat(predic_stats(par),size(conn_lesion_norm_avg,1),1),...
                        repmat(area2test(ar2),size(conn_lesion_norm_avg,1),1),...
                        temp_lesion_mk,...
                        temp_lesion_sess',...
                        'VariableNames',{'conn' 'lesion' 'area' 'mk' 'sess'})];
                end
            end
        end

        models_form = {'conn ~ 1 + lesion + area  + (1|mk) '};
        lme = fitglme(data_tab,models_form{1});

        emm = emmeans(lme,{'lesion'});

        thr = .05;
        %h = emmip(emm,'area');
        clear pval wald
        pval_all=[];
        for cd1 = 1 : length(predic_stats)
            for cd2 = 1 : length(predic_stats)
                if cd1>cd2
                    L_post = strcmp(emm.table.lesion,predic_stats{cd1})' - strcmp(emm.table.lesion,predic_stats{cd2})';
                    H0_post = contrasts_wald(lme,emm,L_post);
                    pval(cd1,cd2) = H0_post.pVal;
                    wald(cd1,cd2) = H0_post.Wald;
                    pval_all = [pval_all ; H0_post.pVal cd1 cd2];
                end
            end
        end
        [h, crit_p, adj_p]=fdr_bh(pval_all(:,1),.05);
        for i = 1 : length(adj_p)
            pval_adj(pval_all(i,2),pval_all(i,3)) = adj_p(i);
        end

        disp(anova(lme));disp(pval_adj);disp(wald);

        subplot(1,5,5)
        for i = 1 : length(emm.table.Estimated_Marginal_Mean)
            plot(i,emm.table.Estimated_Marginal_Mean(i),'.','MarkerSize',30,'Color',colors(i,:));hold on
            plot([i i],emm.table.CI_95_0pct(i,:),'Color',colors(i,:));hold on
        end
        xlim([0 length(predic_stats)+1])
        set(gca,"XTick",1:length(predic_stats),"XTickLabel",predic_stats)
    end

    if strcmp(ana2run,'rewTr')
        %% STAT - comparing between pairs of areas
        p = find(strcmp(periods,'Rew'));
        bin2take = param.bin_name == periods_bin{p} & ( param.time >= periods_time{p}(1)  & param.time <= periods_time{p}(2) );

        data_table = table();
        mks = {'M' 'X'};
        for m = 1 : length(mks)
            for ar = 1 : length(area2test)
                for ar2 = 1 : length(area2test)
                    if ar>ar2 && sum(ismember(LESION{1}.res_CCA.mk{ar,ar2},mks{m}))>=LESION{1}.minReplicates
                        takeme = ismember(LESION{1}.res_CCA.mk{ar,ar2},mks{m});
                        takeme = find(takeme==1);

                        temp_rew = [];temp_norew=[];
                        for j = 1 : size(LESION{1}.res_CCA.cvr{ar,ar2}(takeme,:,:),1)
                            temp_rew(j,:) = diag(squeeze(LESION{1}.res_CCA.cvr{ar,ar2}(takeme(j),:,:)));
                            temp_norew(j,:) = diag(squeeze(LESION{2}.res_CCA.cvr{ar,ar2}(takeme(j),:,:)));
                        end

                        temp_rew = nanmean(temp_rew(:,bin2take),2);
                        temp_norew = nanmean(temp_norew(:,bin2take),2);

                        data_table = [data_table ; table(temp_rew,repmat({[area2test{ar} '-' area2test{ar2}]},size(temp_rew)),repmat(mks(m),size(temp_rew)),repmat({'rew'},size(temp_rew)),'VariableNames',{'corr','pair','mk','reward'}    )];
                        data_table = [data_table ; table(temp_norew,repmat({[area2test{ar} '-' area2test{ar2}]},size(temp_norew)),repmat(mks(m),size(temp_norew)),repmat({'unrew'},size(temp_norew)),'VariableNames',{'corr','pair','mk','reward'}    )];
                    end
                end
            end
        end

        models_form = {'corr ~ 1 + pair*reward  + (1|mk)'};
        lme = fitglme(data_table,models_form{1});
        anova(lme)
    end

%% LESION ANALYSES - CROSSTEMPORAL (FIGURE S5)
elseif strcmp(ana2run,'lesion') && crosstemp

    ar_res = find(ismember(area2test,area2run));

    %- show average correlation for the different lags
    figure;
    bins_lag = -10 : 1 : 10;
    time_lag = bins_lag*(param.time(2)-param.time(1));
    x=0;
    max_crosstemp{1}=NaN(length(area2test),length(area2test));
    max_crosstemp{2}=NaN(length(area2test),length(area2test));
    max_crosstemp{3}=NaN(length(area2test),length(area2test));
    lag_crosstemp = max_crosstemp;

    pp = 2; %- period
    param = LESION{1}.param;
    for par = 1 : length(predic)

        cvr = LESION{par}.res_CCA.cvr ;
        %- the following is to flip maps where ar1>ar2
        for ar1 = 1 : length(area2test)
            for ar2 = 1 : length(area2test)
                if ar1~=ar2 & isempty(cvr{ar1,ar2})
                    for r = 1 : size(cvr{ar2,ar1})
                        cvr{ar1,ar2}(r,:,:) = squeeze(cvr{ar2,ar1}(r,:,:)).';
                    end
                end
            end
        end

        for ar2 = 1 : length(area2test)
            if ar_res~=ar2 && size(cvr{ar_res,ar2},1)>=LESION{1}.minReplicates
                perf = squeeze(nanmean(cvr{ar_res,ar2},1));
                x=x+1;
                subplot(3,3,ar2);
                pp = 2;
                bin2take = param.bin_name == periods_bin{pp} & ( param.time >= periods_time{pp}(1)  & param.time <= periods_time{pp}(2) );
                perf1 = perf(bin2take, bin2take);
                y=0;xx=[];xx_sem=[];
                for d = bins_lag
                    y=y+1;
                    xx(y) = nanmean(diag(perf1,d));
                    xx_sem(y) = nanstd(diag(perf1,d))/sqrt(length(diag(perf1,d)));
                end
                plot(time_lag,xx,'Color',colors(par,:),'LineWidth',2);hold on;
                plot(time_lag,xx+xx_sem,'Color',colors(par,:));hold on;
                plot(time_lag,xx-xx_sem,'Color',colors(par,:));hold on;
                max_crosstemp{par}(ar_res,ar2)=max(xx);
                [~,ii] = max(xx);
                lag_crosstemp{par}(ar_res,ar2)=time_lag(ii);

                title([area2test{ar_res}(2:end) ' -> ' area2test{ar2}(2:end) '  |  ' area2test{ar2}(2:end) ' -> ' area2test{ar_res}(2:end)])
                hold on
                xline(0)
                ylim([0 .2])
            end
        end
    end

    %% directionality ratios
    FF_FB_ratio = cell(length(area2test),length(area2test));
    figure; x = 0;

    for par = 1 : length(predic)
        cvr = LESION{par}.res_CCA.cvr ;
        sess = LESION{par}.res_CCA.sess ;
        for ar1 = 1 : length(area2test)
            for ar2 = 1 : length(area2test)
                if ar1~=ar2 && isempty(cvr{ar1,ar2})
                    for r = 1 : size(cvr{ar2,ar1})
                        cvr{ar1,ar2}(r,:,:) = squeeze(cvr{ar2,ar1}(r,:,:)).';
                    end
                    sess{ar1,ar2} = sess{ar2,ar1};
                end
            end
        end

        for ar2 = 1 : length(area2test)
            if ar_res~=ar2 && size(cvr{ar_res,ar2},1)>=LESION{1}.minReplicates

                subplot(3,3,ar2);
                for p = 1 : size(cvr{ar_res,ar2},1)
                    perf = squeeze(cvr{ar_res,ar2}(p,:,:));
                    bin2take = param.bin_name == periods_bin{pp} & ( param.time >= periods_time{pp}(1)  & param.time <= periods_time{pp}(2) );

                    perf1 = perf(bin2take , bin2take);
                    y=0;xx=[];xx_perm = [];
                    for d = bins_lag
                        y=y+1;
                        xx(y) = nanmean(diag(perf1,d));
                    end
                    minmax_xx = (xx-min(xx))/(max(xx)-min(xx)); %- min max normalization to account for potential negative correlations...
                    FF_FB_ratio{ar_res,ar2}(p,par)=(sum(minmax_xx(time_lag>0)) - sum(minmax_xx(time_lag<0))) / sum(minmax_xx(time_lag~=0));
                end

                boxplot_ind(FF_FB_ratio{ar_res,ar2}(:,par),par,.5,[colors_light(par,:) ; colors(par,:)])
                if par>1
                    [pval,~,sta]=signrank(FF_FB_ratio{ar_res,ar2}(:,par),FF_FB_ratio{ar_res,ar2}(:,1),'method','approximate');
                    if pval<0.001 ;                  text(.9,par,'***','FontSize',20,'Color',colors(par,:))
                    elseif pval>=0.001 & pval<0.01 ; text(.9,par,'**','FontSize',20,'Color',colors(par,:))
                    elseif pval>=0.01 & pval<0.05 ;  text(.9,par,'*','FontSize',20,'Color',colors(par,:))
                    elseif pval>=0.05 & pval<0.1 ;   text(.9,par,'.','FontSize',20,'Color',colors(par,:))
                    end
                end

                xlim([-1 1]);ylim([0 length(predic)+1])
                set(gca,"yTick",1:length(predic),'YTickLabel',predic)
                xline(0);box on;axis ij
                title([area2test{ar_res}(2:end) ' -> ' area2test{ar2}(2:end) '  |  ' area2test{ar2}(2:end) ' -> ' area2test{ar_res}(2:end)])
            end
        end
    end
end



