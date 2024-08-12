%% Canonical Correlation Analysis on triplets of area for POTT dataset
%-
%- Reproduce the cross-area analyses (Figure 4) from:
%- Stoll & Rudebeck (2024) Decision-making shapes dynamic inter-areal communication within macaque ventral frontal cortex
%- https://doi.org/10.1101/2024.07.05.602229
%-
%- Require spk data from all sessions (files like M000000a_SPKpool_counts_1ms.mat)
%-
%- Main difference with this approach compared to original is that we extract the absolute value of the CCA
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Last updated: 2024.08

clear

if isunix
    path2go = '/home/fred/POTTconn/data/'; %- path where SPKpool files are!
elseif ispc
    path2go = 'R:\POTTconn\data\'; %- path where SPKpool files are!
end

area_list = utils_POTT_areas;
area2test = {'a12r' 'a12m' 'a12o' 'a12l' 'a11ml' 'a13l' 'a13m' 'LAI' };
param.minNeurons = 4;
param.smooth = 25; % average FR across bins of XX ms (not sliding window = no overlap, 1 time points every XX ms)
param.sdf = 25; % gaussian smooth, give the time smoothing here, in ms (no smoothing if 0!)
param.SPKthr = .5;
param.nPerm = 100 ; %- number of permutations. If don't need permute, just put [];
param.rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
param.bins4decoding=[2 6]; %- perform CCA on subset on bins (stim and reward here)
param.task = 'INS';
param.nbFold = 10; %- nb of folds for CCA cross-validation
overwrite = false ; %- true or false, if want to recreate matrices!
minReplicates = 3;

%- extract save name to see if should re-run that part!
if ~isempty(param.nPerm)
    savename = [path2go 'FINAL_CCA_crossArea_' param.task '_' num2str(length(area2test)) 'areas_perm_' num2str(param.sdf) 'sdf_' num2str(param.smooth) 'ms.mat'];
else
    savename = [path2go 'FINAL_CCA_crossArea_' param.task '_' num2str(length(area2test)) 'areas_' num2str(param.sdf) 'sdf_' num2str(param.smooth) 'ms.mat'];
end

if exist(savename)~=2 || overwrite

    %- list files to process
    list = dir([path2go '*a_SPKpool_counts_1ms.mat']);

    all_possible_triplets = nchoosek(1:length(area2test),3);
    for n = 1 : length(all_possible_triplets)
        all_possible_triplets_cell{n} = num2str(all_possible_triplets(n,:));
    end

    %- initialize variables
    res_CCA = struct();
    nb = zeros(size(all_possible_triplets_cell));
    xx = 0;
    for s = 1 : length(list)

        %- load SPK data
        disp(list(s).name)
        SPK = load([path2go list(s).name]);

        param.bins2remove = (param.rmv/SPK.subsp)-1;
        SPK.n_evt = length(SPK.times_evts);
        SPK.nTimes = [0 ; (sum(abs(SPK.times(:,:)),2)/SPK.subsp)+1];
        bins = NaN(1,sum(SPK.nTimes));
        time = NaN(1,sum(SPK.nTimes));

        for b = 1 : length(SPK.times_evts)
            time(sum(SPK.nTimes(1:b))+1:sum(SPK.nTimes(1:b+1))) =  (SPK.times(b,1):SPK.subsp:SPK.times(b,2));
            bins(sum(SPK.nTimes(1:b))+1:sum(SPK.nTimes(1:b+1))) = b;
            bins([sum(SPK.nTimes(1:b))+1:sum(SPK.nTimes(1:b))+1+param.bins2remove   sum(SPK.nTimes(1:b+1))-param.bins2remove:sum(SPK.nTimes(1:b+1))    ]) = NaN; %- remove 100 ms each side (avoid overlaps between bins and smoothing problems)
        end
        clear b
        time(~ismember(bins,param.bins4decoding)) = [];
        param.bin_name = bins;
        param.bin_name(~ismember(bins,param.bins4decoding)) = [];
        param.time = time;

        % take the appropriate data (PAV or INS)
        eval(['SPK.SPK_data = SPK.SPK_' param.task ';'])
        eval(['SPK.Tr_Clust_data = SPK.Tr_Clust_' param.task ';'])

        if ~isempty(SPK.SPK_data)

            SPK_data = SPK.SPK_data;
            %- perform smooth and sdf if asked
            if param.sdf~=0
                SPK_data_smo = [];
                parfor tr = 1 : size(SPK_data,1)
                    sdf = single(SPK_data(tr,:))';   % preallocate for speed

                    Gauss_width = max([11 6*param.sdf+1]); % hm, should be an odd number... e.g. 11
                    kernel      = normpdf(-floor(Gauss_width/2):floor(Gauss_width/2),0,param.sdf);
                    dummy       = conv(sdf,kernel);
                    [maxval,maxpos] = max(kernel);
                    kernel      = [-maxpos+1:length(kernel)-maxpos;kernel]';
                    sdf         = dummy(floor(Gauss_width/2)+1:end-floor(Gauss_width/2)); % mode of Gaussian centered on spike -> noncausal
                    SPK_data_smo(tr,:) = sdf';
                end
                SPK_data = SPK_data_smo;
            end

            %- cut the bins you don't consider and the borders of each bin
            SPK_data(:,~ismember(bins,param.bins4decoding))=[];

            if param.smooth~=0
                time_new = []; bin_new = []; SPK_new = [];
                for bi = 1 : length(param.bins4decoding)
                    time_sub = param.time(param.bin_name == param.bins4decoding(bi));
                    curr_bin = find(param.bin_name == param.bins4decoding(bi));

                    timebins = [];
                    timebins(1,:) = 1  : param.smooth : length(time_sub)-param.smooth;
                    timebins(2,:) = param.smooth  : param.smooth : length(time_sub);

                    SPK_temp = [];
                    for t = 1 : size(timebins,2)
                        SPK_temp(:,t) = mean(SPK_data(:,curr_bin(timebins(1,t)) : curr_bin(timebins(2,t))),2);
                    end

                    SPK_new = [SPK_new , SPK_temp];
                    time_new = [time_new , time_sub(timebins(1,:))+(round(param.smooth/2))];
                    bin_new = [bin_new , param.bins4decoding(bi)*ones(1,size(timebins,2))];

                end

                %- update matrices!
                SPK_data = SPK_new;
                param.bin_name = bin_new;
                param.time = time_new;

            end

            %- only run the cross-temporal when in the same alignment!!
            param.mask = param.bin_name.*param.bin_name';
            param.mask = ismember(param.mask,param.bins4decoding.*param.bins4decoding);

            %- look for triplet of simultaneously recorded populations
            goCCA = false(length(area2test),1);
            keep = {};
            for ar = 1 : length(area2test)
                takeit1 = find(ismember(SPK.neurons_area,area_list.(area2test{ar}) ));

                if length(takeit1)>=param.minNeurons

                    % reshape and check the average FR
                    meanFR1 = [];
                    for n = 1 : length(takeit1)
                        meanFR1(n)=1000*mean(mean(SPK_data(ismember(SPK.Tr_Clust_data(:,2),takeit1(n)),:)));
                    end
                    keep{ar} = meanFR1>param.SPKthr;
                    if sum(keep{ar})>=param.minNeurons
                        goCCA(ar) = true;
                    end
                end
            end

            if sum(goCCA)>=3 %- at least 3 areas
                xx = xx + 1;
                disp(sum(goCCA))
                area2take_all = find(goCCA==1);

                all_triplets = nchoosek(area2take_all,3);
                for trip = 1 : length(all_triplets(:,1))

                    area2take = all_triplets(trip,:);
                    which_trip = find(ismember(all_possible_triplets_cell,repmat({num2str(area2take)}, size(all_possible_triplets_cell)  )));
                    nb(which_trip) = nb(which_trip)+1;

                    for t = 1 : length(param.time)
                        clear cca
                        for ar = 1 : length(area2take)
                            takeit = find(ismember(SPK.neurons_area,area_list.(area2test{area2take(ar)}) ));
                            data = 1000*SPK_data(ismember(SPK.Tr_Clust_data(:,2),takeit),:);
                            tr = SPK.Tr_Clust_data(ismember(SPK.Tr_Clust_data(:,2),takeit),:);
                            cca.(['dataNeurons_' num2str(ar)]) = reshape(data(:,t),length(tr(:,1))/length(takeit),length(takeit));
                            cca.(['dataNeurons_' num2str(ar)]) = cca.(['dataNeurons_' num2str(ar)])(:,keep{area2take(ar)});
                        end

                        %- Compute the sample canonical correlation: CROSS-VALIDATED
                        c = cvpartition(length(cca.dataNeurons_1), 'kFold', param.nbFold);
                        cvr = zeros(param.nbFold, 15);
                        for foldIdx = 1:param.nbFold
                            cvr(foldIdx,:) = CanonCorrFitAndPredict_triplet( ...
                                cca.dataNeurons_1(c.training(foldIdx),:), cca.dataNeurons_2(c.training(foldIdx),:), cca.dataNeurons_3(c.training(foldIdx),:),...
                                cca.dataNeurons_1(c.test(foldIdx),:), cca.dataNeurons_2(c.test(foldIdx),:), cca.dataNeurons_3(c.test(foldIdx),:)   );
                        end

                        %- average across cross-validations
                        cvr = nanmean(abs(cvr));

                        res_CCA(which_trip).rep(nb(which_trip)).cvr(:,t) = single(cvr);

                        %- permutations (if needed)
                        if ~isempty(param.nPerm)

                            parfor p = 1 : param.nPerm
                                new_ord_1 = randperm(size(cca.dataNeurons_1,1));
                                new_ord_2 = randperm(size(cca.dataNeurons_2,1));
                                dataNeurons_1_perm = cca.dataNeurons_1(new_ord_1,:); %- change trial order of area 1 and 2, not 3!
                                dataNeurons_2_perm = cca.dataNeurons_2(new_ord_2,:);
                                dataNeurons_3_perm = cca.dataNeurons_3;

                                %- Compute the sample canonical correlation: CROSS-VALIDATED
                                c = cvpartition(length(dataNeurons_1_perm), 'kFold', param.nbFold);
                                cvr_perm = zeros(param.nbFold, 15);
                                for foldIdx = 1:param.nbFold
                                    cvr_perm(foldIdx,:) = CanonCorrFitAndPredict_triplet( ...
                                        dataNeurons_1_perm(c.training(foldIdx),:), dataNeurons_2_perm(c.training(foldIdx),:), dataNeurons_3_perm(c.training(foldIdx),:),...
                                        dataNeurons_1_perm(c.test(foldIdx),:), dataNeurons_2_perm(c.test(foldIdx),:), dataNeurons_3_perm(c.test(foldIdx),:)   );
                                end

                                %- average across cross-validations
                                cvr_perm = nanmean(abs(cvr_perm));
                                cvr_perm_all(p,:) = cvr_perm;
                            end
                            clear cvr_perm
                            res_CCA(which_trip).rep(nb(which_trip)).cvr_perm(t,:,:) = cvr_perm_all;
                        end

                    end

                    res_CCA(which_trip).rep(nb(which_trip)).name =  list(s).name;
                    res_CCA(which_trip).area2test =  area2test(area2take);
                end
            end
        end
    end

    save(savename)
end

%% Post-Processing

load(savename)

%- define time periods (always put one called BL == used as reference for some analyses)
periods =      { 'BL'       'Stim'    'Rew'};
periods_bin =  { 2           2         6};
periods_time = {[-700 -100] [100 700] [100 700]};

%- colors
colors = cbrewer('qual','Set1',9);
colors = colors([9 3 1],:);
colors_light = cbrewer('qual','Pastel1',9);
colors_light = colors_light([9 3 1],:);

%- COLOR ASSIGNMENT
colorsArea = cbrewer('qual', 'Paired', 12);
if length(area2test)==5
    colorsArea = colorsArea([2 4 5 6 12],:);
else
    colorsArea = colorsArea([1:7 12 8:11],:);
end
colorsArea_sub = colorsArea;

%- post hoc
cca_cross_avg = [];
cca_cross_perm = [];
for trip = 1:length(res_CCA)
    if length(res_CCA(trip).rep)>minReplicates
        for sess = 1 : length(res_CCA(trip).rep)
            for p = 1 : length(periods)
                %- bin 2 takes depending on event alignement and time windows
                bin2take = param.bin_name == periods_bin{p} & ( param.time >= periods_time{p}(1)  & param.time <= periods_time{p}(2) );

                posthoc(trip).cca_cross_avg(sess,p,:) = nanmean(res_CCA(trip).rep(sess).cvr(:,bin2take)');
                posthoc(trip).cca_cross_perm(sess,p,:) = nanmean(squeeze(nanmean(res_CCA(trip).rep(sess).cvr_perm(bin2take,:,:),1)));

            end
        end
    end
end

%- correlation of 2 CC from same area (with 2 different areas)
allcorr=[];
allcorr_p=[];
for trip = 1:length(posthoc)
    if length(res_CCA(trip).rep)>minReplicates
        allcorr = [allcorr squeeze(nanmean(posthoc(trip).cca_cross_avg(:,:,13:15)))];
        allcorr_p = [allcorr_p squeeze(nanmean(posthoc(trip).cca_cross_perm(:,:,13:15)))];
    end
end

%% FIGURE 4C
figure;
bins_hist = [0:0.05:1];
for p = 1 : length(periods)
    [a ctr] = histcounts(allcorr(p,:),bins_hist,'Normalization','probability');
    plot(ctr(2:end)-((ctr(2)-ctr(1))/2) ,a,'Color',colors(p,:));hold on

    [b ctr] = histcounts(allcorr_p(p,:),bins_hist,'Normalization','probability');
    plot(ctr(2:end)-((ctr(2)-ctr(1))/2) ,b,'Color',colors_light(p,:))

    [pp,~,sta]=signrank(allcorr(p,:),allcorr_p(p,:),'tail','right','method','approximate');
    % if pp<0.05
    %  text(0.1,.1+(p/20),num2str(round_pval(pp)),'Rotation',0,'Color',colors(p,:))
    text(0.1,.1+(p/20),[periods{p} ' - Z=' num2str(sta.zval) ', p=' num2str(round_pval(pp))],'Rotation',0,'Color',colors(p,:))
    % end
end
set(gca,'FontSize',16);box on;
xlim([0 1]);
xlabel('Average CC rho (btw)')
ylabel('Probability')

allcorr_map=[];
allcorr_sess=[];
allcorr_perm=[];
for trip = 1:length(posthoc)
    if length(res_CCA(trip).rep)>minReplicates
        for a = 1 : 3 %- for each node of a triplet
            curr_corr = squeeze(posthoc(trip).cca_cross_avg(:,:,12+a));
            curr_corr_perm = squeeze(posthoc(trip).cca_cross_perm(:,:,12+a));
            curr_mk = cat(1,res_CCA(trip).rep(:).name); curr_mk = curr_mk(:,1); curr_mk = ismember(curr_mk,'X');
            curr_sess = cat(1,res_CCA(trip).rep(:).name); curr_sess = curr_sess(:,1:8);
            curr_ar = repmat(find(ismember(area2test,res_CCA(trip).area2test(a))),length(curr_mk),1);
            curr_period = sortrows(repmat(1:length(periods),1,length(curr_mk))');
            allcorr_map = [allcorr_map ; curr_corr(:) curr_period repmat(curr_ar,length(periods),1) repmat(curr_mk,length(periods),1)];
            allcorr_sess = [allcorr_sess ; repmat(curr_sess,length(periods),1) ];
            allcorr_perm = [allcorr_perm ; curr_corr_perm(:) curr_period repmat(curr_ar,length(periods),1) repmat(curr_mk,length(periods),1)];
        end
    end
end

%- map it for all areas
for p = 1 : length(periods)
    for ar = 1 : length(area2test)
        avgCorr_mk1(p,ar) = nanmedian(allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar & allcorr_map(:,4)==0 ,1));
        avgCorr_mk2(p,ar) = nanmedian(allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar & allcorr_map(:,4)==1 ,1));
        avgCorr_both(p,ar) = nanmedian(allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar ,1));
        avgCorr_both_std(p,ar) = nanstd(allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar ,1));
        avgCorr_both_perm(p,ar) = nanmedian(allcorr_perm(allcorr_perm(:,2)==p & allcorr_perm(:,3)==ar ,1));
        avgCorr_both_perm_std(p,ar) = nanstd(allcorr_perm(allcorr_perm(:,2)==p & allcorr_perm(:,3)==ar ,1));
    end
end

%% FIGURE 4D

figure;
for p=1:3
    dumm_sess = allcorr_sess(allcorr_map(:,2)==p ,:);
    dumm = allcorr_map(allcorr_map(:,2)==p  ,:);
    modeldata_su = table(double(dumm(:,1)),area2test(dumm(:,3))',categorical(dumm(:,4)),cellstr(dumm_sess), 'VariableNames',{'corr' 'area' 'mk' 'sess'});
    models_form = {'corr ~ 1 + area  + (1|mk) + (1|sess)' ; 'corr ~ 1 + area  + (1|mk)'};
    % [lme,model_final] = model_comparison(modeldata_su,models_form,true);
    lme = fitglme(modeldata_su,models_form{1}); model_final = models_form{1};

    [~,wald,~,pval_adj]  = area_posthoc(lme,area2test,'n');
    disp(anova(lme));disp(pval_adj);disp(wald);
    thr_corr=0.05;

    %-
    subplot(1,3,p);
    xline(0);hold on
    wdth = .6;
    for ar = 1 : length(area2test)
        dumm = allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar ,1);
        dumm1 = allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar & allcorr_map(:,4)==0 ,1);
        dumm2 = allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar & allcorr_map(:,4)==1 ,1);
        boxplot_ind(dumm,ar,wdth,[colorsArea(ar,:) ; colorsArea(ar,:)],false); hold on
        plot(nanmedian(dumm1),ar,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
        plot(nanmedian(dumm2),ar,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
        plot(avgCorr_both_perm(p,:),1:length(area2test),'.','Color','k','MarkerSize',10)
    end
    set(gca,'view',[90 -90],'FontSize',16);box on
    ylim([0 length(area2test)+1])
    xlim([0 1])

    %- plot the significance
    pval2=[pval_adj , zeros(length(area2test),1)]+ [pval_adj' ; zeros(1,length(area2test))];
    for ar1  = 1 : length( area2test)
        Xmax = max(allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar1 ,1));
        Xmin = min(allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar1 ,1));
        Xmean = mean(allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar1 ,1));
        updown=[1.5 2];
        for ar2 = 1 : length(area2test)
            Xmean2 = mean(allcorr_map(allcorr_map(:,2)==p & allcorr_map(:,3)==ar2 ,1));
            if ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean<Xmean2
                text(Xmax+(updown(1)/120),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
                updown(1) = updown(1) + 1;
            elseif ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean>Xmean2
                text(Xmin-(updown(2)/120),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
                updown(2) = updown(2) + 1;
            end
        end
    end
end

%% FIGURE 4B

%- cross-correlation of CCs
allcorr_true=[];
allcorr_perm=[];
for trip = 1:length(posthoc)
    if length(res_CCA(trip).rep)>3
        r_ind_avg = [nanmean(squeeze(nanmean(posthoc(trip).cca_cross_avg(:,:,[4 7 10]))),2) , ...
            nanmean(squeeze(nanmean(posthoc(trip).cca_cross_avg(:,:,[5 8 11]))),2) , ...
            nanmean(squeeze(nanmean(posthoc(trip).cca_cross_avg(:,:,[6 9 12]))),2)];
        allcorr_true = [allcorr_true squeeze(nanmean(posthoc(trip).cca_cross_avg(:,:,1:3)))-r_ind_avg];
        r_ind_avg_perm = [nanmean(squeeze(nanmean(posthoc(trip).cca_cross_perm(:,:,[4 7 10]))),2) , ...
            nanmean(squeeze(nanmean(posthoc(trip).cca_cross_perm(:,:,[5 8 11]))),2) , ...
            nanmean(squeeze(nanmean(posthoc(trip).cca_cross_perm(:,:,[6 9 12]))),2)];
        allcorr_perm = [allcorr_perm squeeze(nanmean(posthoc(trip).cca_cross_perm(:,:,1:3)))-r_ind_avg_perm];
    end
end

figure;
xline(0);hold on
bins_hist = [-.1:0.01:.1];
for p = 1 : length(periods)
    [a ctr] = histcounts(allcorr_true(p,:),bins_hist,'Normalization','probability');
    plot(ctr(2:end)-((ctr(2)-ctr(1))/2) ,a,'Color',colors(p,:))

    [b ctr] = histcounts(allcorr_perm(p,:),bins_hist,'Normalization','probability');
    plot(ctr(2:end)-((ctr(2)-ctr(1))/2) ,b,'Color',colors_light(p,:))

    [pp,~,sta]=signrank(allcorr_true(p,:),allcorr_perm(p,:),'tail','right','method','approximate');
    % if pp<0.05
    text(0.07,.1+(p/20),[periods{p} ' - Z=' num2str(sta.zval) ', p=' num2str(round_pval(pp))],'Rotation',0,'Color',colors(p,:))
    % end
end
set(gca,'FontSize',16);box on;
xlim([-.1 .1]);
xlabel('Average difference in CC rho (dir - ind)')
ylabel('Probability')

%% plot example

expl = find(ismember(all_possible_triplets_cell,{'3  6  8'}));
%expl = find(ismember(all_possible_triplets_cell,{'3  4  5'}));

wdth = .3;

%- example for r_dir/r_ind
y_ax_lim = [0 .5];
figure;
subplot(3,3,[2 5]);
for p = 1 : length(periods)
    ypos = squeeze(posthoc(expl).cca_cross_avg(:,p,1))';
    xpos = boxplot_ind(ypos,p-.2,wdth,[colors_light(p,:) ; colors(p,:)],true);
    ypos2 = nanmean(squeeze(posthoc(expl).cca_cross_avg(:,p,[4 7 10])),2)';
    xpos2 = boxplot_ind(ypos2,p+.2,wdth,[colors_light(p,:) ; colors(p,:)],true);
    for n = 1 : length(xpos)
        line([ypos(n) ypos2(n)],[xpos(n) xpos2(n)],'Color',colors_light(p,:))
    end
    [pp1,~,sta]=signrank(ypos',ypos2','tail','right','method','approximate');
    text(0.01,p+.1,num2str(round_pval(pp1)),'Rotation',90)
end
set(gca,'view',[90 -90],'FontSize',16);box on
ylim([0 length(periods)+1]);xlim(y_ax_lim);set(gca,"YTick",1:length(periods),"YTickLabel",periods);

subplot(3,3,[4 7]);
for p = 1 : length(periods)
    ypos = squeeze(posthoc(expl).cca_cross_avg(:,p,2))';
    xpos = boxplot_ind(ypos,p-.2,wdth,[colors_light(p,:) ; colors(p,:)],true);
    ypos2 = nanmean(squeeze(posthoc(expl).cca_cross_avg(:,p,[5 8 11])),2)';
    xpos2 = boxplot_ind(ypos2,p+.2,wdth,[colors_light(p,:) ; colors(p,:)],true);
    for n = 1 : length(xpos)
        line([ypos(n) ypos2(n)],[xpos(n) xpos2(n)],'Color',colors_light(p,:))
    end
    [pp1,~,sta]=signrank(ypos',ypos2','tail','right','method','approximate');
    text(0.01,p+.1,num2str(round_pval(pp1)),'Rotation',90)
end
set(gca,'view',[90 -90],'FontSize',16);box on
ylim([0 length(periods)+1]);xlim(y_ax_lim);set(gca,"YTick",1:length(periods),"YTickLabel",periods);

subplot(3,3,[6 9]);
for p = 1 : length(periods)
    ypos = squeeze(posthoc(expl).cca_cross_avg(:,p,3))';
    xpos = boxplot_ind(ypos,p-.2,wdth,[colors_light(p,:) ; colors(p,:)],true);
    ypos2 = nanmean(squeeze(posthoc(expl).cca_cross_avg(:,p,[6 9 12])),2)';
    xpos2 = boxplot_ind(ypos2,p+.2,wdth,[colors_light(p,:) ; colors(p,:)],true);
    for n = 1 : length(xpos)
        line([ypos(n) ypos2(n)],[xpos(n) xpos2(n)],'Color',colors_light(p,:))
    end
    [pp1,~,sta]=signrank(ypos',ypos2','tail','right','method','approximate');
    text(0.01,p+.1,num2str(round_pval(pp1)),'Rotation',90)
end
set(gca,'view',[90 -90],'FontSize',16);box on
ylim([0 length(periods)+1]);xlim(y_ax_lim);set(gca,"YTick",1:length(periods),"YTickLabel",periods);

subplot(3,3,1);title(res_CCA(expl).area2test{1});axis off
subplot(3,3,3);title(res_CCA(expl).area2test{2});axis off
subplot(3,3,8);title(res_CCA(expl).area2test{3});axis off

%- example for r_btw
y_ax_lim = [0 1]
% rand_x = ((rand(size(cca_cross_avg,1),1)-.5)/3)';
figure;
subplot(1,3,1);
for p = 1 : length(periods)
    ypos = squeeze(posthoc(expl).cca_cross_avg(:,p,13))';
    xpos = boxplot_ind(ypos,p-.2,wdth,[colors_light(p,:) ; colors(p,:)]);
    ypos2 = squeeze(posthoc(expl).cca_cross_perm(:,p,13))';
    xpos2 = boxplot_ind(ypos2,p+.2,wdth,[colors_light(p,:) ; colors(p,:)]);
    % for n = 1 : length(xpos)
    %     line([ypos(n) ypos2(n)],[xpos(n) xpos2(n)],'Color',colors_light(p,:))
    % end
    [pp,~,sta]=signrank(ypos',ypos2','tail','right','method','approximate');
    text(.1,p,num2str(round_pval(pp)),'Rotation',90)
end
set(gca,'view',[90 -90],'FontSize',16);box on
ylim([0 length(periods)+1]);xlim(y_ax_lim);set(gca,"YTick",1:length(periods),"YTickLabel",periods);axis
title(res_CCA(expl).area2test{1})
subplot(1,3,2);
for p = 1 : length(periods)
    ypos = squeeze(posthoc(expl).cca_cross_avg(:,p,14))';
    xpos = boxplot_ind(ypos,p-.2,wdth,[colors_light(p,:) ; colors(p,:)]);
    ypos2 = squeeze(posthoc(expl).cca_cross_perm(:,p,14))';
    xpos2 = boxplot_ind(ypos2,p+.2,wdth,[colors_light(p,:) ; colors(p,:)]);
    % for n = 1 : length(xpos)
    %     line([ypos(n) ypos2(n)],[xpos(n) xpos2(n)],'Color',colors_light(p,:))
    % end
    [pp,~,sta]=signrank(ypos',ypos2','tail','right','method','approximate');
    text(.1,p,num2str(round_pval(pp)),'Rotation',90)
end
set(gca,'view',[90 -90],'FontSize',16);box on
ylim([0 length(periods)+1]);xlim(y_ax_lim);set(gca,"YTick",1:length(periods),"YTickLabel",periods);axis
title(res_CCA(expl).area2test{2})
subplot(1,3,3);
for p = 1 : length(periods)
    ypos = squeeze(posthoc(expl).cca_cross_avg(:,p,15))';
    xpos = boxplot_ind(ypos,p-.2,wdth,[colors_light(p,:) ; colors(p,:)]);
    ypos2 = squeeze(posthoc(expl).cca_cross_perm(:,p,15))';
    xpos2 = boxplot_ind(ypos2,p+.2,wdth,[colors_light(p,:) ; colors(p,:)]);
    % for n = 1 : length(xpos)
    %     line([ypos(n) ypos2(n)],[xpos(n) xpos2(n)],'Color',colors_light(p,:))
    % end
    [pp,~,sta]=signrank(ypos',ypos2','tail','right','method','approximate');
    text(.1,p,num2str(round_pval(pp)),'Rotation',90)
end
set(gca,'view',[90 -90],'FontSize',16);box on
ylim([0 length(periods)+1]);xlim(y_ax_lim);set(gca,"YTick",1:length(periods),"YTickLabel",periods);axis
title(res_CCA(expl).area2test{3})




%% FOR REVISION !!! SHOULD BE ABLE TO REMOVE SOON!!!!!!!

allcorr_map=[];
allcorr_sess=[];
allcorr_perm=[];
pair = [1 2 ; 1 3 ; 2 3];

for trip = 1:length(posthoc)
    if length(res_CCA(trip).rep)>minReplicates
        for a = 1 : 3 %- for each node of a triplet
            curr_corr = squeeze(posthoc(trip).cca_cross_avg(:,:,12+a));
            curr_corr_perm = squeeze(posthoc(trip).cca_cross_perm(:,:,12+a));

            curr_corr_ind = nanmean(squeeze(posthoc(trip).cca_cross_avg(:,:,[3+a 6+a 9+a])),3);%-squeeze(posthoc(trip).cca_cross_avg(:,:,3+a));
            curr_corr_dir = squeeze(posthoc(trip).cca_cross_avg(:,:,a)) ;
            curr_corr_ratio = curr_corr_dir./curr_corr_ind;

            curr_mk = cat(1,res_CCA(trip).rep(:).name); curr_mk = curr_mk(:,1); curr_mk = ismember(curr_mk,'X');
            curr_sess = cat(1,res_CCA(trip).rep(:).name); curr_sess = curr_sess(:,1:8);
            curr_ar = repmat(find(ismember(area2test,res_CCA(trip).area2test(pair(a,:)))),length(curr_mk),1);
            curr_alt_ar = repmat(find(ismember(area2test,res_CCA(trip).area2test(find(~ismember([1 2 3],pair(a,:)))))),length(curr_mk),1);
            curr_period = sortrows(repmat(1:length(periods),1,length(curr_mk))');
            allcorr_map = [allcorr_map ; curr_corr_dir(:) curr_corr_ratio(:) curr_period repmat(curr_ar,length(periods),1) repmat(curr_alt_ar,length(periods),1) repmat(curr_mk,length(periods),1)];
            allcorr_sess = [allcorr_sess ; repmat(curr_sess,length(periods),1) ];
            allcorr_perm = [allcorr_perm ; curr_corr_perm(:) curr_period repmat(curr_ar,length(periods),1) repmat(curr_alt_ar,length(periods),1) repmat(curr_mk,length(periods),1)];
        end
    end
end

% btw dir-ind period ar1 ar2 arALT mk
figure;x=0;
for ar1 = 1 : length(area2test)
    for ar2 =  1 : length(area2test)
        if ar1<ar2
            for p = 1 : length(periods)
                temp = allcorr_map(allcorr_map(:,3)==p & allcorr_map(:,4)==ar1 & allcorr_map(:,5)==ar2 ,[1 2 6]);
                ratio_dir_ind = [grpstats(temp(:,1:2),temp(:,3),{'mean'}) , grpstats(temp(:,1),temp(:,3),{'numel'}) ,unique(temp(:,3))];
                if ~isempty(ratio_dir_ind)
                    if p==1
                        x = x + 1 ;
                    end
                    nullratio = NaN(8,1);
                    nullratio(ratio_dir_ind(:,4),1) = ratio_dir_ind(:,2);
                    subplot(5,5,x)
                    bar(p:3:length(area2test)*length(periods),nullratio,'FaceColor',colors(p,:),'BarWidth',0.2);hold on
                    for t = 1 : length(ratio_dir_ind(:,1))
                        text((ratio_dir_ind(t,end)*length(periods))-2,.8,num2str(ratio_dir_ind(t,3)))
                    end

                    xlim([0 1+length(area2test)*length(periods)]);ylim([.8 1.5]);% ylim([-0.01 0.1])
                    xloc = 2:3:length(area2test)*length(periods);
                    set(gca,"XTick",xloc(ratio_dir_ind(:,end)),'XTickLabel',area2test(ratio_dir_ind(:,end)))
                    title([area2test{ar1} ' - ' area2test{ar2} ' (dir=' num2str(nanmean(ratio_dir_ind(:,1))) ')'] )
                end
            end
        end
    end
end

% different plot
% btw dir-ind period ar1 ar2 arALT mk
figure;x=0;
for ar1 = 1 : length(area2test)
    for ar2 =  1 : length(area2test)
        if ar1<ar2
            takeme = allcorr_map(:,4)==ar1 & allcorr_map(:,5)==ar2;
            data4stat = allcorr_map(takeme,[2 3 6]);
            mk4stat = allcorr_sess(takeme,1);
            nbtrip = length(unique(data4stat(:,end)));

            % stats
            if nbtrip>=4
                data4stat_tab = table(double(data4stat(:,1)),periods(data4stat(:,2))',area2test(data4stat(:,3))',mk4stat,'VariableNames',{'conn' 'period' 'area' 'mk'});

                model_form = 'conn ~ 1 + area*period  + (1|mk)';
                lme = fitglme(data4stat_tab,model_form);
                disp(anova(lme));
                % [~,wald,~,pval_adj]  = area_posthoc(lme,area2test,'n');
                % disp(pval_adj);disp(wald);
                x = x + 1 ;

                res(x).anov = anova(lme);
                res(x).emm = emmeans(lme,{'area' 'period'});
                [pval,wald,pval_adj] = cca_triplet_posthoc(lme,periods,area2test(unique(data4stat(:,end))),true);

            end



            clear avg_dir
            for p = 1 : length(periods)
                temp = allcorr_map(allcorr_map(:,3)==p & allcorr_map(:,4)==ar1 & allcorr_map(:,5)==ar2 ,[1 2 6]);
                nbtrip = [grpstats(temp(:,1),temp(:,3),{'numel'}) , unique(temp(:,3))];
                wdth = .6;
                if length(nbtrip)>=4
                    ar3=unique(temp(:,3));
                    subplot(3,5,x)
                    for ot = 1 : length(ar3)
                        boxplot_ind(temp(temp(:,3)==ar3(ot),2),(ar3(ot)*length(periods))-2+p,wdth,[colors_light(p,:) ; colors(p,:)],true); hold on
                    end
                    for t = 1 : length(nbtrip(:,1))
                        text(.6,(nbtrip(t,2)*length(periods))-2,num2str(nbtrip(t,1)))
                    end
                    ylim([0 1.5+length(area2test)*length(periods)]); xlim([.5 3]);% ylim([-0.01 0.1])
                    yloc = 3:3:length(area2test)*length(periods);
                    set(gca,"YTick",yloc(ar3),'YTickLabel',area2test(ar3))
                    set(gca,'view',[90 -90]);box on

                    avg_dir(p,1) = mean(temp(:,1));

                    if p == length(periods)
                        title([area2test{ar1} ' - ' area2test{ar2} ' / BL=' num2str(round(avg_dir(1,1),3)) ',Stim=' num2str(round(avg_dir(2,1),3)) ',Rew=' num2str(round(avg_dir(3,1),3))])
                    end

                end
            end
        end
    end
end





