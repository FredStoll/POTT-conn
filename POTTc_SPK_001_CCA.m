%% Canonical Correlation Analysis on POTT dataset - Generate CCA matrices
%-
%- Reproduce the analyses from:
%- Stoll & Rudebeck (2024) Decision-making shapes dynamic inter-areal communication within macaque ventral frontal cortex
%- https://doi.org/10.1101/2024.07.05.602229
%-
%- Require spk data from all sessions (files like M000000a_SPKpool_counts_1ms.mat)
%- Should run a few times with different parameters to perform all possible analyses: 
%- 
%- MAIN ANALYSES: 
%-      FIG 2/5/6/S2/S3: param.crosstemp = true; param.resid = false;
%-      FIG S1: param.crosstemp = false; param.resid = false; param.bins4decoding=[1 2 4 5 6]; %- that one takes many hours!
%-      FIG S4D: param.crosstemp = false; param.resid = false; param.rewardedOrNot = 0;
%-      FIG S4D: param.crosstemp = false; param.resid = false; param.rewardedOrNot = 1;
%- LESION ANALYSES: 
%-      FIG 3B/S4B: param.crosstemp = false; param.resid = true; param.predic='allstim' ; param.area4resid = 3; param.bothresid = false;
%-      FIG S4B/C: param.crosstemp = false; param.resid = true; param.predic='chosenproba' ; param.area4resid = 3; param.bothresid = false;
%-      FIG S4B: param.crosstemp = false; param.resid = true; param.predic='chosenjuice' ; param.area4resid = 3; param.bothresid = false;
%-      FIG S4B/C: param.crosstemp = false; param.resid = true; param.predic='chosenside' ; param.area4resid = 3; param.bothresid = false;
%-      FIG 3C/S4C: param.crosstemp = false; param.resid = true; param.predic='allrew' ; param.area4resid = 3; param.bothresid = false;
%-      FIG S4C: param.crosstemp = false; param.resid = true; param.predic='rew' ; param.area4resid = 3; param.bothresid = false;
%-      FIG S4C: param.crosstemp = false; param.resid = true; param.predic='chosenjuicerew' ; param.area4resid = 3; param.bothresid = false;
%-      FIG S4A: param.crosstemp = false; param.resid = true; param.predic='allstim' ; param.area4resid = 3; param.bothresid = true;
%-      FIG S4A: param.crosstemp = false; param.resid = true; param.predic='allrew' ; param.area4resid = 3; param.bothresid = true;
%-      FIG S5A: param.crosstemp = true; param.resid = true; param.predic='allstim' ; param.area4resid = 5; param.bothresid = false;
%-      FIG S5A: param.crosstemp = true; param.resid = true; param.predic='chosenproba' ; param.area4resid = 5; param.bothresid = false;
%-      FIG S5A: param.crosstemp = true; param.resid = true; param.predic='chosenjuice' ; param.area4resid = 5; param.bothresid = false;
%-      FIG S5A: param.crosstemp = true; param.resid = true; param.predic='chosenside' ; param.area4resid = 5; param.bothresid = false;
%-      FIG S5B: param.crosstemp = true; param.resid = true; param.predic='allrew' ; param.area4resid = 8; param.bothresid = false;
%-      FIG S5B: param.crosstemp = true; param.resid = true; param.predic='rew' ; param.area4resid = 8; param.bothresid = false;
%-      FIG S5B: param.crosstemp = true; param.resid = true; param.predic='chosenjuicerew' ; param.area4resid = 8; param.bothresid = false;
%-      FIG S5B: param.crosstemp = true; param.resid = true; param.predic='chosenproba' ; param.area4resid = 8; param.bothresid = false;
%-      FIG S5B: param.crosstemp = true; param.resid = true; param.predic='chosenside' ; param.area4resid = 8; param.bothresid = false;
%- CROSS-AREA ANALYSES:
%-      FIG 4: run 'POTTc_SPK_003_CCA_triplets' instead
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
area2test = {'a12r' 'a12m' 'a12o' 'a12l' 'a11ml' 'a13l' 'a13m' 'LAI' }; % areas considered
param.minNeurons = 4; % number of minimum neurons per population
param.smooth = 25; % average FR across bins of XX ms (not sliding window = no overlap, 1 time points every XX ms)
param.sdf = 25; % gaussian smooth, give the time smoothing here, in ms (no smoothing if 0!)
param.crosstemp = false; % crosstemporal CCA (correlation at 0 lag + lags if true, only lag 0 otheriwse

param.resid = false; % apply ANOVA to remove influence of selected decision related parameters
param.predic = 'allstim' ; % parameters to remove: allstim / chosenproba / chosenjuice / chosenside / allrew / rew / chosenjuicerew / chosenproba / chosenside
param.area4resid = 5; % which area (position in the area2test matrix) to use the residuals (e.g. area4resid = 3 means a12o)
param.bothresid = false; % whether to also use residuals on the other area (true for Figure S4A or false for the others)

param.rewardedOrNot = 0; %- empty = ignore // 0 = unrewarded trials only // 1 = rewarded trials only

param.SPKthr = .5 ; %- FR threshold (0.5Hz)
param.nPerm = 1 ; %- number of permutations. 1 ids enough when cross-temporal! If don't need permute, just put [];
param.rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
param.bins4decoding=[2 6]; %- [2 6] perform CCA on subset on bins (stim and reward alignments here)
param.task = 'INS'; % instrumental task
param.nbFold = 10; %- nb of folds for CCA cross-validation
overwrite = false ; %- true or false, if want to recreate matrices!
minReplicates = 3; %- for postprocessing, only keep pairs of areas which have X sessions at least

%- extract save name to see if should re-run that part!
if param.resid && ~param.bothresid
    dumm = ['resid' area2test{param.area4resid}(2:end) '_' param.predic] ;
    param.nPerm = [] ; %- remove permutations if residual models
elseif param.resid && param.bothresid
    dumm = ['resid' area2test{param.area4resid}(2:end) '_both_' param.predic] ;
    param.nPerm = [] ; %- remove permutations if residual models
else
    dumm = 'raw';
    param.predic = [];
    param.area4resid = [];
end
if param.rewardedOrNot==0
    dumm = 'raw_unrewarded';
    param.predic = [];
    param.area4resid = [];
elseif param.rewardedOrNot==1
    dumm = 'raw_rewarded';
    param.predic = [];
    param.area4resid = [];
end

if ~isempty(param.nPerm)
    savename = [path2go 'FINAL_CCA_' param.task '_' num2str(length(area2test)) 'areas_perm_' num2str(param.sdf) 'sdf_' num2str(param.smooth) 'ms_' dumm '.mat'];
else
    savename = [path2go 'FINAL_CCA_' param.task '_' num2str(length(area2test)) 'areas_' num2str(param.sdf) 'sdf_' num2str(param.smooth) 'ms_' dumm '.mat'];
end


if exist(savename)~=2 || overwrite %- only if CCA file doesn't exist or if overwrite is asked

    %- start a parallel pool (8 for local / 36 for servers)
    [~,name] = system('hostname');
    if isempty(gcp('nocreate')) && strcmp(name(1:5),'DESKT')
        parpool(8);
    elseif isempty(gcp('nocreate')) && ~strcmp(name(1:5),'DESKT')
        parpool(36);
    end

    %- list files to process
    list = dir([path2go '*a_SPKpool_counts_1ms.mat']);

    %- initialize variables
    matrix2create = {'cvr' 'cvr_perm' 'cvr2' 'cvr2_perm' 'coeffvar' 'coeffvar2' 'cva' 'cvb' 'mk' 'sess'};
    for i = 1 : length(matrix2create)
        res_CCA.(matrix2create{i}) = cell(length(area2test),length(area2test));
    end
    nb = zeros(length(area2test),length(area2test));

    for s = 1 : length(list) %- for every session

        %- load SPK data
        disp(list(s).name)
        SPK = load([path2go list(s).name]);

        %- recreate time vectors from the multiple alignments
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

        % take the appropriate task (PAV or INS)
        eval(['SPK.SPK_data = SPK.SPK_' param.task ';'])
        eval(['SPK.Tr_Clust_data = SPK.Tr_Clust_' param.task ';'])

        if ~isempty(param.area4resid) %- early check that enough neuron in area used for residuals so doesn't spend time for nothing! 
            tryit = (length(find(ismember(SPK.neurons_area,area_list.(area2test{param.area4resid}) ))) >= param.minNeurons );
        else
            tryit = true;
        end

        if ~isempty(SPK.SPK_data) & tryit

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

            %- smooth data if required
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

            %- only run the cross-temporal when in the same event alignment!!
            if param.crosstemp
                param.mask = param.bin_name.*param.bin_name';
                param.mask = ismember(param.mask,param.bins4decoding.*param.bins4decoding);
            else
                param.mask = eye(length(param.bin_name));
            end

            for ar = 1 : length(area2test)
                takeit1 = find(ismember(SPK.neurons_area,area_list.(area2test{ar}) )); %- find neurons from considered area
                for ar2 = 1 : length(area2test)
                    takeit2 = find(ismember(SPK.neurons_area,area_list.(area2test{ar2}) )); %- find neurons from considered area

                    %- only run if one of the area is the one for residual
                    if param.resid
                        ok = ar==param.area4resid | ar2==param.area4resid ;
                    else
                        ok = true;
                    end

                    if length(takeit1)>=param.minNeurons && length(takeit2)>=param.minNeurons && ar>ar2 && ok

                        % reshape and check the average FR
                        meanFR1 = []; meanFR2 = [];
                        for n = 1 : length(takeit1)
                            meanFR1(n)=1000*mean(mean(SPK_data(ismember(SPK.Tr_Clust_data(:,2),takeit1(n)),:)));
                        end
                        for n = 1 : length(takeit2)
                            meanFR2(n)=1000*mean(mean(SPK_data(ismember(SPK.Tr_Clust_data(:,2),takeit2(n)),:)));
                        end
                        keep1 = meanFR1>param.SPKthr;
                        keep2 = meanFR2>param.SPKthr;

                        %- extract data for the considered neurons
                        data1 = SPK_data(ismember(SPK.Tr_Clust_data(:,2),takeit1),:)*1000;
                        data2 = SPK_data(ismember(SPK.Tr_Clust_data(:,2),takeit2),:)*1000;
                        tr1 = SPK.Tr_Clust_data(ismember(SPK.Tr_Clust_data(:,2),takeit1),:);
                        tr2 = SPK.Tr_Clust_data(ismember(SPK.Tr_Clust_data(:,2),takeit2),:);

                        if param.resid | ~isempty(param.rewardedOrNot)
                            %- extract behavioral variables
                            predic_full = {'chosenproba' ;'chosenjuice' ; 'chosenside' ; 'rew'; 'chosenjuicerew'};
                            for pa = 1:length(predic_full)
                                if strcmp(predic_full(pa),'chosenjuicerew')
                                    factors_full(:,pa) = SPK.TrialType_INS(ismember(SPK.TrialType_header,'I_chosenjuice'),:)';
                                elseif strcmp(predic_full(pa),'rew')
                                    factors_full(:,pa) = SPK.TrialType_INS(ismember(SPK.TrialType_header,predic_full(pa)),:)';
                                else
                                    factors_full(:,pa) = SPK.TrialType_INS(ismember(SPK.TrialType_header,{['I_' predic_full{pa}]}),:)';
                                end
                            end

                            %- trials to consider for instrum
                            when = @(arg,cd) SPK.TrialType_INS(strcmp(SPK.TrialType_header,arg),:)==cd ;
                            norm_pb = @(data) -1+((data-0.1)*2)/(0.9-0.1) ;

                            diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
                            same_juice = (when('I_juiceL',1) & when('I_juiceR',1)) | (when('I_juiceL',2) & when('I_juiceR',2)); %- take only the different juice trials

                            %- normalize predictors
                            factors_full(factors_full(:,2)==1,2)=-1; %- norm chosen juice (-1 vs 1)
                            factors_full(factors_full(:,2)==2,2)=1;
                            factors_full(factors_full(:,3)==1,3)=-1; %- norm chosen side (-1 vs 1)
                            factors_full(factors_full(:,3)==2,3)=1;
                            factors_full(factors_full(:,5)==1,5)=-1; %- %- norm chosen juice (-1 vs 1)
                            factors_full(factors_full(:,5)==2,5)=1;

                            factors_full_nested = factors_full;
                            factors_full_nested(factors_full_nested(:,4)==-1,5) = 0;

                            if param.resid
                                if strcmp(param.predic,'chosenproba') | strcmp(param.predic,'allstim') | strcmp(param.predic,'allrew')
                                    tr2consider = find(diff_juice' &  ismember(factors_full(:,1),[10:20:90])); %- limit to certain proba
                                elseif strcmp(param.predic,'chosenjuicerew')
                                    tr2consider = find(diff_juice' & factors_full(:,4)==1);
                                else
                                    tr2consider = find(diff_juice');
                                end
                            else
                                if param.rewardedOrNot==0
                                    tr2consider = find(factors_full(:,4)==-1);
                                elseif param.rewardedOrNot==1
                                    tr2consider = find(factors_full(:,4)==1);
                                end
                            end

                            factors_full(:,1) = norm_pb(factors_full(:,1)/100); %- normalize proba

                        else
                            tr2consider = unique(tr1(:,1)) ; % if no residual activity, take all trials
                            factors_full=[];
                        end

                        if sum(keep1)>=param.minNeurons && sum(keep2)>=param.minNeurons % CHECK THAT ENOUGH NEURONS

                            nb(ar,ar2) = nb(ar,ar2) + 1; %- update where to save the data
                            disp(['    --> Computing ' area2test{ar} ' vs ' area2test{ar2} ])

                            clear dataNeurons1 dataNeurons2
                            for t1 = 1  : length(data1(1,:))

                                dataNeurons1 = reshape(data1(:,t1),length(tr1(:,1))/length(takeit1),length(takeit1));

                                if (param.resid && ar == param.area4resid) || (param.resid && param.bothresid)
                                    dataNeurons1 = residualFR(dataNeurons1(tr2consider,:),param,factors_full(tr2consider,:),predic_full);
                                else
                                    dataNeurons1 = dataNeurons1(tr2consider,:);
                                end

                                % initialize some variables...
                                cvr_all = [];cvr_all_perm =[];
                                cvr2_all = [];cvr2_all_perm =[];
                                coeffvar_all=[];coeffvar2_all=[];
                                cva_all=NaN(sum(keep1),length(data1(1,:)));
                                cvb_all=NaN(sum(keep2),length(data1(1,:)));

                                parfor t2 = 1 : length(data1(1,:))
                                    if param.mask(t1,t2)

                                        dataNeurons2 = reshape(data2(:,t2),length(tr2(:,1))/length(takeit2),length(takeit2));

                                        if (param.resid && ar2 == param.area4resid) || (param.resid && param.bothresid)
                                            dataNeurons2 = residualFR(dataNeurons2(tr2consider,:),param,factors_full(tr2consider,:),predic_full);
                                        else
                                            dataNeurons2 = dataNeurons2(tr2consider,:);
                                        end

                                        %- Compute the sample canonical correlation: CROSS-VALIDATED
                                        c = cvpartition(length(dataNeurons1), 'kFold', param.nbFold);
                                        numPairs = min( size(dataNeurons1(:,keep1), 2), size(dataNeurons2(:,keep2), 2) );
                                        cvr = zeros(param.nbFold, numPairs);
                                        A = zeros(param.nbFold, size(dataNeurons1(:,keep1), 2));
                                        B = zeros(param.nbFold, size(dataNeurons2(:,keep2), 2));

                                        for foldIdx = 1:param.nbFold
                                            [cvr(foldIdx,:),r,A(foldIdx,:),B(foldIdx,:)] = CanonCorrFitAndPredict( ...
                                                dataNeurons1(c.training(foldIdx),keep1), dataNeurons2(c.training(foldIdx),keep2), ...
                                                dataNeurons1(c.test(foldIdx),keep1)    , dataNeurons2(c.test(foldIdx),keep2)   );
                                        end

                                        %- average across cross-validations
                                        coeffvar = nanstd(cvr)./nanmean(cvr);
                                        cvr = nanmean(cvr);
                                        cva = nanmean(abs(A),1);
                                        cvb = nanmean(abs(B),1);

                                        %- compute the full CCA (not cross-validated), to check. Pretty well correlated!
                                        % [~, ~, r] = canoncorr( dataNeurons1(:,keep1), dataNeurons2(:,keep2) );

                                        %- keep values!
                                        cvr_all(1,t2) = single(cvr(1));
                                        cvr2_all(1,t2) = single(cvr(2));
                                        coeffvar_all(1,t2) = single(coeffvar(1));
                                        coeffvar2_all(1,t2) = single(coeffvar(2));
                                        cva_all(:,t2) = single(cva);
                                        cvb_all(:,t2) = single(cvb);

                                        %- perm
                                        if ~isempty(param.nPerm)

                                            new_ord = randperm(size(dataNeurons1,1));
                                            dataNeurons1_perm = dataNeurons1;
                                            dataNeurons2_perm = dataNeurons2(new_ord,:);

                                            %- Compute the sample canonical correlation: CROSS-VALIDATED
                                            numPairs = min( size(dataNeurons1(:,keep1), 2), size(dataNeurons2(:,keep2), 2) );
                                            cvr_perm = zeros(param.nbFold, numPairs);
                                            for foldIdx = 1:param.nbFold
                                                cvr_perm(foldIdx,:) = CanonCorrFitAndPredict( ...
                                                    dataNeurons1_perm(c.training(foldIdx),keep1), dataNeurons2_perm(c.training(foldIdx),keep2), ...
                                                    dataNeurons1_perm(c.test(foldIdx),keep1)    , dataNeurons2_perm(c.test(foldIdx),keep2)   );
                                            end

                                            %- average across cross-validations
                                            cvr_perm = nanmean(cvr_perm);

                                            %- keep values!
                                            cvr_all_perm(1,t2) = single(cvr_perm(1));
                                            cvr2_all_perm(1,t2) = single(cvr_perm(2));
                                        else
                                            cvr_all_perm(1,t2) = NaN;
                                            cvr2_all_perm(1,t2) = NaN;

                                        end
                                    else
                                        cvr_all(1,t2) = NaN;
                                        cvr_all_perm(1,t2) = NaN;
                                        cvr2_all(1,t2) = NaN;
                                        cvr2_all_perm(1,t2) = NaN;
                                        coeffvar_all(1,t2) = NaN;
                                        coeffvar2_all(1,t2) = NaN;
                                        cva_all(:,t2) = NaN(size(cva_all(:,t2)));
                                        cvb_all(:,t2) = NaN(size(cvb_all(:,t2)));
                                    end
                                end
                                res_CCA.cvr{ar,ar2}(nb(ar,ar2),t1,:) = single(cvr_all);
                                res_CCA.cvr_perm{ar,ar2}(nb(ar,ar2),t1,:) = single(cvr_all_perm);
                                res_CCA.cvr2{ar,ar2}(nb(ar,ar2),t1,:) = single(cvr2_all);
                                res_CCA.cvr2_perm{ar,ar2}(nb(ar,ar2),t1,:) = single(cvr2_all_perm);
                                res_CCA.coeffvar{ar,ar2}(nb(ar,ar2),t1,:) = single(coeffvar_all);
                                res_CCA.coeffvar2{ar,ar2}(nb(ar,ar2),t1,:) = single(coeffvar2_all);
                                res_CCA.cva{ar,ar2}{nb(ar,ar2)}(t1,:,:) = single(cva_all)';
                                res_CCA.cvb{ar,ar2}{nb(ar,ar2)}(t1,:,:) = single(cvb_all)';
                            end

                            res_CCA.mk{ar,ar2}(nb(ar,ar2),1) = list(s).name(1);
                            res_CCA.sess{ar,ar2}(nb(ar,ar2)) = {list(s).name(1:8)};
                            res_CCA.nbunits{ar,ar2}(nb(ar,ar2),:) = [sum(keep1) sum(keep2)];
                        end
                    end
                end
            end
        end
        clearvars -except s list res_CCA nb param area2test path2go area_list minReplicates savename
    end
    param.list = list;
    clear s list area_list

    save(savename)
end

