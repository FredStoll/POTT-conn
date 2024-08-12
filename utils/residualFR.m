function dataNeurons = residualFR(dataNeurons,param,factors_full,predic_full)


parfor j = 1 : size(dataNeurons,2) %- for each time bin

    if strcmp(param.predic,'allstim')
        params2keep = [1 2 3];
        [Ap,At,Amdl,~] = anovan(dataNeurons(:,j),factors_full(:,params2keep),'continuous',1,'varname',predic_full(params2keep),'model',eye(length(params2keep)),'display','off');
    elseif strcmp(param.predic,'allrew')
        params2keep = [4 5 1 3];
        [Ap,At,Amdl,~] = anovan(dataNeurons(:,j),factors_full(:,params2keep),'continuous',3,'varname',predic_full(params2keep),'nested',[0 0 0 0; 1 0 0 0; 0 0 0 0 ; 0 0 0 0],'model',eye(length(params2keep)),'display','off');
    elseif strcmp(param.predic,'chosenproba')
        params2keep = 1;
        [Ap,At,Amdl,~] = anovan(dataNeurons(:,j),factors_full(:,params2keep),'continuous',1,'varname',predic_full(params2keep),'model',eye(length(params2keep)),'display','off');
    else
        params2keep = find(ismember(predic_full,param.predic));
        [Ap,At,Amdl,~] = anovan(dataNeurons(:,j),factors_full(:,params2keep),'continuous',[],'varname',predic_full(params2keep),'model',eye(length(params2keep)),'display','off');
    end

    dataNeurons(:,j) = Amdl.resid; % replace with the residuals

end

