function r = CanonCorrFitAndPredict_triplet(Atrain, Btrain, Ctrain, Atest, Btest, Ctest)
warning('off','stats:canoncorr:NotFullRank');
lastwarn('')

[betaA_AB, betaB_AB, ~] = canoncorr( Atrain, Btrain );
[~, id1] = lastwarn;
[betaA_AC, betaC_AC, ~] = canoncorr( Atrain, Ctrain );
[~, id2] = lastwarn;
[betaB_BC, betaC_BC, ~] = canoncorr( Btrain, Ctrain );
[~, id3] = lastwarn;

if ~isempty(id1) && ~isempty(id2)  && ~isempty(id3) && strcmp(id1, 'stats:canoncorr:NotFullRank')  && strcmp(id2, 'stats:canoncorr:NotFullRank') && strcmp(id3, 'stats:canoncorr:NotFullRank')
    r = nan(1,15);
else
    %- r(1) = AB dir / r(2) = AC dir / r(3) = BC dir 
    %- r(4) = AB ind(A*C-B*C) / r(5) = AC ind(A*B-BC*) / r(6) = BC ind(AB*-AC*) 
    %- r(7) = AB ind(A*B-B*C) / r(8) = AC ind(A*B-AC*) / r(9) = BC ind(B*C-AC*) 
    %- r(10)= AB ind(A*C-AB*) / r(11)= AC ind(A*C-BC*) / r(12)= BC ind(AB*-BC*) 
    %- r(13)= AA btw(A*B-A*C) / r(14)= BB btw(AB*-B*C) / r(15)= CC btw(AC*-BC*) 

    r(1) = corr( Atest*betaA_AB(:,1), Btest*betaB_AB(:,1) );
    r(2) = corr( Atest*betaA_AC(:,1), Ctest*betaC_AC(:,1) );
    r(3) = corr( Btest*betaB_BC(:,1), Ctest*betaC_BC(:,1) );
    r(4) = corr( Atest*betaA_AC(:,1), Btest*betaB_BC(:,1) );
    r(5) = corr( Atest*betaA_AB(:,1), Ctest*betaC_BC(:,1) );
    r(6) = corr( Btest*betaB_AB(:,1), Ctest*betaC_AC(:,1) );
    r(7) = corr( Atest*betaA_AB(:,1), Btest*betaB_BC(:,1) );
    r(8) = corr( Atest*betaA_AB(:,1), Ctest*betaC_AC(:,1) );
    r(9) = corr( Btest*betaB_BC(:,1), Ctest*betaC_AC(:,1) );
    r(10) = corr( Atest*betaA_AC(:,1), Btest*betaB_AB(:,1) );
    r(11) = corr( Atest*betaA_AC(:,1), Ctest*betaC_BC(:,1) );
    r(12) = corr( Btest*betaB_AB(:,1), Ctest*betaC_BC(:,1) );
    r(13) = corr( Atest*betaA_AB(:,1), Atest*betaA_AC(:,1) );
    r(14) = corr( Btest*betaB_AB(:,1), Btest*betaB_BC(:,1) );
    r(15) = corr( Ctest*betaC_AC(:,1), Ctest*betaC_BC(:,1) );
    
end

warning('on','stats:canoncorr:NotFullRank');
end
