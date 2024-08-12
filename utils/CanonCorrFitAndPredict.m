function [r,r_train,A,B] = CanonCorrFitAndPredict(Xtrain, Ytrain, Xtest, Ytest)
%function [r,r_train,cc] = CanonCorrFitAndPredict(Xtrain, Ytrain, Xtest, Ytest)
warning('off','stats:canoncorr:NotFullRank');
lastwarn('')

[A, B, r_train] = canoncorr( Xtrain, Ytrain );
[~, id] = lastwarn;

numPairs = min( size(Xtrain, 2), size(Ytrain, 2) );
if ~isempty(id) && strcmp(id, 'stats:canoncorr:NotFullRank')
	r = nan(1, numPairs);
else
	r = zeros(1, numPairs);
	for pairIdx = 1:numPairs
		r(pairIdx) = corr( Xtest*A(:,pairIdx), Ytest*B(:,pairIdx) );
	end
end

warning('on','stats:canoncorr:NotFullRank');
%cc.A = A;
%cc.B = B;
A = A(:,1); %- keep 1st comp only
B = B(:,1);

end
