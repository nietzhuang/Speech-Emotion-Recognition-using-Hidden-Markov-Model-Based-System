function [TRAN,EMISGmm,logliks] = GmmHmmtrain(seqs,TRAN_ini,EMISGmm_ini,varargin)
% A hybrid model GMM/HMM developed from function hmmtrain.
%
%
%
%


if nargin > 0
    if isstring(seqs)
        seqs = cellstr(seqs);
    end
end
if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

tol = 1e-6;
trtol = tol;
etol = tol;
maxiter = 500;
verbose = false;
[numStates, checkTr] = size(TRAN_ini);
if checkTr ~= numStates
    error(message('stats:hmmtrain:BadTransitions'));
end

% number of rows of EMIS must be same as number of states
[checkE, checkEparam] = size(EMISGmm_ini);
if checkE ~= numStates
    error(message('stats:hmmtrain:InputSizeMismatch'));
end
if checkEparam ~= 3
	error(message('Insufficient GMM Parameters.'))
end
%if (numStates ==0 || numEmissions == 0)
%    guessTR = [];
%    guessE = [];
%    return
%end
numComp = size(EMISGmm_ini{1, 1}, 2);
if size(EMISGmm_ini{1, 2}, 1) ~= numComp || size(EMISGmm_ini{1, 3}, 3) ~= numComp
	error(message('stats:hmmtrain:InputSizeMismatch'));
end



baumwelch = true;

if nargin > 3
    if rem(nargin,2)== 0
        error(message('stats:hmmtrain:WrongNumberArgs', mfilename));
    end
    okargs = {'tolerance','maxiterations','verbose','algorithm','trtol','etol'};    
    dflts  = {[]           maxiter         verbose   ''           []      []};
    [tol,maxiter,verbose,alg,trtol,etol] = ...    
        internal.stats.parseArgs(okargs, dflts, varargin{:}); 
     
    % deal with a single sequence first
    if ~iscell(seqs) || ischar(seqs{1})
        [~, seqs]  = ismember(seqs,symbols);
        if any(seqs(:)==0)
            error(message('stats:hmmtrain:MissingSymbol'));
        end
    else  % now deal with a cell array of sequences
        %numSeqs = numel(seqs);
        %newSeqs = cell(numSeqs,1);
        %for count = 1:numSeqs
        %    [~, newSeqs{count}] = ismember(seqs{count},symbols);
        %    if any(newSeqs{count}(:)==0)
        %        error(message('stats:hmmtrain:MissingSymbol'));
        %    end
        %end
        %seqs = newSeqs;
    end
    if ischar(verbose)
        verbose = any(strcmpi(verbose,{'on','true','yes'}));
    end
     
    if ~isempty(alg)
        alg = internal.stats.getParamVal(alg,{'baumwelch','viterbi'},'Algorithm');
        baumwelch = strcmpi(alg,'baumwelch');
    end
end

if isempty(tol)
    tol = 1e-6;
end
if isempty(trtol)
    trtol = tol;
end
if isempty(etol)
    etol = tol;
end

if isnumeric(seqs)
	error(message('stats:hmmtrain:BadSequence'));  % Not support single seqence. Update afterward.
elseif iscell(seqs)
    numSeqs = numel(seqs);
    cellflag = true;
else
    error(message('stats:hmmtrain:BadSequence'));
end

% initialize the counters
TR = zeros(size(TRAN_ini));
TRAN = TRAN_ini;
EMISGmm = EMISGmm_ini;

%if ~pseudoTRcounts
%    pseudoTR = TR;
%end
%E = zeros(numStates,numEmissions);

%if ~pseudoEcounts
%    pseudoE = E;
%end

converged = false;
loglik = 1; % loglik is the log likelihood of all sequences given the TR and E
logliks = zeros(1,maxiter);
for iteration = 1:maxiter
    oldLL = loglik;
    loglik = 0;
    oldEMISGmm = EMISGmm;
    oldTRAN = TRAN;
    for Seqcount = 1:numSeqs
        if cellflag
            seq = seqs{Seqcount};
            seqLength = length(seq);
			dim = size(seq, 1);
        else
			error(message('stats:hmmtrain:BadSequence'));  % Not support single seqence. Update afterward.
            %seq = seqs(Seqcount,:);
        end
		
		% Construct emission matrix with GMM
		E = zeros(numStates, seqLength);
		for count = 1:seqLength
			for state = 1:numStates
				bj = 0;
				for comp = 1:numComp
					bj = bj + EMISGmm{state, 1}(comp) * mvnpdf(seq(:, count)', EMISGmm{state, 2}(comp, :), EMISGmm{state, 3}(:, :, comp)); 			
				end
				E(state, count) = bj;
			end
        end        
	    
        if baumwelch   % Baum-Welch training
            % get the scaled forward and backward probabilities
            [~,logPseq,fs,bs,scale] = GmmHmmdecode(seq,TRAN,EMISGmm);
 				           
            loglik = loglik + logPseq;
            logfs = log(fs);
            logbs = log(bs);
            logE = log(E);
            logTR = log(TRAN);
            % f and b start at 0 so offset seq by one
            seq = [zeros(dim, 1) seq];
			% Calculate rho
			rho = zeros(numStates, numComp, seqLength+1);	 % From now on, seq is shifted.		
			for count = 2:seqLength+1
				for statej = 1:numStates
					for comp = 1:numComp					
						if count == 2
							rho(statej, comp, count) = EMISGmm{statej, 1}(comp) * mvnpdf(seq(:, count)', EMISGmm{statej, 2}(comp, :), EMISGmm{statej, 3}(:, :, comp));
						else
							rho(statej, comp, count) = sum(exp(logfs(:, count-1) + logTR(:, statej))) * ...
													   EMISGmm{statej, 1}(comp) * mvnpdf(seq(:, count)', EMISGmm{statej, 2}(comp, :), EMISGmm{statej, 3}(:, :, comp));
						end
					end							
				end			
            end
            % Update transition matrix
			TR = zeros(numStates, numStates);
            for statei = 1:numStates
                for statej = 1:numStates
					numerator = 0;
					denominator = 0;
                    for count = 2:seqLength  % from t=1 to t=T-1.
						numerator = numerator + exp(logfs(statei, count) + logTR(statei, statej) + logE(statej,count) + logbs(statej, count+1)); % the seq in emission matrix is not shifted.
						denominator = denominator + exp(logfs(statei, count) + logbs(statei, count)) .* scale(count);  % elimate scaled effect.  % scale(count+1) should be check!!!!
                    end
					TR(statei, statej) = numerator/denominator;
                end
            end	
			% Update weight values
			weight = zeros(numStates, numComp);
			for statej = 1:numStates
				for comp = 1:numComp
					numerator = 0;
					denominator = 0;
					for count = 2:seqLength+1
						numerator = numerator + rho(statej, comp, count) .* exp(logbs(statej, count));  % rho is shifed.
						denominator = denominator + exp(logfs(statej, count) + logbs(statej, count)) .* scale(count);
					end				
					weight(statej, comp) = numerator/denominator;
				end
			end			
			% Update mean
			mean = zeros(numComp, dim, numStates);
			for statej = 1:numStates
				for comp = 1:numComp
					numerator = zeros(dim, 1);
					denominator = zeros(dim, 1);
					for count = 2:seqLength+1
						numerator = numerator + repmat(rho(statej, comp, count), [dim 1]).*repmat(exp(logbs(statej, count)), [dim 1]).*seq(:, count);
						denominator = denominator + repmat(rho(statej, comp, count), [dim 1]).*repmat(exp(logbs(statej, count)), [dim 1]);							
					end
					mean(comp, :, statej) = reshape(numerator./denominator, [1 dim]);
				end
			end
			% Update covariance
			cov = zeros(dim, dim, numComp);
			for statej = 1:numStates
				for comp = 1:numComp
					numerator = zeros(dim, 1);
					denominator = zeros(dim, 1);
					for count = 2:seqLength+1
						numerator = numerator + repmat(rho(statej, comp, count), [dim dim]).*repmat(exp(logbs(statej, count)), [dim dim]).*((seq(:, count)-EMISGmm{statej, 2}(comp, :)')*(seq(:, count)-EMISGmm{statej, 2}(comp, :)')');
						denominator = denominator + repmat(rho(statej, comp, count), [dim dim]).*repmat(exp(logbs(statej, count)), [dim dim]);
					end
					cov(:, :, comp, statej) = numerator./denominator;
				end
			end
			% Combine parameters			
			TRAN = TR;
			for state = 1:numStates
				EMISGmm(state, :) = {weight(state, :) mean(:, :, state) cov(:, :, :, state)};
            end            
			
			
			% Check invalidation of paramerters
			if sum(sum(isnan(TRAN))) % Check transition matrix
				TRAN = oldTRAN;
				EMISGmm = oldEMISGmm;
			end										
			for state = 1:numStates
				for comp = 1:numComp
				    [R, p] = cholcov(EMISGmm{state, 3}(:, :, comp), 0);  % used to check positive defined.
					if p ~= 0
						TRAN = oldTRAN;
						EMISGmm{state, 3}(:, :, comp) = oldEMISGmm{state, 3}(:, :, comp);
					end
					if sum(isnan(EMISGmm{state, 1})) || sum(isnan(EMISGmm{state, 2}(comp, :))) || isinf(sum(log(diag(R))))
					
					end																		
				end						
            end
			% Check singularity
			for state = 1:numStates
				for comp = 1:numComp
					if det(EMISGmm{state, 3}(:, :, comp)) < 1e-150  % singularity constraint.                        
                        fprintf("%d-state, %d-th component covariance singularity caused! at %d iteration times.\n", state, comp, iteration);
						TRAN = oldTRAN;  % ignore the update this time.
					    EMISGmm{state, 3}(:, :, comp) = oldEMISGmm{state, 3}(:, :, comp);
					end
				end
			end
			
			% Check whether new model produce impossible likelihood
            seq = seq(:, 2:end);  % eliminate zero padding.
			[~,logPseq] = GmmHmmdecode(seq,TRAN,EMISGmm);
			if isinf(logPseq) || isnan(logPseq)  % Check invalidation of P(O|£f), ignore updated if invalid. 
                fprintf("%d-th sequence produces impossible likelihood at iteration times: %d\n", Seqcount, iteration-1);
                TRAN = oldTRAN;
		    	EMISGmm = oldEMISGmm;
            end
			
			
			
			%for state = 1:numStates
			%	for comp = 1:numComp
			%		if isnan(EMISGmm{state, 1}) % Check weight parameter
			%			EMISGmm{state, 1} = oldEMISGmm{state, 1};
			%		end	
			%		if isnan(EMISGmm{state, 2}(comp, :)) % Check mean parameter
			%			EMISGmm{state, 2}(comp, :) = oldEMISGmm{state, 2}(comp, :);
			%		end					
			%		[R, p] = cholcov(EMISGmm{state, 3}(:, :, comp), 0); % Check covariance parameter
			%		if p~=0
			%			%fprintf("Not positive defined Covariance, Set as last one\n");
			%			EMISGmm{state, 3}(:, :, comp) = oldEMISGmm{state, 3}(:, :, comp);
			%		end
			%		if isinf(sum(log(diag(R))))
			%			%fprintf("Bad Covariance, Set as last one.\n");
			%			EMISGmm{state, 3}(:, :, comp) = oldEMISGmm{state, 3}(:, :, comp);
			%		end
			%	end
			%end			
			
        end
    end
    %totalEmissions = sum(E,2);
    %totalTransitions = sum(TR,2);
    
    % avoid divide by zero warnings
    %guessE = E./(repmat(totalEmissions,1,numEmissions));
    %guessTR  = TR./(repmat(totalTransitions,1,numStates));
    % if any rows have zero transitions then assume that there are no
    % transitions out of the state.
    %if any(totalTransitions == 0)
    %    noTransitionRows = find(totalTransitions == 0);
    %    guessTR(noTransitionRows,:) = 0;
    %    guessTR(sub2ind(size(guessTR),noTransitionRows,noTransitionRows)) = 1;
    %end
    % clean up any remaining Nans
    %guessTR(isnan(guessTR)) = 0;
    %guessE(isnan(guessE)) = 0;
    
    %if verbose
    %    if iteration == 1
    %        fprintf('%s\n',getString(message('stats:hmmtrain:RelativeChanges')));
    %        fprintf('   Iteration       Log Lik    Transition     Emmission\n');
    %    else 
    %        fprintf('  %6d      %12g  %12g  %12g\n', iteration, ...
    %            (abs(loglik-oldLL)./(1+abs(oldLL))), ...
    %            norm(guessTR - oldGuessTR,inf)./numStates, ...
    %            norm(guessE - oldGuessE,inf)./numEmissions);
    %    end
    %end
	
    % Durbin et al recommend loglik as the convergence criteria  -- we also
    % use change in TR and E. Use (undocumented) option trtol and
    % etol to set the convergence tolerance for these independently.
    %
	

	
	% Convergence constrain
    logliks(iteration) = loglik;
    if (abs(loglik-oldLL)/(1+abs(oldLL))) < tol
        if norm(TRAN - oldTRAN,inf)/numStates < trtol
			for statej = 1:numStates
				delta(statej) = sum(sum(sum((EMISGmm{statej, 3} - oldEMISGmm{statej, 3}).^2)));
			end
			if sum(delta)/numStates < etol
				for statej = 1:numStates
					delta(statej) = sum(sum((EMISGmm{statej, 2} - oldEMISGmm{statej, 2}).^2));
				end				
				if sum(delta)/numStates < etol
					for statej = 1:numStates
						delta(statej) = sum((EMISGmm{statej, 1} - oldEMISGmm{statej, 1}).^2);
					end
					if sum(delta)/numStates < etol
						converged = true;
						break;
					end
				end
			end
        end
    end
    %E =  pseudoE;
    %TR = pseudoTR;
end
if ~converged
    warning(message('stats:hmmtrain:NoConvergence', num2str( tol ), maxiter));
else
	fprintf("Converged at iteration times: %d\n\n", iteration-1);
end
logliks(logliks ==0) = [];
