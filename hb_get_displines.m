function [cSP_v11,cSP_v11_sum,cSP_v11_fullRange,S_nonNorm_fullRange_v11,S_nonNorm_sum_v11]=hb_get_displines(sOrder,sz,xx,intRange,KsPrimaryDomain,C,lFact,logiCenter,ee,varargin)
%
%
%
%
%
% Hamid Behjat
% March 2017.


control_params = {
    'sensitivity','minimal',...
    'weightDSF',1,...
    'checkUsefullnessOfDomainInfo',0,...
    };
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

kAdjust =0;
if length(intRange)~=length(KsPrimaryDomain) 
    error('ooopps..')
else
    Ns=length(intRange);
    dummy=zeros(1,intRange(end));
    dummy(intRange)=KsPrimaryDomain;
    KsPrimaryDomain=dummy;
end
L_end=xx(1);
U_end=xx(end); % xx(end) should be an integer, i.e. sample index. 

linearSP=hb_get_spline(sOrder,sz);
% Domain Similarity Matrix for v4 .........................................
delta=(sOrder+1)/2;
dummy = (-delta:sz:delta)';
D_v4 = zeros(numel(dummy),U_end); % keeping it large, i.e. U_end columns rather than Ns, for compatability with cSpline function, etc.
%fprintf('\n. computing domain similarity lookup table..\n')
for n=1:Ns
    k=intRange(n);
    if weightDSF
        numer = zeros(size(dummy));
        for j=1:numel(C)
            numer = numer + abs(C{j}(k+dummy)-C{j}(k));
        end
        if any(isnan(numer))
            %fprintf('\nsample %d skipped\n',k);
            continue
        end
        numer = 1-(numer/numel(C));
    else
        numer = ones(size(dummy));
    end
    %logiCenter=1/numel(C);
    if KsPrimaryDomain(k+kAdjust)>0
        numer = myLogistic(C{KsPrimaryDomain(k+kAdjust)}(k+dummy)-logiCenter,lFact) .* numer;
    else % if there exists two or more subdomains with equal probability
        % The below code is a loose fix for now. Not thouroughly tested.
        % i.e. its ok if number of constraints is 2, & intersection of the
        % constraints at integer k, where constraint values are equal,
        % is cross like. (15 Nov. 2016)
        
        for l=1:numel(dummy)
            Cs = zeros(1,numel(C));
            for j=1:numel(C)
                Cs(j) = C{j}(k+dummy(l));
            end
            [minC,~] = min(Cs,[],2);
            numer(l) = myLogistic(minC-logiCenter,lFact) * numer(l);
        end
    end
    D_v4(:,k+kAdjust) = numer;
end
%fprintf(' done.\n')

% initialize stuff ........................................................
SS_v11 = zeros(numel(xx),1);
S_sum_v11 = zeros(numel(xx),1);
S_nonNorm_sum_v11 = zeros(numel(xx),1);
cSP_v11_sum = zeros(numel(xx),1);
S_fullRange_v11 = cell(Ns,1);
S_nonNorm_fullRange_v11 = cell(Ns,1);
cSP_v11_fullRange = cell(Ns,1);
cSP_fullRange_std = cell(Ns,1);
for n = 1:Ns
    S_fullRange_v11{n}=zeros(numel(xx),1);
    S_nonNorm_fullRange_v11{n}=zeros(numel(xx),1);
    cSP_v11_fullRange{n}=zeros(numel(xx),1);
    cSP_fullRange_std{n}=zeros(numel(xx),1);
end
S_v11 = cell(Ns,1);
S_nonNorm_v11 = cell(Ns,1);
cSP_v11 = cell(Ns,1); % (uses cSpline_v4.m) -- cSP which has the dsf at non-homogenious supports (slightly modified compared to cSpline_v1.m -- consider homogenious domain if equal at 4 or 5 locations, for non-integer and interger points respectively.) + weighting factor

% compute DI splines ......................................................
count=1;
nRange=1+3:Ns-3;%L_end+3:U_end-3;
%fprintf('\n. computing DI splines..')
for n=nRange
    k=intRange(n);
    count = count+1;
    fastDSF=1;
    
    [cSP_v11{n},x, S_v11{n},~,S_nonNorm_v11{n}]...
        =cSpline_v4(k,C,sOrder,delta,sz,ee,D_v4,KsPrimaryDomain,kAdjust,...
        lFact,fastDSF,...
        'sensitivity',sensitivity,...
        'weightDSF',weightDSF,...
        'logiCenter',logiCenter,...
        'checkUsefullnessOfDomainInfo',checkUsefullnessOfDomainInfo);
    
    indCent = (k-L_end)*(1/sz)+1;
    temp = double((numel(cSP_v11{n})-1)/2);
    if rem(numel(x),2)
        SS_v11(indCent-temp:indCent+temp)...
            =SS_v11(indCent-temp:indCent+temp)+S_v11{n};
        S_fullRange_v11{n}(indCent-temp:indCent+temp)...
            =S_fullRange_v11{n}(indCent-temp:indCent+temp)+S_v11{n};
        S_nonNorm_fullRange_v11{n}(indCent-temp:indCent+temp)...
            =S_nonNorm_fullRange_v11{n}(indCent-temp:indCent+temp)+S_nonNorm_v11{n};
        S_sum_v11...
            =S_sum_v11+S_fullRange_v11{n};
        S_nonNorm_sum_v11...
            =S_nonNorm_sum_v11+S_nonNorm_fullRange_v11{n};
        cSP_v11_fullRange{n}(indCent-temp:indCent+temp)...
            =cSP_v11_fullRange{n}(indCent-temp:indCent+temp)+cSP_v11{n};
        cSP_fullRange_std{n}(indCent-temp:indCent+temp)...
            =cSP_fullRange_std{n}(indCent-temp:indCent+temp)+linearSP;
        cSP_v11_sum...
            =cSP_v11_sum+cSP_v11_fullRange{n};
    else
        fprintf('ooops')
    end
    %hb_progress(n,nRange(end),'init',n==nRange(1));
end
%fprintf(' done.\n')

