function D = hb_get_samplesdomain(I,varargin)
% HB_GET_SAMPLESDOMAIN determines primary subdomain associated to samples.
% 
% Inputs:
%   I: sample indices - (integer values).
%   C: a cell array of Nc function handles, one specifying each subdomain.
%
% Outputs:
%   D: a vector of length(I), specifying the index of the subdomain 
%      associated to each sample; i.e. a value in the range 1:length(C). 
%      If a sample is equally, maximally associated to more than one 
%      subdomain, say k subdomain, the value -k will be returned. 
%      [The function currently suppors upto 5 subdomins, but can be easily 
%       extended if more subdomains are used.]    
%
% 
% Hamid Behjat 
% March 2017.

control_params = {
    'C_funchandle',[],... % Ncx1 cell arary 
    'C_vect',[]};         % NsxNc matrix
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

Ns=numel(I);

if ~isempty(C_funchandle)
    C=C_funchandle;
    Nc=numel(C);
    CCs = zeros(Ns,Nc);
    for n=1:Ns
        for i=1:Nc
            CCs(n,i) = C{i}(I(n));
        end
    end
elseif ~isempty(C_vect)
    if size(C_vect,2)>size(C_vect,1)
        error('each subdomain description should be in one column, not row.');
        % I am assumng number of samples > number of subdomains.  
    end
    Nc=size(C_vect,2);
    CCs=C_vect;
else
    error('ooopps..')
end

[c_indCmax,indCmax] = max(CCs,[],2);

D = zeros(1,Ns);
for n=1:Ns
    Cs=CCs(n,:);
    if nnz(Cs)==0
        fprintf('watch out...! loose fix\n')
        D(n)=3; % csf - a lose fix for now.
        continue
    elseif any(isnan(Cs)) && n~=Ns
        error('fishy..')
        %Cs(indNan) = 0;
    end
    dummy1 = heaviside(c_indCmax(n) - 1/Nc);
    if dummy1>=0 % i.e. 0.5 or 1
        if dummy1==0 %it means that partition of unity was not satisfied
            fprintf('[!] partion of unity not satisfied on sample: %d\n',n);
        end
        if numel(uniquetol(Cs))==numel(sort(Cs)) && all(uniquetol(Cs)==sort(Cs))
            D(n) = indCmax(n);
        else
            dummy2 = sort(Cs);
            if round(dummy2(end)*1e3)==round(dummy2(end-1)*1e3)
                if numel(dummy2)>2 && round(dummy2(end-1)*1e3)==round(dummy2(end-2)*1e3)
                    if numel(dummy2)>3 && round(dummy2(end-2)*1e3)==round(dummy2(end-3)*1e3)
                        if numel(dummy2)>4 && round(dummy2(end-3)*1e3)==round(dummy2(end-4)*1e3)
                            D(n) = -5; 
                        else
                            D(n) = -4; 
                        end
                    else
                        D(n) = -3; 
                    end
                else
                    D(n) = -2; % i.e. Integer k has equal probability at 2 subdomains
                end
            else
                D(n) = indCmax(n); % i.e. Integer k has highest probability at subdomain indCmax(k)
            end
        end
        
    elseif isnan(dummy1)
        D(n) = D(n-1); % a loos fix
    end
end

