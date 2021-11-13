function [spl,cc,xx,x,ic,im,ixr,l_segs,l_half,p]=hb_get_spline(sOrder,sz,varargin)
%
% ic: ind center
% im: ind mid integers (integers: sample points. If T=1, then its integers)
% ixr: ind range excluding center & ends.
%
% sOrder: B-spline order. Using the 'analytic' approach, B-spline of any 
% order >=0 can be specified, but for the 'piecewise' scheme, I have only 
% coded for orders 0,1,3.   
% 
% Hamid Behjat
% March/Aug 2017.

cntr = {
    'T',1,... % sampling resolution
    'approach','analytic'};
argselectAssign(cntr);
argselectCheck(cntr,varargin);
argselectAssign(varargin);

szT = sz/T;

switch approach
    case 'piecewise'
        switch sOrder
            case 0
                cc{1} = 1;
                xx{1} = -0.5+szT:szT:0.5-szT;
                d = polyval(cc{1},xx{1});
                spl = [0.5,d,0.5];
            case 1 % first order B-splines
                cc{1} = [1,1];
                cc{2} = [-1,1];
                xx{1} =    -1:szT:0;
                xx{2} = 0+szT:szT:1;
                d = [polyval(cc{1},xx{1}),polyval(cc{2},xx{2})];
                spl = d;%[0,d,0];
                x = [xx{1},xx{2}];
                ic = ceil(numel(x)/2);
                im = [];
            case 2
                
            case 3 % Cubic B-splines
                cc{1} = [1/6,1,2,4/3];
                cc{2} = [-1/2,-1,0,2/3];
                cc{3} = [1/2,-1,0,2/3];
                cc{4} = [-1/6,1,-2,4/3];
                xx{1} = -2:szT:-1;
                xx{2} = -1+szT:szT:0;
                xx{3} =  0+szT:szT:1;
                xx{4} =  1+szT:szT:2;
                p{1} = polyval(cc{1},xx{1});
                p{2} = polyval(cc{2},xx{2});
                p{3} = polyval(cc{3},xx{3});
                p{4} = polyval(cc{4},xx{4});
                spl = [p{1},p{2},p{3},p{4}];%[0,d,0];
                x = [xx{1},xx{2},xx{3},xx{4}];
                ic = ceil(numel(x)/2);
                d = round(T/sz);
                im = [1+d:d:ic-1,ic+d:d:numel(x)-1];
            otherwise
                error('oops..')
        end
        ixr = [2:ic-1,ic+1:numel(x)-1];
        spl = spl(:);
        x = x(:)*T;
        
        if isempty(im)
            l_segs = [];
        else
            l_segs = zeros(1,length(im)/2); % NOTE 1
            for iL=1:length(l_segs)
                l_segs(iL) = im(iL)-1;
            end
        end
        l_half = ic-1;
        
    case 'analytic'
        % cf. Unser1999,p.24,Box1,Eq.(10)

        n = sOrder;
        switch n
            case 0
                x = -0.5:szT:0.5;
                spl = [0.5,ones(1,numel(x)-2),0.5];
            otherwise
                delta = (n+1)/2;
                x = -delta:szT:delta;
                
                d1 = zeros(size(x));
                for k = 0:n+1
                    d2 = x-k+(n+1)/2;
                    d1 = d1+nchoosek(n+1,k)*(-1)^k*((d2>=0).*d2).^n;
                end
                spl = 1/factorial(n)*d1;
        end
        
        ic = ceil(numel(x)/2);
        if n<=1
            im = [];
        else
            rTsz = round(T/sz);
            if rem(n,2) % odd order B-spline
                d = rTsz;
            else        % even order B-spline
                d = rTsz/2; 
            end
            im = [1+d:rTsz:ic-1,ic+rTsz:rTsz:numel(x)-1];
            
            %if any(rem(round(im*100),100)), error; %sanity check  
            %else, im = round(im);
            %end
            
        end
        ixr = [2:ic-1,ic+1:numel(x)-1];
        
        spl = spl(:);
        x = x(:)*T;
        
        cc = [];
        xx = [];
        p = [];
        
        if isempty(im)
            l_segs = [];
        else
            l_segs = zeros(1,length(im)/2); % NOTE 1
            for iL=1:length(l_segs)
                l_segs(iL) = im(iL)-1;
            end
        end
        l_half = ic-1;
end

% NOTES
%
% NOTE 1
% l_segs(1): length of the first segment of the B-spline.
% l_segs(2): length of the first+second segments of the B-spline.
