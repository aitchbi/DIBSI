function [spl,cc,xx]=hb_plot_splines(L_end,U_end,sOrder,sz,varargin)
%
%
%
%
%
% Hamid Behjat
% March 2017.

control_params = {
    'hFig',[],...
    'subplot_pos',[],...
    };
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

if any([isempty(hFig),isempty(subplot_pos)])%#ok<NODEF>
    hFig=figure;
    set(hFig,'Position',[50,2000,2500,100])
else
    figure(hFig);
    subplot(subplot_pos(1),subplot_pos(2),subplot_pos(3))
end

if 1
    [spl,~,~,x] = hb_get_spline(sOrder,sz);
    for shift=L_end+1:U_end-1
        plot(x+shift,spl,'linewidth',2);
        if shift==L_end+1
            hold on
        end
    end
    set(gca, 'XLim',[L_end,U_end],'XTick',(L_end:U_end),'Box','off','YLim',[-0.05,1.2],'XGrid','on')
else
    switch sOrder
        case 1 % first order B-splines
            [spl,cc,xx,x] = hb_get_spline(sOrder,sz);
            for shift=L_end+1:U_end-1
                if 0 % HB 7 oct 2017
                    plot([xx{1},xx{2}]+shift,...
                        [polyval(cc{1},xx{1}),polyval(cc{2},xx{2})],'linewidth',2);
                else
                    plot(x+shift,spl,'linewidth',2);
                end
                if shift==L_end+1
                    hold on
                end
            end
            set(gca, 'XLim',[L_end,U_end],'XTick',(L_end:U_end),'Box','off','YLim',[-0.05,1.2],'XGrid','on')
            
        case 3 % Cubic B-splines
            [spl,cc,xx,x]=hb_get_splines(sOrder,sz);
            for shift = L_end+2:U_end-2
                if 0 % HB 7 oct 2017
                    plot([xx{1},xx{2},xx{3},xx{4}]+shift,...
                        [polyval(cc{1},xx{1}),polyval(cc{2},xx{2}),...
                        polyval(cc{3},xx{3}),polyval(cc{4},xx{4})],'linewidth',2)
                else
                    plot(x+shift,spl,'linewidth',2)
                end
                if shift==L_end+2
                    hold on
                end
            end
            set(gca, 'XTick',(L_end:U_end),'Box','off','YLim',[-0.05 2],'XGrid','on')
    end
end

