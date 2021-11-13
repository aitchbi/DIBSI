function bjt=hb_get_vects(bjt)

x=bjt.x_lowres;
y=bjt.y_lowres;
bjt.vect.f_lowres=bjt.slice.f_lowres(x,y);
bjt.vect.gm_lowres=bjt.slice.gm_lowres(x,y); % for pve correction
bjt.vect.wm_lowres=bjt.slice.wm_lowres(x,y); % for pve correction
bjt.vect.csf_lowres=bjt.slice.csf_lowres(x,y); % for pve correction
bjt.vect.sum_lowres=bjt.slice.sum_lowres(x,y);

x=bjt.x_highres;
y=bjt.y_highres;
bjt.vect.f_lowres_from_highres=bjt.slice.f_highres(x,y); 
bjt.vect.gm_lowres_from_highres=bjt.slice.gm_highres(x,y); %lazeme?
bjt.vect.wm_lowres_from_highres=bjt.slice.wm_highres(x,y); %lazeme?
bjt.vect.csf_lowres_from_highres=bjt.slice.csf_highres(x,y); %lazeme?
bjt.vect.sum_lowres_from_highres=bjt.slice.sum_highres(x,y); %lazeme?
if bjt.phase_pvecorr
    bjt.vect.f_highres_pvecorrected=bjt.slice.f_highres_pvecorrected(x,y);
    bjt.vect.f_highres_pvecorrected_std=bjt.slice.f_highres_pvecorrected_std(x,y);
end