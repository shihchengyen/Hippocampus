function r = plus(p,q,varargin)
%@dirfiles/plus Overloaded plus function for dirfiles objects.
%   R = plus(P,Q) combines dirfiles objects P and Q and returns the
%   dirfiles object R.

% For vmcorr

% get name of class
classname = mfilename('class');

% check if first input is the right kind of object
if(~isa(p,classname))
	% check if second input is the right kind of object
	if(~isa(q,classname))
		% both inputs are not the right kind of object so create empty
		% object and return it
		r = feval(classname);
	else
		% second input is the right kind of object so return that
		r = q;
	end
else
	if(~isa(q,classname))
		% p is the right kind of object but q is not so just return p
		r = p;
    elseif(isempty(p))
        % p is right object but is empty so return q, which should be
        % right object
        r = q;
    elseif(isempty(q))
        % p are q are both right objects but q is empty while p is not
        % so return p
        r = p;
	else
		% both p and q are the right kind of objects so add them 
		% together
		% assign p to r so that we can be sure we are returning the right
		% object
		r = p;
		% useful fields for most objects
		r.data.numSets = p.data.numSets + q.data.numSets;
        r.data.origin = [p.data.origin; q.data.origin];

		% object specific fields %% Full session
        r.data.pv = [p.data.pv; q.data.pv];
        r.data.ph = [p.data.ph; q.data.ph];
        r.data.gridSteps = [p.data.gridSteps; q.data.gridSteps]; % Is one enough?
        r.data.binDepths = [p.data.binDepths; q.data.binDepths];
        r.data.placebins = [p.data.placebins; q.data.placebins];
        r.data.viewbins = [p.data.viewbins; q.data.viewbins];
        r.data.headdirectionbins = [p.data.headdirectionbins; q.data.headdirectionbins];
        r.data.Args = p.data.Args;
        
        % add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
        
        
%         r.data.maps_dist_p = [p.data.maps_dist_p; q.data.maps_dist_p];
%         r.data.maps_dist_p_adsm = [p.data.maps_dist_p_adsm; q.data.maps_dist_p_adsm];
%         r.data.maps_dist_sv = [p.data.maps_dist_sv; q.data.maps_dist_sv];
%         r.data.maps_dist_sv_adsm = [p.data.maps_dist_sv_adsm; q.data.maps_dist_sv_adsm];
%         r.data.distratio_p = [p.data.distratio_p; q.data.distratio_p];
%         r.data.distratio_sv = [p.data.distratio_sv; q.data.distratio_sv];
%         r.data.distcorr_p = [p.data.distcorr_p; q.data.distcorr_p];
%         r.data.distcorr_sv = [p.data.distcorr_sv; q.data.distcorr_sv];
%         
%         r.data.maps_raw_corrp = [p.data.maps_raw_corrp; q.data.maps_raw_corrp];
%         r.data.maps_raw_corrpset = [p.data.maps_raw_corrpset; {q.data.maps_raw_corrpset}];
%         r.data.maps_adsm_corrp = [p.data.maps_adsm_corrp; q.data.maps_adsm_corrp];
%         r.data.maps_adsm_corrpset = [p.data.maps_adsm_corrpset; {q.data.maps_adsm_corrpset}];
%         r.data.maps_bcsm_corrp = [p.data.maps_bcsm_corrp; q.data.maps_bcsm_corrp];
%         r.data.maps_bcsm_corrpset = [p.data.maps_bcsm_corrpset; {q.data.maps_bcsm_corrpset}];
%         r.data.maps_dksm_corrp = [p.data.maps_dksm_corrp; q.data.maps_dksm_corrp];
%         r.data.maps_dksm_corrpset = [p.data.maps_dksm_corrpset; {q.data.maps_dksm_corrpset}];
%         r.data.dur_adsm_corrp = [p.data.dur_adsm_corrp; q.data.dur_adsm_corrp];
%         r.data.SIC_adsm_corrp = [p.data.SIC_adsm_corrp; q.data.SIC_adsm_corrp];
%         r.data.SIC_adsm_corrpset = [p.data.SIC_adsm_corrpset; {q.data.SIC_adsm_corrpset}];
%         r.data.SIC_bcsm_corrp = [p.data.SIC_bcsm_corrp; q.data.SIC_bcsm_corrp];
%         r.data.SIC_bcsm_corrpset = [p.data.SIC_bcsm_corrpset; {q.data.SIC_bcsm_corrpset}];
%         r.data.SIC_dksm_corrp = [p.data.SIC_dksm_corrp; q.data.SIC_dksm_corrp];
%         r.data.SIC_dksm_corrpset = [p.data.SIC_dksm_corrpset; {q.data.SIC_dksm_corrpset}];
%         r.data.ISE_adsm_corrp = [p.data.ISE_adsm_corrp; q.data.ISE_adsm_corrp];
%         r.data.ISE_adsm_corrpset = [p.data.ISE_adsm_corrpset; {q.data.ISE_adsm_corrpset}];
% %         r.data.ratesplaceMethod2 = [p.data.ratesplaceMethod2; q.data.ratesplaceMethod2];
% %         r.data.ratesviewMethod2 = [p.data.ratesviewMethod2; q.data.ratesviewMethod2];
%         
%         r.data.maps_raw_corrsv = [p.data.maps_raw_corrsv; q.data.maps_raw_corrsv];
%         r.data.maps_raw_corrsvset = [p.data.maps_raw_corrsvset; {q.data.maps_raw_corrsvset}];
%         r.data.maps_adsm_corrsv = [p.data.maps_adsm_corrsv; q.data.maps_adsm_corrsv];
%         r.data.maps_adsm_corrsvset = [p.data.maps_adsm_corrsvset; {q.data.maps_adsm_corrsvset}];
%         r.data.maps_bcsm_corrsv = [p.data.maps_bcsm_corrsv; q.data.maps_bcsm_corrsv];
%         r.data.maps_bcsm_corrsvset = [p.data.maps_bcsm_corrsvset; {q.data.maps_bcsm_corrsvset}];
%         r.data.maps_dksm_corrsv = [p.data.maps_dksm_corrsv; q.data.maps_dksm_corrsv];
%         r.data.maps_dksm_corrsvset = [p.data.maps_dksm_corrsvset; {q.data.maps_dksm_corrsvset}];
%         r.data.dur_adsm_corrsv = [p.data.dur_adsm_corrsv; q.data.dur_adsm_corrsv];
%         r.data.SIC_adsm_corrsv = [p.data.SIC_adsm_corrsv; q.data.SIC_adsm_corrsv];
%         r.data.SIC_adsm_corrsvset = [p.data.SIC_adsm_corrsvset; {q.data.SIC_adsm_corrsvset}];
%         r.data.SIC_bcsm_corrsv = [p.data.SIC_bcsm_corrsv; q.data.SIC_bcsm_corrsv];
%         r.data.SIC_bcsm_corrsvset = [p.data.SIC_bcsm_corrsvset; {q.data.SIC_bcsm_corrsvset}];
%         r.data.SIC_dksm_corrsv = [p.data.SIC_dksm_corrsv; q.data.SIC_dksm_corrsv];
%         r.data.SIC_dksm_corrsvset = [p.data.SIC_dksm_corrsvset; {q.data.SIC_dksm_corrsvset}];
%         r.data.ISE_adsm_corrsv = [p.data.ISE_adsm_corrsv; q.data.ISE_adsm_corrsv];
%         r.data.ISE_adsm_corrsvset = [p.data.ISE_adsm_corrsvset; {q.data.ISE_adsm_corrsvset}];
%         r.data.covmat = [p.data.covmat; {q.data.covmat}];
%         r.data.covmat_norm = [p.data.covmat_norm; {q.data.covmat_norm}];
%         r.data.l1norm = [p.data.l1norm; q.data.l1norm];
%         r.data.l2norm = [p.data.l2norm; q.data.l2norm];
%         
%         r.data.convergewithsv = [p.data.convergewithsv; q.data.convergewithsv];
%         r.data.convergewithp = [p.data.convergewithp; q.data.convergewithp];
%         r.data.smoothpick = [p.data.smoothpick; q.data.smoothpick];
%         r.data.llhpick = [p.data.llhpick; q.data.llhpick];
%         r.data.NumIterLlh = [p.data.NumIterLlh; q.data.NumIterLlh];
%         r.data.llhpicklabel = [p.data.llhpicklabel; {q.data.llhpicklabel}];
%         r.data.llh = [p.data.llh; {q.data.llh}];
%         
%         % 1st half
%         r.data.maps_dist_p1 = [p.data.maps_dist_p1; q.data.maps_dist_p1];
%         r.data.maps_dist_sv1 = [p.data.maps_dist_sv1; q.data.maps_dist_sv1];
%         r.data.distratio_p1 = [p.data.distratio_p1; q.data.distratio_p1];
%         r.data.distratio_sv1 = [p.data.distratio_sv1; q.data.distratio_sv1];
%         r.data.distcorr_p1 = [p.data.distcorr_p1; q.data.distcorr_p1];
%         r.data.distcorr_sv1 = [p.data.distcorr_sv1; q.data.distcorr_sv1];
%         
%         r.data.maps_adsm_corrp1 = [p.data.maps_adsm_corrp1; q.data.maps_adsm_corrp1];
%         r.data.maps_adsm_corrpset1 = [p.data.maps_adsm_corrpset1; {q.data.maps_adsm_corrpset1}];
%         r.data.maps_bcsm_corrp1 = [p.data.maps_bcsm_corrp1; q.data.maps_bcsm_corrp1];
%         r.data.maps_bcsm_corrpset1 = [p.data.maps_bcsm_corrpset1; {q.data.maps_bcsm_corrpset1}];
%         r.data.maps_dksm_corrp1 = [p.data.maps_dksm_corrp1; q.data.maps_dksm_corrp1];
%         r.data.maps_dksm_corrpset1 = [p.data.maps_dksm_corrpset1; {q.data.maps_dksm_corrpset1}];
%         r.data.SIC_adsm_corrp1 = [p.data.SIC_adsm_corrp1; q.data.SIC_adsm_corrp1];
%         r.data.SIC_adsm_corrpset1 = [p.data.SIC_adsm_corrpset1; {q.data.SIC_adsm_corrpset1}];
%         r.data.SIC_bcsm_corrp1 = [p.data.SIC_bcsm_corrp1; q.data.SIC_bcsm_corrp1];
%         r.data.SIC_bcsm_corrpset1 = [p.data.SIC_bcsm_corrpset1; {q.data.SIC_bcsm_corrpset1}];
%         r.data.SIC_dksm_corrp1 = [p.data.SIC_dksm_corrp1; q.data.SIC_dksm_corrp1];
%         r.data.SIC_dksm_corrpset1 = [p.data.SIC_dksm_corrpset1; {q.data.SIC_dksm_corrpset1}];
%         r.data.ISE_adsm_corrp1 = [p.data.ISE_adsm_corrp1; q.data.ISE_adsm_corrp1];
%         r.data.ISE_adsm_corrpset1 = [p.data.ISE_adsm_corrpset1; {q.data.ISE_adsm_corrpset1}];
% %         r.data.ratesplaceMethod2 = [p.data.ratesplaceMethod2; q.data.ratesplaceMethod2];
% %         r.data.ratesviewMethod2 = [p.data.ratesviewMethod2; q.data.ratesviewMethod2];
%         
%         r.data.maps_adsm_corrsv1 = [p.data.maps_adsm_corrsv1; q.data.maps_adsm_corrsv1];
%         r.data.maps_adsm_corrsvset1 = [p.data.maps_adsm_corrsvset1; {q.data.maps_adsm_corrsvset1}];
%         r.data.maps_bcsm_corrsv1 = [p.data.maps_bcsm_corrsv1; q.data.maps_bcsm_corrsv1];
%         r.data.maps_bcsm_corrsvset1 = [p.data.maps_bcsm_corrsvset1; {q.data.maps_bcsm_corrsvset1}];
%         r.data.maps_dksm_corrsv1 = [p.data.maps_dksm_corrsv1; q.data.maps_dksm_corrsv1];
%         r.data.maps_dksm_corrsvset1 = [p.data.maps_dksm_corrsvset1; {q.data.maps_dksm_corrsvset1}];
%         r.data.SIC_adsm_corrsv1 = [p.data.SIC_adsm_corrsv1; q.data.SIC_adsm_corrsv1];
%         r.data.SIC_adsm_corrsvset1 = [p.data.SIC_adsm_corrsvset1; {q.data.SIC_adsm_corrsvset1}];
%         r.data.SIC_bcsm_corrsv1 = [p.data.SIC_bcsm_corrsv1; q.data.SIC_bcsm_corrsv1];
%         r.data.SIC_bcsm_corrsvset1 = [p.data.SIC_bcsm_corrsvset1; {q.data.SIC_bcsm_corrsvset1}];
%         r.data.SIC_dksm_corrsv1 = [p.data.SIC_dksm_corrsv1; q.data.SIC_dksm_corrsv1];
%         r.data.SIC_dksm_corrsvset1 = [p.data.SIC_dksm_corrsvset1; {q.data.SIC_dksm_corrsvset1}];
%         r.data.ISE_adsm_corrsv1 = [p.data.ISE_adsm_corrsv1; q.data.ISE_adsm_corrsv1];
%         r.data.ISE_adsm_corrsvset1 = [p.data.ISE_adsm_corrsvset1; {q.data.ISE_adsm_corrsvset1}];
% %         r.data.covmat1 = [p.data.covmat1; {q.data.covmat1}];
% %         r.data.covmat_norm1 = [p.data.covmat_norm1; {q.data.covmat_norm1}];
% %         r.data.l1norm1 = [p.data.l1norm1; q.data.l1norm1];
% %         r.data.l2norm1 = [p.data.l2norm1; q.data.l2norm1];
%         
%         % 2nd half
%         r.data.maps_dist_p2 = [p.data.maps_dist_p2; q.data.maps_dist_p2];
%         r.data.maps_dist_sv2 = [p.data.maps_dist_sv2; q.data.maps_dist_sv2];
%         r.data.distratio_p2 = [p.data.distratio_p2; q.data.distratio_p2];
%         r.data.distratio_sv2 = [p.data.distratio_sv2; q.data.distratio_sv2];
%         r.data.distcorr_p2 = [p.data.distcorr_p2; q.data.distcorr_p2];
%         r.data.distcorr_sv2 = [p.data.distcorr_sv2; q.data.distcorr_sv2];
%         
%         r.data.maps_adsm_corrp2 = [p.data.maps_adsm_corrp2; q.data.maps_adsm_corrp2];
%         r.data.maps_adsm_corrpset2 = [p.data.maps_adsm_corrpset2; {q.data.maps_adsm_corrpset2}];
%         r.data.maps_bcsm_corrp2 = [p.data.maps_bcsm_corrp2; q.data.maps_bcsm_corrp2];
%         r.data.maps_bcsm_corrpset2 = [p.data.maps_bcsm_corrpset2; {q.data.maps_bcsm_corrpset2}];
%         r.data.maps_dksm_corrp2 = [p.data.maps_dksm_corrp2; q.data.maps_dksm_corrp2];
%         r.data.maps_dksm_corrpset2 = [p.data.maps_dksm_corrpset2; {q.data.maps_dksm_corrpset2}];
%         r.data.SIC_adsm_corrp2 = [p.data.SIC_adsm_corrp2; q.data.SIC_adsm_corrp2];
%         r.data.SIC_adsm_corrpset2 = [p.data.SIC_adsm_corrpset2; {q.data.SIC_adsm_corrpset2}];
%         r.data.SIC_bcsm_corrp2 = [p.data.SIC_bcsm_corrp2; q.data.SIC_bcsm_corrp2];
%         r.data.SIC_bcsm_corrpset2 = [p.data.SIC_bcsm_corrpset2; {q.data.SIC_bcsm_corrpset2}];
%         r.data.SIC_dksm_corrp2 = [p.data.SIC_dksm_corrp2; q.data.SIC_dksm_corrp2];
%         r.data.SIC_dksm_corrpset2 = [p.data.SIC_dksm_corrpset2; {q.data.SIC_dksm_corrpset2}];
%         r.data.ISE_adsm_corrp2 = [p.data.ISE_adsm_corrp2; q.data.ISE_adsm_corrp2];
%         r.data.ISE_adsm_corrpset2 = [p.data.ISE_adsm_corrpset2; {q.data.ISE_adsm_corrpset2}];
% %         r.data.ratesplaceMethod2 = [p.data.ratesplaceMethod2; q.data.ratesplaceMethod2];
% %         r.data.ratesviewMethod2 = [p.data.ratesviewMethod2; q.data.ratesviewMethod2];
%         
%         r.data.maps_adsm_corrsv2 = [p.data.maps_adsm_corrsv2; q.data.maps_adsm_corrsv2];
%         r.data.maps_adsm_corrsvset2 = [p.data.maps_adsm_corrsvset2; {q.data.maps_adsm_corrsvset2}];
%         r.data.maps_bcsm_corrsv2 = [p.data.maps_bcsm_corrsv2; q.data.maps_bcsm_corrsv2];
%         r.data.maps_bcsm_corrsvset2 = [p.data.maps_bcsm_corrsvset2; {q.data.maps_bcsm_corrsvset2}];
%         r.data.maps_dksm_corrsv2 = [p.data.maps_dksm_corrsv2; q.data.maps_dksm_corrsv2];
%         r.data.maps_dksm_corrsvset2 = [p.data.maps_dksm_corrsvset2; {q.data.maps_dksm_corrsvset2}];
%         r.data.SIC_adsm_corrsv2 = [p.data.SIC_adsm_corrsv2; q.data.SIC_adsm_corrsv2];
%         r.data.SIC_adsm_corrsvset2 = [p.data.SIC_adsm_corrsvset2; {q.data.SIC_adsm_corrsvset2}];
%         r.data.SIC_bcsm_corrsv2 = [p.data.SIC_bcsm_corrsv2; q.data.SIC_bcsm_corrsv2];
%         r.data.SIC_bcsm_corrsvset2 = [p.data.SIC_bcsm_corrsvset2; {q.data.SIC_bcsm_corrsvset2}];
%         r.data.SIC_dksm_corrsv2 = [p.data.SIC_dksm_corrsv2; q.data.SIC_dksm_corrsv2];
%         r.data.SIC_dksm_corrsvset2 = [p.data.SIC_dksm_corrsvset2; {q.data.SIC_dksm_corrsvset2}];
%         r.data.ISE_adsm_corrsv2 = [p.data.ISE_adsm_corrsv2; q.data.ISE_adsm_corrsv2];
%         r.data.ISE_adsm_corrsvset2 = [p.data.ISE_adsm_corrsvset2; {q.data.ISE_adsm_corrsvset2}];
% %         r.data.covmat2 = [p.data.covmat2; {q.data.covmat2}];
% %         r.data.covmat_norm2 = [p.data.covmat_norm2; {q.data.covmat_norm2}];
% %         r.data.l1norm2 = [p.data.l1norm2; q.data.l1norm2];
% %         r.data.l2norm2 = [p.data.l2norm2; q.data.l2norm2];
        
		
	end
end
