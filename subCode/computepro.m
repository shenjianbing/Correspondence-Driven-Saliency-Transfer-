 function [bp_R,fp_R]= computepro(BSsaliency,UMsaliency_voting,bias_mask,superpixels)
    back_mask = double(UMsaliency_voting<4);
    back_mask(bias_mask)=0;
    BSsaliency_f = BSsaliency;
    BSsaliency_f(~bias_mask) =0;
    R = getSuperpixelMeanValues(BSsaliency_f,superpixels.Label,double(max(superpixels.Label(:))));     
    fp_R = normalize(R);   
    fp = fp_R(superpixels.Label);
    
    R = getSuperpixelMeanValues(double(back_mask),superpixels.Label,double(max(superpixels.Label(:))));       
    bp_R = normalize(R);   
    
    bdIds=[superpixels.Label(1,:)';superpixels.Label(end,:)';superpixels.Label(:,1);superpixels.Label(:,end)];
    bdIds = unique(bdIds);
    bdIds_bias = bdIds;
    bp_R(bdIds_bias) = max(bp_R(bdIds_bias),1-fp(bdIds_bias));
    fp_R(bdIds_bias) = 1-bp_R(bdIds_bias);
