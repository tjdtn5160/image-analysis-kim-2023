function [POI,badtraces]=GetPointsOfInterest_6(signal,traces,stats)
switch signal
    case 'Cdk2_R'
        %[POI,badtraces]=getCdk2features_1(traces,stats,minlength);
        %[POI,badtraces]=getCdk2features_steve(traces,stats,minlength);
        %[POI,badtraces]=getCdk2features_2(traces,stats,minlength);
        [POI,badtraces]=getCdk2features_R_serumrelease(traces,stats); %to call R even if CDK2 is added afterwards (doesn't work that well)
    case 'Cdk2'
         [POI,badtraces]=getCdk2features_cycling(traces,stats);
    case 'Cdk2_G1'
        [POI,badtraces]=getCdk2features_Sstart_1(traces,stats);
    case 'Cdk2_SG2'
        [POI,badtraces]=getCdk2features_SG2_1(traces,stats);
    case 'Geminin'
        %[POI,badtraces]=getGemininfeatures_1(traces,stats);
        %[POI,badtraces]=getGemininfeatures_2(traces,stats); %norm by 1st frame of daughter
        [POI,badtraces]=getGemininfeatures_3(traces,stats); %norm by total max
    case 'Geminin_no_norm'
        [POI,badtraces]=getGemininfeatures_3_noNorm(traces,stats);
    case 'Geminin_no_norm_reversible'
        [POI,badtraces]=getGemininfeatures_4_noNorm_reversible(traces,stats);
    case 'Geminin_G1S'
        [POI,badtraces]=getGemininfeatures_G1S(traces,stats); %norm by total max    
    case 'Geminin_SG2'
        [POI,badtraces]=getGemininfeatures_SG2_1(traces,stats); %norm by total max
    case 'Geminin_SR'
%         [POI,badtraces]=getGemininfeatures_serumrelease_chad(traces,stats);
        [POI,badtraces]=getGemininfeatures_serumrelease(traces,stats);
    case 'Mitosis_SR'
        [POI,badtraces]=getMfeatures_serumrelease(traces,stats);
    case 'PipdegFall'
        %[POI,badtraces]=getPipdegFall_NormByFirst(traces,stats,minlength);
        %[POI,badtraces]=getPipdegFall_NormByMother(traces,stats,minlength);
        %[POI,badtraces]=getPipdegFall_Dropstart(traces,stats,minlength);
        [POI,badtraces]=getPipdegFall_Dropstart_1(traces,stats);
    case 'PipdegRise'
        %[POI,badtraces]=getPipdegRise(traces,stats);
        [POI,badtraces]=getPipdegRise_1(traces,stats);
end