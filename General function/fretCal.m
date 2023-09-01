function sig = fretCal(cellid,numgood,img,thresh,numElement,cellMask_info)
if isempty(cellid)
    cellid=1:numgood;
end
sig=ones(numgood,1)*NaN;
for n=1:numgood
    cc=cellid(n);
    
    if cc>numel(cellMask_info)
        break;
    end
    
    sig_all=img(cellMask_info(cc).PixelIdxList);
    sig_all(sig_all>prctile(sig_all,98))=[];
    sigForeground=sig_all(sig_all>thresh);
    if numel(sigForeground)>numElement
        sig(n)=nansum(sigForeground);
    else
        sig(n)=nan;
    end
end
end

