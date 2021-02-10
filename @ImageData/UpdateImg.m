function out = UpdateImg(this, img)
%% update ImageData object with new img data while keeping other parameters and masks
%
% -------------------------------------------------------------------------------------------------------------------------

    % check dimension
    if any(size(img)~=size(this.img)), error('%s: The sizes of new and old images do not match',mfilename); end
    
    % update image data
    out = this ; 
    out.img = img ; 
    
end