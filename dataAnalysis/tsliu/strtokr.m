%find the last token after delimiter in string
% by Taosheng Liu
% usage: y=strtokr(pwd, filesep): return the last part of a whole path
function tok=strtokr(str,delim)

idx=strfind(str,delim);
tok=str(idx(end)+1:end);