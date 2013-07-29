function [coords b]=rspGetClick(dis)
% Get a click or keypress and return the coordinates relative to the
% original micrograph, zero-based.

[x y b]=Myginput(1,'square');
% [x y b]=ginput(1);
%             rawCoords=[x y]

coords=zeros(1,3);
if numel(x)<1 || numel(b)<1
    b=0;
    return
end;
coords(1)=(x+dis.org(1))*dis.ds;
coords(2)=(dis.size(2)-y+dis.org(2)-1)*dis.ds;
% coords