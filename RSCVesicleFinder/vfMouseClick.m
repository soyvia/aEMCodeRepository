function [h, doUpdateDisplay]=vfMouseClick(h)
% Called from a mouse click on the display, to change add or remove a
% vesicle.  What happens depends on whether the click is near the center of
% an existing vesicle.
% Right click: change to or add a good vesicle
% Left (ctrl) click: change to or add a bad vesicle
% Center (shift) click: delete the vesicle
% the function updates these variables
%   h.markedVesicleIndex (=0 if nothing changed)
%   h.goodVesImage, badVesImage
%   h.mi.vesicle (element markedVesicleIndex, if changed).


nrgn=25;  % size of box, in display pixels, allowed for clicking on
          % existing vesicle
p=get(h.axes1,'CurrentPoint');
n=size(h.rawImage);
p1=([p(1,1) n(2)-p(1,2)]-1)*h.ds0;
b=get(gcf,'SelectionType'); % normal, alt, extend are the 3 buttons.
h.changedVesicleIndex=0;   % default
doUpdateDisplay=false;     % default
badVes=h.mi.vesicle.ok(:,1) & ~h.mi.vesicle.ok(:,2);

% Assign the nature of the desired new vesicle newVes
switch b
    case 'normal'  % replace or add a new 'good' vesicle
        newVesFlag=1;
    case 'alt'     % replace or add a new 'bad' vesicle
        newVesFlag=2;
    case 'extend'  % Delete a vesicle
        newVesFlag=0;
end;

% find the index vind of the existing vesicle oldVes, and set the value of
% the newVesFlag.
% 
vesToSearch=h.mi.vesicle.ok(:,1);
ys=h.mi.vesicle.y;
ys(~vesToSearch)=inf;  % exclude non-vesicles from search
dists=sqrt((h.mi.vesicle.x-p1(1)).^2 + (ys-p1(2)).^2);
oldVesFlag=0;  % default is, no vesicle found
vind=0;

if numel(dists)>0
    [mnv, tempInd]=min(dists);
disp([mnv tempInd]);
    if mnv<nrgn*h.ds0/2  % we found something within the nrgn box.
        vind=tempInd;
        oldVesFlag=1+badVes(vind);  % 1: good vesicle; 2: bad vesicle.
    end;
end;

% Now consider the new vesicle
if newVesFlag==oldVesFlag  % nothing to do
    h.markedVesicleIndex=vind;  % mark the vesicle if found
    return

elseif oldVesFlag==0  % We're creating a new one
    if max(h.ccValsScaled(:))==0  % No cross correlation available
        h.markedVesicleIndex=0;  % no new vesicle
        return
    end;
    rctr=ceil((nrgn+1)/2);
    rgn=ExtractImage(h.ccValsScaled,round(p1/h.ds0)+1,nrgn);
    [mxv, i, j]=max2di(rgn);
    disp([mxv i j]);
    if (all([i j]>1) && all([i j]<nrgn) && mxv>0)  % we have a valid maximum
% disp('Valid maximum');
        newcoords=h.ds0*([i j]-rctr)+h.ds0*round(p1/h.ds0);  % replace the position
        vind=numel(h.mi.vesicle.x)+1;  % default: add a new vesicle
        h.mi.vesicle.x(vind)=newcoords(1);
        h.mi.vesicle.y(vind)=newcoords(2);
        h.mi.vesicle.s(vind)=mxv;
        p2=round([i j]+p1/h.ds0-rctr+1);
        p2=max(1,min(p2,n));
        h.mi.vesicle.r(vind)=h.ccRadii(p2(1),p2(2));
        h.mi.vesicle.ok(vind,1)=true;  % extend the array
        h.mi.vesicle.ok(vind,2)=(newVesFlag==1);
% disp('Adding a vesicle: ');
% disp(p2);
% disp(h.mi.vesicle.s(vind))
% disp(h.mi.vesicle.r(vind))
        v0=meMakeModelVesicles(h.mi,n,vind,0);
        v=real(ifftn(fftn(v0).*ifftshift(h.ctf))); % filter by ctf
    else
        v=0;
    end;
else  % an old vesicle exists, remove it
    v0=meMakeModelVesicles(h.mi,n,vind,0);
    v=real(ifftn(fftn(v0).*ifftshift(h.ctf))); % filter by ctf

    h.mi.vesicle.ok(vind,1:2)=false;  % mark it empty
    if oldVesFlag==1 % remove from the good image
        h.goodVesImage=h.goodVesImage-v;
    elseif oldVesFlag==2
        h.badVesImage=h.badVesImage-v;
    end;
end;

if newVesFlag>0 && vind>0 % We're to add something
    if newVesFlag==1  % add to good ves
        h.goodVesImage=h.goodVesImage+v;
        h.mi.vesicle.ok(vind,1:2)=true;
    else
        h.badVesImage=h.badVesImage+v;
        h.mi.vesicle.ok(vind,1:2)=[true false];
    end;        
end;

h.markedVesicleIndex=vind; 
doUpdateDisplay=true;
disp('mc'); 
