function [boxesX, boxesY]=rspMakeBoxes(mi,dis,boxCoords)
% dis is the rsp display structure
% boxRadius is half-size of boxes in original pixel size.
% boxCoords is an n x 3 array (x,y,ok) in original-micrograph coordinates.
% If ok<=0 the box at that location isn't plotted.
% Creates 6n element arrays of coordinates which create boxes e.g. as
% plot(boxesX,boxesY,'color',[.8 .8 0]);

rFrac=0.8;  % boxes lying farther than rFrac x radius from the center of a
% vesicle get the radial vector shown.

if ~ismatrix(boxCoords)  % handle the case that we've copied from a 3D array
    boxCoords=shiftdim(boxCoords,1);
end;
n=mi.imageSize/dis.ds;  % size of the output image before adjusting coordinates
[numB ne]=size(boxCoords); % Number of boxes to construct, no. of columns
if numB<1
    boxesX=NaN;
    boxesY=NaN;
    return;
end;
if ne<4
    boxCoords(1,4)=0;  % expand the array with zeros.
end;

bx=boxCoords(:,1)/dis.ds-dis.org(1)+1;
by=dis.size(2)-boxCoords(:,2)/dis.ds+dis.org(2);
% rb=boxRadius/dis.ds;  % 'radius' of a box
rb=dis.pars(20)/dis.ds;  % 'radius' of a box
sqx0=[-rb rb rb -rb -rb NaN 0 0 NaN]';
sqy0=[-rb -rb rb rb -rb NaN 0 0 NaN]';
nL=numel(sqx0);

boxesX=single(ones(nL,numB))*NaN;
boxesY=single(ones(nL,numB))*NaN;


for i=1:numB
    if boxCoords(i,3)>0
        boxesX(:,i)=bx(i)+sqx0;
        boxesY(:,i)=by(i)+sqy0;
        j=boxCoords(i,4);
        if j>numel(mi.vesicle.x)
            warning('Inconsistent vesicles and picks, can''t draw boxes');
            return;
        end;
        if j>0 % a vesicle is marked, draw the vector toward the center
            % Get the full-length vector
            vx=mi.vesicle.x(j)/dis.ds-dis.org(1)+1-bx(i);
            vy=dis.size(2)-mi.vesicle.y(j)/dis.ds+dis.org(2)-by(i);
            dist=hypot(vx,vy);
            if dist>rFrac*mi.vesicle.r(j)/dis.ds % if the vector is long enough, draw it.
                rbe=rb/max(abs(vx),abs(vy));
                s1=rbe;
                s2=1;  % draw the vector all the way to the center
                boxesX(7:8,i)=[s1 s2]'*vx+bx(i);
                boxesY(7:8,i)=[s1 s2]'*vy+by(i);
            end;
        end;
    end;
end;
boxesX=boxesX(:);
boxesY=boxesY(:);
