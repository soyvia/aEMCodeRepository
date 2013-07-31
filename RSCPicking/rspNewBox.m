function [picks, ptrs, doUpdate]=rspNewBox(mi,rscc,dis,picks,ptrs,coords,clickType) % insert a manual pick
% The doUpdate flag is set to 1 when a bad vesicle is marked, so any
% deleted particles are also marked.

% flags        pick array index and color
delFlag=1;      % color=1
vesFlag=2;      % color=5
badVesFlag=3;   % color=7;
manFlag=16;     % color=2
manRawFlag=17;  % color=6;
% autoFlag=32;  % color=3 (not used here)
bkdFlag=48;     % color=4

nsets=numel(ptrs);

dsc=mi.imageSize(1)/size(rscc.mxCC,1);  % Downsample factor of mxCC

doUpdate=0;

coords(3)=0;
if numel(coords)<8
    coords(8)=0;   % The written coords are: [x y type vesInd RSO alpha beta gamma]
end;

allInds=zeros(nsets,1);
minDistances=zeros(nsets,1);
minDistances(1)=inf;
allInds(1)=1;
for i=2:nsets  % we ignore i=1 which are deleted particles
    dists=hypot(picks(i,1:ptrs(i),1)-coords(1),...
        picks(i,1:ptrs(i),2)-coords(2));
    if numel(dists)<1
        dists=inf;
    else
        dists(picks(i,:,3)==0)=inf;  % Don't count bad vesicle coords
    end;
    [minDistances(i), allInds(i)]=min(dists);
end;

[globalMin, setInd]=min(minDistances);

switch clickType
    case {1,'.'}  % simple click, should be a manual pick, goes into picks(2,..)
        %         We'll select it if we're not at an old location
        if globalMin>=dis.minDist
            %                 || picks(globalInd+1,allInds(globalInd),3)==0
            coords(3)=manRawFlag;  % mark as manually picked
            if clickType==1
                %             Check against the CC map
                nrgn=round(0.6*dis.pars(20)/dis.ds)*2+1;  % compute from box size
                rctr=ceil((nrgn+1)/2);
                rgn=ExtractImage(rscc.mxCC,round(coords(1:2)/dsc)+1,nrgn);
                [mxv, i, j]=max2di(rgn);
                if ~(any([i j]==0) || any([i j]==nrgn) || mxv<=0)  % we have a valid maximum
                    newcoords=dsc*([i j]-rctr)+dsc*round(coords(1:2)/dsc);  % replace the position
                    %                 disp(['corrected coords: ' num2str([coords(1:2) newcoords])]);
                    coords(1:2)=newcoords;
                    coords(3)=manFlag;
                end;
            end;
            %             Insert the vesicle index
            if all(coords>=0)
                dsCoords=round(coords(1:2)/dsc)+1;
                coords(4)=single(rscc.mxVesInds(dsCoords(1),dsCoords(2)));
                coords(5)=single(rscc.mxCC(dsCoords(1),dsCoords(2)));
                coords(6)=single(rscc.mxTemplInds(dsCoords(1),dsCoords(2)));
                %                 coords(5)=single(rscc.mxRsos(dsCoords(1),dsCoords(2)));
            else
                error('Negative coordinates');
            end;
            if coords(3)==manFlag  % This has corrected coordinates
                boxColor=2;
                ptrs(2)=ptrs(2)+1;
                picks(2,ptrs(2),:)=coords;
            else
                boxColor=6;
                ptrs(6)=ptrs(6)+1;
                picks(6,ptrs(6),:)=coords;
            end;
            ShowInfo(coords,'New:  ');
            [bX, bY]=rspMakeBoxes(mi,dis,coords);
            plot(bX,bY,'color',dis.boxColors(boxColor,:));  % draw a box.
        end;
    case 2  % shift-click or center button erases
        if globalMin<dis.minDist  % We clicked on an existing marker
            %             allInds(setInd)
            partInd=allInds(setInd);
            picks(setInd,partInd,3)=0;  % Mark it deleted.
            boxLoc=picks(setInd,partInd,:);
            ShowInfo(boxLoc,'Deleted: ');
            boxLoc(3)=1;  % mark it erased
            [bX, bY]=rspMakeBoxes(mi,dis,boxLoc);
            plot(bX,bY,'color',dis.boxColors(1,:));  % draw a black box.
            if setInd==3  % special case:  it was an auto-picked particle.
                %                     We need to store the location.
                ptrs(1)=ptrs(1)+1;  % we'll make an entry in the deleted stack
                picks(1,ptrs(1),:)=picks(setInd,partInd,:);
                picks(1,ptrs(1),3)=delFlag;  % Mark the entry as deleted.
            end;
        end;
        
    case 3  % ctrl-click marks a non-particle.
        coords(3)=bkdFlag;  % mark as a background particle
        ptrs(4)=ptrs(4)+1;
        picks(4,ptrs(4),:)=coords;
        [bX, bY]=rspMakeBoxes(mi,dis,coords);
        plot(bX,bY,'color',dis.boxColors(4,:));  % draw a box.
        
    case 'l'  % mark the vesicle (Liposome) center
        vdist=inf;
        if numel(mi.vesicle.x)>0 % vesicles exist
            vdist=hypot(mi.vesicle.x(:)-coords(1),mi.vesicle.y(:)-coords(2));
        end;
        [minVes, indVes]=min(vdist);
        if minVes<dis.minDist    % We are close to a vesicle center
            ptrs(5)=ptrs(5)+1;
            boxColor=5;                    %             Insert the vesicle index
            coords(1)=mi.vesicle.x(indVes);
            coords(2)=mi.vesicle.y(indVes);
            coords(3)=vesFlag;
            coords(4)=indVes;
            picks(5,ptrs(5),:)=coords;
            [bX, bY]=rspMakeBoxes(mi,dis,coords);
            plot(bX,bY,'color',dis.boxColors(boxColor,:));  % draw a box.
        end;
    case 'x'  % mark a bad vesicle
        vdist=inf;
        if numel(mi.vesicle.x)>0 % vesicles exist
            vdist=hypot(mi.vesicle.x(:)-coords(1),mi.vesicle.y(:)-coords(2));
        end;
        [minVes, indVes]=min(vdist);
        %         [minVes indVes]
        if minVes<2*dis.minDist    % We are close to a vesicle center
            ptrs(7)=ptrs(7)+1;   % make a vesicle entry
            coords(1)=mi.vesicle.x(indVes);
            coords(2)=mi.vesicle.y(indVes);
            coords(3)=badVesFlag;
            coords(4)=indVes;
            picks(7,ptrs(7),:)=coords;
            [picks, ptrs]=rspDeleteBadVesiclePicks(picks,ptrs);
            doUpdate=1;  % We'll let the calling program update the display.
        else % draw a miniature box to show that we've responded.
            origBoxSize=dis.pars(20);
            dis.pars(20)=dis.pars(20)/2;
            [bX, bY]=rspMakeBoxes(mi,dis,coords);
            dis.pars(20)=origBoxSize;
            plot(bX,bY,'color',dis.boxColors(7,:));  % draw a box.
        end;
        
    case 'i'  % get info about a particle        if globalMin<dis.minDist  % We clicked on an existing marker
        if globalMin<dis.minDist  % We clicked on an existing marker
            partInd=allInds(setInd);  % Look up the particle info
            coords=squeeze(picks(setInd,partInd,:))';  % Mark it deleted.
            if all(coords>=0)
                dsCoords=round(coords(1:2)/dsc)+1;
                coords(4)=single(rscc.mxVesInds(dsCoords(1),dsCoords(2)));
                coords(5)=single(rscc.mxCC(dsCoords(1),dsCoords(2)));
                coords(6)=single(rscc.mxTemplInds(dsCoords(1),dsCoords(2)));
            end;
            ShowInfo(coords,'Info: ');
        end;
end;  % switch

function ShowInfo(coords,str)
if dis.listParticleInfo  % if we've turned this on
    rcoords=round(coords);
    fprintf(1,[str 'x= %d  y= %d flag= %d ves= %d cc=%5.2f ref= %d\n'],...
        rcoords(1),rcoords(2),rcoords(3),rcoords(4),coords(5),rcoords(6));
    c1=round(coords(1:2)/dsc+1);
    disp(['        Amplitude ' num2str(rscc.mxCC(c1(1),c1(2)))]);
    disp(['        Variance  ' num2str(rscc.mxVars(c1(1),c1(2)))]);
    vind=coords(4);
    if vind>0  % vesicle number
        vesr=mi.vesicle.r(vind);
        vesc=[mi.vesicle.x(vind) mi.vesicle.y(vind)];
        r=sqrt(sum(([coords(1) coords(2)]-vesc).^2)); % dist to ves ctr
        dists=sqrt((coords(1)-mi.vesicle.x(:)).^2+(coords(2)-mi.vesicle.y(:)).^2);
        dists(vind)=inf;
        [minv, vind2]=min(dists);
        disp(['        Bobbing dist ' num2str(round((r-vesr)*mi.pixA))]);
        disp(['        Next vesicle dist '...
            num2str(round((minv-mi.vesicle.r(vind2))*mi.pixA))]);
    end;
end;

end % showInfo
end % main function