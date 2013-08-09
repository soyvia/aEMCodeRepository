% SimpleRSPicker.m
% F. Sigworth, Jan '13
%  See rspLoadPicksFromMi.m to see assignment of ptrs
%  See rspNewBox to see assignment of flag values

version=130;

% Retrieve parameters from a file in the program directory
pa=fileparts(which('SimpleRSPicker'));
datName=[AddSlash(pa) 'SimpleRSPickerDat.mat'];
disOk=0;
if exist(datName,'file')>0
    dis=load(datName);
    dis=dis.dis;
    dis.datName=datName;
    if isfield(dis,'version')
        disOk=dis.version==version && dis.miValid;
        disp('Loading settings');
    end;
end;

newDis=0;
if ~disOk
    disp('Loading defaults');
    dis=struct;
    dis.version=version;
    dis.datName=datName;
    % Box colors
    dis.infoPath='./';
    bColorErased=[.3 .3 .3];
    bColorMan=[1 .9 0];  % yellow
    bColorManNotRefined=[1 .75 .1]; %orangish
    bColorAuto=[.8 1 0];  % greenish
    bColorBkd=[.4 1 .8]; % cyan
    bColorVesicle=[.2 0 1]; % blue vesicle
    bColorBadVesicle=[1 .1 .3]; % red bad vesicle
    dis.boxColors=[bColorErased; bColorMan; bColorAuto; bColorBkd;
        bColorVesicle; bColorManNotRefined; bColorBadVesicle];
    
    % overlay colors
    dis.maskColor=[.9 .75 .8];
    dis.blankColor=[1 .85 .8];
    dis.overlapColor=[.95 .8 1];
    dis.ghostColor=[.7 .7 1];
    dis.ghostColorBad=[1 .7 .6];
    dis.ghostAmp=.6;
    
    % initial display
    dis.ndis=1024; %%
    dis.size=[1200 800];
    dis.size=min(dis.size,dis.ndis); % Force the display to be no bigger than the image
    dis.org=[0 0]; %pixels run from org+1 to org+size
    
    %     default parameters
    dis.mode=2;  % vesicles subtracted
    dis.showMask=1;
    dis.showGhosts=1;
    dis.showBoxes=1;
    dis.listParticleInfo=1;
    dis.contrast=[3e-3 3e-3];
    dis.varThresh=40;
    dis.pars=[.55 1 40 35 100 120 120 200];  % min, max amp; max var; rso offset;
%     particle blank radius, vesicle blank radius, maxBob, border.
    dis.pars(20)=32;  % box size.
    dis.minDist=dis.pars(20);  % distance in original pixels, based on box size
    dis.filter=[1000 20];  % inverse frequency in A
    dis.infoName='';
    dis.defaultPixA=2.9;
    
    dis.finished=0;
    dis.miValid=0;
    newDis=1;
end;
%%
% Try to load the previous mi file

if exist(dis.infoPath,'dir')
    cd(dis.infoPath);
end;

% Set the first command
if (dis.finished || newDis || ~dis.miValid) % We need to open a new file
    b='o';  % ask for a new file
else
    b='v';  % 'revert' to latest file.
end;
oldB=0;

% disp('Make figure');
screenSize=get(0,'screensize');
xsiz=min(screenSize(3),dis.size(1)+3);
ysiz=min(screenSize(4)-50,dis.size(2)+3);
figure(1);
clf;
% Put the window near the top middle of the screen
set(gcf,'position',[(screenSize(3)-dis.size(1))/2 ((screenSize(4)-50)-ysiz)*0.9 xsiz ysiz],'toolbar','none','resize','off');
% Main display
axsiz=xsiz-3;
aysiz=ysiz-3;
ax1=axes('units','pixels','position',[2 3 axsiz aysiz]); %,'ticklength',[0 0]);
ax2=axes('position',[.8 0 .2 .2]);
ax3=axes('outerposition',[.8 0 .2 .2]);
axes(ax1);

refreshReconstruct=0;  % flag to update the reconstruction display
%%
% % interactive loop
while (b~='q') && (b~='Q') % q = quit; Q = quit but return to this image later
    switch b

        case {1 2 3 '.' 'i' 'l' 'x'}  % simple click or erase, info, vesicle
            [picks, ptrs, doUpdate]=rspNewBox(mi,rscc,dis,picks,ptrs,coords,b); % insert a manual pick
            refreshReconstruct=1;
            if doUpdate
                rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            end;
            if b==1 || b==2 || b=='.' || b=='x' % changed a particle
                refreshReconstruct=1;
            end;
        case 'a'    % autopicking
            [qx, qy, b]=Myginput(1,'circle'); % Put up a circle while we enter data
            ok=0;
            switch b
                case 'a'  % aa: Go ahead and do the auto-picking
                    disp(['Auto-picker parameters: ' num2str(dis.pars)]);
                    ok=1;
                case 'p'  % ap: change parameters, then auto-pick
                    disp('Auto-picker parameters:');
                    dis.pars(1)=MyInput(' Min amplitude  ',dis.pars(1));
                    dis.pars(2)=MyInput(' Max amplitude  ',dis.pars(2));
                    dis.pars(3)=MyInput(' Max variance   ',dis.pars(3));
                    dis.pars(4)=MyInput(' RSO offset, A  ',dis.pars(4));
                    dis.pars(5)=MyInput(' Blank radius, A',dis.pars(5));
                    dis.pars(6)=MyInput(' Overlap blanking radius, A',dis.pars(6));
                    dis.pars(7)=MyInput(' Max bob, A',dis.pars(7));
                    dis.pars(8)=MyInput(' Border, A', dis.pars(8));
                    ok=1;
                    
                    % optional online parameter entry
                case 'v'  % av: set the variance threshold then go
                    [dis.pars(3) ok]=rspGetGValue(dis.pars(3));
                    disp(['Max variance = ' num2str(dis.pars(3))]);
                otherwise
                    beep;
            end;
            if ok
                if max(rscc.mxCC(:))==0  % No cross correlation data
                    disp('No picking for lack of rscc preprocessor data.');
                else
                    disp('Auto finding particles...');
                    % Show status in the window title bar
                    set(gcf,'name',[dis.infoName ' - autopicking...']);
                    drawnow;
                    % ----- Do the autopicking here -----------------------
                    [coords, ovMask, endCC]=rscAutoParticleFinder(mi,rscc,dis.pars);
                    imgs(:,:,7)=150*(ovMask-.5);
                    masks(:,:,5)=ovMask>1.5;
                    imgs(:,:,8)=imscale(endCC,256);
                    [ptrs(3), ncf]=size(coords);
                    picks(3,1:ptrs(3),1:ncf)=coords;
                    [picks, ptrs]=rspDeleteBadAutoPicks(dis,picks,ptrs);
                    [picks, ptrs]=rspDeleteBadVesiclePicks(picks,ptrs);
                end;
%                 update the found particles display
                refreshReconstruct=1;
                % update the mask display
                masks(:,:,3)=((rs.mVars>dis.pars(3))|rawMask);
                rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
                disp([num2str(ptrs(3)) ' particles found.']);
                set(gcf,'name',dis.infoName);
            end;
            
        case 'b'  % toggle box display
            dis.showBoxes=~dis.showBoxes;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            if dis.showBoxes
                disp('Box display on');
            else
                disp('Box display off');
            end;
        case 'B'
            dis.pars(20)=MyInput(' Box size ',dis.pars(20));
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
        case 'c'  % set contrast
            disp('Setting contrast');
            dis.contrast(1)=MyInput('Black contrast',dis.contrast(1));
            dis.contrast(2)=MyInput('White contrast',dis.contrast(2));
            imgs(:,:,1:2)=rspFilterAndScaleImages(mi,dis,m0,mVes);
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'C'  % set overlay colors
             dis.maskColor=MyInput('Mask color',dis.maskColor);
             dis.blankColor=MyInput('Blank color',dis.blankColor);
             dis.overlapColor=MyInput('Overlap color',dis.overlapColor);
             dis.ghostColor=MyInput('Ghost color',dis.ghostColor);
             dis.ghostColorBad=MyInput('Bad ghost color',dis.ghostColorBad); 
             rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case {'d' 'e' 'D'}  % display mode (toggle 1-2, 2-3, 1-end)
            dis.mode=dis.mode+1;
            if b=='d' && dis.mode>2 % lower case toggles norm / ves-subtr
                dis.mode=1;
            elseif b=='e' && dis.mode>3
                dis.mode=2;
            end;
            if dis.mode > size(imgs,3) % upper case cycle through all imgs
                dis.mode=1;
            end;
            if dis.mode==3 && refreshReconstruct
                imgs(:,:,3)=rspReconstructParticles(dis,mi,picks,ptrs,rscc);
                refreshReconstruct=0;
            end;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'f'  % set filtering
            disp('Setting the filter:');
            dis.filter(1)=MyInput('Highpass filter, A',dis.filter(1));
            dis.filter(2)=MyInput('Lowpass filter, A',dis.filter(2));
            imgs(:,:,1:2)=rspFilterAndScaleImages(mi,dis,m0,mVes);
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'g'  % toggle "ghost" vesicles
            dis.showGhosts=~dis.showGhosts;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
        case 'G'  % set ghost amplitude
            disp('Setting ghost amplitude');
            dis.ghostAmp=MyInput('Ghost amplitude',dis.ghostAmp);
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
        case 'h'  % Show autopicker histogram
            if ptrs(3)>0
                    axes(ax2);
                    cla;
                    axes(ax3);
                    midAmp=median(mi.vesicle.s(~isnan(mi.vesicle.s)));
                    amps=[];
                    for i=1:ptrs(3)
                        c=squeeze(picks(3,i,:));
                        if c(3)>15 && c(4)>0 % a pick
%                             amps=[amps; c(5)*midAmp/mi.vesicle.s(c(4))];
                            amps=[amps; c(5)];
                        end;
                    end;
%                     set(ax2,'color','w');
                    hist(amps,30);
                   drawnow; 
            end;

        case 'm'  % toggle mask display
            dis.showMask=~dis.showMask;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case {'n' 'N' 'o' 'v'}  % open a file
            %             First, store the present results
            if b~='v'
                if numel(dis.infoName)>3 && dis.miValid && exist('mi','var')
                    mi=rspStorePicksInMi(mi,picks,ptrs);
                    save([dis.infoPath dis.infoName],'mi');  % store the mi structure
                    disp([dis.infoName ' written.']);
                    disp([' with ' num2str(sum(rspCountParticles(picks))) ' particles']);
                    
                end;
            end;
            save(dis.datName,'dis');
            disp('Loading mi file');
            if b=='o'  % get a new file
                [dis.infoName,dis.infoPath] = uigetfile({'*.mat';'*.jpg'},'Load Info File');
                if ~isa(dis.infoName,'char')
                    disp('No file selected');
                    return
                end;
                dis.infoPath=AddSlash(dis.infoPath);
                dis.basePath=ParsePath(dis.infoPath);
                cd(dis.basePath);
                
            else % open next or previous or current
                [pa, nm, ex]=fileparts(dis.infoName);
                fileTypes=strcmp(ex,{'.mat'; '.jpg'});
                if fileTypes(1)
                    names=FindFilenames(dis.infoPath,'.+mi\.mat');
                else
                    names=FindFilenames(dis.infoPath,'.+\.jpg');
                end;
                currentIndex=find(strcmp(dis.infoName, names));
                if numel(currentIndex)<1
                    disp(['can''t find the file ' dis.infoName]);
                    break;
                end;
                if b=='n'
                    disp('Get next file');
                    ind=currentIndex+1;
                    if ind>numel(names)
                        beep;
                    else
                        dis.infoName=names{ind};
                        
                    end;
                    
                elseif b=='N'
                    disp('Get previous file');
                    ind=currentIndex-1;
                    if ind<1
                        beep;
                        
                    else
                        dis.infoName=names{ind};
                    end;
                    
                else % 'v': reload the current file.
                    cd(dis.basePath);
                end;
            end;
            if exist([dis.infoPath dis.infoName],'file')  % a valid file
                [dis, mi, rscc, rs, imgs, masks, m0, mVes, rawMask]=rspLoadFiles(dis);
                mi.basePath=AddSlash(pwd);
                % Load the previous picks from the mi file.
                mi=rspUpdateMiStructure(mi);  % Change the mi.particle fields
                [picks, ptrs]=rspLoadPicksFromMi(mi);
                [picks, ptrs]=rspDeleteBadAutoPicks(dis,picks,ptrs);
                refreshReconstruct=1;
                disp([num2str(rspCountParticles(picks)) ' particles loaded']);
                rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
                set(gcf,'name',mi.baseFilename);
            else
                dis.miValid=0;
            end;
        case 's'  % save mi file
            if numel(dis.infoName)>3
                disp('Saving mi file');
                mi=rspStorePicksInMi(mi,picks,ptrs);
                save([dis.infoPath dis.infoName],'mi');  % store the mi structure
                disp(dis.infoName);
                save(dis.datName,'dis');
            end;
            
        case 'r'  % shift right
            disp('shift right');
            if dis.org(1)>=dis.ndis-dis.size(1)
                dis.org(1)=0;  % wrap around
            else
                dis.org(1)=round(min(dis.org(1)+dis.size(1)/4, ...
                    dis.ndis-dis.size(1)));
            end;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
        case 'u'  % shift up
            disp('shift up');
            if dis.org(2)>=dis.ndis-dis.size(2)
                dis.org(2)=0;  % wrap around
            else
                dis.org(2)=round(min(dis.org(2)+dis.size(2)/3, ...
                    dis.ndis-dis.size(2)));
            end;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);

        case 'Z'  % zero out all the boxes
            disp('Zeroing out all boxes.');
            ptrs=ptrs*0;
            picks=picks*0;
            refreshReconstruct=1;
            rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);
            
    end;  %switch
    oldB=b;  % store the previous key
    [coords b]=rspGetClick(dis);
    axes(ax1);
end;  % while b~='q'
%%
hold off;


mi=rspStorePicksInMi(mi,picks,ptrs);
mi.particle.autopickPars=dis.pars;
if numel(dis.infoName)>3 && dis.miValid
    save([dis.infoPath dis.infoName],'mi');  % store the mi structure
    disp([dis.infoName ' written']);
    count=rspCountParticles(picks);
    disp([' with ' num2str(count) ' particles.']);
end;

dis.finished=(b=='q');  % lower-case q causes final exit.
set(gcf,'name','Done');
if dis.finished
    dis.miValid=0;
    disp('Exiting');
else
    disp('Exiting, to continue this micrograph.');
end;
save(datName,'dis');
