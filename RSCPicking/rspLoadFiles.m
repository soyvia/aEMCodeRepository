function [dis, mi, rscc, rs, imgs, masks, m0, mVes, rawMask]=rspLoadFiles(dis)
% Given the info file name and path, load the mi, merged image and rscc
% files.  If a vesicle file is also present load it, otherwise compute the
% vesicles.

dis.infoPath=AddSlash(dis.infoPath);
dis.basePath=AddSlash(pwd);
disp('Loading the files:');
[~, ~, ex]=fileparts(dis.infoName);
% Possibly a jpeg file...
if strcmp(ex,'.jpg')  % we are picking a simple image
    mfileName=[dis.infoPath dis.infoName];
    disp(mfileName);
    mc=ReadEMFile(mfileName); % get the jpeg image
    
    %      Make a fake info structure
    dis.miValid=0;  % This mi structure shouldn't be stored
    mi=meCreateMicrographInfoStruct13; % don't create mi.mask
    [mi.procPath, mi.basePath]=ParsePath(dis.infoPath);
    len=numel(dis.infoName);
    mi.baseFilename=dis.infoName(1:len-5);
    mi.pixA=mi.defaultPixA;
    mi.imageSize=size(mc);
    mi.kV=200;
    dis.infoName='';
else  % normal case
    disp([dis.infoPath dis.infoName]);
    mi=load([dis.infoPath dis.infoName]);  % get the mi structure
    mi=mi.mi;
    mi.basePath=dis.basePath;
    set(gcf,'name',[dis.infoName ' - loading...']);
    dis.miValid=1;
    drawnow;
    mc=meReadMergedImage(mi);
    if numel(mc)<2  % unsuccessful read
        mc=single(zeros(1024));
        disp('No merged image found.');
    end;
end;
rsccName=[mi.procPath mi.baseFilename 'rscc.mat'];
disp(rsccName);
if exist(rsccName,'file')
    rscc=load(rsccName);
else
    disp('...not found');
    rscc.mxVars=single(zeros(mi.imageSize/2));
    rscc.mxCC=rscc.mxVars;
    rscc.mxVesInds=rscc.mxVars;
end;
dis.ds=mi.imageSize(1)/dis.ndis;
if isfield(mi.particle,'autopickPars') && numel(mi.particle.autopickPars)>2 ...
    && mi.particle.autopickPars(1)>0 % Has valid parameters
    npars=numel(mi.particle.autopickPars);
    dis.pars(1:npars)=mi.particle.autopickPars;
end;
save(dis.datName,'dis');  % store the preliminary settings for next time.

% % Pick up the list of names in the info directory
% dis.infoNames=FindFilenames(dis.infoPath,'.+mi\.mat');
% if numel(dis.infoNames)<1
%     error(['No mi.mat files found in ' dis.infoPath]);
% end;

% Set up parameters that depend on the image
dis.ds=mi.imageSize(1)/dis.ndis;
% pixA=mi.pixA*dis.ds;

disp('Downsample');
rs.mVars=DownsampleGeneral(rscc.mxVars,dis.ndis);
rs.blanks=DownsampleGeneral(rscc.mxVars==0,dis.ndis);
rs.mCC=DownsampleGeneral(rscc.mxCC,dis.ndis);
rs.mVesInds=DownsampleGeneral(rscc.mxVesInds<1,dis.ndis);
m0=DownsampleGeneral(mc,dis.ndis);
if numel(mi.vesicle.ok)<1
    disp('No vesicles found.');
    mVes=single(zeros(dis.ndis));
    mVesGood=mVes;
    mVesBad=mVes;
end;
% Get the vesicle model mVes
disp('Make Vesicles');
if size(mi.vesicle.ok,2)<4  % handle old mi files; all ves are good.
    mi.vesicle.ok=repmat(mi.vesicle.ok,1,4);
else
    goodVes=all(mi.vesicle.ok(:,2:3),2);  % vesicle in range and refined
badVes=mi.vesicle.ok(:,1) & ~goodVes; % refined but not in range
    mi2=mi;
    mi2.vesicle.s=0*mi2.vesicle.s+median(mi2.vesicle.s);
    mVesGood=meMakeModelVesicles(mi2,dis.ndis,find(goodVes));
    mVesBad=meMakeModelVesicles(mi2,dis.ndis,find(badVes));
    mVes=meMakeModelVesicles(mi,dis.ndis,find(mi.vesicle.ok(:,1)));  % Actual vesicle model
    % %%%%% experimental nonlinear code.
    % th1=1.5;
    % wi1=1;
    % q=mVes>th1;
    % mVes(q)=mVes(q).*exp(-(mVes(q)-th1).^2/(2*wi1^2));
    
    
    %     if dis.miValid  % let's store the model vesicles
    %         if ~isfield(mi,'tempPath')
    %             mi.tempPath='Temp/';
    %         end;
    %         if ~exist([dis.basePath 'Temp'],'dir') || ~isfield(mi,'tempPath')
    %             mkdir([dis.basePath 'Temp/']);
    %             mi.tempPath='Temp/';
    %         end;
    %         outVesName=[dis.basePath mi.tempPath mi.baseFilename 'v.mrc'];
    %         WriteMRC(mVes,mi.pixA*mi.imageSize(1)/size(mc,1),outVesName);
    %         disp([outVesName ' saved']);
    %     end;
end;
rawMask=~meGetMask(mi,dis.ndis);

disp('FilterAndScale');
imgs=zeros(dis.ndis,dis.ndis,7,'single');  % Images for display
imgs(:,:,1:2)=rspFilterAndScaleImages(mi,dis,m0,mVes);
imgs(:,:,3)=0;  % to receive model picked particles
imgs(:,:,4)=imscale(rs.mCC,256,1e-3);
imgs(:,:,5)=imscale(rs.mVars,256,1e-2);
imgs(:,:,6)=imscale(mVes,256,1e-4);
imgs(:,:,7)=0;  % to receive overlap mask
% masks are:  1: good vesicle ghosts;  2: bad vesicle ghosts;
%             3: variance above thresh 4: blank region
%             5: vesicle overlap
masks=zeros(dis.ndis,dis.ndis,5,'single');
masks(:,:,1)=imscale(max(-mVesGood,0),1,1e-3);
if any(badVes)
    masks(:,:,2)=imscale(max(-mVesBad,0),1,1e-3);
end;
masks(:,:,3)=((rs.mVars>dis.pars(3))|rawMask);
masks(:,:,4)=rs.blanks;

dis.org=[0 0];
dis.mode=min(2,dis.mode);
disp('Done.');