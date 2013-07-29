function [picks, ptrs, mi]=rspLoadPicksFromMi(mi)
% Load the arrays from the mi file
% particle types are the following
% 0-15 not displayed
%  0 = null
%  1 = deleted from auto
% 16-31 manual picking (we use only 16 at present)
% 32-47 auto picking (we use only 32 at present)
% 48-63 background (we use only 48 at present)

% Make a lookup table to map particle types to entries in the local picks
% array.
maxLut=63;
lut=zeros(1,maxLut);
lut(1)=1;
lut(16)=2;
lut(17)=6;
lut(32)=3;
lut(48)=4;
lut(2)=5;  % vesicle marker
lut(3)=7;  % bad vesicle marker
nflags=max(lut);
% The new field mi.particle.picks is an np x 8 array.  Each row is
% [x y type vesicleIndex RSO alpha beta gamma]

np=size(mi.particle.picks,1);  % number of total entries
% np1=max(np,100);
np1=np;
picks=single(zeros(nflags,np1,8));
ptrs=zeros(1,nflags);
for i=1:np
    c=mi.particle.picks(i,:);
    ncf=numel(c);  % Number of c fields
    q=round(c(3));  % particle type
    if q>0 && q<=maxLut  % type is in range
        ind=lut(q);
        if ind>0
            ptrs(ind)=ptrs(ind)+1;
            picks(ind,ptrs(ind),1:ncf)=c;
        else
            disp(['rspLoadPicks: ignored particle type: ' num2str(c(3))]);
        end;
    elseif q>0
        disp(['rspLoadPicks: invalid particle type: ' num2str(c(3))]);
    end;
end;
