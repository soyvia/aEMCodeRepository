function [picks ptrs]=rspDeleteBadVesiclePicks(picks,ptrs)
% scan each auto-picked particle and delete it if it's associated with a bad
% vesicle.


if ptrs(3)<1 || ptrs(7)<1  % no auto-picked particles or vesicles to scan
    return
end;

vesPicks=picks(7,1:ptrs(7),4);  % vesicle indices
vesBad=vesPicks(picks(7,1:ptrs(7),3)>0);  % type == 0 is deleted bad vesicle.

for i=1:ptrs(3)  % scan all picked particles
    vesInd=picks(3,i,4);
    if any(vesBad==vesInd)
        picks(3,i,3)=0;  % a non-particle.
    end;
end;
