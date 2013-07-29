function rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs)
% Draw the initial image
cimg=rspMakeRGB(dis,imgs(:,:,dis.mode),masks);
hold off;
image(cimg);
hold on;
if dis.showBoxes
    for i=2:numel(ptrs)
       [bX, bY]=rspMakeBoxes(mi,dis,picks(i,1:ptrs(i),:));
       plot(bX,bY,'color',dis.boxColors(i,:));
    end;
end;
