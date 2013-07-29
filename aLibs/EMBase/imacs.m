function h=imacs(x,y,m)% function h=imacs(x,y,m);% function h=imacs(m);% Draw an image using cartesian coordinates and autoscaling,% with m(1,1) being in the standard position at lower left.% if x and y arguments are given, they set the scaling of ticks.% We assume a 256-level grayscale, as is set by the function SetGrayscale.% If an output argument h is given, the function returns the handle to the% image object.if nargin==3    if isa(m,'integer')        m=single(m);    end;    m=squeeze(m);    mn=min(min(m));    mx=max(max(m));    if (mx-mn)<eps*max(1,(mx+mn))  % Nothing there: show gray        mx=mx+1;        mn=mn-1;    end;    h=image(x,y,255*(m'-mn)/(mx-mn));    axis xy;    else % no scale arguments, just display the x variable.        if isa(x,'integer')        x=single(x);    end;    x=squeeze(x);    mn=min(min(x));    mx=max(max(x));    if (mx-mn)<eps*max(1,(mx+mn))  % Nothing there: show gray        mx=mx+1;        mn=mn-1;    end;    h=image(255*(x'-mn)/(mx-mn));    axis xy;end;if nargout<1    clear hend;