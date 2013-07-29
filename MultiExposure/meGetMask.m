function msk=meGetMask(mi,n,indices)
% function msk=meGetMask(mi,n,indices)
% From the mi.mask array of structures, generate a boolean image of size n.
% If the optional argument indices is given, use only those elements of the
% mask stack.
% The masked points have the value logical false.
% If no masks are present, return the default msk = true(n).
% If n=0 then we return only one mask at its native size.
% 
% mi.mask has the fields
%   merge  - a text string AND, OR, OVER, OFF
%   encoding - a text string RLE, beam
%       where RLE is run-length encoding of a binary image, and beam is
%       simply a fuzzy disc of given radius and position in original image
%       units.
%   data - data encoding the mask.

if n==0
    indices=indices(1);
    msk=true;
else
msk=true(n);  % default
end;
if isfield(mi,'mask')
    nim=numel(mi.mask);
    if nargin<3
        indices=1:nim;
    else
        if ~any(indices<=nim)
            return
        end;
        indices=indices(indices<=nim);
    end;
    for i=indices
        if numel(mi.mask(i).merge)>0  % something there
            switch mi.mask(i).encoding
                case 'RLE'
                    m1=RunLengthDecode(mi.mask(i).data);
%                 case 'beam'
                    % for now, do nothing.
                otherwise
                    error(['Unexpected mi.mask.encoding: ' mi.mask(i).encoding]);
            end;
            if n==0
                msk=m1;
                return
            end;

            n1=size(m1);
            if n1(1) ~= n(1)  % need to change size
                lcf=LeastCommonFactor(n(1),n1(1));
                if n(1)/lcf < 8  % don't allow huge expansions
                    m1=BinImage(ExpandImage(m1,n(1)/lcf),n1(1)/lcf);
                else
                    m1=DownsampleGeneral(m1,n,n(1)/n1(1))>0.5;
                end;
            end;
            switch mi.mask(i).merge
                case 'AND'
                    msk=msk & m1;
                case 'OR'
                    msk=msk | m1;
                case 'OVER' % Overwrite the earlier masks.
                    msk=m1;
                case 'OFF'
                    % do nothing
                otherwise
                    error(['Unexpected mi.mask.merge value: ' mi.mask(i).merge]);
            end;
        end;
    end;
else
    msk=single(ones(n));  % default is no mask.
end;