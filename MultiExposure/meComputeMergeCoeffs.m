function [coeffs, mergedCTF, dctfs]=meComputeMergeCoeffs( freqs, CTPars, doses, nzeros, weights)
% [coeffs mergedCTF dctfs]=meComputeMergeCoeffs(freqs,CTPars,doses,nzeros,weights)
% or: [coeffs mergedCTF dctfs]=meComputeMergeCoeffs( freqs, mi, nzeros)
% Given frequency values, an array of nims CTPars structures and a vector
% of nims doses, in e/A^2, compute the merging coefficients, which are the
% coefficients by which the fourier images are weighted before summing.
% Freqs (in A^-1) can be of any dimension. Suppose freqs is nx x ny.  Then
% the returned variable coeffs is nx x ny x nims.
% Image 1 is taken to be the lowest defocus.  The coeffs also effect phase
% flipping. The CTPars structure is the same as used for the
% ContrastTransfer() function.
% - doses is the dose in e/A^2. The scalar or vector nzeros gives the number of zeros to preserve in a
% given image (image 1 always has the entire ctf preserved).  This prevents
% magnification of low-defocus images to fill zeros in high-defocus images.
% The default is 1.
% - weights is an optional array of weights given the various images.  To merge the
% first 2 of 3 images for example, you would set weights = [q r 0] where q
% and r are nonzero, but otherwise their values don't affect the
% coefficients.  If the mi structure is given, weights=mi.weights*mi.ctf.ampFactor
% - mergedCTF is the effective CTF computed according to the merging, but
% also includes the effect of the weights (here the values of the weights
% are additional scale factors for each image's ctf).  The dctfs is an
% array of the ctfs of each image, including radiation damage and truncated
% after nzeros zeros, but *not* multiplied by weights.
%
% Example 1, to graph the effect of merging,
% f=0:.001:.1;
% [coeffs mergedCTF]=ComputeMergeCoeffs(f,CTPars,[10 20]);
% plot(f, mergedCTF);
%
% Example 2, to do merging of images m1 and m2, already aligned:
% f1=fftn(m1);
% f2=fftn(m2);
% freq=fftshift(Radius(n)/(n*pixA));  % image size is n x n, and the
%                                     % pixel size is pixA angstroms.
% [coeffs ceff]=ComputeMergeCoeffs( freq, CTPars, [10 20]);
% fmerge=f1.*coeffs(:,:,1)+f2.*coeffs(:,:,2);
% mmerge=real(ifftn(fmerge));
%
% memory use for a 4k x 4k frequency array is about (6+3*nim)*68MB.
% fixed modeling of decay (using dctfs in computing merged CTF) fs Apr 2013

% % test code:
% nargin=5;
% freqs=(0:.0001:.3)';
% % freqs=Radius(4096)/1e4;
% % disp('start');
% %
% P.lambda=.025;
% P.Cs=2;
% P.alpha=.07;
%
% defvals=[.6 3 10];
% doses=[10 10 20];
% nzeros=1;
% weights=[1 1 1];
% for i=1:3
%     P.defocus=defvals(i);
%     P.B=100*defvals(i);
%     CTPars(i)=P;
% end;

% -----------
epsi=1e-20;  % smallest noise variance to preserve.
% underflow of single.

% Handle the two options for arguments.
if isfield(CTPars,'ctf') % This is actually an mi structure
    mi=CTPars;
    CTPars=mi.ctf;
    nzeros=doses;
    doses=mi.doses;
    weights=mi.weights;
    if isfield(mi.ctf(1),'ampFactor')
        for i=1:numel(CTPars)
            weights(i)=mi.weights(i)*mi.ctf(i).ampFactor;
        end;
    end;
    if nargin<3
        nzeros=1;
    end;
else
    if nargin<4
        nzeros=1;
    end;
    if nargin<5
        weights=ones(1,numel(CTPars));
    end;
end;


nim=numel(CTPars);

if numel(nzeros)<nim
    nzeros=ones(nim,1)*nzeros;
end;

weights=weights(:);

siz=size(freqs);

% create the dimension of the output variables
sizx=siz;
ndims=numel(sizx);
if sizx(ndims)==1  % remove a trailing singleton dimension
    sizx=siz(1:ndims-1);
end;
sizx=[sizx nim];

nel=numel(freqs);  % 1d variable index

f=abs(single(freqs(:)));  % take the absolute value!
doses=single(doses(:));
cumdose=cumsum([0; doses]);

% find the first exposure with a nonzero dose, and let that be the standard
q=find(doses,1);
dose1=doses(q);

% Compute the decay constant due to radiation, as a function of frequency
% --a fit to the data in the Rubenstein paper.
%%%%%% n0=2*0.16./(abs(f)+.004)+5;  % twice the critical dose at 200 kV
% modified to better match vesicle decay.
n0=.32./(abs(f.^2*10)+.002)+5;

iFirst=find(weights,1); % mark the first nonzero weight
nzeros(iFirst)=inf;  % by default, no masking of 1st exposure.

% Get the coefficients coeff = ctf*signal
% Coefficients should be proportional to signal amp / rms noise
% where signal is the relative signal after radiation damage.
dctfs=single(zeros(nel,nim));
coeffs=single(zeros(nel,nim));

if ~isfield(CTPars(nim),'ampFactor') || numel(CTPars(nim).ampFactor)<1
    for i=1:nim
        CTPars(i).ampFactor=1;
    end;
end;

for i=find(doses(:)'>0)
    % compute the normalized signal amplitude after decay (=1
    %     given no decay)
    signal=(exp(-cumdose(i)./n0).*n0.*(1-exp(-doses(i)./n0)))/doses(i);
    [ctf,chi]=ContrastTransfer(freqs,CTPars(i));
    fmask=abs(chi(:))<nzeros(i);  % frequency mask beyond zeros.
    %     ctf including decay and mask
    dctfs(:,i)=ctf(:).*signal.*fmask;  % no ampFactor is applied.
    %     but we check whether the weight is nonzero.
    coeffs(:,i)=dctfs(:,i)*doses(i)/dose1*(weights(i)>0);
end;

% compute the frequency-dependent denominator
D=single(zeros(nel,1));
for i=1:nim
    D=D+coeffs(:,i).^2;  % sum the noise variance
end;
D(D<epsi)=1;
D=sqrt(D);  % noise standard deviation given non-normalized coefficients
%%
mergedCTF=single(zeros(nel,1));
for i=1:nim
    coeffs(:,i)=coeffs(:,i)./D;  % normalize to signal at dose=1.
    %     Note that in computing the MergedCTF we *do* use the weights.
    mergedCTF=mergedCTF+coeffs(:,i).*dctfs(:,i)*weights(i)*doses(i)/dose1;
end;

coeffs=reshape(coeffs,sizx);
dctfs=reshape(dctfs,sizx);
mergedCTF=reshape(mergedCTF,siz);

% freqs=reshape(freqs,siz);

% % more test code...
% subplot(2,1,1);
% semilogx(freqs,[abs(signals.*ctfs)/10 mergedCTF]);
%
% subplot(2,1,2);
% semilogx(freqs,abs(coeffs));
%
% % imacs(mergedCTF);
