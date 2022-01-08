function [Y,xmcv,optjmp,stime] = psp(model,x0,xbound,options,varargin)

    %PSP Finds data patterns a computational model can possibly generate.
    %   The PSP (Parameter Space Partitioning) algorithm finds discrete data patterns
    %   the given model could produce by identifying their corresponding regions in
    %   the unknown partition of the model's parameter space.
    %
    %   Y=PSP(MODEL,X0,XBOUND) searches the parameter space defined by XBOUND
    %   for the data patterns of MODEL, with starting point(s), X0, and
    %   returns discovered data patterns. The function MODEL accepts parameter
    %   values in a column vector and returns a pattern of data in a row vector.
    %   MODEL must always return a row vector of the same length for any given
    %   point ranged within XBOUND. X0 must be an NDIM-by-M matrix where NDIM
    %   is the number of model parameters and M is the number of starting
    %   points--the user can supply multiple starting points.
    %   XBOUND is an NDIM-by-2 matrix whose first column sets the lower bounds
    %   of parameters, and second column, the upper bounds. XBOUND must be
    %   finite (i.e., Inf or -Inf cannot be used). In case some unconstrained
    %   parameters are unavoidable (i.e., plausible data patterns could be
    %   generated from their extreme values), it is recommended to
    %   reparameterize MODEL, making use of log, inverse logistic, or other
    %   transformation function. Y stores the discovered data patterns in rows.
    %
    %   Y=PSP(MODEL,X0,XBOUND,OPTIONS) runs with the default search parameters
    %   replaced by values in the structure OPTIONS, an argument created with the
    %   OPTSETPSP function. See OPTSETPSP for details. Used options are MaxPSP,
    %   IniJump, SmpSize, VolEst, VolSmpSize, Display, DispInterval.
    %   Use OPTIONS=[] as a place holder if no options are set.
    %
    %   Y=PSP(MODEL,X0,XBOUND,OPTIONS,P1,P2,...) passes the problem-dependent
    %   parameters P1,P2,... directly to the function MODEL. In PSP, MODEL is
    %   evaluated in the form: feval(MODEL,X,P1,P2,...).
    %
    %   [Y,XMCV]=PSP(MODEL,X0,XBOUND,...) returns information about the partition of
    %   parameter space. XMCV is a 1-by-4 cell array whose first element is a matrix
    %   containing parameter vectors that can generate the discovered data patterns,
    %   and second element is a matrix containing mean parameter vectors of the
    %   regions identified in the partition, and third element is a 3-dimensional
    %   matrix in which the covariance matrix of parameters over each region stacks
    %   up along the 3rd dimension, and, finally fourth element is a row vector of
    %   approximated volumes of the regions on a log scale.
    %
    %   [Y,XMCV,OPTJMP]=PSP(MODEL,X0,XBOUND,...) returns the results of adapting the
    %   Markov chain in each discovered region. The size of adapted jumping
    %   distribution (radius of a sphere) for the i-th region can be retrieved by
    %   IniJump*2^OPTJMP(i).
    %
    %   [Y,XMCV,OPTJMP,STIME]=PSP(MODEL,X0,XBOUND,...) returns the information of
    %   (cumulative) search time elapsed before finding each data pattern. STIME is a
    %   2-by-(N+1) matrix whose first row lists time points at which each of N data
    %   patterns is discovered, and the time point of terminating the algorithm in the
    %   last column. The second row lists the same information except that it is in
    %   the total number of search trials made before finding each new pattern.
    
    %   Authorship: Woojae Kim, Department of Psychology, Ohio State University 
    %   $Revision: 3.0 $  $Date: 2005/07/19 $
    
    format compact
    rand('state', sum(100*clock)), randn('state', sum(100*clock))
    
    if nargin < 3, error('Not enough input arguments')
    elseif nargin == 3, options = [];
    end
    
    xrange = xbound(:,2)-xbound(:,1);
    if any(xrange<=0), error('Invalid XBOUND specification'), end
    
    if isempty(x0) | any(any(x0<repmat(xbound(:,1),1,size(x0,2)) | x0>repmat(xbound(:,2),1,size(x0,2))))
       error('Invalid X0')
    end
    
    if isempty(options), options=optsetpsp('Joke'); end
    
    maxpsp = options.MaxPSP;
    inijmp = options.IniJump;
    smpsz = options.SmpSize;
    volest = options.VolEst;
    vsmpsz = options.VolSmpSize;
    dsp = options.Display;
    dispintv = options.DispInterval;
    
    ndim = size(xbound,1);
    
    % === Default values of options ===
    if isempty(maxpsp), maxpsp = 6; end
    if isempty(inijmp), inijmp = .1; end
    if isempty(smpsz), smpsz = ceil([100 200]*1.2^ndim); end
    if isempty(vsmpsz), vsmpsz = ceil(500*1.2^ndim); end
    if isempty(dsp) | (~strncmpi(dsp,'off',1)&~strncmpi(dsp,'min',1)&~strncmpi(dsp,'final',1))
       dsp = 'all';
    end
    if isempty(dispintv), dispintv = 30000; end
    
    if ~isstr(model), mdstr = func2str(model); else, mdstr = model; end
    
    % ********** MCMC-based Parameter Space Partitioning Algorithm ***************
    
    t0 = deal(cputime);
    [x, ptnfound, stime, smpcnt] = deal([]);
    [ntr, icnt1, psz] = deal(0);
    
    if strncmpi(dsp,'all',1) | strncmpi(dsp,'min',1)
       disp('=================================================================')
       disp(['PSP SEARCH FOR ' mdstr ' STARTS...']), disp(' ')
    end
    
    for h = 1:size(x0,2)
       y = x0(:,h);
       ntr = ntr + 1;
       currptn = feval(model,y,varargin{:});
       if isempty(ptnfound) | all(any(ptnfound ~= repmat(currptn,psz,1),2))
          x = [x, y];
          ptnfound = [ptnfound; currptn];
          psz = size(ptnfound,1);
          stime = [stime, [cputime-t0; ntr]];
          if strncmpi(dsp,'all',1) | strncmpi(dsp,'min',1)
             disp(['New data pattern found: Region #' num2str(psz)])
          end
          if strncmpi(dsp,'all',1)
             disp(['w/ supplied starting point(s), Total elapsed time: ' ...
                   num2str(round(stime(1,end))) ' secs (' num2str(ntr) 'trials)'])
          end
       end
    end
    
    [smpcnt, optjmp, levels, alps] = deal( zeros(1, psz) );
    [icnt1, icnt2] = deal(0);
    [cnt1, cnt2] = deal(cputime);
    xsum = zeros(ndim,psz);
    xcsum = zeros(ndim,ndim,psz);
    
    maxpspp = maxpsp*smpsz(2);
    
    while any(levels<2) | min(smpcnt) <= maxpspp
       
       tmp = min(levels);
       ind1 = find(levels==tmp);
       [tmp, ind2] = min(smpcnt(ind1));
       rgind = ind1(ind2);
       smpcnt(rgind) = smpcnt(rgind) + 1;
       
       tmp = randn(ndim,1); tmp = rand^(1/ndim)* tmp /sqrt(tmp'*tmp);
       y = x(:,rgind) + xrange.*( inijmp*(2^optjmp(rgind))*tmp );
       ntr = ntr + 1;
       
       if y>=xbound(:,1) & y<=xbound(:,2)
          
          currptn = feval(model,y,varargin{:});
          
          if all(currptn == ptnfound(rgind,:))
             x(:,rgind) = y;
             alps(rgind) = alps(rgind) + 1;
          else
             
             if all(any(ptnfound ~= repmat(currptn,psz,1),2))
                x = [x, y];
                ptnfound = [ptnfound; currptn];
                psz = size(ptnfound,1);
                stime = [stime, [cputime-t0, ntr]'];
                smpcnt = [smpcnt, 0]; optjmp = [optjmp, 0];
                levels = [levels, 0]; alps = [alps,0];
                xsum = [xsum, zeros(ndim,1)]; xcsum = cat(3, xcsum,zeros(ndim));
                [icnt1, icnt2] = deal(0); [cnt1, cnt2] = deal(cputime);
                if strncmpi(dsp,'all',1) | strncmpi(dsp,'min',1)
                   disp(' '), disp(['New data pattern found: Region #' num2str(psz)])
                end
                if strncmpi(dsp,'all',1)
                   disp(['PSP, Total elapsed time: ' ...
                         num2str(round(stime(1,end))) ' secs (' num2str(ntr) ' trials)'])
                end
             end
             
          end
       end
        
        if levels(rgind) == 0
            tmp = smpcnt(rgind) / smpsz(1);
            if tmp == ceil(tmp)
                acrate = alps(rgind) / smpsz(1);
                alps(rgind) = 0;
             if strncmpi(dsp,'all',1)
                disp(' '), disp(['Level 1 adaptation of MCMC in Region #' num2str(rgind)])
                disp(['Cycle #' num2str(tmp) ', Acceptance rate: ' num2str(acrate,4)])
             end
                if acrate < .12
                    if optjmp(rgind) > 0
                        optjmp(rgind) = optjmp(rgind)-.5;
                   levels(rgind) = 1; smpcnt(rgind)=0;
                    else
                        optjmp(rgind) = optjmp(rgind)-1;
                    end
                elseif acrate >= .12 & acrate < .36
                    levels(rgind) = 1; smpcnt(rgind)=0;
                elseif acrate >= .36
                    if optjmp(rgind) < 0
                        optjmp(rgind) = optjmp(rgind)+.5;
                   levels(rgind) = 1; smpcnt(rgind)=0;
                    else
                        optjmp(rgind) = optjmp(rgind)+1;
                    end
                end
            end
            
        elseif levels(rgind) == 1
            tmp = smpcnt(rgind) / smpsz(2);
            if tmp == ceil(tmp)
                acrate = alps(rgind)/ smpsz(2);
                alps(rgind) = 0;
             if strncmpi(dsp,'all',1)
                disp(' '), disp(['Level 2 adaptation of MCMC in Region #' num2str(rgind)])
                disp(['Cycle #' num2str(tmp) ', Acceptance rate: ' num2str(acrate,4)])
             end
                if acrate < .15
                    optjmp(rgind)=optjmp(rgind)-.25/ceil(tmp/2);
                    if tmp==4, levels(rgind)=2; smpcnt(rgind)=0; end
             elseif acrate >= .15 & acrate < .19
                optjmp(rgind)=optjmp(rgind)-.125;
                levels(rgind)=2; smpcnt(rgind)=0;
             elseif acrate >= .19 & acrate < .24
                    levels(rgind)=2; smpcnt(rgind)=0;
             elseif acrate >= .24 & acrate < .3
                optjmp(rgind)=optjmp(rgind)+.125;
                levels(rgind)=2; smpcnt(rgind)=0;
                elseif acrate >= .3
                    optjmp(rgind)=optjmp(rgind)+.25/ceil(tmp/2);
                    if tmp==4, levels(rgind)=2; smpcnt(rgind)=0; end
                end
            end
            
        elseif levels(rgind) == 2
          if strncmpi(dsp,'all',1)
             tmp = smpcnt(rgind) / smpsz(2);
             if smpcnt(rgind)==1
                disp(['Adaptation of MCMC in Region #' num2str(rgind) ' finished'])
             elseif tmp == ceil(tmp)
                acrate = alps(rgind)/ smpcnt(rgind);
                disp(' '), disp(['Monitoring after adaptation in Region #' num2str(rgind)])
                disp(['Cycle #' num2str(tmp) ', Acceptance rate (cumulative): ' num2str(acrate,4)])
             end
          end
            xsum(:,rgind) = xsum(:,rgind) + x(:,rgind);
            xcsum(:,:,rgind) = xcsum(:,:,rgind) +  x(:,rgind)*x(:,rgind)';
            
        end
        
        icnt1 = icnt1 + 1;
        if any(levels<2) | max(smpcnt)-min(smpcnt)>1
            icnt2=0; cnt2=cputime;
        else
            icnt2 = icnt2 + 1;
        end
        
       if strncmpi(dsp,'all',1) | strncmpi(dsp,'min',1)
          tmp = icnt1/dispintv;
          if tmp == ceil(tmp)
             disp([num2str(round(cnt2-cnt1)) ' + ' num2str(round(cputime-cnt2)) ' secs (' ...
                   num2str(icnt1-icnt2) ' + ' num2str(round(icnt2/psz)) 'x' num2str(psz) ...
                   ' trials) elapsed since last discovery'])
          end
       end
       
    end
    
    Y = ptnfound;
    
    xmean = xsum ./ repmat(smpcnt,ndim,1);
    
    xcovmat = zeros(ndim,ndim,psz);
    for m = 1:psz
        xcovmat(:,:,m) = xcsum(:,:,m)/smpcnt(m) - xsum(:,m)*xsum(:,m)'/smpcnt(m)^2;
    end
    
    logvol = zeros(1,psz);
    nhalf = ndim/2; nfloor = floor(nhalf);
    if nhalf == nfloor
       tmp = nhalf*log(pi)-gammaln(nhalf+1);
    else
       tmp = ndim*log(2)+gammaln(nfloor+1)-gammaln(ndim+1)+nfloor*log(pi);
    end
    for k=1:psz
       logvol(k)= tmp + .5*sum( log(ndim+2) + log( eig(xcovmat(:,:,k)) ) );
    end
    
    if strncmpi(volest,'hitmiss',1)
       if strncmpi(dsp,'all',1) | strncmpi(dsp,'min',1)
          disp(' '), disp('Volume estimation by hit-or-miss method begins...')
       end
       for k=1:psz
          nhit = 0;
            if strncmpi(dsp,'all',1)
                disp(['Estimating the volume of Region #' num2str(k) ' / ' num2str(psz) ' ...'])
            end
          for r=1:vsmpsz
             tmp = randn(ndim,1); tmp = rand^(1/ndim)* tmp /sqrt(tmp'*tmp);
             y = xmean(:,k) + sqrtm((ndim+2)*xcovmat(:,:,k))*tmp;
                if y>=xbound(:,1) & y<=xbound(:,2)
                    cptn = feval(model,y,varargin{:});
                    if cptn == ptnfound(k,:), nhit = nhit + 1; end
                end
          end
          logvol(k) = logvol(k) + log(nhit) - log(vsmpsz);
       end
       if strncmpi(dsp,'all',1) | strncmpi(dsp,'min',1)
          disp(['...Volume estimation terminated for all ' num2str(psz) ' regions.'])
       end
    end
    
    xmcv = {x, xmean, xcovmat, logvol};
    
    stime = [stime, [cputime-t0, ntr]'];
    
    tmp = stime(1,end);
    hrs = floor(tmp/3600);
    mins = floor(tmp/60 - hrs*60);
    secs = floor(tmp - mins*60 -hrs*3600);
    
    if ~strncmpi(dsp,'off',1)
       disp(' ')
       disp(['PSP SEARCH FOR ' mdstr ' TERMINATED.'])
       disp(['TOTAL ' num2str(psz) ' DATA PATTERNS FOUND.'])
       disp(['TOTAL ' num2str(hrs) ' hours ' num2str(mins) ' mins ' ...
             num2str(secs) ' secs (' num2str(ntr) ' trials) ELAPSED.']);
       disp('=================================================================')
    end