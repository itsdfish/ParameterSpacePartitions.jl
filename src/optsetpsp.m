function options = optsetpsp(varargin)

    %OPTIONS = OPTSETPSP('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates a PSP search
    %   options structure OPTIONS in which the named parameters have the specified
    %   values. Any unspecified parameters are set to [] (parameters with value []
    %   indicate to use the default value for that parameter when OPTIONS is passed
    %   to PSP). It is sufficient to type only the leading characters that uniquely
    %   identify the parameter. Case is ignored for parameter names.
    %     
    %   OPTSETPSP PARAMETERS:
    %
    %   MaxPSP - Maximum number of search cycles allowed before termination. 
    %   The PSP search stops when the number of search cycles of every Markov chain
    %   (in every discovered region) reaches MaxPSP. A cycle consists of a certain
    %   number of search trials (i.e., MCMC sample points). The cycle size is set
    %   by the 2nd element of SmpSize (see below). The couting of cycles does not
    %   begin until the adaptive stage of every Markov chain enters the stage of a
    %   regular MCMC sampling. The default value is 6. Setting MaxPSP to a low
    %   value may cause an incomplete search or unreliable information on the
    %   partition of parameter space.
    %   Note: PSP has an internal mechanism that makes the total number of MCMC
    %   trials in the search history eventually equalized for all discovred
    %   parameter regions. It is due to this mechanism that everytime a new region
    %   is identified, the algorithm concentrates for a while on the search
    %   near it.
    %
    %   IniJump - Initial size of jumping distribution with which the Markov chain of
    %   a newly found region starts. The jumping distribution of the MCMC process is a
    %   uniform density over a hyper-sphere. IniJump specifies the initial size of the
    %   sphere's radius. The regular MCMC sampling process in each discovered region
    %   does not begin until the size of jumping distribution is adapted so that it
    %   meets a certain rate of sample acceptance. IniJump does not depend on the scale
    %   of parameters because the range of every parameter is rescaled to [0 1] in PSP.
    %   Even so, IniJump is problem-dependent and should affect the efficiency of PSP.
    %   The default value is 0.1.
    %
    %   SmpSize - Sizes of MCMC samples used to adapt the Markov chain of each
    %   discovered region. The adaptation is performed in two stages each of which
    %   consists of a few sampling cycles. SmpSize specifies the size of MCMC sample
    %   in one cycle. SmpSize is a row vector of two integers, e.g., [SMPSZ1 SMPSZ2],
    %   where SMPSZ1 defines the size of one cycle in the Level 1 stage, and SMPSZ2 in
    %   the Level 2 stage. The former must be smaller than the latter (e.g., SMPSZ1 is
    %   half of SMPSZ2) since the Level 1 adaptation is designed to be fast but coarse-
    %   grained whereas the Level 2 adaptation to be slow but fine-grained. SMPSZ2
    %   also serves as the size of regular MCMC cycles following the adaptation, that
    %   are counted to determine the termination of PSP (see MaxPSP above). The default
    %   setting in PSP is CEIL([100 200]*1.2^NDIM).
    %
    %   VolEst - Method of volume estimation of dicovered regions. 'Ellipsoid'
    %   (default) estimates the volume through ellipsoid approximation using eigen
    %   decomposition of the estimated covariance matrix from each region's MCMC
    %   sample. Although this method does not require additional computing time, it
    %   causes a bias toward overestimation because an ellipsoid fills up the, if any,
    %   concavity of the region. 'HitMiss' uses a method inspired by the hit-or-miss
    %   Monte Carlo integration, which corrects the bias due to ellipsoid
    %   approximation. This method is the same as the regular hit-or-miss method
    %   except that the estimated ellipsoid is employed as a base sampling domain,
    %   which may exclude some portion of the region and may thus result in
    %   underestimation. This method, however, will give better estimates in general.
    %   Additional computing time is required after the completion of search process.
    %
    %   VolSmpSize - Size of a sample used to estimate a region's volume when VolEst
    %   is set as 'HitMiss.' Default is CEIL(500*1.2^NDIM).
    %
    %   Display - Level of display. 'off' displays no output; 'final' displays just the
    %   final output (number of patterns found & elapsed time); 'min' displays the
    %   region index each time a new data pattern is discovered; 'all' (DEFAULT)
    %   displays all the information regarding the status of search and MCMC adaptation.
    %
    %   DispInterval - Length of intervals (in search trials) at which PSP displays
    %   elapsed time while no new pattern is found. By DEFAULT, PSP does so every
    %   30,000 search trials.
    
    %   Authorship: Woojae Kim, Department of Psychology, Ohio State University 
    %   $Revision: 3.0 $  $Date: 2005/07/19 $
    
    [maxpsp,inijump,smpsize,volest,volsmpsize,display,dispinterval] = deal([]);
    
    nopt = floor(nargin/2);
    
    for k = 1:nopt
        par = varargin{2*k-1};
        val = varargin{2*k};
       if strncmpi(par,'MaxPSP',4)
            maxpsp = val;
        elseif strncmpi(par,'IniJump',1)
            inijump = val;
        elseif strncmpi(par,'SmpSize',1)
          smpsize = val;
       elseif strncmpi(par,'VolEst',4)
          volest = val;
       elseif strncmpi(par,'VolSmpSize',4)
          volsmpsize = val;
       elseif strncmpi(par,'Display',5)
            display = val;
       elseif strncmpi(par,'DispInterval',5)
          dispinterval = val;
       end
    end
    
    options = struct('MaxPSP',maxpsp,'IniJump',inijump,'SmpSize',smpsize,...
       'VolEst',volest,'VolSmpSize',volsmpsize,'Display',display,'DispInterval',dispinterval);