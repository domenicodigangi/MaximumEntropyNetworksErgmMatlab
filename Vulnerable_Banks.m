function [AV,SYS,VUL] = Vulnerable_Banks(mode,input_data,equity,shock,varargin)....liq,sizes,leverage_cap,n_sample,precision)
    
% 
% This function computes the systemic risk quantifiers from the Vulnerable 
% Banks framework described in 
%   [1] Robin Greenwood, Augustin Landier, David Thesmar, Vulnerable banks, Journal of Financial Economics, Volume 115, Issue 3, March 2015, Pages 471-485, ISSN 0304-405X, http://dx.doi.org/10.1016/j.jfineco.2014.11.006.
%       (http://www.sciencedirect.com/science/article/pii/S0304405X14002529)
%
% and also employed in 
%   [2] Duarte, Fernando and Eisenbach, Thomas M., Fire-Sale Spillovers 
%       and Systemic Risk (February 1, 2015). FRB of New York Staff Report 
%       No. 645. Available at SSRN: http://ssrn.com/abstract=2340669 

% The most important feature is that it can ESTIMATE the Vulnerable Banks'
% quantities also when the BIPARTITE WEIGHTED NETWORK of investors and 
% assets is UNKNOWN BUT ONLY PARTIAL INFORMATION IS AVAILABLE
% It does so with the application of different network ensembles described
% in 
%   [3] Di Gangi, Domenico and Lillo, Fabrizio and Pirino, Davide,
%       Assessing Systemic Risk Due to Fire Sales Spillover Through 
%       Maximum Entropy Network Reconstruction (August 3, 2015). 
%       Available at SSRN: http://ssrn.com/abstract=2639178 
%
% The estimation requires the Max_entr_nets_ensembles package of witch this
% function is only one component. In particular the Max_Entr_Nets function
% is needed, as weel as  at least one BIPARTITE network model, that needs
% to be in the models folder.

% INPUT:
%       mode -> Indicates whether the function is required to compute the
%               real statistics or estimate them from partial information.
%           
%            mode = 'REAL' -->          Return the systemic risk measures
%                                       as computed by DUARTE-EISENBACH 
%                                       based on the knowledge of the 
%                                       weighted adjacency matrix X, 
%                                       expected as input_data. As a
%                                       convention each row of the matrix
%                                       is associated to an investor, while
%                                       each column to an asset.
%       mode = 'REAL-GREEN' -->          Return the systemic risk measures
%                                       as computed by GREENWOOD.
%            mode = 'ESTIMATE'->        Estimate the vulnerable banks
%                                       statistics from partial information
%                                       . Using the MECAPM model described
%                                       in [3]
%           mode = 'ESTIMATE-LIST' ->   Print a list of the models
%                                       available for the estimation from 
%                                       partial informaiton
%           mode = 'ESTIMATE-BIP***'->  Substituting *** with the name of a
%                                       model available, i.e.
%                                       'ESTIMATE-BIPWCM' returnns the
%                                       estimates based on the
%                                       corresponding network ensemble. The
%                                       input_data needs to be exactly that
%                                       required from the particular model
%                                       BIP***. 
%       
%       input_data -->  Contains the available information about the
%                       configuration of the network structure. This info
%                       can be complete or partial. 
%
%                       When mode == 'REAL' then input_data = X; where X is 
%                       the weighted adjacency matrix, that describes the 
%                       bipartite weighted network of investors and assets.
%
%                       When mode == 'ESTIMATE' , input_data needs to be a
%                       cell array with 2 columns. the fist cell element
%                       needs to contain a column array with the investors' 
%                       capitalizaiton, the second cell element a column 
%                       of the assets capitalization.
%
%                       When mode == 'ESTIMATE-BIP***' the input data needs
%                       to contain the information required by the model
%                       BIP*** . Type help BIP*** for the requirements of a
%                       specific model BIP***.

%
%       equity -->      Contains a vector of the equities of each investor.
%                       The order needs to be the same of input_data, i.e.
%                       equity(1) is the equity of the investor associated
%                       to the first row of the adjacency matrix.
%
%       shock-->         shock needs to be a ROW vector
%                       of real numbers between 0 and 1 that indicates the entity of
%                       the initial shock that triggers the vulnerable
%                       banks dynamics described in [1]. If not provided a
%                       uniform
%                       shock of 0.01 is assumed, i.e. a 1% depreciation
%                       of al the assets present in the system.
%
%       liq-->          The vector of illiquidities associated to each
%                       asset. If not provided the first asset is assumed
%                       to be cash, while all other assets have an assumed
%                       liquidity of 10^(-13), means that "€10 billion $ of 
%                       trading imbalnces lead to a price change of 10 
%                       basis points."[1]
%                       Note that if the unit measure is not 1$ the
%                       liquidity needs to be rescaled.
%               
%      sizes-->        A vector with the sizes of alll the banks if
%                       sizes
%                       are different from the sum of all the assets owned 
%                       by each bank. For example Greenwood et.al. use such
%                       an approach in their paper. 
%
%       leverage_cap--> The maximum value allowed for institutions 
%                       leverages. 



% Output:
%           AV --->     The Aggregate Vulnerability of the system, as defined 
%                       in [1].
%
%           
%           SYS -->     A vector of the systemicnes of each investor.
%
%
%           VUL -->     A vector of the Indirect Vulnerabilities of each
%                       investor.

%%

        
if strcmpi(mode,'REAL')
    % Vulnerable banks when the network is available
        
        %In order to compute the "real" quantities wee need to know the 
        % weighted adjacency matrix a.k.a. matrix of expositions
        X = input_data;
        n_assets = length(X(1,:));
        n_invest = length(X(:,1));
        
    
        
        % set assets' liquidity as suggested in [1] equal for all assets
        liq = ones(n_assets,1).*10^(-13);
        liq(1) = 0; %the first asset is assumed to have zero liquidity, e.g. cash
        
        % use user provided liquidities, if given
        if nargin>=5 && ~isempty(varargin{1}) 
            liq = varargin{1};            
        end
        
        
        sr = sum(X,2);
        % check for size zero banks
        ind_banks_to_remove = (sr == 0);
        
        sizes = sr;
         % use user provided liquidities, if given
        if nargin>=6 && ~isempty(varargin{2}) 
            sizes = varargin{2};            
        end      
        % compute book leverage as in [1] and [2]
        lev = (sizes-equity)./equity;   

        Lost_amount = (ones(n_invest,1)*shock).*X;
        r =  sum(Lost_amount,2)./sr;    
        r(ind_banks_to_remove) = 0;
        % System Equity 
        e_tot =sum(equity);        
        % compute the vector of banks' returns that follow the shock
        SYS =  (lev.*r).*(  sum( repmat( sum(X.*repmat(liq',n_invest,1)),n_invest,1 ).*X ,2) )/e_tot  *100;
        AV = sum(SYS);
        
        VUL = sum(X.*repmat(liq'.*  sum(X.*repmat(lev.*r,1,n_assets)),n_invest,1) ,2)./equity;

     
elseif strcmpi(mode,'REAL_THESIS')
    % Vulnerable banks when the network is available. 
    % in this version we compute the W matrix from X_input and then obtain
    % X as W*sizes_input
        
        %In order to compute the "real" quantities wee need to know the 
        % weighted adjacency matrix a.k.a. matrix of expositions
        X = input_data;
        n_assets = length(X(1,:));
        n_invest = length(X(:,1));
        
    
        
        % set assets' liquidity as suggested in [1] equal for all assets
        liq = ones(n_assets,1).*10^(-13);
        liq(1) = 0; %the first asset is assumed to have zero liquidity, e.g. cash
        
        % use user provided liquidities, if given
        if nargin>=5 && ~isempty(varargin{1}) 
            liq = varargin{1};            
        end
        
        
        sr = sum(X,2);
        % check for size zero banks
        ind_banks_to_remove = (sr == 0);
        
 
         % use user provided sizes, if given, for both leverage and matrix
         % computation
        if nargin>=6 && ~isempty(varargin{2}) 
            sizes = varargin{2};  
            W = X./repmat(sr,1,n_assets);
            % null banks have null weights (null banks might result from
            % ensemble sampling)
            W(ind_banks_to_remove,:) = 0;
            X = W.*repmat(sizes,1,n_assets);
           
        else
            error('REAL2 Method requires the sizes to be a separate input wrt X')
        end      
        % compute book leverage as in [1] and [2]
        lev = (sizes-equity)./equity;   

        Lost_amount = (ones(n_invest,1)*shock).*X;
        r =  sum(Lost_amount,2)./sizes;    
        r(ind_banks_to_remove) = 0;
        % System Equity 
        e_tot =sum(equity);        
        % compute the vector of banks' returns that follow the shock
        SYS =  (lev.*r).*(  sum( repmat( sum(X.*repmat(liq',n_invest,1)),n_invest,1 ).*X ,2) )/e_tot  *100;
        AV = sum(SYS);
        
        VUL = sum(X.*repmat(liq'.*  sum(X.*repmat(lev.*r,1,n_assets)),n_invest,1) ,2)./equity;
      
elseif strcmpi(mode,'REAL-GREEN')
    % Vulnerable banks when the network is available
        
        %In order to compute the "real" quantities wee need to know the 
        % weighted adjacency matrix a.k.a. matrix of expositions
        X = input_data;
        n_assets = length(X(1,:));
        n_invest = length(X(:,1));
         sizes = varargin{2};
      W = X./repmat(sizes,1,n_assets);
       
        % set assets' liquidity as suggested in [1] equal for all assets
        liq = ones(n_assets,1).*10^(-13);
        liq(1) = 0; %the first asset is assumed to have zero liquidity, e.g. cash
        
        % use user provided liquidities, if given
        if nargin>=5 && ~isempty(varargin{1}) 
            liq = varargin{1};            
        end
       
        sr = sum(X,2);
        % compute book leverage as in [1] and [2]
        lev = (sizes-equity)./equity;   
   
        Lost_frac = (ones(n_invest,1)*shock).*W;
        r =  sum(Lost_frac,2);           
        % check for size zero banks
        ind_banks_to_remove = (sr == 0); 
        r(ind_banks_to_remove) = 0;       
        
        
        lev_cap_greenwood = 30;
         lev_cap = lev_cap_greenwood;
        % use user provided leverage cap if any
        if nargin>=7 && ~isempty(varargin{3})            
            lev_cap = varargin{3};
        end 
        
        ind_lev_sup_cap = lev>lev_cap;           
         
        if sum(ind_lev_sup_cap)>0
            lev(ind_lev_sup_cap) = lev_cap;
            %   warning( [num2str(sum(ind_lev_sup_cap)) ' Investors Excedded cap leverage'])
        end            
    
        
        % System Equity 
        e_tot =sum(equity); 
        %fraction of each portfolio to be sold in order to target leverage
        sale_part= (lev.*r);
        %or as consequence of default
        ind_def = (((r.*sizes)./equity)>=1);
%       sum(ind_def)
        sale_part(ind_def) = 1-r(ind_def);
             
        SYS =  sale_part.*(  sum( repmat( sum(X.*repmat(liq',length(sr),1)),length(sr),1 ).*X ,2) )/e_tot  *100;
        AV = sum(SYS);
        VUL = sum(X.*repmat(liq'.*  sum(X.*repmat(sale_part,1,n_assets)),n_invest,1) ,2)./equity;
        
   
            

        
        
elseif strcmpi(mode(1:8),'ESTIMATE')
    % Estimate the statistics if required
         
        n_invest = length(input_data{1});
        n_assets = length(input_data{2});
        n_sample = 1000; 
        prec = 10^-5;        
        
      
        

        % set assets' liquidity as suggested in [1] equal for all assets
        liq = ones(n_assets,1).*10^(-13);
        liq(1) = 0; %the first asset is assumed to have zero liquidity, e.g. cash
        
        % use user provided liquidities, if given
        if nargin>=5 && ~isempty(varargin{1}) 
            liq = varargin{1};
            
        end
        
         if nargin>=6 && ~isempty(varargin{2}) 
            sizes = varargin{2};
         else
            sizes=[];
         end
        
        
      
         if nargin>=7 && ~isempty(varargin{3})            
            lev_cap = varargin{3};
         else
             lev_cap = [];
         end
        
          if nargin>=8 && ~isempty(varargin{4})            
            n_sample = varargin{4};
         
          end
          
          
          if nargin>=9 && ~isempty(varargin{5})            
            prec = varargin{5};
         
          end
          
          
          
          if strcmpi(mode,'ESTIMATE') || strcmpi(mode,'ESTIMATE-MECAPM')
            % if the model is not specified use MECAPM
            model_ens = Max_Entr_Nets('MECAPM',input_data,10^-3);
            else
            method = mode(10:end);
            model_ens = Max_Entr_Nets(method,input_data,prec);
            if strcmpi(method,'EMECAPM');
                %n_sample = 1000;
            end
          end
 
        %sample the ensemble
        model_sample = model_ens.sample(n_sample);
 
        %compute the statistics for each element of the ensemble
        AV_sample = zeros(1,n_sample);
        SYS_sample = zeros(n_invest,n_sample);
        VUL_sample = zeros(n_invest,n_sample);
        
       %toc
        for i = 1:n_sample
           X_i = model_sample{i} ;
            [AV_tmp,SYS_tmp,VUL_tmp] = Vulnerable_Banks_duarte('REAL',X_i,equity,shock,liq,sizes);
            AV_sample(i) = AV_tmp;
            SYS_sample(:,i) = SYS_tmp;
            VUL_sample(:,i) = VUL_tmp;      
        
        end
       % vul_time = toc
        AV = mean(AV_sample);
        SYS = mean(SYS_sample,2);
        VUL = mean(VUL_sample,2);
        
else
    error('invalid mode')
end



end
















