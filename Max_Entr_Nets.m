% Code written by Domenico Di Gangi, domenico.digangi@sns.it   

 function [Out_MENE, Out_details] =...
                        Max_Entr_Nets(model,in_data,varargin)
                                                  %, prec)   
                     
% This function allows the estimation and sampling of maximum entropy 
% network ensembles, according to specific structures (e.g. unipartite
% directed, bipartite weighted, etc.), and observed values of the 
% constraining statistics that define the model.
% In order to estimate and sample a specific model the corresponding
% file need to be present in the models folder
%
 
% DESCRIPTION OF THE INPUTS:
%       model ---->  fixes the structure common to all the networks in the ensemble and
%                   the constrained statistics. Can be set equal to each one of the 
%                   models present in the models folder. e.g.:
%                   
%                   meth   == BIPCM      ->     Bipartite Configuration 
%                                               Model: fixed degree 
%                                               sequences of a binary
%                                               network. 
%                                               Data needed: degree sequences.


%                           == BIPWCM      ->   Bipartite Weigthed Configuration 
%                                               Model: fixed strength
%                                               sequences. Data
%                                               needed:strength sequences.
%
%
%                           == BIPECM      ->   Bipartite Enhanced Configuration 
%                                               Model: fixed strength and
%                                               degree sequences. Data
%                                               needed:strength sequences,
%                                               and degree sequences.

%                                              
%                        

%                           == BIPMECAPM   ->   Bipartite Maximum Entropy CAPM:
%                                               Fixed expected matrix equal to the
%                                               one predicted by the CAPM.
%                                               Data needed: strength
%                                               sequences

%                           For a complete list of the models check the
%                           model folder. Each file should contains a brief
%                           description.
%                   .
%                   .
%
%       in_data  ---->  The input data is a  cell array of columns. In
%                       particular, the first column represents the
%                       strenghts(or degrees) of the rows (of the adjacenct matrix), the second
%                       the strengths of the columns and any additional
%                       sequence of constraints is added as a new entry in
%                       the cell array

%                       Example: for BIPBCM  given matrix   [1,0;
%                                                            1,0] 
%
%                       it should be  in_data_list = {[1;1],[2;0]};
%    
%       prec    --->    When specified it sets the required maximum of the
%                       relative difference between in_data values and
%                       ensemble averages of the constraints, in order to
%                       stop the optimization
%                       
%                     
%                      
%  
% DESCRIPTION OF THE OUTPUT:
%       Out_MENE --->   A struct variable that is associated to the Maximum
%                       Entropy Network Ensemble (MENE) required.
%                       
%       Out_MENE.parameters --> Field that contains the vector of the
%                               numerical values for the parameters of the
%                               model
%
%       Out_MENE.sample -->     A function that allows sampling from the
%                               ensemble. It takes only one argument, i.e.
%                               the sample size. e.g.  Out_MENE.sample(10)
%                               returns a cell array composed of 10
%                               adjacency matrices sampled from the
%                               ensemble.
%
%       Out_MEME.prec.errlike-> The obtained precision, i.e. the maximum,
%                               among the constrained statistics,relative
%                               of difference between their input value and
%                               the analytical average over the ensemble,
%                               obtained from the parameters
%                               Out_MENE.parameters.





%%

    % add the needed folders to matlab path 
    tmp_path = mfilename('fullpath');
    addpath(genpath(tmp_path(1,1:end-13)));
    
    %The user can ask for a list of the models available
    if strcmpi(model,'LIST-BIP') % list all bipartite models
       ls([tmp_path(1,1:end-13), 'models/BIP*'])
       Out_MENE = [];
       return
    elseif strcmpi(model,'LIST')% list all models available
        ls([tmp_path(1,1:end-13), 'models/*M.m'])
        Out_MENE = [];
        return
        
    end
        


% The required minimum relative precision for the constraints can be set 
if nargin >=3
     precision = varargin{1};
else
    precision = 10^(-2);
end
%  prec fixes the required precision, i.e. the maximum value of the expected
%  relative error  for each constraint
%  e.g. 
%  prec == 10^(-3); (the expected error over the ensembles is less than 0.001 for each constraint)
%     


%Select model among those available
eval(['max_ent_model = ' upper(model) '(in_data,precision);'])

%%  Estimate the Model 
% Look for a maximum of the appropriate log-likelihood function.
% For the optimization we need Matlab's optimization Toolboox, in
% particular fzero is repeatedly used, for most of the models
if nargin >=4 &&  ~isempty(varargin{2})
    % if, for some reasons the parameters of the model are already known,
    % then do not estimate them
    in_parameters = varargin{2};
    max_ent_model.estimation_time= 0;
    max_ent_model.parameters = in_parameters;
    max_ent_model.iterations_n = 0;
else
    tic
    % run the optimization
    tmp_opt =max_ent_model.optim_run();
    max_ent_model.estimation_time= toc;
    max_ent_model.parameters = tmp_opt.par;
    max_ent_model.iterations_n = tmp_opt.iter_num;
end
%Compute and store relative errors  
switch upper(max_ent_model.estimation)
    case 'ANALYTICAL'
        max_ent_model.precision.errors = 0; 
        max_ent_model.precision.maxerr  = 0;
        max_ent_model.precision.meanerr  = 0;    
    case 'NUMERICAL'              
        max_ent_model.precision.errors = max_ent_model.check_sys(max_ent_model.parameters ,1); 
        max_ent_model.precision.maxerr  = max(max(max_ent_model.precision.errors ));
        max_ent_model.precision.meanerr  = mean(mean(max_ent_model.precision.errors) );    
end
 
%If the errors obtained are too large signal it with a warning
if (max_ent_model.precision.maxerr  > precision)
    if (max_ent_model.precision.maxerr  > 10^(-2))
        warning('Obtained precision is not sufficient, i.e. errors are larger than 1%')
    else
        warning(['The required precision has not been reached. Instead ' ...
                    num2str(max_ent_model.precision.maxerr) ])
    end    
end
      
%% Define the output structure
% field containing a function that samples n matrices from the estimated 
% model
Out_MENE.sample = @(n)max_ent_model.sampling_fun(max_ent_model.parameters,n);
% Function that computes the expected matrix
Out_MENE.exp_matrix = @()max_ent_model.check_sys(max_ent_model.parameters,31);
% The parameters are available in the main output structure
Out_MENE.parameters = max_ent_model.parameters;
% Errors and estimation time are available
Out_MENE.precision.errors = max_ent_model.precision.errors;
Out_MENE.estimation_time = max_ent_model.estimation_time;

% Define the detailed output structure, to be returned only if required
Out_details = max_ent_model;
  
end

