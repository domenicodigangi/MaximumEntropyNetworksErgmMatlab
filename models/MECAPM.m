function model_struct = MECAPM(in_data,precision)
      
%   This function contains all the functions, optimization options, and
%   description of the input data, required to estimate and sample 
%   the network ensemble known as:

%----------------------------------------------------------------------
%-------------------------MECAPM------------------------------------
%----------------------------------------------------------------------
%-------BIPARTITE FIXED STRENGTH SEQUENCE + CAPM EXPECTED MATRIX--------------
%
%   INPUT: in_data  is a cell array that contains the strength sequences.
%                   The first elementis a column with rows strengths, the
%                   second one a column with the columns strengths
%           
%   OUTPUT: model_struct    is a structure with a number of fields 
%                           described in the following. Each field refers
%                           to informations, funcitons or options specific
%                           for the particular model. 

%   To be called only at the beginning of the main script Max_Entr_Nets
%%
%Variables used by different sub functions defined in the following 
sr = in_data{1};       sc = in_data{2};
n_row = length(sr);    n_col = length(sc);
n_tot = n_row+n_col;

% Model Name
model_struct.name = 'MECAPM';           
model_struct.in_data = in_data; 
model_struct.estimation = 'Analytical';
model_struct.is_unipartite = false;
model_struct.is_bipartite = true;
model_struct.is_weighted = true;
%% Field for the values of model's parameters
% the parameters for BIPMECAPM can be obtained analytically
         
 

%% Optimization run
     function opt_out = optim_r()           
        exp_X_mecapm = sr*sc'./sum(sr);
        % for mecapm the parameters are compited from the expected matrix
        z_mecapm = 1./(1+exp_X_mecapm);
       
       
        
        
                
        opt_out.par =  z_mecapm;   
        opt_out.iter_num = 0;  
         
    end

model_struct.optim_run = @optim_r;
 %% Constraints' System Check

 % The following function evaluates the system of "soft" constraints 
 % imposed in the definition of a given model and returns the relative
 % errors between the observed values of the constraint statistics and
 % their Analytical expectations over the ensemble
 
function F = system_check(z,arg)
   
    % INPUT:    z       ---> Vector of model's parameters
    %           arg     ---> Flag that selects the desired form of the
    %                        output:
    %                       "arg" == 0 produce the difference between the 
    %                                  observed and the expected values of 
    %                                  the rows and colums sums.
    %                       "arg" == 1 Relative difference 
    %                       "arg" == 31 The Expected matrix 
    

    
    %Compute the expected values for the constraints 
    EXP_mat_s =  sr*sc'./sum(sr);
    Exp=[sum(EXP_mat_s,2);sum(EXP_mat_s)'];
    Obs=[sr;sc];
    
    if (arg == 0 )
        diff = Obs - Exp;
        F = diff;
    elseif (arg == 1 )
        diff = Obs - Exp;
        F = abs(diff)./Obs;
        if ~isempty(find(Obs==0,1))        
            F(Obs==0) = diff(Obs==0);
        end
     elseif arg == 31
          F = EXP_mat_s;
               
    end
end

model_struct.check_sys = @(z,arg)system_check(z,arg);


%% Sampling
% This function enables a sampling of the model at hand
function Ensemble = sampling(z_sample,n_ensemble)
    Ensemble = cell(n_ensemble,1);     
    
    for i  = 1:n_ensemble          
        W = (geornd(z_sample));
        Ensemble{i} = W;
    end
    if i==1
        Ensemble = W;
    end
end
model_struct.sampling_fun = @(x,n)sampling(x,n);
end
    
 






















   
  
