function model_struct =  BIPCM(in_data,precision)      
%%%%%%%%%% Progress: funzionante. si potrebbe migliorare scegliendo
%%%%%%%%%% starting  point in modo piu' sensato .
      
%   This function contains all the functions, optimization options, and
%   description of the input data, required to estimate and sample 
%   the network ensemble known as:

%--------------------------------------------------------------------
%-------------------------- BIPCM------------------------------------
%--------------------------------------------------------------------
%---------------BIPARTITE FIXED  DEGREE SEQUENCES------------------------- 
%%
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
kr = in_data{1};       kc = in_data{2};
n_row = length(kr);    n_col = length(kc);
 

% Model Name
model_struct.name = 'BIPCM';           
model_struct.in_data = in_data; 
model_struct.estimation = 'Numerical';
model_struct.is_unipartite = false;
model_struct.is_bipartite = true;
model_struct.is_weighted = false;
%% Field for the values of model's parameters
% the parameters need to be numerically estimated
model_struct.parameters =  [];    

%% Optimization Starting Point
% Choose a point in the available parameters space to start the
% loglikelihood optimization

%Starting point that fixes the columns degree sequence and a uniform
%row deg sequence (because is analytical)
L = sum(kr); 
psi_0 = 1000;
tmp_0 = kc/n_row;
ind_0 = round(tmp_0) == 1;
gamma_0 =   ( (tmp_0) ./((1- (tmp_0) )) )./psi_0 ;
gamma_0(ind_0) = 100;
z_start = [psi_0* ones(n_row,1);gamma_0 ] ;
model_struct.opt.start_point = z_start;
model_struct.opt.par_store = z_start;

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
    


    %decompose the parameters vector
    psi_sc = z(1:n_row);
    gamma_sc = z(n_row +1 : n_row + n_col);
        
    phi_tms_phi_sc = psi_sc*gamma_sc';
    %Compute the expected values for the constraints 
    EXP_mat_s = phi_tms_phi_sc./(1+phi_tms_phi_sc);
    Exp=[sum(EXP_mat_s,2); sum(EXP_mat_s,1)'];   
    Obs= [kr;kc];    
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
    
%% Run Optimization 
% for CM the iterative solution of and approximate version of the system
% yelds the correct parameters
    function opt_out = optim_r()  
       
        n_iter = 10000;               
        psi_opt = model_struct.opt.start_point(1:n_row);
        gamma_opt  = model_struct.opt.start_point(n_row+1:n_row+n_col);
             
        for t_opt = 1:n_iter-1                     
            %approximate solution
            new_psi = kr./( sum( (ones(n_row,1)*gamma_opt' )./(1 + (psi_opt*gamma_opt') ),2) )  ;               
            psi_opt = new_psi;       
                  
            new_gamma = kc./( sum( (psi_opt*ones(1,n_col) )./(1 + (psi_opt*gamma_opt') ) )' )  ; 
            gamma_opt = new_gamma;        
            
            exp_X = (psi_opt*gamma_opt')./(1+psi_opt*gamma_opt');
            %abs_max(t_opt) = max(abs([sum(exp_X)'-kc; sum(exp_X,2)-kr ]));
            rel_max = max(abs([(sum(exp_X)'-kc)./kc; (sum(exp_X,2)-kr)./kr ]));
            if rel_max< precision
                break
            end           
        end     
        
        opt_out.par = [psi_opt;gamma_opt]; 
        opt_out.iter_num = t_opt;  
         
    end


model_struct.optim_run = @optim_r;


        

%% Sampling
% This function enables a sampling of the model at hand

function Ensemble = sampling(z,n_ensemble)
    Ensemble = cell(n_ensemble,1);       
    psi_samp = z(1:n_row);    
    gamma_samp = z(n_row+1:n_row + n_col);
     
        W_init = zeros(n_row,n_col,'single');
       
        par_prod = psi_samp * gamma_samp';
        bin_par = par_prod ./(1+par_prod) ;
        for i = 1:n_ensemble
            W=W_init;
           %Bernoulli trial 
            W_bin = logical(binornd(ones(n_row,n_col),bin_par));
            Ensemble{i} = W_bin;

        end
end

model_struct.sampling_fun = @(x,n)sampling(x,n);
    
 

end




















   
  
