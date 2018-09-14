function model_struct =  DCMECAPM(in_data,precision)
%%%%%%%%%% Progress: manca solo il sampling che non sembra funzionare
      
%   This function contains all the functions, optimization options, and
%   description of the input data, required to estimate and sample 
%   the network ensemble known as:

%--------------------------------------------------------------------
%-------------------------- DCMECAPM------------------------------------
%--------------------------------------------------------------------
%----------------BIPARTITE FIXED DENSITY AND CAPM EXPECTED WEIGHTS------------------------- 
%%
%   INPUT: in_data  is a cell array that contains the strength sequences, and desity.
%                   The first elementis a column with rows strengths, the
%                   second one a column with the columns strengths, the
%                   third is the network desity
%           
%   OUTPUT: model_struct    is a structure with a number of fields 
%                           described in the following. Each field refers
%                           to informations, funcitons or options specific
%                           for the particular model. 

%   To be called only at the beginning of the main script Max_Entr_Nets

%%
%Variables used by different sub functions defined in the following 
sr = in_data{1};       sc = in_data{2}; L = in_data{3};
n_row = length(sr);    n_col = length(sc);


% Model Name and details
model_struct.name = 'DCMECAPM';           
model_struct.in_data = in_data; 
model_struct.estimation = 'Numerical';
model_struct.is_unipartite = false;
model_struct.is_bipartite = true;
model_struct.is_weighted = true;

%% Field for the values of model's parameters
% the parameters need to be numerically estimated
model_struct.parameters =  [];          
 
  

%% Optimization Starting Point
 % no need for starting points

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
    phi_sc = z(1:n_row,1:n_col);
    delta_sc = z(n_row+1,1) ;

    %X CAPM
    X_c = sr * sc'./(sum(sr));   
 
    %Compute the expected values for the constraints, the capm matrix with
    %an additional row in the end for the density
    Exp_mat_logic =  (phi_sc*delta_sc)./(1-phi_sc*(1-delta_sc)) ;  
    Exp_mat_s =phi_sc./(1-phi_sc) - (phi_sc*(1-delta_sc) )./(1-phi_sc*(1-delta_sc));
    Exp=[Exp_mat_s; [sum(sum(Exp_mat_logic)),zeros(1,n_col-1)]] ;
    
     
    Obs=[X_c; [L,zeros(1,n_col-1)]];
   
   
    if (arg == 0 )
        diff = Obs - Exp;
        F = diff;
    elseif (arg == 1 )
        diff = Obs - Exp;
        F = abs(diff)./Obs;
        if ~isempty(find(Obs==0,1))        
            F(Obs==0) = diff(Obs==0);
        end
      % max(max(F));
       
     elseif arg == 31
          F = Exp_mat_s;
               
    end
end

model_struct.check_sys = @(z,arg)system_check(z,arg);




    
    
%% Run Optimization 
% parameters' estimation for DCMECAPM is almost completely analytical

    function opt_out = optim_r()        
        
        
        X_c = sr * sc'./(sum(sr)); %X CAPM
        %  the relation between X delta and L
        est_delta = @(del)(  (L -.....
                            sum(sum(....
                            (del/(1-del))* (2*X_c./(del*X_c -del +sqrt(del^2*(X_c-1).^2 + 4*del*X_c) )  - 1 )....
                            )))  /L  );
                          
        bnd_ch_sign = [ 10^(-10),10^7];
        delta_opt =  fzero( @(x)(est_delta(x)) ,bnd_ch_sign);
        %phi is a matrix of parameters, that, toghether with delta(1 par), defines
        %the model
        phi_num = (X_c*(2-delta_opt)+delta_opt - sqrt(delta_opt^2*(X_c -1).^2 +4*delta_opt*X_c) );
        phi_den = (2*X_c*(1 -delta_opt));
        phi_opt = phi_num./phi_den;
      
       % par is the matrix of phis with a last row just for delta
        opt_out.par = [phi_opt; [delta_opt, zeros(1,n_col-1)]] ; 
        opt_out.iter_num = 1; 
         
    end

model_struct.optim_run = @optim_r;

%% Sampling
% This function enables a sampling of the model at hand

function Ensemble = sampling(z,n_ensemble)
    Ensemble = cell(n_ensemble,1);       
    phi_samp = z(1:n_row,1:n_col);    
    del_samp = z(n_row+1,1);
     
        W_init = zeros(n_row,n_col,'single');
       
        bin_par = (del_samp*phi_samp)./(1- phi_samp + phi_samp*del_samp)  ;
        for i = 1:n_ensemble
            W=W_init;
           %Bernoulli trial preliminare
            W_bin = logical(binornd(ones(n_row,n_col),bin_par));

            %per gli elementi di matrice estratti W_bin~=0 estrai
            %da distribuzione geometrica con probabilità 
            geom_par = 1-  phi_samp;
           
            W(W_bin) =single( 1 + geornd(geom_par(W_bin)));
            Ensemble{i} = W;
            %sum_ens = sum_ens + W;
        end
end
model_struct.sampling_fun = @(x,n)sampling(x,n);
end
    
 






















   
  
