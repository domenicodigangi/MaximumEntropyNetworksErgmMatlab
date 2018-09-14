function model_struct = DCBIPWCM(in_data,precision)

%%%%%%%%%% Progress: funziona bene. Rimane da capire l'effettiva efficienza
%%%%%%%%%% dello starting point utilizzato

      
%   This function contains all the functions, optimization options, and
%   description of the input data, required to estimate and sample 
%   the network ensemble known as:

%--------------------------------------------------------------------
%--------------------------DCBIPWCM------------------------------------
%--------------------------------------------------------------------
%---------------BIPARTITE FIXED STRENGTH SEQUENCES and Density------------------------- 
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

%% Field for the values of model's parameters
%Variables used by different sub functions defined in the following 
sr = in_data{1};       sc = in_data{2};
L = in_data{3};
n_row = length(sr);    n_col = length(sc);
n_tot = n_row+n_col;

% Model Name
model_struct.name = 'DCBIPWCM';           
model_struct.in_data = in_data; 
model_struct.estimation = 'Numerical';
model_struct.is_unipartite = false;
model_struct.is_bipartite = true;
model_struct.is_weighted = true;


model_struct.parameters =  [];          

%% Optimization Starting Point
% Choose a point in the available parameters space to start the
% loglikelihood optimization

%-----Starting point from BIPWCM------------------
indata_bipwcm = in_data; indata_bipwcm{3} = [];
tmp_BIPWCM = Max_Entr_Nets('BIPWCM',indata_bipwcm,10^-3);

phi_0 = tmp_BIPWCM.parameters(1:n_row);
xi_0 = tmp_BIPWCM.parameters(1+n_row:n_row + n_col);
 
% parameters vector. the last number is the starting value for delta (does not matter)
z_start = [phi_0 ; xi_0;1];

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
    phi_sc = z(1:n_row);
    xi_sc = z(n_row+1:n_row+n_col);    
    delta_sc = z(n_row +n_col+1);
    phi_tms_xi_sc = phi_sc*xi_sc';
    %Compute the expected values for the constraints 
    EXP_mat_s = (phi_tms_xi_sc)./(1- (phi_tms_xi_sc)) - ((phi_tms_xi_sc)*(1-delta_sc))./( 1 - (phi_tms_xi_sc)*(1-delta_sc)) ;
    Exp_L = sum(sum( (delta_sc *phi_tms_xi_sc)./(1+phi_tms_xi_sc*(delta_sc -1)) ));
    Exp=[sum(EXP_mat_s,2);sum(EXP_mat_s)';Exp_L];
    Obs=[sr;sc;L];
    
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

    function opt_out = optim_r()  
       
        n_iter = 10000;         
        phi_opt = model_struct.opt.start_point(1:n_row);
        xi_opt  = model_struct.opt.start_point(n_row+1:n_row+n_col);
        delta_opt = model_struct.opt.start_point(n_row + n_col + 1);
         
        % Functions that define the constraints system
        Exp_L = @(p,x,d)(sum(sum( d*(p*x')./(1 - (p*x')*(1-d)) )) );
        Exp_X = @(p,x,d)( (p*x')./(1- (p*x')) - ((p*x')*(1-d))./( 1 - (p*x')*(1-d))  );
       
        
        
           counter = 0;    
           rel_max_preced = 100;
         for t_opt = 1:n_iter-1
            %t_opt
            % estimate delta given phi and xi
           
%             
            bnd_ch_sign = [10^-10,10^9];
%             L - Exp_L(phi_opt,xi_opt,bnd_ch_sign(1))
%             L - Exp_L(phi_opt,xi_opt,bnd_ch_sign(2))
            delta_opt =  fzero( @(x)(L - Exp_L(phi_opt,xi_opt,x) ) ,bnd_ch_sign);            
            % find the largest xi parameter
            xi_max_t = max(xi_opt);
            for r = 1:n_row
                %r
                tmp = (sr(r)/n_col)/(1+(sr(r)/n_col));
                %il 10^-9 va sostituito da  qyualcosa che segua la scala
                %della rete
                bnd_ch_sign = [tmp/xi_max_t, 1/(xi_max_t+10^-12)];                             
                new_phi_r = fzero( @(x)(sr(r) - sum(Exp_X(x,xi_opt,delta_opt)) ) ,bnd_ch_sign);
                phi_opt(r) = new_phi_r;
            end
            phi_max_t = max(phi_opt);            
            for c=1:n_col
                %c
                tmp = (sc(c)/n_row)/(1+sc(c)/n_row);
                bnd_ch_sign = [tmp/phi_max_t, 1/(phi_max_t+10^-12)];
               
                
                new_xi_c = fzero( @(x)(sc(c) - sum(Exp_X(phi_opt,x,delta_opt)) ) ,bnd_ch_sign);
                xi_opt(c) = new_xi_c;

            end

            %xi_store(:,t_opt+1) = xi_opt;

           exp_X = Exp_X(phi_opt,xi_opt,delta_opt);
           exp_L = Exp_L(phi_opt,xi_opt,delta_opt);
          
           rel_max  = max(abs([(sum(exp_X)'-sc)./sc; (sum(exp_X,2)-sr)./sr; (exp_L -L)/L ]));

         
            if (abs(rel_max_preced- rel_max)<10^-10)
                counter =counter+1;
                if counter>50;
                    break
                end
            end

           if rel_max < precision
              break
           end
           rel_max_preced = rel_max;
         end
     
        opt_out.par =[phi_opt;xi_opt;delta_opt]; 
        opt_out.iter_num = t_opt; 
   
    end


model_struct.optim_run = @optim_r;

%% Sampling
% This function enables a sampling of the model at hand

function Ensemble = sampling(z,n_ensemble)
    Ensemble = cell(n_ensemble,1);       
    phi_samp = z(1:n_row)*z(n_row+1:n_row+n_col)';    
    del_samp = z(n_row+n_col+1);
     
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

%% Analytical Expected Matrix
 function out_exp_m = exp_m(z)  
       
 
        phi_exp = z(1:n_row);
        xi_exp  = z(n_row+1:n_row+n_col);
        delta_exp = z(n_row + n_col + 1);
         
        % Functions that define the constraints system
        Exp_L = @(p,x,d)(sum(sum( d*(p*x')./(1 - (p*x')*(1-d)) )) );
        Exp_X = @(p,x,d)( (p*x')./(1- (p*x')) - ((p*x')*(1-d))./( 1 - (p*x')*(1-d))  );
       
        out_exp_m = Exp_X(phi_exp,xi_exp,delta_exp);
             
 end
end
    
 






















   
  
