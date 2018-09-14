function model_struct = BIPECM(in_data,precision)

%%%%%%%%%% Progress: tutto da scrivere. Starting point: confrontare BIPWCM
%%%%%%%%%% con DCBIPWCM( l'idea sarebbe di prendere òle phi dauno di questi
%%%%%%%%%% modelli e le psi dal bipcm(binario)). DA VERIFICARE IN FONDO!!
%%%%%%%%%% Optimization: provare versione con doppio
%%%%%%%%%% approccio(esatto sulle phi e approx sulle psi), confrontare con
%%%%%%%%%% esatto su entrambe (in tal caso tenere conto di limiti su
%%%%%%%%%% entrambi i gruppi di variabili).
      
%   This function contains all the functions, optimization options, and
%   description of the input data, required to estimate and sample 
%   the network ensemble known as:

%--------------------------------------------------------------------
%-------------------------BIPECM------------------------------------
%--------------------------------------------------------------------
%---------BIPARTITE  FIXED STRENGTH and DEGREE SEQUENCES ------------------------- 
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
sr = in_data{1};       sc = in_data{2};
kr = in_data{3};       kc = in_data{4};
n_row = length(sr);    n_col = length(sc);

% Model Name
model_struct.name = 'BIPECM';           
model_struct.in_data = in_data; 
model_struct.estimation = 'Numerical';
model_struct.is_unipartite = false;
model_struct.is_bipartite = true;
model_struct.is_weighted = true;
%% Field for the values of model's parameters
% the parameters need to be numerically estimated
model_struct.parameters =  [];          
  
 
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
    psi_sc = z(n_row+n_col+1:2*n_row+n_col);
    gamma_sc = z(2*n_row +n_col+1:2*n_row+2*n_col);
    psi_gamma_sc = psi_sc*gamma_sc';
    phi_tms_xi_sc = phi_sc*xi_sc';
    
    %Compute the expected values for the constraints 
    EXP_mat_s = (phi_tms_xi_sc)./(1- (phi_tms_xi_sc)) - ((phi_tms_xi_sc).*(1-psi_gamma_sc))./( 1 - (phi_tms_xi_sc).*(1-psi_gamma_sc)) ;
    EXP_mat_bin = (psi_gamma_sc .*phi_tms_xi_sc)./(1+phi_tms_xi_sc.*(psi_gamma_sc -1)) ;
    Exp=[sum(EXP_mat_s,2);sum(EXP_mat_s)';sum(EXP_mat_bin,2);sum(EXP_mat_bin)';];
    Obs=[sr;sc;kr;kc];
    
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

    function opt_out= optim_r()  
       
        
        %% Optimization Starting Point
% Choose a point in the available parameters space to start the
% loglikelihood optimization

%-----Starting point from BIPWCM------------------
indata_dcbipwcm = in_data; indata_dcbipwcm{3} = sum(in_data{3}); indata_dcbipwcm{4}=[];
tmp_DCBIPWCM = Max_Entr_Nets('DCBIPWCM',indata_dcbipwcm,10^-3);

phi_0 = tmp_DCBIPWCM.parameters(1:n_row);
xi_0 = tmp_DCBIPWCM.parameters(1+n_row:n_row + n_col);
 
uniform_psi_0 = sqrt( tmp_DCBIPWCM.parameters(1+n_row + n_col) );
psi_0 = ones(n_row,1)*uniform_psi_0;
gamma_0 = ones(n_col,1)*uniform_psi_0;

% 
% indata_bipcm{1} = in_data{3}; indata_bipcm{2} = in_data{4};
% 
% tmp_BIPCM = Max_Entr_Nets('BIPCM',indata_bipcm,10^-3);
% 
% psi_0 = tmp_BIPCM.parameters(1:n_row);
% gamma_0 = tmp_BIPCM.parameters(1+n_row:n_row + n_col);
%  

%-----Starting point from DCBIPWCM and uniform degrees------------------

%-----Starting point from DCBIPWCM and semiuniform (only one group) degrees-

% parameters vector.
z_start = [phi_0 ; xi_0;psi_0;gamma_0];

model_struct.opt.start_point = z_start;
model_struct.opt.par_store = z_start;

        n_iter = 1000;
        n_subiter = 1000;
        phi_opt = model_struct.opt.start_point(1:n_row);
        xi_opt  = model_struct.opt.start_point(n_row+1:n_row+n_col);
        psi_opt =  model_struct.opt.start_point(n_row+n_col+1:2*n_row+n_col);
        gamma_opt =  model_struct.opt.start_point(2*n_row +n_col+1:2*n_row+2*n_col);
         
        % Functions that define the constraints system
        Exp_X_bin = @(ph,xi,ps,ga)( (ps*ga').*(ph*xi')./(1 - (ph*xi').*(1-(ps*ga')))  );
        Exp_X = @(ph,xi,ps,ga)( (ph*xi')./(1- (ph*xi')) - ((ph*xi').*(1-(ps*ga')))./( 1 - (ph*xi').*(1-(ps*ga')))  );
       
        %    exp_X_bin = Exp_X_bin(phi_opt,xi_opt,psi_opt,gamma_opt)
%           phi_opt*xi_opt'
%          Exp_X(phi_opt,xi_opt,psi_opt,gamma_opt)
              
counter = 0;
rel_max_preced = 100;
rel_max_1 = 100;
         for t_opt = 1:n_iter-1
            t_opt           ;
           
            % Phi Xi part
            % calcola il piu' grande dei parametri dell altro gruppo
            xi_max_t = max(xi_opt);  
            for r = 1:n_row
                
                tmp = (sr(r)/n_col)/(1+(sr(r)/n_col));
                %il 10^-9 va sostituito da  qyualcosa che segua la scala
                %della rete
                bnd_ch_sign = [tmp/xi_max_t, 1/(xi_max_t+10^-12)];            
                new_phi_r = fzero( @(x)(sr(r) - sum(Exp_X(x,xi_opt,psi_opt(r),gamma_opt)) ) ,bnd_ch_sign);
                phi_opt(r) = new_phi_r;
            end
            phi_max_t = max(phi_opt);
            for c=1:n_col
                %c
                tmp = (sc(c)/n_row)/(1+sc(c)/n_row);
                bnd_ch_sign = [tmp/phi_max_t, 1/(phi_max_t+10^-12)];              
                new_xi_c = fzero( @(x)(sc(c) - sum(Exp_X(phi_opt,x,psi_opt,gamma_opt(c))) ) ,bnd_ch_sign);
                xi_opt(c) = new_xi_c;
            end
            
            %Approximate Psi gamma part
            for t_opt_sub=1:n_subiter
              
                t_opt_sub;
            for r = 1:n_row
                %r
                new_psi_r = kr(r)/ (sum( (phi_opt(r)*(xi_opt.*gamma_opt))./(phi_opt(r)*psi_opt(r)*(xi_opt.*gamma_opt) +1 - phi_opt(r)*xi_opt  )   ));
                psi_opt(r) = new_psi_r;
            end
            for c=1:n_col
                %c
                new_gamma_c = kc(c)/(sum( (xi_opt(c)*(phi_opt.*psi_opt))./(xi_opt(c)*gamma_opt(c)*(phi_opt.*psi_opt)+1 - phi_opt*xi_opt(c)  )   ));
                gamma_opt(c) = new_gamma_c;      
            end
            exp_X_bin = Exp_X_bin(phi_opt,xi_opt,psi_opt,gamma_opt);
            rel_max_2 = max(abs([(sum(exp_X_bin)'-kc)./kc; (sum(exp_X_bin,2)-kr)./kr; ] ));
            
            if ((rel_max_1<precision))
                if rel_max_2 < precision
                    break
                end
            else
                if rel_max_2 <10^-3
                    break
                end       
            end
            end
           exp_X = Exp_X(phi_opt,xi_opt,psi_opt,gamma_opt);
           
           %abs_max(t_opt) = max(abs([sum(exp_X)'-sc; sum(exp_X,2)-sr ; exp_L - L]));
           rel_max_1 = max(abs([(sum(exp_X)'-sc)./sc; (sum(exp_X,2)-sr)./sr;  ] ));         
           rel_max = max([rel_max_1;rel_max_2]);
          
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
        
      
     
        opt_out.par = [phi_opt;xi_opt;psi_opt;gamma_opt]; 
        opt_out.iter_num = t_opt; 

         
    end


model_struct.optim_run = @optim_r;

%% Sampling
% This function enables a sampling of the model at hand

function Ensemble = sampling(z,n_ensemble)
    Ensemble = cell(n_ensemble,1);       
    phi_samp = z(1:n_row)*z(n_row+1:n_row+n_col)';
    
    del_samp = z(n_row+n_col+1:2*n_row + n_col)*....
            z(2 * n_row + n_col +1: 2*n_row + 2*n_col)';
     
        W_init = zeros(n_row,n_col,'single');
       
        bin_par = (del_samp.*phi_samp)./(1- phi_samp + phi_samp.*del_samp)  ;
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
    
 






















   
  
