function model_struct =  EMECAPM(in_data,precision)
%%%%%%%%%% Progress: Funziona bene.Sampling
%%%%%%%%%% implementato, da chiarire la lentezza della convergenza sui
%%%%%%%%%% vincoli dalle media campionarie


%   This function contains all the functions, optimization options, and
%   description of the input data, required to estimate and sample 
%   the network ensemble known as:

%--------------------------------------------------------------------
%-------------------------- EMECAPM------------------------------------
%--------------------------------------------------------------------
%-----BIPARTITE FIXED DEGREE SEQUENCES AND CAPM EXPECTED WEIGHTS------------------ 
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
sr = in_data{1};        sc = in_data{2}; 
kr = in_data{3};        kc = in_data{4};
n_row = length(sr);    n_col = length(sc);


% Model Name and details
model_struct.name = 'EMECAPM';           
model_struct.in_data = in_data; 
model_struct.estimation = 'Numerical';
model_struct.is_unipartite = false;
model_struct.is_bipartite = true;
model_struct.is_weighted = true;

%% Field for the values of model's parameters
% the parameters need to be numerically estimated
model_struct.parameters =  [];          

%
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
    psi_sc = z(1:n_row,n_col+1) ;
    gamma_sc = z(n_row+1,1:n_col);

    %X CAPM
    X_c = sr * sc'./(sum(sr));   
    
    %Compute the expected values for the constraints, the capm matrix with
    %an additional row in the end for the density
    psi_gam_sc = (psi_sc*gamma_sc);
    Exp_mat_logic =  (psi_gam_sc.*phi_sc)./(1-phi_sc + phi_sc.*psi_gam_sc);
    Exp_mat_s =  phi_sc./(1-phi_sc) + (phi_sc.*(psi_gam_sc-1))./(1 + phi_sc.*(psi_gam_sc-1) );
    
    Exp= [[Exp_mat_s, sum(Exp_mat_logic,2)] ; sum(Exp_mat_logic),0]   ;
    Obs=[[X_c, kr] ; [kc',0]];
      
    if (arg == 0 )
        diff = Obs - Exp;
        F = diff;
    elseif (arg == 1 )
        diff = Obs - Exp;
        F = abs(diff)./Obs;
        if ~isempty(find(Obs==0,1))        
            F(Obs==0) = diff(Obs==0);
        end
       %max(max(F))       
     elseif arg == 31
          F = Exp_mat_s;               
    end
end

model_struct.check_sys = @(z,arg)system_check(z,arg);    
    
%% Run Optimization 
% parameters' estimation for EMECAPM is Partially analytical, partially
% numerical


   


    function opt_out = optim_r() 
        
        
        % Optimization Starting Point
%As a starting point for psi gamma take the sqr root of delta in DCMECAPM
%for L =sum(kr)


% appropriately convert the in data
 indata_dcmecapm = cell(1,3);
 indata_dcmecapm{1,1} = in_data{1,1};
  indata_dcmecapm{1,2} = in_data{1,2};
 indata_dcmecapm{1,3} = sum(kr);
 %As starting point I take an approx solution of the system for psi gamma
 %with numerical values for phi obtained from DCMECAPM 
 tmp_dcmecapm = Max_Entr_Nets('DCMECAPM',indata_dcmecapm,10^-3);
 phi_dcmecapm = tmp_dcmecapm.parameters(1:n_row,1:n_col);
 
    
 n_iter_start = 10000;
 rel_max_start = zeros(1,n_iter_start);
 start_value_common = 10; 
 psi_start_tmp = start_value_common* ones(n_row,1);
 gamma_start_tmp  = start_value_common* ones(n_col,1);   

    
    
% iterations of approximate solutions for the psi and gammas given phis            
        for t_opt_start = 1:n_iter_start-1            
            new_psi_s = kr ./( sum( phi_dcmecapm.*(ones(n_row,1)*gamma_start_tmp')./(1-phi_dcmecapm + phi_dcmecapm.*(psi_start_tmp*gamma_start_tmp')  ),2) )  ;               
            psi_start_tmp = new_psi_s;                           
            new_gamma_c = kc ./( sum( phi_dcmecapm.*(psi_start_tmp*ones(1,n_col))./(1-phi_dcmecapm + phi_dcmecapm.*(psi_start_tmp*gamma_start_tmp')  ),1)' )  ;         
            gamma_start_tmp  = new_gamma_c;                           
            exp_X_start = ((psi_start_tmp*gamma_start_tmp').*phi_dcmecapm)./(1-phi_dcmecapm + phi_dcmecapm.*( psi_start_tmp*gamma_start_tmp'));
            rel_max_start(t_opt_start) = max(abs([(sum(exp_X_start)'-kc)./kc; (sum(exp_X_start,2)-kr)./kr ]));            
            if rel_max_start(t_opt_start)< precision
                break
            end           
        end
%start_iter = t_opt_start
%start_time = toc


 z_start  = [psi_start_tmp;gamma_start_tmp];
 
model_struct.opt.start_point = z_start;
model_struct.opt.par_store = z_start;
     
        %Here is the numerical part: estimate psi and gamma from iterative
        %solution of approximated equations
        n_iter = 10000; 
        rel_max = zeros(1,n_iter);
        X_c = sr * sc'./(sum(sr)); %X CAPM   
        
        %a part of the function that generates the expected binary matrix 
        fun_den = @(x,y)((x*y').*X_c - (x*y') +  sqrt( ((x*y').^2).* (X_c -1).^2 + 4*X_c.*(x*y')  )  ); 
        %Expected binary matrix from a psi gamma and X_c
       % fun_exp_bin = @(x,y)( (x*y')./(1-x*y') .* (2*X_c./fun_den(x,y) -1) ); 
        %phi parameter as funciton of psi gamma and X_c
        fun_phi = @(x,y)( (2*X_c - fun_den(x,y))./(2*X_c .*(1 -(x*y') ))   );
         
        % import the starting point
        psi_opt = model_struct.opt.start_point(1:n_row);
        gamma_opt  = model_struct.opt.start_point(n_row+1:n_row+n_col);
        %solve the system by iteratively solving each equation   
        
        
counter = 0;
rel_max_preced = 100;
        for t_opt = 1:n_iter-1
           % t_opt
            vec_phi = fun_phi(psi_opt,gamma_opt);
                  
            %approximate solution
            new_psi = kr./(sum( (vec_phi.*(ones(n_row,1)*gamma_opt') )./(1-vec_phi + vec_phi.*(psi_opt*gamma_opt')),2) );        
            psi_opt  = new_psi;               
            
            vec_phi = fun_phi(psi_opt,gamma_opt);
            new_gamma = kc./(sum( (vec_phi.*(psi_opt*ones(1,n_col)) )./(1-vec_phi + vec_phi.*(psi_opt*gamma_opt'))  ) )'; 
            gamma_opt = new_gamma;
                 
            
            vec_phi = fun_phi(psi_opt,gamma_opt);
            exp_X = ((psi_opt*gamma_opt').*vec_phi)./(1-vec_phi + vec_phi.*( psi_opt*gamma_opt'));
            
            rel_max  = max(abs([(sum(exp_X)'-kc)./kc; (sum(exp_X,2)-kr)./kr ]));
           % rel_max(t_opt)
%             if t_opt>1
%                 improv = abs((rel_max(t_opt)-rel_max(t_opt-1))) ;
%             end

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
        

        
        % analytically obtain the matrix of phis as function of psi, gamma
        % and X_c 
        phi_opt = vec_phi;
        
     
        
        % the parameters are grouped in one matrix. The first n_row rows
        % and n_col_columns are the phi parameters, the last row contains
        % the psi and the last column contains  gamma
       
        
        opt_out.par =  [[phi_opt, psi_opt]; [gamma_opt',0] ];   
        opt_out.iter_num = t_opt;  
         
    end

model_struct.optim_run = @optim_r;
%% Sampling
% This function enables a sampling of the model at hand

function Ensemble = sampling(z,n_ensemble)
    Ensemble = cell(n_ensemble,1);       
    phi_samp = z(1:n_row,1:n_col);    
    gamma_samp = z(n_row+1,1:n_col);
    psi_samp  = z(1:n_row,n_col+1);
        W_init = zeros(n_row,n_col,'single');
    psi_gamma_samp = psi_samp*gamma_samp;   
%     size(psi_gamma_samp)
%     size(phi_samp)
        bin_par = (psi_gamma_samp.*phi_samp)./(1-phi_samp + psi_gamma_samp.*phi_samp )  ;
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
    
 






















   
  
