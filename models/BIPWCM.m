function model_struct = BIPWCM(in_data,precision)

%%%%%% Progress: Funziona bene. Resta da capire come scala su reti molto
%%%%%% grandi e se è possibile migliorare la scalabilità con miglioramento
%%%%%% dei bounds o addirittura cambiamento dell algoritmo in favore di uno
%%%%%% che utilizzi le derivate(che noi possiamo calcolare analiticamente).
%%%%%% Un cambio in WCM porterebbe miglioramenti a cascata in tutti gli
%%%%%% altri modelli.
      
%   This function contains all the functions, optimization options, and
%   description of the input data, required to estimate and sample 
%   the network ensemble known as:

%--------------------------------------------------------------------
%--------------------------BIPWCM------------------------------------
%--------------------------------------------------------------------
%---------------BIPARTITE FIXED STRENGTH SEQUENCES------------------------- 
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
n_row = length(sr);    n_col = length(sc);
n_tot = n_row+n_col;

% Model Name
model_struct.name = 'BIPWCM';           
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
    xi_sc = z(n_row+1:n_tot);    
    phi_tms_xi_sc = phi_sc*xi_sc';
    %Compute the expected values for the constraints 
    EXP_mat_s = phi_tms_xi_sc./(1-phi_tms_xi_sc);
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

 
    
    
%% Run Optimization 

    function opt_out = optim_r()  
        
        % Optimization Starting Point
% Choose a point in the available parameters space to start the
% loglikelihood optimization

%-----Starting point that produces a uniform mean matrix-------------------
S = sum(sr);
Sl  = S/(n_col*n_row);
%  use the phi xi version of the parameters
phi_or_xi0 =  sqrt((Sl/((1+Sl))));

%uniform phi
phi_0 = phi_or_xi0*ones(n_row,1) ;
% the exact ones that follow from phi_0 for xi_0
xi_0  = (sc./(n_row*phi_0(1))) ./ (1 + sc/n_row)  ;

%diff = (sum((sum((phi_0*xi_0')./(1-(phi_0*xi_0')),2))) - sum(sr))./sum(sr)

z_start = [phi_0 ; xi_0];

model_struct.opt.start_point = z_start;
model_struct.opt.par_store = z_start;
       
        n_iter = 1000;
         
        phi_opt = model_struct.opt.start_point(1:n_row);
        xi_opt  = model_struct.opt.start_point(n_row+1:n_row+n_col);
        
        %phi_store = zeros(n_row,n_iter);
        %xi_store = zeros(n_col,n_iter);
        

        %The constraint associated to a single parameter
        con = @(x,y)( sum( x*y./((1)- x*y)) );
        %con_v = @(x,y)(x*y./((1)- x*y));
%         sr
%         sc
        
   rel_max_preced = 100;
   counter = 0;
        for t_opt = 1:n_iter-1
            %t_opt
            
            % calcola il piu' grande dei parametri dell altro gruppo
            xi_max_t = max(xi_opt);
            %xi_min_t = min(xi_opt);
           % t_opt
            for r = 1:n_row
                
                %r
                tmp = (sr(r)/n_col)/(1+(sr(r)/n_col));
                %il 10^-9 va sostituito da  qyualcosa che segua la scala
                %della rete
                 bnd_ch_sign = [tmp/xi_max_t, 1/(xi_max_t+10^-12)];
               
                new_phi_r = fzero( @(x)(sr(r) -con(x,xi_opt)) ,bnd_ch_sign);
                phi_opt(r) = new_phi_r;

            end

            %phi_store(:,t_opt+1) = phi_opt;

            phi_max_t = max(phi_opt);
            %phi_min_t = min(phi_opt);

            for c=1:n_col
                %c
                tmp = (sc(c)/n_row)/(1+sc(c)/n_row);
                bnd_ch_sign = [tmp/phi_max_t, 1/(phi_max_t+10^-12)];
                new_xi_c = fzero( @(x)(sc(c) -con(x,phi_opt)) ,bnd_ch_sign);
                xi_opt(c) = new_xi_c;

            end
            
            exp_X = (phi_opt*xi_opt')./(1-phi_opt*xi_opt');
           rel_max  = max(abs([(sum(exp_X)'-sc)./sc; (sum(exp_X,2)-sr)./sr ]));

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
        opt_out.par = [phi_opt;xi_opt]; 
        opt_out.iter_num = t_opt; 
               
    end
model_struct.optim_run = @optim_r;

%% Sampling
% This function enables a sampling of the model at hand
function Ensemble = sampling(z,n_ensemble)
    Ensemble = cell(n_ensemble,1);       
    phi_samp = z(1:n_row);
    xi_samp = z(n_row+1:n_tot);
    phi_tms_xi = phi_samp*xi_samp';
    % P is the matrix of parameters of the geometric distribution
    P =(1 -phi_tms_xi);
    for l  = 1:n_ensemble
        W = single(geornd(P));
        Ensemble{l} = W;
    end   
end
model_struct.sampling_fun = @(x,n)sampling(x,n);
end
    
 






















   
  
