%Function that computes the required statistics from a given adjacency
%matrix (could be unipartite, bipartire, weighted), as required for the
%given method (e.g. BIPCM Binary Configuration model requires the degree
%sequences)
function [indata] = indata_from_matrix_Nets(method,X)

% Returns the appropriate input data for a given network ensemble (specified by "method")
switch upper(method)
    %Return the degree sequences
     case 'BIPCM' % Bipartite configuration model
     
         X_bin = logical(X);
         k_row = sum(X_bin,2);
         k_col = sum(X_bin)';
         
         indata = cell(1,2);
         indata{1,1} = k_row;
         indata{1,2} = k_col;
         
       
     %Return the strenght sequences
     case {'BIPWNB', 'BIPWCM','MECAPM','PPCAPM','BIPWMB'}  
     
         s_row = sum(X,2);
         s_col = sum(X)';
         
         indata = cell(1,2);
         indata{1,1} = s_row;
         indata{1,2} = s_col;
    
    % Return the strength sequence for unipartite networks         
    case {'WCM','WNB'} 
     
         s_row = full(sum(X,2));
         
         s_row(s_row == 0) = []; 
         indata = cell(1,2);
         indata{1,1} = s_row;
         
             
    %return strenth sequences and and density
    case {'DCMECAPM','DCBIPWCM'}
         s_row = sum(X,2);
         s_col = sum(X)';
         
         L = sum(sum(logical(X)));
         
         indata = cell(1,2);
         indata{1,1} = s_row;
         indata{1,2} = s_col;
         indata{1,3} = L;
         
    % return strenght and degree sequences
    case {'EMECAPM','BIPECM'}
         X_bin = logical(X);
         s_row = sum(X,2);
         s_col = sum(X)';
         k_row = sum(X_bin,2);
         k_col = sum(X_bin)';
         L = sum(sum(logical(X)));
         
         indata = cell(1,2);
         indata{1,1} = s_row;
         indata{1,2} = s_col;
         indata{1,3} = k_row;
         indata{1,4} = k_col;
end

    
end













