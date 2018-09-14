function mat = list_to_mat_uni_und(list)
%converti lista in matrice

L = length(list(:,1));
N = max(max(list(:,1)),max(list(:,2)));



if N>10^3 % large matrix
    if L<N^1.5  %check sparcity
    mat = spalloc(N,N,round(N^1.5));

    else error
        
    end
else  
    
    mat = zeros(N);
end
for i=1:L
    
    mat(list(i,1),list(i,2)) = mat(list(i,1),list(i,2)) + list(i,3);

    mat(list(i,2),list(i,1)) = mat(list(i,2),list(i,1)) + list(i,3);
end

end