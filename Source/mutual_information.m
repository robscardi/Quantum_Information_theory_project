function [I_x_y] = mutual_information(prob_x,prob_matrix) 


%% Setup matrix

%Esempio:
%p = 1
%err = 0.4
%prob_x = [p 1-p]
%prob_matrix = [(1-err) err  0; 0 err (1-err)]

prob_x = prob_x;
prob_matrix = prob_matrix;
prob_y = prob_x*prob_matrix;

%% Calcolo H(y):

H_y = 0;
[c,r] = size(prob_y);

for y = 1:r
    a = -prob_y(y)*log2(prob_y(y));
    if isnan(a) == 1
            a = 0;
    end 
       H_y = H_y + a; %H_y = -∑P(y)log2(P(y))
end


%% Calcolo H(y|x):

H_y_x = 0;
[c,r] = size(prob_matrix);

for x = 1:c
    temp = 0;
    for y = 1:r 
        t = - prob_matrix(x,y)*log2(prob_matrix(x,y));
        if isnan(t)
            t = 0;
        end 
        temp = temp + t;
    end 
    temp = temp*prob_x(x);
    H_y_x = H_y_x + temp;
end 
    
%% Calcolo I(x,y) = H(y) - H(y|x):

I_x_y = H_y - H_y_x;
disp(['The mutual information is I(x,y)=', num2str(I_x_y)])

end 