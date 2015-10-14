data=xlsread('C:/Users/caped_crusader/Bayesian Data

Analyses/Assignments/binary.xlsx');

%Hyperparameters

V_pr=[0.25 0 0 0;0 100 0 0;0 0 1.5  0; 0 0 0 1]; %%Values assigned based 

on sd of data of individual var

b_pr=[0.5;1/600;1/3.4;1/2.5];       %Prior beta from reciprocal of average

N=400;

k=4;        %No. of parameters

y=data(1:400,1);   %Dependent Variable

x=ones(400,4);      %Covariate X1

x(1:400,2:4)=data(1:400,2:4); %Covariates

b_ols=(((transpose(x))*x)^-1)*(transpose(x)*y);   %OLS beta

xsq=(transpose(x)*x);

G=10000;

ystar=ones(400,1);          %Latent variable is y*

sum_b_exp=0;

for i=1:G

   V_posterior=(V_pr^-1 + xsq)^-1;

   b_posterior=V_posterior*(V_pr^-1*b_pr + transpose(x)*ystar);     

%posterior beta

   for j=1:N            %Latent variable generation from random draw

       x_i(1,1:4)=x(j,1:4);

       if y(j,1)==1     %Right truncated

           while ystar>=0

           ystar(i,1)=normrnd((x_i)*b_posterior,1);

           end

       end

       if y(j,1)==0     %Left truncated

           while ystar<0

            ystar(i,1)=normrnd((x_i)*b_posterior,1);

           end

       end         

   end

   b=transpose(mvnrnd(b_posterior,V_posterior));    %Beta draw

   if i>=0.1*G

   sum_b_exp=sum_b_exp+b;           %Acceptance after 1000 burns

   end

end

b_exp=sum_b_exp/(0.9*G)  %Average of beta

%Inference

%%Since variables cannot be interpreted based on absolute value but signs,

%%hence the concluding remark is "All variables have positive impact on

%%admission"