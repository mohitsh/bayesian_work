%%Mohit Shukla

%Q1

%Generate Y-Data series

y=zeros(100,1);

w=zeros(98,2);

y(1,1)=normrnd(0,1);

y(2,1)=normrnd(0,1);

y2=zeros(2,1);

phi1=1;

phi2=-0.5;

sum_w=zeros(2,2);

n=98;

%Generate lag y

for i=3:100

   eps=normrnd(0,1);

   y(i,1)=phi1*y(i-1,1) + phi2*y(i-2,1) + eps; 

   w(i,1)=y(i-1,1);

   w(i,2)=y(i-2,1);

   sum_w=sum_w+w(i,:)'*w(i,:);

end

y2(1,1)=y(1,1);

y2(2,1)=y(2,1);

%Calculate G,Phi hat etc

G=sum_w;

phi_hat=zeros(1,2);

error=zeros(98,1);

sum_err_sq=0;

for i=3:100

    phi_hat=phi_hat+y(i,:)*w(i,:);

    error(i-2,1)=y(i,1)-w(i,1)*phi1-w(i,2)*phi2;

    sum_err_sq=sum_err_sq+(error(i-2,1)^2);

end

%Calculate Posterior parameters

    phi_hat=phi_hat*(G^-1);

    G_Inv=G^-1;

    V_Inv=[1-phi2^2,-phi1*(1+phi2);-phi1*(1+phi2),1-phi2^2];

    V=V_Inv^-1;

    V_post_Inv=(y2'*V_Inv*y2 + sum_err_sq);

    V_post=V_post_Inv^-1;

    

    s=50000;

    s0=1000;

    df=n/2;

    phi_prev=[1,-0.5];

    iteration=1;

    sigma_prev=1;

    sum=phi_hat;

    

%Gibbs Sampling

    for i=1:s

       sigma_sq_inv=wishrnd(V_post,df);

       sigma_sq=sigma_sq_inv^-1;

       phi_draw=mvnrnd(phi_hat,sigma_sq*G_Inv);

       sum=sum+phi_draw;

       %Acceptance region of phi1 and phi2

       if (phi_draw(1,2)+phi_draw(1,1)<1) &  (-

phi_draw(1,1)+phi_draw(1,2)<1) & (phi_draw(1,2)>-1)

        V2_inv=[1-phi_draw(1,2)^2,-phi_draw(1,1)*(1+phi_draw(1,2));-

phi_draw(1,1)*(1+phi_draw(1,2)),1-phi_draw(1,2)^2];

        V1_inv=[1-phi_prev(1,2)^2,-phi_prev(1,1)*(1+phi_prev(1,2));-

phi_prev(1,1)*(1+phi_prev(1,2)),1-phi_prev(1,2)^2];

        ratio=((det(V2_inv)^0.5)*exp(-

(y2'*V2_inv*y2)/(2*sigma_prev)))/((det(V1_inv)^0.5)*exp(-

(y2'*V1_inv*y2)/(2*sigma_prev)));

       if ratio>=1

            phi_prev=phi_draw;

       end

       u=rand;

       if ratio<1 & u<ratio

            phi_prev=phi_draw;

       end

       if i>s0   %Gibbs sampling after s0 burns

           gibbs_phi(iteration,1)=phi_prev(1,1);

           gibbs_phi(iteration,2)=phi_prev(1,2);

           gibbs_sigma(iteration,1)=sigma_sq;

           sigma_prev=sigma_sq;

           iteration=iteration+1;

       end

       end

    end

%Summary Phi

mean(gibbs_phi)

median(gibbs_phi)

min(gibbs_phi)

max(gibbs_phi)

std(gibbs_phi)

%Summary Sigma

mean(gibbs_sigma)

median(gibbs_sigma)

min(gibbs_sigma)

max(gibbs_sigma)

std(gibbs_sigma)