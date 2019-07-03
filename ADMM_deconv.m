function X=ADMM_deconv(Y,H0,c,lambda)



mu=10;
N=length(c);

Nr=64;
Nc=64;
X=zeros(Nr*Nc,1);
W=X;
D=X;
tol = 1e-4;
tol1 = sqrt(N)*tol;
tol2 = sqrt(N)*tol;
res_p = inf;
res_d = inf;
B=inv(H0'*H0+mu*eye(N));
% B=1/(N/2+mu)* (eye(N) - ones(N)*(N/4)/(1+N^2/4) ); 
Hty=H0'*Y;
t=1;
Niter=100;
while (t<Niter) && ((abs (res_p) > tol1) || (abs (res_d) > tol2)) 
% t
    X0=X;
    % Update X
%     1/2*norm(Y-H0*X)^2 + mu/2*norm(W-X-D)^2
    X=B*(mu*(W-D)+Hty);
%     1/2*norm(Y-H0*X)^2 + mu/2*norm(W-X-D)^2
    % Update W
%     1/2*norm(Y-H0*X)^2 + mu/2*norm(W-X-D)^2 + lambda/2*norm(real(ifft(c.*fft(W)))).^2
    W=real(ifft( (1./(mu +(lambda)*abs(c).^2)).*fft(mu*(X+D))));
%     1/2*norm(Y-H0*X)^2 + mu/2*norm(W-X-D)^2 + lambda/2*norm(real(ifft(c.*fft(W)))).^2
    % Update D
    D=D-(W-X);
    
        % update mu
    if mod(t,10) == 1
    % primal residue
    res_p = sqrt(norm(W-X,'fro')^2);
        % dual residue
        res_d = mu*sqrt(norm(X-X0,'fro')^2 ) ;
                % update mu
        if res_p > 5*res_d
            mu = mu*4;
            D = D/4;
            B=inv(H0'*H0+mu*eye(N));
        elseif res_d > 5*res_p
            mu = mu/4;
            D = D*4;
            B=inv(H0'*H0+mu*eye(N)); 
        end
        fprintf(' t = %f, res_p = %f, res_d = %f \n',t,res_p,res_d)
    end
    t=t+1;
end

X=reshape(X,Nr,Nc);