function df=numerical_jacobian(f,x) 
n=length(x); 
E=speye(n); 
e=1e-4;
df = zeros(n,n);
for i=1:n 
 df(:,i)=(f(x+e*E(:,i))-f(x-e*E(:,i)))/(2*e); % zentraler Differenzenquotient 
end 
end