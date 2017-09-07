function PA=funcesp(A,Z,Fs,Tn)

t=(0:length(A)-1)*1/Fs;

for n=1:length(Tn)

    T=Tn(n);
    wn=2*pi/T;


    H=tf([0 0 -1],[1 2*Z*wn wn^2]);
    H2=tf([-1 0 0],[1 2*Z*wn wn^2]);
    
    De=lsim(H,A,t);
    Ac=lsim(H2,A,t);

    D(n)=max(abs(De));
    EA(n)=max(abs(Ac));

    PA(n)=D(n)*(2*pi/Tn(n))^2/9.81;
end

