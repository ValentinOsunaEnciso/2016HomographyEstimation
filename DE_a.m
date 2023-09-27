function Qt  = DE_a(Pt, M, d, x_low, x_high,X)
F=0.25;         %Factor de escalamiento
Cr=0.8;         %Probabilidad de cruza
Qt=[]; Np=size(Pt,1);
for in1=1:Np
       r1=randi(Np);r2=randi(Np); r3=randi(Np);
       while(r1==r2==r3==in1)
          r1=randi(Np);r2=randi(Np);r3=randi(Np);%Generados sean diferentes
       end
       v(1,:)=Pt(r3,1:d)+F*(Pt(r1,1:d)-Pt(r2,1:d));%Mutacion
       id1=find(v<x_low);id2=find(v>x_high);  
       id=union(id1,id2);
       v(1,id)=x_low(1,id)+(x_high(1,id)-x_low(1,id)).*rand(1,length(id));  
       u(1,:)=Pt(in1,1:d);
       j_rand=randi(d);
       for in2=1:d                         %Genero vector de prueba
          if(rand()<Cr || in2==j_rand)     %Cruza
             u(1,j_rand)=v(1,j_rand);
          end
       end 
       Qt(in1,1:d)=u;
       %Qt(in1,d+1: M+d) =functionObjective(u, M, d);
       Qt(in1,d+1: M+d) = fH_a(round(u(1,1:d-1)), X,d-1,u(1,d));
end    
end 