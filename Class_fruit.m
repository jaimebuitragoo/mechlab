
clear all;
close all;
load('data_base.mat')


%% Datos para Clasificar
%p=[a b];
p=Fuerza_LAB(:,2:3);
c1=1;c2=1;c3=1;c4=1;


for i=1:length(p(:,1))
    
state(1) = sort_ellipse2([p(i,1) p(i,2)],[Categ(1,:)]);
state(2) = sort_ellipse2([p(i,1) p(i,2)],[Categ(2,:)]);
state(3) = sort_ellipse2([p(i,1) p(i,2)],[Categ(3,:)]);
state(4) = sort_ellipse2([p(i,1) p(i,2)],[Categ(4,:)]);

j=sum([1 2 3 4].*state);
 
if (j>=1)&(j<=4)&(state(1)*state(2)==0)

switch j
    
    case 1
    Verde(c1,1)=[i];
    c1=c1+1;
    case 2
   SemiMaduro(c2,1)=[i];
    c2=c2+1;
    case 3
     Maduro(c3,1)=[i];
    c3=c3+1;
    case 4
      SobreMaduro(c4,1)=[i];
    c4=c4+1;
       
end


else
   
 D(1)=sqrt(sum([(p(i,1)-Categ(1,1))^2 (p(i,2)-Categ(1,2))^2]));  
 D(2)=sqrt(sum([(p(i,1)-Categ(2,1))^2 (p(i,2)-Categ(2,2))^2]));  
 D(3)=sqrt(sum([(p(i,1)-Categ(3,1))^2 (p(i,2)-Categ(3,2))^2]));  
 D(4)=sqrt(sum([(p(i,1)-Categ(4,1))^2 (p(i,2)-Categ(4,2))^2]));  

[v vmin]=min(D);           
 
switch vmin
    
   
    case 1
    Verde(c1,1)=[i];
    c1=c1+1;
    case 2
   SemiMaduro(c2,1)=[i];
    c2=c2+1;
    case 3
     Maduro(c3,1)=[i];
    c3=c3+1;
    case 4
      SobreMaduro(c4,1)=[i];
    c4=c4+1;
       
end

end   
  
end   

kk=1;

% Clasificación de la masa.
Mv=Masa_Fuerza(Verde);
Msm=Masa_Fuerza(SemiMaduro);
Msbm=Masa_Fuerza(SobreMaduro);
Mm=Masa_Fuerza(Maduro);

% Clasificación del volumen.
Vv=Volumen_Fuerza(Verde);
Vsm=Volumen_Fuerza(SemiMaduro);
Vm=Volumen_Fuerza(Maduro);
Vsbm=Volumen_Fuerza(SobreMaduro);

% Clasificación de la densidad.
Dv=Densidad_Fuerza(Verde);
Dsm=Densidad_Fuerza(SemiMaduro);
Dm=Densidad_Fuerza(Maduro);
Dsbm=Densidad_Fuerza(SobreMaduro);

%Clasificación de los frutos desprendidos por torque.
Torque_LAB = Db_Torque(:, 3:5);

p = Torque_LAB(:, 2:3);
c1=1; c2=1; c3=1; c4=1;


for i=1:length(p(:,1))
    
    state(1) = sort_ellipse2([p(i,1) p(i,2)],[Categ(1,:)]);
    state(2) = sort_ellipse2([p(i,1) p(i,2)],[Categ(2,:)]);
    state(3) = sort_ellipse2([p(i,1) p(i,2)],[Categ(3,:)]);
    state(4) = sort_ellipse2([p(i,1) p(i,2)],[Categ(4,:)]);

    j=sum([1 2 3 4].*state);
 
    if (j<=4)&(state(1)*state(2)==0)

        switch j

            case 1
                TVerde(c1,1)=[i];
                c1=c1+1;
            case 2
                TSemiMaduro(c2,1)=[i];
                c2=c2+1;
            case 3
                TMaduro(c3,1)=[i];
                c3=c3+1;
            case 4
                TSobreMaduro(c4,1)=[i];
                c4=c4+1;

        end


    else

        D(1)=sqrt(sum([(p(i,1)-Categ(1,1))^2 (p(i,2)-Categ(1,2))^2]));  
        D(2)=sqrt(sum([(p(i,1)-Categ(2,1))^2 (p(i,2)-Categ(2,2))^2]));  
        D(3)=sqrt(sum([(p(i,1)-Categ(3,1))^2 (p(i,2)-Categ(3,2))^2]));  
        D(4)=sqrt(sum([(p(i,1)-Categ(4,1))^2 (p(i,2)-Categ(4,2))^2]));  

        [v vmin]=min(D);           

        switch vmin

            case 1
                TVerde(c1,1)=[i];
                c1=c1+1;
            case 2
                TSemiMaduro(c2,1)=[i];
                c2=c2+1;
            case 3
                TMaduro(c3,1)=[i];
                c3=c3+1;
            case 4
                TSobreMaduro(c4,1)=[i];
                c4=c4+1;

        end

    end   
  
end   

% kk=1;

Masa_Torque = Db_Torque(:, 2);

% Separación de los datos de color por estados de maduración.
% Datos de color por estados de maduración.
TLAB_Verde = Torque_LAB(TVerde, :);
TLAB_SemiMaduro = Torque_LAB(TSemiMaduro, :);
TLAB_Maduro = Torque_LAB(TMaduro, :);
TLAB_SobreMaduro = Torque_LAB(TSobreMaduro, :);

% Clasificación de la masa.
TMv=Masa_Torque(TVerde);
TMsm=Masa_Torque(TSemiMaduro);
TMm=Masa_Torque(TMaduro);
TMsbm=Masa_Torque(TSobreMaduro);
