function [Iz,Iy,Cz,Cy, A] = TrapezoidInertia_Top(Beta,t5,t7,chord,h)  
    %inercia trapecio acimetrico 
    % con ejes de referencia rotados beta grados
    % eje de X pocitivo asia la izquierda
 
    beta=Beta;
    t1=t5;
    t6=t7;
    c=chord;
    h=h;
    L=(0.70*c-.45*c-.125*c);

    %data de la geometria interna del trapecio
    gama1=90-beta;
    gama1=(gama1*pi)/180;
    beta=beta*pi/180;

    a1=t1*sin(gama1); %a=altura b=base
    b1=t1*cos(gama1);
    a4=t6*sin(beta);
    b4=t6*cos(beta);
    a2=a4;
    b2=(L/cos(beta))-b1-b4;
    a3=a1-a2;
    b3=b2;


    % centroides en eje y2
    y1=(1/3)*a1;
    y2=(1/2)*a2;
    y3=((1/3)*a3+a2);
    y4=(1/3)*a4;

    %centroides en eje x
    x1=(2/3)*b1;
    x2=((1/2)*b2+b1);
    x3=((1/3)*b3+b1);
    x4=((1/3)*b4+b2+b1);

    %areas
    A1=(1/2)*b1*a1;
    A2=a2*b2;
    A3=(1/2)*b3*a3;
    A4=(1/2)*b4*a4;

    A=A1+A2+A3+A4;

    %centroide del trapecio
    X=(A1*x1+A2*x2+A3*x3+A4*x4)/(A1+A2+A3+A4);
    Cz=(X*cos(beta)+0.45*c)*-1;
    Y=(A1*y1+A2*y2+A3*y3+A4*y4)/(A1+A2+A3+A4);
    Cy=-(Y*cos(beta)-(h/2));

    %distancia al eje neutro para x
    r1=Y-y1;
    r2=Y-y2;
    r3=y3-Y; %se encuentra sobre el eje neutro
    r4=Y-y4;
    %para z
    r1y=X-x1;
    r2y=x2-X;
    r3y=X-x3;
    r4y=x4-X;

    %inercia en x
    Ix1=((1/36)*b1*a1^3)+A1*r1^2;
    Ix2=((1/12)*b2*a2^3)+A2*r2^2;
    Ix3=((1/36)*b3*a3^3)+A3*r3^2;
    Ix4=((1/36)*b4*a4^3)+A4*r4^2;

    Iz=Ix1+Ix2+Ix3+Ix4;

    %inercia en y
    Iy1=((1/36)*a1*b1^3)+A1*r1y^2;
    Iy2=((1/12)*a2*b2^3)+A2*r2y^2;
    Iy3=((1/36)*a3*b3^3)+A3*r3y^2;
    Iy4=((1/36)*a4*b4^3)+A4*r4y^2;

    Iy=Iy1+Iy2+Iy3+Iy4;
end