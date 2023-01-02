fprintf ("please enter the WG resonator dimensions in meters (b<a<d) ");
a=input("a=");
b=input("b=");
d=input("d=");
epsr=input("please enter the relative permittivity ");
losstan=input("please enter the loss tangent ");
cond=input("please enter the conductivity ");
Qd=1/losstan;
c=3*10^8;
muo=4*pi*(10^-7);
epso=8.854187817*10^(-12);
TE101(1,0,1,c,epso,muo,epsr,a,b,d,Qd,cond);
if( ((1/d)^2 + (1/b)^2) > ((1/a)^2 + (2/d)^2) )
TE102(1,0,2,c,epso,muo,epsr,a,b,d,Qd,cond);
end
if ( ((1/d)^2 + (1/b)^2) < ((1/a)^2 + (2/d)^2) )
TE011(0,1,1,c,epso,muo,epsr,a,b,d,Qd,cond);
else
TM110(1,1,0,c,epso,muo,epsr,a,b,d,Qd,cond);
end
function TE101(m,n,l,c,epso,muo,epsr,a,b,d,Qd,cond)
f=c*((m*pi/a)^2+(n*pi/b)^2+(l*pi/d)^2)^(1/2)/(2*pi*epsr^(1/2));
beta = l*pi/d;
Rs=(pi*f*muo/cond)^(1/2);
We = (epso*epsr*a*b*d)*((2*pi*f*muo*a/pi)^(2))/4;
Pcond = Rs* (((beta*a/pi)^2)*(2*a*b+a*d) + 2*b*d+ a*d);
Qcond = 2*2*pi*f*We/Pcond;
BW = 1/Qcond + 1/Qd ;
Qtot = 1/BW;
fprintf ("resonance frequency of TE101 is: %fGHz\n", f*10^(-9) );
fprintf ("resonance fractional bandwidth of TE101 in percentage is: %f \n ", BW *100);
fprintf ("resonance absolute bandwidth of TE101 is: %f Hz\n ", BW *f);
fprintf ("quality factor of TE101 is: %f\n", Qtot );
end
function TE102(m,n,l,c,epso,muo,epsr,a,b,d,Qd,cond)
f=c*((m*pi/a)^2+(n*pi/b)^2+(l*pi/d)^2)^(1/2)/(2*pi*epsr^(1/2));
beta = l*pi/d;
Rs=(pi*f*muo/cond)^(1/2);
We = (epso*epsr*a*b*d)*((2*pi*f*muo*a/pi)^(2))/4;
Pcond = Rs* (((beta*a/pi)^2)*(2*a*b+a*d) + 2*b*d+ a*d);
Qcond = 2*2*pi*f*We/Pcond;
BW = 1/Qcond + 1/Qd ;
Qtot = 1/BW;
fprintf ("resonance frequency of TE102 is: %fGHz\n", f*10^(-9) );
fprintf ("resonance fractional bandwidth of TE102 in percentage is: %f \n ", BW*100 );
fprintf ("resonance absolute bandwidth of TE101 is: %f Hz\n ", BW *f);
fprintf ("quality factor of TE102 is: %f\n", Qtot );
end
function TE011(m,n,l,c,epso,muo,epsr,a,b,d,Qd,cond)
f=c*((m*pi/a)^2+(n*pi/b)^2+(l*pi/d)^2)^(1/2)/(2*pi*epsr^(1/2));
beta = l*pi/d;
Rs=(pi*f*muo/cond)^(1/2);
We = (epso*epsr*a*b*d)*((2*pi*f*muo*b/pi)^(2))/4;
Pcond = Rs* (((beta*b/pi)^2)*(2*a*b+b*d) + 2*a*d+ b*d);
Qcond = 2*2*pi*f*We/Pcond;
BW = 1/Qcond + 1/Qd ;
Qtot = 1/BW;
fprintf ("resonance frequency of TE011 is: %fGHz\n", f*10^(-9) );
fprintf ("resonance fractional bandwidth of TE011 in percentage is: %f \n ", BW*100 );
fprintf ("resonance absolute bandwidth of TE011 is: %f Hz\n ", BW *f);
fprintf ("quality factor of TE011 is: %f\n", Qtot );
end
function TM110(m,n,l,c,epso,muo,epsr,a,b,d,Qd,cond)
f=c*((m*pi/a)^2+(n*pi/b)^2+(l*pi/d)^2)^(1/2)/(2*pi*epsr^(1/2));
Rs = (pi*f*muo/cond)^(1/2);
Wm=muo*a*b*d/16;
Pcond = 0.5* Rs *((1/2)*(2*pi*f)^2 *muo*epso*epsr*a*b + ((pi^2)/(b^2))*a*d+
((pi^2)/(a^2))*b*d);
Qcond = 2*2*pi*f*Wm/Pcond;
BW = 1/Qcond + 1/Qd ;
Qtot = 1/BW;
fprintf ("resonance frequency of TM110 is: %fGHz\n", f*10^(-9) );
fprintf ("resonance fractional bandwidth of TM110 in percentage is: %f \n ", BW*100 );
fprintf ("resonance absolute bandwidth of TE101 is: %f Hz\n ", BW *f);
fprintf ("quality factor of TM110 is: %f\n", Qtot );
end
