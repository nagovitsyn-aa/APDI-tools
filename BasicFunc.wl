(* ::Package:: *)

BeginPackage["BasicFunc`"];
\[Omega]2f::usage="\[Omega]2f@\[Omega] Transform angular frequency(rad/s) to ordinary frequency in GHz ";
f2\[Omega]::usage="f2\[Omega]@f Transform ordinary frequency in GHz to angular frequency(rad/s) ";

vTi::usage="vTi[x] Ions thermal speed";
vTe::usage="vTe[x] Electrons thermal speed";
\[Omega]ce::usage="\[Omega]ce[x] \[Omega]ce>0";
\[Omega]ci::usage="\[Omega]ci[x]";
\[Omega]pe::usage="\[Omega]pe[x]";
\[Omega]pi::usage="\[Omega]pi[x]";
\[Omega]uh::usage="\[Omega]uh[x] Frequency of upper-hybryd resonance";
\[Omega]lh::usage="\[Omega]uh[x] Frequency of lower-hybryd resonance";
\[CurlyEpsilon]::usage="\[CurlyEpsilon][\[Omega],x] Transverse diagonal component of dielectric permittivity tensor";
g::usage="g[\[Omega],x] Non-diagonal component of dielectric permittivity tensor";
\[Eta]::usage="\[Eta][\[Omega],x] Parallel diagonal component of dielectric permittivity tensor";
kX::usahe="kX[\[Omega],x[ Wave number of X-wave";
lt::usage="lt[\[Omega],x] The Bohm-Gross amendment \[CurlyEpsilon]=\[CurlyEpsilon]0(1+(lt k)^2)) \[Omega]ce<\[Omega]<2\[Omega]ce";
qplus::usage="qplus[\[Omega],x] Slow branch of the upper hybrid wave";
qminus::usage="qminus[\[Omega],x] Fast branch of the upper hybrid wave";
Duh::usage="Duh[q,\[Omega],x] Dispersion relation of the upper hybrid wave";

me=9.10938215`*10^-28;
mi=2 me*1.8362`*10^3;
c=3*10^10;
el=4.80320451 10^-10;
\[Omega]0=f2\[Omega]@82.5;

Begin["`Private`"];
Needs["Profiles`"];
\[Omega]2f[\[Omega]_]=\[Omega]/(2. \[Pi] 10^9);
f2\[Omega][f_]=f 2. \[Pi] 10^9;

eVtoerg=1.604 10^-12;
vTi[x_]=(2 Ti[x] eVtoerg/mi)^(1/2);
vTe[x_]=(2 Te[x] eVtoerg/me)^(1/2);
\[Omega]ce[x_]=(el B[x])/(me c);
\[Omega]ci[x_]=(el B[x])/(mi c);
\[Omega]pe[x_]=Sqrt[4 \[Pi] den[x] el^2/me];
\[Omega]pi[x_]=Sqrt[4 \[Pi] den[x] el^2/mi];

(*\[Omega]uh[x_]=Sqrt[\[Omega]pe[x]^2+\[Omega]ce[x]^2];*)
\[Omega]lh[x_]=Sqrt[(\[Omega]ce[x]^2+\[Omega]pe[x]^2+\[Omega]ci[x]^2+\[Omega]pi[x]^2)/2-1/2 Sqrt[(\[Omega]ce[x]^2+\[Omega]pe[x]^2+\[Omega]ci[x]^2+\[Omega]pi[x]^2)^2-4(\[Omega]ce[x]^2 \[Omega]ci[x]^2+\[Omega]pe[x]^2 \[Omega]ci[x]^2+\[Omega]pi[x]^2 \[Omega]ce[x]^2)]];
\[Omega]uh[x_]=Sqrt[(\[Omega]ce[x]^2+\[Omega]pe[x]^2+\[Omega]ci[x]^2+\[Omega]pi[x]^2)/2+1/2 Sqrt[(\[Omega]ce[x]^2+\[Omega]pe[x]^2+\[Omega]ci[x]^2+\[Omega]pi[x]^2)^2-4(\[Omega]ce[x]^2 \[Omega]ci[x]^2+\[Omega]pe[x]^2 \[Omega]ci[x]^2+\[Omega]pi[x]^2 \[Omega]ce[x]^2)]];

\[CurlyEpsilon][\[Omega]_,x_]=1-\[Omega]pe[x]^2/(\[Omega]^2-\[Omega]ce[x]^2)-\[Omega]pi[x]^2/(\[Omega]^2-\[Omega]ci[x]^2);
\[Eta][\[Omega]_,x_]=1-\[Omega]pe[x]^2/(\[Omega]^2);
g[\[Omega]_,x_]=-((\[Omega]ce[x] \[Omega]pe[x]^2)/(\[Omega] (\[Omega]^2-\[Omega]ce[x]^2)))+(\[Omega]ci[x] \[Omega]pi[x]^2)/(\[Omega] (\[Omega]^2-\[Omega]ci[x]^2));
kX[\[Omega]_,x_]=\[Omega]/c Sqrt[\[CurlyEpsilon][\[Omega],x]-g[\[Omega],x]^2/\[CurlyEpsilon][\[Omega],x]];

lt[\[Omega]_,x_]= vTe[x] \[Omega]pe[x]((3)/(2(\[Omega]^2-\[Omega]ce[x]^2) (4 \[Omega]ce[x]^2-\[Omega]^2)))^(1/2);
qplus[\[Omega]_,x_]=((-\[CurlyEpsilon][\[Omega],x]+(\[CurlyEpsilon][\[Omega],x]^2-(4 g[\[Omega],x]^2 \[Omega]^2 lt[\[Omega],x]^2)/c^2)^(1/2))/(2 lt[\[Omega],x]^2))^(1/2);
qminus[\[Omega]_,x_]=((-\[CurlyEpsilon][\[Omega],x]-(\[CurlyEpsilon][\[Omega],x]^2-(4 g[\[Omega],x]^2 \[Omega]^2 lt[\[Omega],x]^2)/c^2)^(1/2))/(2 lt[\[Omega],x]^2))^(1/2);

Duh[k_,\[Omega]_,x_]=lt[\[Omega],x]^2 k^4+\[CurlyEpsilon][\[Omega],x] k^2+g[\[Omega],x]^2 \[Omega]^2/c^2;

(*SetAttributes[{"\[Eta]","c","Duh","el","g","kX","lt","me","mi","qminus","qplus","vTe","vTi","\[CurlyEpsilon]","\[Omega]ce","\[Omega]ci","\[Omega]pe","\[Omega]pi"},{Protected}]*)

End[];
EndPackage[];
