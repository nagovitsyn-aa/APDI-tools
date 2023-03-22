(* ::Package:: *)

(* ::Text:: *)
(*1 dimension only. Parallel component of wave vector is 0. Z-function expands at infinite argument.*)
(**)


BeginPackage["HotDisp`",{"BasicFunc`"}];
\[CurlyEpsilon]11Hot::usage="\[CurlyEpsilon]11Hot[q,\[Omega],x] First component of warm permittivity tensor in case of perpendicular propagation";
\[Lambda]i;
DIBq::usage="DIBq[q,\[Omega],x] Group velocity of IB wave. Derivative of IB wave dispersion relation by q wave number";
\[CurlyEpsilon]Exact;
\[CurlyEpsilon]Approx;
DIBqExact;
DIBqApprox;
DIBqPW;


Begin["Private`"];


With[{nApprox=20, \[Lambda]b=650,ntermExact=100},
BesselSeries[n_,\[Lambda]_]=Plus@@((List@@(Series[BesselI[n,\[Lambda]],{\[Lambda],\[Infinity],nApprox}]//Normal))[[2;;]]);
nmax=n/.FindRoot[Abs[BesselSeries[n,\[Lambda]b]/BesselI[n,\[Lambda]b]]==0.5,{n,80}];

termSeries[n_,\[Lambda]_]=Plus@@((List@@(Series[E^-\[Lambda]/\[Lambda] BesselI[n,\[Lambda]],{\[Lambda],\[Infinity],nApprox}]//Normal))[[2;;]]);
sumSeries[nres_,\[Lambda]_]=Sum[n^2/((nres)^2-n^2) termSeries[n,\[Lambda]],{n,1,Ceiling@nmax}];
sum[nres_,\[Lambda]_]=Sum[n^2/((nres)^2-n^2) BesselI[n,\[Lambda]] E^-\[Lambda]/\[Lambda],{n,1,100}];


termSeries2[n_,\[Lambda]_]=Plus@@((List@@(Series[E^-\[Lambda]/\[Lambda] (\[Lambda] BesselI[n-1,\[Lambda]]-(\[Lambda]+n)BesselI[n,\[Lambda]]),{\[Lambda],\[Infinity],nApprox+20}]//Normal//N))[[3;;]]);
sumSeries2[nres_,\[Lambda]_]=Sum[n^2/((nres)^2-n^2) termSeries2[n,\[Lambda]],{n,1,20+Ceiling@nmax}];
sumD[nres_,\[Lambda]_]=Sum[n^2/((nres)^2-n^2) E^-\[Lambda]/\[Lambda](\[Lambda] BesselI[n-1,\[Lambda]]-(\[Lambda]+n)BesselI[n,\[Lambda]]),{n,1,ntermExact}];




\[Lambda]i[q_,x_]=((q vTi[x])/\[Omega]ci[x])^2;

\[CurlyEpsilon]11Hot[q_,\[Omega]_,x_]=1-\[Omega]pe[x]^2/(\[Omega]^2-\[Omega]ce[x]^2)-2 \[Omega]pi[x]^2/\[Omega]ci[x]^2 Piecewise[{{sumSeries[\[Omega]/\[Omega]ci[x],\[Lambda]i[q,x]],\[Lambda]i[q,x]>\[Lambda]b},{sum[\[Omega]/\[Omega]ci[x],\[Lambda]i[q,x]],\[Lambda]i[q,x]<=\[Lambda]b}}];
\[CurlyEpsilon]Exact[q_,\[Omega]_,x_]=1-\[Omega]pe[x]^2/(\[Omega]^2-\[Omega]ce[x]^2)-2 \[Omega]pi[x]^2/\[Omega]ci[x]^2 sum[\[Omega]/\[Omega]ci[x],\[Lambda]i[q,x]];
\[CurlyEpsilon]Approx[q_,\[Omega]_,x_]=1-\[Omega]pe[x]^2/(\[Omega]^2-\[Omega]ce[x]^2)-2 \[Omega]pi[x]^2/\[Omega]ci[x]^2 sumSeries[\[Omega]/\[Omega]ci[x],\[Lambda]i[q,x]];

DIBqExact[q_,\[Omega]_,x_]=1-\[Omega]pe[x]^2/(\[Omega]^2-\[Omega]ce[x]^2)-2 \[Omega]pi[x]^2/\[Omega]ci[x]^2 sumD[\[Omega]/\[Omega]ci[x],\[Lambda]i[q,x]];
DIBqApprox[q_,\[Omega]_,x_]=1-\[Omega]pe[x]^2/(\[Omega]^2-\[Omega]ce[x]^2)-2 \[Omega]pi[x]^2/\[Omega]ci[x]^2 sumSeries2[\[Omega]/\[Omega]ci[x],\[Lambda]i[q,x]];
DIBqPW[q_,\[Omega]_,x_]=1-\[Omega]pe[x]^2/(\[Omega]^2-\[Omega]ce[x]^2)-2 \[Omega]pi[x]^2/\[Omega]ci[x]^2 Piecewise[{{sumSeries2[\[Omega]/\[Omega]ci[x],\[Lambda]i[q,x]],\[Lambda]i[q,x]>\[Lambda]b},{sumD[\[Omega]/\[Omega]ci[x],\[Lambda]i[q,x]],\[Lambda]i[q,x]<=\[Lambda]b}}];
DIBq[q_,\[Omega]_,x_]=DIBqExact[q,\[Omega],x];
]



End[];
EndPackage[];
