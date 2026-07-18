(* ::Package:: *)

(* PolyDet` -- single scalar section.  Loaded by Kernel/PolyDet.wl from inside
   PolyDet`Private` (Get'd by the loader); NOT a standalone package -- it has no
   BeginPackage/Begin of its own and relies on the loader's context.
   Transfer matrix Ts (det T^s=1), PolyDetLnR/Ip/Sigma/Ren, optimal truncation
   (g_p Laurent + zeta add-back), PolyDetInputs, scalar zero removal PolyDetR1Prime. *)


(* per-matching transfer matrix T^s (det T^s = 1) *)
Ts[nu_, Rs_, ms_, ms1_] := Module[{I0, I1, K0, K1, Ie, Iff, Ke, Kf},
  I0 = BesselI[nu, ms Rs];   I1 = BesselI[nu - 1, ms Rs];
  K0 = BesselK[nu, ms Rs];   K1 = BesselK[nu - 1, ms Rs];
  Ie = BesselI[nu, ms1 Rs];  Iff = BesselI[nu - 1, ms1 Rs];
  Ke = BesselK[nu, ms1 Rs];  Kf = BesselK[nu - 1, ms1 Rs];
  Rs {{ms1 I0 Kf + ms I1 Ke, ms1 K0 Kf - ms K1 Ke},
      {ms1 I0 Iff - ms I1 Ie, ms1 K0 Iff + ms K1 Ie}}];

(* zero-removal source  d T^s / d(m^2)  (offset m^2 -> m^2 + mu^2);
   slot 1 = d/d(m_s^2) [inner], slot 2 = d/d(m_{s+1}^2) [outer].  Building block
     d/d(m^2) I_nu(m R) =  R/(2m) I_{nu-1} - nu/(2m^2) I_nu,
     d/d(m^2) K_nu(m R) = -R/(2m) K_{nu-1} - nu/(2m^2) K_nu,
   reduced to orders {nu,nu-1} -> reuses the same 8 Bessels as Ts; det T^s=1
   preserved (tr(T^-1 dT)=0). Verified vs central differences to ~1e-22. *)
dTs[nu_, Rs_, ms_, ms1_, slot_] := Module[
  {I0, I1, K0, K1, Ie, Iff, Ke, Kf},
  I0 = BesselI[nu, ms Rs];   I1 = BesselI[nu - 1, ms Rs];
  K0 = BesselK[nu, ms Rs];   K1 = BesselK[nu - 1, ms Rs];
  Ie = BesselI[nu, ms1 Rs];  Iff = BesselI[nu - 1, ms1 Rs];
  Ke = BesselK[nu, ms1 Rs];  Kf = BesselK[nu - 1, ms1 Rs];
  If[slot === 1,
   Module[{dI0 =  Rs/(2 ms) I1 - nu/(2 ms^2) I0, dI1 =  Rs/(2 ms) I0 + (nu - 1)/(2 ms^2) I1,
           dK0 = -Rs/(2 ms) K1 - nu/(2 ms^2) K0, dK1 = -Rs/(2 ms) K0 + (nu - 1)/(2 ms^2) K1},
    Rs {{ms1 dI0 Kf  + I1 Ke/(2 ms) + ms dI1 Ke,   ms1 dK0 Kf  - K1 Ke/(2 ms) - ms dK1 Ke},
        {ms1 dI0 Iff - I1 Ie/(2 ms) - ms dI1 Ie,   ms1 dK0 Iff + K1 Ie/(2 ms) + ms dK1 Ie}}],
   Module[{dIe =  Rs/(2 ms1) Iff - nu/(2 ms1^2) Ie, dIff =  Rs/(2 ms1) Ie + (nu - 1)/(2 ms1^2) Iff,
           dKe = -Rs/(2 ms1) Kf - nu/(2 ms1^2) Ke, dKf = -Rs/(2 ms1) Ke + (nu - 1)/(2 ms1^2) Kf},
    Rs {{I0 Kf/(2 ms1) + ms1 I0 dKf + ms I1 dKe,    K0 Kf/(2 ms1) + ms1 K0 dKf - ms K1 dKe},
        {I0 Iff/(2 ms1) + ms1 I0 dIff - ms I1 dIe,  K0 Iff/(2 ms1) + ms1 K0 dIff + ms K1 dIe}}]]];

PolyDetXMax[m2_, radii_] := Max@Table[
   Max[Abs[Sqrt[m2[[s]]] radii[[s]]], Abs[Sqrt[m2[[s + 1]]] radii[[s]]]],
   {s, Length[radii]}];

(* MASSLESS false vacuum (m2[[-1]] == 0, symmetric FV): the final-segment solutions
   degenerate from Bessels to the power laws rho^(+-nu), so the ratio is normalized
   against the regular solution rho^nu of the massless free operator,
     ln R_nu = LogGamma[nu+1] + nu Log[2/m_1] + Log[A_N] ,
   A_N the rho^nu coefficient after the chain.  Applying the power-law boundary
   inverse to the last matching gives the closed-form final row
     A_N = (m_s R^(1-nu)/(2 nu)) (I_{nu-1}(m_s R) v_1 - K_{nu-1}(m_s R) v_2) ,
   (m_s the last WALL mass, R the last radius, (v_1, v_2) the first column of the
   interior chain product), evaluated in log form so no explicit R^(+-nu) appears.
   Equals the m_FV -> 0 limit of the massive chain; single matching:
   R_nu = Gamma[nu] (2/z)^(nu-1) I_{nu-1}(z), z = m_1 R_1. *)
lnRMLfin[nu_, Rs_, ms_, v1_, v2_, m1_] := Module[
  {I1 = BesselI[nu - 1, ms Rs], K1 = BesselK[nu - 1, ms Rs]},
  Re[LogGamma[nu + 1] + nu Log[2/m1] + (1 - nu) Log[Rs] +
     Log[ms (I1 v1 - K1 v2)/(2 nu)]]];

PolyDetLnR[nu_, m2_, radii_, prec_: 60] := Module[{mm, rr, Nm, Tlist, rest},
  mm = Sqrt[SetPrecision[m2, prec]];  rr = SetPrecision[radii, prec];  Nm = Length[radii];
  If[TrueQ[m2[[-1]] == 0],
    rest = If[Nm > 1,
       Dot @@ Reverse[Table[Ts[nu, rr[[s]], mm[[s]], mm[[s + 1]]], {s, 1, Nm - 1}]],
       IdentityMatrix[2]];
    lnRMLfin[nu, rr[[Nm]], mm[[Nm]], rest[[1, 1]], rest[[2, 1]], mm[[1]]],
    Tlist = Table[Ts[nu, rr[[s]], mm[[s]], mm[[s + 1]]], {s, 1, Nm}];
    Re[Log[(mm[[-1]]/mm[[1]])^nu (Dot @@ Reverse[Tlist])[[1, 1]]]]]];

(* l=1 translational zero-mode removal.  R_1' replaces the (artifact) l=1 ratio
   R_{l=1}=0 in the determinant (degeneracy d_1=D).  Offset the interior shell
   masses m_s^2 -> m_s^2+mu^2 (FV m_N fixed); since R_{l=1}(inf)=0 for the exact
   bounce, R_1' is the linear coefficient
     R_1' = mhat^2 d/d(mu^2)[ (m_N/m_1)^nu T^tot_11 ]_{nu=D/2}
          = mhat^2 (m_N/m_1)^{D/2} ( dT^tot_11 - (D/2)/(2 m_1^2) T^tot_11 ),
   with dT^tot from the chain product rule using dTs.  (The prefactor term
   vanishes as the polygon refines, where T^tot_11 -> 0.)  -> e^{D-1}/12 thin-wall. *)
Options[PolyDetR1Prime] = {"Dimension" -> 4, "Multipole" -> 1, WorkingPrecision -> 80};
PolyDetR1Prime[m2_, radii_, OptionsPattern[]] := Module[
  {dim, prec, l, nu, mm, M, rr, Tlist, Ttot, dTtot},
  dim = OptionValue["Dimension"];  prec = OptionValue[WorkingPrecision];  l = OptionValue["Multipole"];
  nu = l + dim/2 - 1;  M = Length[radii];               (* l=1 -> nu=D/2 (translational); l=0 -> nu=D/2-1 (Goldstone) *)
  mm = Sqrt[SetPrecision[m2, prec]];  rr = SetPrecision[radii, prec];
  Tlist = Table[Ts[nu, rr[[s]], mm[[s]], mm[[s + 1]]], {s, M}];
  Ttot = Dot @@ Reverse[Tlist];
  dTtot = Sum[                                           (* d T^tot / d(mu^2) *)
     Dot @@ Reverse[ReplacePart[Tlist, s ->
        dTs[nu, rr[[s]], mm[[s]], mm[[s + 1]], 1] +
        If[s < M, dTs[nu, rr[[s]], mm[[s]], mm[[s + 1]], 2], 0]]], {s, M}];
  Re[ m2[[-1]] (mm[[-1]]/mm[[1]])^nu (dTtot[[1, 1]] - nu/(2 m2[[1]]) Ttot[[1, 1]]) ]];

(* vectorised inner loop: ALL lnR_nu for l=0..Ns-1 (nu=l+dim/2-1) in two Bessel
   calls.  The chain needs Bessels only on the 2(N-1) arguments {m_s R_s, m_{s+1} R_s}
   at orders {nu-1, nu}; tabulate the whole (orders x arguments) grid at once
   (BesselI/K are Listable in both slots), then assemble each T^s from lookups.
   Reproduces PolyDetLnR exactly; ~1.7x faster than the per-nu loop. *)
lnRall[m2_, radii_, dim_, Ns_, prec_, lmin_: 0] := Module[
  {mlfv, mm, rr, Nm, Nb, av, bv, nus, orders, allargs, om, am, Itab, Ktab, No},
  mlfv = TrueQ[m2[[-1]] == 0];                     (* massless FV: power-law final segment *)
  mm = Sqrt[SetPrecision[m2, prec]];  rr = SetPrecision[radii, prec];  Nm = Length[radii];
  Nb = If[mlfv, Nm - 1, Nm];                       (* skip the vanishing FV Bessel argument *)
  av = Table[mm[[s]] rr[[s]], {s, Nm}];  bv = Table[mm[[s + 1]] rr[[s]], {s, Nb}];
  nus = Table[l + dim/2 - 1, {l, lmin, Ns - 1}];   (* lmin>0: bottom-cut, skip low multipoles *)
  orders = Join[{nus[[1]] - 1}, nus];  No = Length[orders];   (* nu-1 sits one below each nu *)
  allargs = Join[av, bv];                                     (* a_1..a_Nm, then b_1..b_Nb *)
  om = Transpose[ConstantArray[orders, Length@allargs]];      (* conforming arrays *)
  am = ConstantArray[allargs, No];
  Itab = BesselI[om, am];  Ktab = BesselK[om, am];            (* two vectorised calls total *)
  (* high-nu modes can underflow machine numbers (negligible terms); silence munfl *)
  Quiet[Table[Module[{nu = nus[[i]], Tlist, rest},               (* row i+1 = nu, row i = nu-1 *)
    Tlist = Table[Module[
       {ms = mm[[s]], ms1 = mm[[s + 1]], Rs = rr[[s]],
        I0 = Itab[[i + 1, s]], I1 = Itab[[i, s]], K0 = Ktab[[i + 1, s]], K1 = Ktab[[i, s]],
        Ie = Itab[[i + 1, Nm + s]], Iff = Itab[[i, Nm + s]], Ke = Ktab[[i + 1, Nm + s]], Kf = Ktab[[i, Nm + s]]},
      Rs {{ms1 I0 Kf + ms I1 Ke, ms1 K0 Kf - ms K1 Ke},
          {ms1 I0 Iff - ms I1 Ie, ms1 K0 Iff + ms K1 Ie}}], {s, Nb}];
    If[mlfv,
      rest = If[Nb > 0, Dot @@ Reverse[Tlist], IdentityMatrix[2]];
      (* final massless matching in log form (no explicit R^(+-nu)); Bessels at the
         a-side of matching Nm are grid column Nm: I1/K1 = order nu-1 = row i *)
      Re[LogGamma[nu + 1] + nu Log[2/mm[[1]]] + (1 - nu) Log[rr[[Nm]]] +
         Log[mm[[Nm]] (Itab[[i, Nm]] rest[[1, 1]] - Ktab[[i, Nm]] rest[[2, 1]])/(2 nu)]],
      Re[Log[(mm[[-1]]/mm[[1]])^nu (Dot @@ Reverse[Tlist])[[1, 1]]]]]], {i, Length@nus}],
   {General::munfl}]];

(* ---- FrozenWall thin-wall fast band ----
   Case-b thin wall: exact Bessels ONLY for the TV-plateau log-derivative and the
   FV decomposition (4 per mode, N-independent, vectorised over the band); across
   the wall shells the centrifugal term is FROZEN per sub-interval (at the
   log-mean value, see subgrid below) and the piecewise-constant m^2 problem is
   solved EXACTLY by cosh/sinh (cos/sin for kappa^2 < 0 -- no turning-point
   breakdown) on nsub frozen sub-intervals per shell.  Variables: u = Sqrt[rho] psi,
   so u'' = [m^2 + (nu^2-1/4)/rho^2] u;  Wronskian[uI, uK] = -1 gives the growing
   amplitude A_N = u' uK - u uK' = uK (u' - u uK'/uK).  The wall product carries
   only e^{O(kappa * wall width)} dynamic range -> MACHINE precision suffices;
   the huge plateau/FV exponentials enter additively as logs at precB digits.
   Band error falls ~1/x_max^3 and ~nsub^-1.7. *)
lnRallFrozen[m2_, radii_, dim_, Ns_, lminF_, nsub_, precB_: 40] := Module[
  {m2n, rr, m1, mhat, R0, RN, nus, x0, xN, I0v, I1v, K0v, K1v,
   sig0, lnu0, sigK, lnuK, lnpref, m2w, subgrid},
  m2n = N[m2];  rr = N[radii];
  m1 = Sqrt[SetPrecision[m2[[1]], precB]];  mhat = Sqrt[SetPrecision[m2[[-1]], precB]];
  R0 = SetPrecision[radii[[1]], precB];     RN = SetPrecision[radii[[-1]], precB];
  nus = Table[l + dim/2 - 1, {l, lminF, Ns - 1}];
  x0 = m1 R0;  xN = mhat RN;
  I0v = BesselI[nus, x0];  I1v = BesselI[nus - 1, x0];   (* Listable in the order *)
  K0v = BesselK[nus, xN];  K1v = BesselK[nus - 1, xN];
  sig0 = N[1/(2 R0) + m1 (I1v/I0v - nus/x0)];            (* u'/u at R0: O(1) *)
  sigK = N[1/(2 RN) + mhat (-K1v/K0v - nus/xN)];         (* uK'/uK at RN: O(1) *)
  lnu0 = N[Log[Sqrt[R0]] + Log[I0v]];                    (* additive logs, precB *)
  lnuK = N[Log[Sqrt[RN]] + Log[K0v]];
  lnpref = N[nus Log[mhat/m1]];
  m2w = m2n[[2 ;; -2]];                                   (* wall-shell m^2 (may be < 0) *)
  (* per frozen sub-interval [a,b] keep {h, cinv2}: the centrifugal 1/rho^2 is frozen
     at the LOG-MEAN value cinv2 = (Log[b/a]/h)^2, NOT 1/(ab) -- this makes the
     large-nu exponent kappa h -> nu Log[b/a] the EXACT eikonal, so the frozen
     chain's UV tail matches the exact Laurent subtraction (a 1/(ab) freeze leaves
     an O(nu u^3) exponent mismatch whose d_l-weighted sum accumulates). *)
  subgrid = Table[Module[{a, b, h},
      a = rr[[s]] + (rr[[s + 1]] - rr[[s]]) (j - 1)/nsub;
      b = rr[[s]] + (rr[[s + 1]] - rr[[s]]) j/nsub;  h = b - a;
      {h, (Log[b/a]/h)^2}],
    {s, 1, Length[rr] - 1}, {j, nsub}];                   (* shell s: [rr[s], rr[s+1]], mass m2w[[s]] *)
  Table[Module[{nu = N[nus[[i]]], w = 1. + 0. I, wp, kap, h, ch, sh},
     wp = sig0[[i]] + 0. I;
     Do[Do[
        {h, kap} = subgrid[[s, j]];
        kap = Sqrt[m2w[[s]] + (nu^2 - 0.25) kap + 0. I];
        ch = Cosh[kap h];  sh = Sinh[kap h];
        {w, wp} = {ch w + sh/kap wp, kap sh w + ch wp},
        {j, nsub}], {s, Length[m2w]}];
     Re[lnu0[[i]] + lnuK[[i]] + Log[wp - w sigK[[i]]] + lnpref[[i]]]],
   {i, Length[nus]}]];

(* Universal per-matching Laurent coefficients g_p(ux,uy) of ln G_s, hardcoded as
   exact rationals (g_p = sum (c) ux^i uy^j), derived once by a Cauchy-product /
   Bell recursion of the large-order Bessel remainder series.  Model-independent. *)
GPMAX = 14;
GP = <|
  1 -> {{{1,0},1/1}, {{0,1},-1/1}},
  2 -> {},
  3 -> {{{2,0},-1/2}, {{0,2},1/2}},
  4 -> {{{2,0},1/2}, {{1,1},-1/1}, {{0,2},1/2}},
  5 -> {{{3,0},2/3}, {{2,0},-1/2}, {{0,3},-2/3}, {{0,2},1/2}},
  6 -> {{{3,0},-2/1}, {{2,1},2/1}, {{2,0},1/2}, {{1,2},2/1}, {{1,1},-1/1}, {{0,3},-2/1}, {{0,2},1/2}},
  7 -> {{{4,0},-5/4}, {{3,0},14/3}, {{2,1},-4/1}, {{2,0},-1/2}, {{1,2},4/1}, {{0,4},5/4}, {{0,3},-14/3}, {{0,2},1/2}},
  8 -> {{{4,0},29/4}, {{3,1},-5/1}, {{3,0},-10/1}, {{2,2},-9/2}, {{2,1},10/1}, {{2,0},1/2}, {{1,3},-5/1}, {{1,2},10/1}, {{1,1},-1/1}, {{0,4},29/4}, {{0,3},-10/1}, {{0,2},1/2}},
  9 -> {{{5,0},14/5}, {{4,0},-59/2}, {{3,1},24/1}, {{3,0},62/3}, {{2,1},-20/1}, {{2,0},-1/2}, {{1,3},-24/1}, {{1,2},20/1}, {{0,5},-14/5}, {{0,4},59/2}, {{0,3},-62/3}, {{0,2},1/2}},
  10 -> {{{5,0},-26/1}, {{4,1},14/1}, {{4,0},209/2}, {{3,2},12/1}, {{3,1},-94/1}, {{3,0},-42/1}, {{2,3},12/1}, {{2,2},-21/1}, {{2,1},42/1}, {{2,0},1/2}, {{1,4},14/1}, {{1,3},-94/1}, {{1,2},42/1}, {{1,1},-1/1}, {{0,5},-26/1}, {{0,4},209/2}, {{0,3},-42/1}, {{0,2},1/2}},
  11 -> {{{6,0},-7/1}, {{5,0},160/1}, {{4,1},-116/1}, {{4,0},-1383/4}, {{3,2},-32/1}, {{3,1},324/1}, {{3,0},254/3}, {{2,3},32/1}, {{2,1},-84/1}, {{2,0},-1/2}, {{1,4},116/1}, {{1,3},-324/1}, {{1,2},84/1}, {{0,6},7/1}, {{0,5},-160/1}, {{0,4},1383/4}, {{0,3},-254/3}, {{0,2},1/2}},
  12 -> {{{6,0},281/3}, {{5,1},-42/1}, {{5,0},-828/1}, {{4,2},-35/1}, {{4,1},684/1}, {{4,0},4407/4}, {{3,3},-100/3}, {{3,2},144/1}, {{3,1},-1059/1}, {{3,0},-170/1}, {{2,4},-35/1}, {{2,3},144/1}, {{2,2},-171/2}, {{2,1},170/1}, {{2,0},1/2}, {{1,5},-42/1}, {{1,4},684/1}, {{1,3},-1059/1}, {{1,2},170/1}, {{1,1},-1/1}, {{0,6},281/3}, {{0,5},-828/1}, {{0,4},4407/4}, {{0,3},-170/1}, {{0,2},1/2}},
  13 -> {{{7,0},132/7}, {{6,0},-801/1}, {{5,1},520/1}, {{5,0},19554/5}, {{4,2},208/1}, {{4,1},-3456/1}, {{4,0},-3434/1}, {{3,2},-408/1}, {{3,1},3348/1}, {{3,0},1022/3}, {{2,4},-208/1}, {{2,3},408/1}, {{2,1},-340/1}, {{2,0},-1/2}, {{1,5},-520/1}, {{1,4},3456/1}, {{1,3},-3348/1}, {{1,2},340/1}, {{0,7},-132/7}, {{0,6},801/1}, {{0,5},-19554/5}, {{0,4},3434/1}, {{0,3},-1022/3}, {{0,2},1/2}},
  14 -> {{{7,0},-340/1}, {{6,1},132/1}, {{6,0},5623/1}, {{5,2},108/1}, {{5,1},-4286/1}, {{5,0},-17502/1}, {{4,3},100/1}, {{4,2},-1143/1}, {{4,1},16098/1}, {{4,0},10559/1}, {{3,4},100/1}, {{3,3},-388/1}, {{3,2},1404/1}, {{3,1},-10388/1}, {{3,0},-682/1}, {{2,5},108/1}, {{2,4},-1143/1}, {{2,3},1404/1}, {{2,2},-342/1}, {{2,1},682/1}, {{2,0},1/2}, {{1,6},132/1}, {{1,5},-4286/1}, {{1,4},16098/1}, {{1,3},-10388/1}, {{1,2},682/1}, {{1,1},-1/1}, {{0,7},-340/1}, {{0,6},5623/1}, {{0,5},-17502/1}, {{0,4},10559/1}, {{0,3},-682/1}, {{0,2},1/2}}
|>;

pw0[x_, 0] = 1;  pw0[x_, n_] := x^n;   (* x^0 -> 1 also at x = 0 (massless FV argument) *)
Ip[p_, ux_, uy_] := Total[(#[[2]] pw0[ux, #[[1, 1]]] pw0[uy, #[[1, 2]]]) & /@ GP[p]];  (* per-matching
   I_p contribution (evaluates the GP table above); I_2 = 0 *)

PolyDetIp[m2_, radii_, pmax_] := Table[
   Sum[Ip[p, m2[[s]] radii[[s]]^2/4, m2[[s + 1]] radii[[s]]^2/4], {s, Length[radii]}],
   {p, 1, pmax}];

(* shifted index, degeneracy d_l on S^{D-1}, convergent zeta add-back Z_D(p) *)
nuOf[l_, dim_] := l + dim/2 - 1;
dDeg[l_, dim_] := (2 l + dim - 2) Gamma[l + dim - 2]/(Gamma[l + 1] Gamma[dim - 1]);
ZetaD[p_, 4] := Zeta[p - 2];                 (* Sum_l (l+1)^(2-p) *)
ZetaD[p_, 3] := 2 Zeta[p - 1, 1/2];          (* Sum_l (2l+1)/(l+1/2)^p *)
(* offset (Hurwitz) add-backs for the "SubtractFrom" split: the tail sums
   Sum_{l >= lstar} d_l/nu^p, i.e. the zeta identity zeta(s) = H_{nu*}^(s) +
   zeta(s, nu*+1) read from the Hurwitz side *)
ZetaDOff[p_, 4, lstar_] := HurwitzZeta[p - 2, lstar + 1];
ZetaDOff[p_, 3, lstar_] := 2 HurwitzZeta[p - 1, lstar + 1/2];

(* extreme-thin-wall bottom cut: the multipole sum is dominated by nu ~ x_max; the low ones
   carry weight ~ (nu/x_max)^(D-1)  [Delta^2 in D=3].  "NuMin"->Automatic drops nu below
   nu_min = c_D x_max target^(1/(D-1)) with conservative c_D [errs safe: the dropped weight
   stays well below the target; D=4 cubic and clean, D=3 quadratic with a larger
   prefactor; c_D calibrated empirically]. *)
cNuMin = <|3 -> 0.45, 4 -> 0.45|>;   (* conservative: keeps the RENORMALISED error < target despite
   the |Sigma_D|/|det_ren| ~ 5 counterterm cancellation in D=4 (the cut errs to the sum magnitude) *)
numinAuto[xs_, tgt_, dim_, nm_] := Replace[nm, {
   Automatic :> If[NumericQ[tgt] && KeyExistsQ[cNuMin, dim], Floor[cNuMin[dim] xs tgt^(1/(dim - 1))], 0],
   n_?NumericQ :> n, _ :> 0}];
lminOf[numin_, dim_] := Max[0, Ceiling[numin - (dim/2 - 1)]];   (* smallest l with nuOf[l] >= numin *)

Options[PolyDetSigma] = {"Dimension" -> 4, WorkingPrecision -> Automatic, "Nsum" -> Automatic, "Pmax" -> 11, "Target" -> Automatic, "NuMin" -> 0, "Method" -> "Exact", "NSub" -> 8, "SubtractFrom" -> Automatic};
PolyDetSigma::frozen = "\"Method\" -> \"FrozenWall\" needs a case-b-like chain (positive plateau and FV m^2, ordered radii); falling back to \"Exact\".";
PolyDetSigma[m2_, radii_, OptionsPattern[]] := Module[
  {dim, prec, pmax, tgt, xs, Ns, lmin, lstar, meth, lminF, precLow, Ip, sub, lnRv, direct, addback, lowcorr, nu},
  dim = OptionValue["Dimension"];
  pmax = OptionValue["Pmax"];  tgt = OptionValue["Target"];
  xs = PolyDetXMax[m2, radii];
  (* "Method": "Exact" (default) | "FrozenWall" | Automatic.  FrozenWall = the
     thin-wall fast band (lnRallFrozen above): l >= 2 at MACHINE precision, l = 0,1
     kept on the exact chain.  The full subtraction machinery (Pmax, zeta add-back,
     lowcorr) is unchanged and evaluated at the INPUT precision -- high-precision
     inputs remain the standing requirement, as for the exact branch.  Automatic
     switches on for a case-b-like chain with x_max >= 30 and a numeric
     "Target" >= 1e-3 (there the frozen band error is well below that, falling as
     ~1/x_max^3 and ~NSub^-1.7). *)
  meth = Replace[OptionValue["Method"], Automatic :>
     If[NumericQ[tgt] && tgt >= 1/1000 && xs >= 30 && Length[radii] >= 3 &&
        N[m2[[1]]] > 0 && N[m2[[-1]]] > 0 && N[radii[[1]]/radii[[-1]]] >= 1/2,
        "FrozenWall", "Exact"]];
  If[meth === "FrozenWall" && ! (N[m2[[1]]] > 0 && N[m2[[-1]]] > 0),
     Message[PolyDetSigma::frozen]; meth = "Exact"];
  (* "Nsum"->n overrides; else "Target"->eps sets the cutoff
     nu_max ~ x_max eps^(-1/(pmax+2-dim)) (relative form; conservative, errs safe),
     floored just above the turning point x_max; Automatic target -> ~5 x_max (~1e-8). *)
  Ns = Replace[OptionValue["Nsum"], Automatic :>
       If[NumericQ[tgt],
          Max[Ceiling[xs] + 2, Ceiling[xs tgt^(-1/(pmax + 2 - dim))]],
          Max[100, Ceiling[5 xs]]]];
  (* WorkingPrecision->Automatic: MachinePrecision is ~4x faster and keeps full
     accuracy ONLY when (i) nu stays out of the nu>>xs regime (Ns<=3.5 xs, xs<=35)
     AND (ii) the ~xs^(pmax-dim+1) cancellation of the I_p sum against the zeta
     add-back still fits machine precision.  Without (ii) the high-pmax subtraction
     is silently corrupted for moderate xs (e.g. D=3 thin-wall at small Delta); then
     scale arbitrary precision with both xs and that cancellation depth. *)
  prec = Replace[OptionValue[WorkingPrecision], Automatic :>
         With[{cancel = (pmax - dim + 1) Log10[Max[xs, 2]]},
           If[Ns <= 7/2 xs && xs <= 35 && cancel <= 10,
              MachinePrecision, Ceiling[Max[40, 11/5 xs, cancel + 25]]]]];
  Ip = PolyDetIp[m2, radii, pmax];                          (* {I_1..I_pmax}, I_2=0 *)
  sub[v_] := Sum[Ip[[p]]/v^p, {p, 1, pmax}];
  lmin = lminOf[numinAuto[xs, tgt, dim, OptionValue["NuMin"]], dim];  (* bottom cut: drop l<lmin *)
  (* "SubtractFrom" -> nu* : apply the Laurent subtraction only for nu >= nu*
     (l >= lstar) and restore with HURWITZ zetas at offset lstar; below lstar the
     divergent orders are compensated by their EXPLICIT finite sums (I_1 Sum nu +
     I_3 H_{nu*} in D=4 -- the harmonic numbers of the thin-wall low sum).  The
     returned Sigma_D is IDENTICAL for every nu* (exact zeta-split identity); what
     changes is the ARITHMETIC: full-range subtraction evaluates the Laurent far
     below its radius of convergence (|sub| ~ x_max^pmax at nu = O(1)), routing the
     answer through a ~x_max^(pmax-dim+1) cancellation against the add-back, which
     is what forced high-precision inputs/accumulation; with nu* ~ x_max no
     intermediate exceeds ~x_max^dim and the assembly is machine-safe.
     Automatic: nu* = x_max when xs >= 30, else 0 (the proven moderate-x path). *)
  lstar = Replace[OptionValue["SubtractFrom"], {
     Automatic :> If[xs >= 30, lminOf[Ceiling[xs], dim], 0],
     n_?NumericQ /; n > 0 :> lminOf[n, dim], _ :> 0}];
  If[lstar > 0, lstar = Min[Max[lstar, lmin], Ns]];
  lnRv = If[meth === "FrozenWall",
     (* l = 0,1 (negative/zero modes) on the exact chain at chain-conditioned
        precision; the l >= 2 band on the frozen-wall evaluator at machine.
        With full-range subtraction (lstar = 0) the machine VALUES must be PADDED
        to high precision for the assembly arithmetic (the sub-turning-point
        summand reaches ~x_max^(pmax-dim+1) before cancelling against the zeta
        add-back); with the offset split (lstar > 0) a flat 40-digit pad is ample. *)
     precLow = Max[50, Ceiling[11/5 xs]];
     lminF = Max[lmin, 2];
     SetPrecision[
       Join[If[lmin < 2, lnRall[m2, radii, dim, Min[2, Ns], precLow, lmin], {}],
            If[Ns > lminF, lnRallFrozen[m2, radii, dim, Ns, lminF, OptionValue["NSub"]], {}]],
       If[lstar > 0, 40, Max[40, (pmax - dim + 1) Log10[Max[xs, 2]] + 25]]],
     lnRall[m2, radii, dim, Ns, prec, lmin]];               (* lnR_nu for l=lmin..Ns-1, two Bessel calls *)
  direct = Sum[nu = nuOf[l, dim];
     dDeg[l, dim] (lnRv[[l - lmin + 1]] - If[l >= lstar, sub[nu], 0]), {l, lmin, Ns - 1}];
  addback = If[lstar == 0,
     Sum[Ip[[p]] ZetaD[p, dim], {p, dim, pmax}],            (* p>=D convergent; p<D are subtracted (renorm) *)
     Sum[Ip[[p]] ZetaDOff[p, dim, lstar], {p, dim, pmax}] -
       (* divergent orders p<D: below lstar they were NOT subtracted, so their formal
          full-range subtraction (restored by ln R^r in PolyDetRen) is completed by
          the explicit finite low sums -- the harmonic numbers *)
       Sum[Ip[[p]] Sum[dDeg[l, dim]/nuOf[l, dim]^p, {l, 0, lstar - 1}], {p, 1, dim - 1}]];
  (* restore the Laurent that the skipped low modes would have subtracted (lmin=0 or
     lstar>0 -> 0; below lstar nothing is subtracted); the PHYSICAL low-nu sum
     Sum_{l<lmin} d_l ln R_l is DROPPED either way -- the controlled approximation *)
  lowcorr = If[lstar == 0,
     Sum[Ip[[p]] Sum[dDeg[l, dim]/nuOf[l, dim]^p, {l, 0, lmin - 1}], {p, 1, pmax}],
     0];
  direct + addback - lowcorr];

(* renormalisation counterterm moments on the polygon *)
i3Log[m2_, radii_] := -(1/32) Sum[                          (* I_3^ln: log-weighted I_3 moment (D=4) *)
   radii[[s]]^4 (m2[[s]]^2 - m2[[s + 1]]^2) Log[radii[[s]]], {s, Length[radii]}];
momD3[m2_, radii_] := (1/3) Sum[                            (* (1/3) Sum R_s^3 (m_s^2-m_{s+1}^2)  (D=3) *)
   radii[[s]]^3 (m2[[s]] - m2[[s + 1]]), {s, Length[radii]}];

PolyDetRen::dim = "Renormalised add-back implemented for Dimension 3 and 4 only (got `1`).";
PolyDetRen::mlfv = "\"ZeroRemoval\" is not defined for a massless false vacuum (masses2[[-1]] == 0); computing with \"ZeroRemoval\" -> False.";
Options[PolyDetRen] = Join[{"Scale" -> 1, "ZeroRemoval" -> True}, Options[PolyDetSigma]];
PolyDetRen[m2_, radii_, opts : OptionsPattern[]] := Module[
  {dim, mu, prec, sig, i3, muT, ct, zr, zrPrec, lmin},
  dim = OptionValue["Dimension"];  mu = OptionValue["Scale"];
  prec = Replace[OptionValue[WorkingPrecision], Automatic -> MachinePrecision];  (* ct is a cheap scalar *)
  lmin = lminOf[numinAuto[PolyDetXMax[m2, radii], OptionValue["Target"], dim, OptionValue["NuMin"]], dim];
  sig = PolyDetSigma[m2, radii, FilterRules[{opts}, Options[PolyDetSigma]]];
  (* l=1 translational zero modes: PolyDetSigma includes the (artifact) term
     d_1 ln|R_{l=1}|, which -> -inf as the polygon refines (R_{l=1}->0).  Replace it
     with the removed-zero d_1 ln|delta R_1(inf)| = d_1 ln|R_1'/mhat^2|, d_1 = D.
     The mhat^2 (= m_N^2) carries the dimensionful m^D of the rate prefactor. *)
  If[TrueQ[OptionValue["ZeroRemoval"]] && TrueQ[m2[[-1]] == 0], Message[PolyDetRen::mlfv]];
  zr = If[TrueQ[OptionValue["ZeroRemoval"]] && !TrueQ[m2[[-1]] == 0],
     zrPrec = Replace[OptionValue[WorkingPrecision],
        {Automatic | MachinePrecision :> Max[50, Ceiling[2 PolyDetXMax[m2, radii]]]}];
     (* subtract the bulk-sum l=1 artifact ONLY if it is still in the sum (lmin<=1); a bottom
        cut with lmin>1 already dropped it, so we just add the removed-zero R_1' contribution *)
     dDeg[1, dim] (Log[Abs[PolyDetR1Prime[m2, radii, "Dimension" -> dim, WorkingPrecision -> zrPrec]/m2[[-1]]]]
                   - If[lmin <= 1, PolyDetLnR[dim/2, m2, radii, zrPrec], 0]),
     0];
  ct = Which[
    dim == 4,                                                (* MSbar; 1/eps coeff = I_3 (cancels vs action) *)
      i3 = PolyDetIp[m2, radii, 3][[3]];
      muT = Sqrt[4 Pi] Exp[-EulerGamma/2] SetPrecision[mu, prec];
      i3 (3/4 + EulerGamma + Log[muT/2]) + i3Log[m2, radii],
    dim == 3,                                                (* super-renormalisable: the only divergent order
        (p=1) is a POWER (linear) divergence -- scaleless, hence 0 in MSbar dim reg -- so there is no
        counterterm and ln det_ren = Sigma_3 directly. *)
      0,
    True, Message[PolyDetRen::dim, dim]; Return[$Failed]];
  sig + ct + zr];

(* {masses2, radii} from a FindBounce BounceFunction; Vpp = V''(phi).
   case b (R_1>0, field waits at TV): {V''(TV), ddVL..., V''(FV)}, radii.
   case a (R_1=0, field rolls from rho=0): {ddVL..., V''(FV)}, Rest[radii]. *)
PolyDetInputs[Vpp_, bf_] := Module[{nodes, radii, ddVL},
  nodes = Flatten[bf["Path"]];  radii = bf["Radii"];
  (* polygonal (extended-polygonal) wall masses: the per-segment AVERAGE of V''(phi),
     m^2_s = [V'(phi_{s+1}) - V'(phi_s)]/(phi_{s+1} - phi_s), evaluated by Simpson's rule on
     Vpp = V'' (exact for a quartic, O(h^4) in general).  We do NOT use
     bf["CoefficientsExtension"], a different FindBounce parametrisation that does not equal
     the segment average and converges in N more slowly (it cost a ~5% bias in the
     thin-wall determinant). *)
  ddVL = Table[(Vpp[nodes[[i]]] + 4 Vpp[(nodes[[i]] + nodes[[i + 1]])/2] + Vpp[nodes[[i + 1]]])/6,
     {i, Length[nodes] - 1}];
  If[Abs[N@radii[[1]]] < 10^-10,
    {Join[ddVL, {Vpp[nodes[[-1]]]}], Rest[radii]},
    {Join[{Vpp[nodes[[1]]]}, ddVL, {Vpp[nodes[[-1]]]}], radii}]];

(* one-call RENORMALISED determinant from a BounceFunction + the potential;
   the spacetime dimension is read from the bounce (bf["Dimension"]). *)
PolyDet[Vexpr_, field_, bf_, opts : OptionsPattern[PolyDetRen]] := Module[{d2, inp},
  d2 = D[Vexpr, {field, 2}];
  inp = PolyDetInputs[(d2 /. field -> # &), bf];
  If[inp === $Failed, $Failed,
     PolyDetRen[inp[[1]], inp[[2]], opts, "Dimension" -> bf["Dimension"]]]];

(* ======================================================================
   COUPLED MULTI-FIELD (matrix m^2) -- no internal symmetry.
   N_f coupled scalars: each segment mass masses2[[s]] is an N_f x N_f matrix.
   Diagonal masses reduce to a sum of scalar PolyDet results; the chain is
   basis-independent and keeps det T^s = 1.  Implementation in Multiscalar.wl.
   ====================================================================== *)

