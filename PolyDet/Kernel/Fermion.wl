(* ::Package:: *)

(* PolyDet` -- Dirac fermion (Yukawa) section.  Loaded by Kernel/PolyDet.wl from inside
   PolyDet`Private` (Get'd by the loader); NOT a standalone package -- it has no
   BeginPackage/Begin of its own and relies on the loader's context.
   Signed shell masses M_s=m+y phi_s, first-order adjacent-order Dirac chain,
   PolyDet*Psi (GPF Laurent through p=15, Z^psi add-back). Self-contained. *)

(* ======================================================================
   DIRAC FERMION (Yukawa), D=4 -- SIGNED shell masses M_s = m + y phi_s.
   First-order radial Dirac transfer chain (ADJACENT Bessel orders nu, nu+1);
   ln R_nu^psi = 2 ln|r_nu|, degeneracy d_nu^psi = nu(nu+1),
   zeta add-back Z^psi(p) = zeta(p-2)+zeta(p-1).
   ====================================================================== *)

(* hardcoded universal per-matching Laurent coefficients of ln R_nu^psi (= 2 ln G^psi),
   GPF[p] = {{{i,j}, num/den}, ...}, I_p^psi(xs,xp) = sum (num/den) xs^i xp^j with
   xs = M_s R_s, xp = M_{s+1} R_s (SIGNED real); exact rationals, p<=15. *)
GPMAXPsi = 15;
GPF = <|
  1 -> {{{2,0},1/2}, {{0,2},-1/2}},
  2 -> {{{2,0},-1/2}, {{1,1},1/2}},
  3 -> {{{4,0},-1/16}, {{0,4},1/16}, {{2,0},1/2}, {{1,1},-1/2}},
  4 -> {{{4,0},1/4}, {{3,1},-1/8}, {{2,2},-1/16}, {{1,3},-1/8}, {{0,4},1/16}, {{2,0},-1/2}, {{1,1},1/2}},
  5 -> {{{6,0},1/48}, {{0,6},-1/48}, {{4,0},-11/16}, {{3,1},1/2}, {{2,2},1/8}, {{0,4},1/16}, {{2,0},1/2}, {{1,1},-1/2}},
  6 -> {{{6,0},-1/6}, {{5,1},1/16}, {{4,2},1/32}, {{3,3},1/24}, {{2,4},1/32}, {{1,5},1/16}, {{0,6},-1/16}, {{4,0},13/8}, {{3,1},-11/8}, {{2,2},-3/16}, {{1,3},-1/8}, {{0,4},1/16}, {{2,0},-1/2}, {{1,1},1/2}},
  7 -> {{{8,0},-5/512}, {{0,8},5/512}, {{6,0},5/6}, {{5,1},-1/2}, {{4,2},-5/32}, {{3,3},-1/8}, {{2,4},-1/32}, {{1,5},1/8}, {{0,6},-7/48}, {{4,0},-57/16}, {{3,1},13/4}, {{2,2},1/4}, {{0,4},1/16}, {{2,0},1/2}, {{1,1},-1/2}},
  8 -> {{{8,0},1/8}, {{7,1},-5/128}, {{6,2},-5/256}, {{5,3},-3/128}, {{4,4},-9/512}, {{3,5},-3/128}, {{2,6},-5/256}, {{1,7},-5/128}, {{0,8},29/512}, {{6,0},-27/8}, {{5,1},5/2}, {{4,2},1/2}, {{3,3},5/16}, {{2,4},1/16}, {{1,5},5/16}, {{0,6},-5/16}, {{4,0},15/2}, {{3,1},-57/8}, {{2,2},-5/16}, {{1,3},-1/8}, {{0,4},1/16}, {{2,0},-1/2}, {{1,1},1/2}},
  9 -> {{{10,0},7/1280}, {{0,10},-7/1280}, {{8,0},-245/256}, {{7,1},1/2}, {{6,2},11/64}, {{5,3},5/32}, {{4,4},9/128}, {{3,5},1/32}, {{2,6},-1/64}, {{1,7},-3/16}, {{0,8},59/256}, {{6,0},195/16}, {{5,1},-81/8}, {{4,2},-21/16}, {{3,3},-2/3}, {{2,4},-1/16}, {{1,5},5/8}, {{0,6},-31/48}, {{4,0},-247/16}, {{3,1},15/1}, {{2,2},3/8}, {{0,4},1/16}, {{2,0},1/2}, {{1,1},-1/2}},
  10 -> {{{10,0},-1/10}, {{9,1},7/256}, {{8,2},7/512}, {{7,3},1/64}, {{6,4},3/256}, {{5,5},9/640}, {{4,6},3/256}, {{3,7},1/64}, {{2,8},7/512}, {{1,9},7/256}, {{0,10},-13/256}, {{8,0},23/4}, {{7,1},-245/64}, {{6,2},-117/128}, {{5,3},-11/16}, {{4,4},-53/256}, {{3,5},-1/8}, {{2,6},-9/128}, {{1,7},-47/64}, {{0,8},209/256}, {{6,0},-329/8}, {{5,1},585/16}, {{4,2},99/32}, {{3,3},11/8}, {{2,4},3/32}, {{1,5},21/16}, {{0,6},-21/16}, {{4,0},251/8}, {{3,1},-247/8}, {{2,2},-7/16}, {{1,3},-1/8}, {{0,4},1/16}, {{2,0},-1/2}, {{1,1},1/2}},
  11 -> {{{12,0},-7/2048}, {{0,12},7/2048}, {{10,0},273/256}, {{9,1},-1/2}, {{8,2},-93/512}, {{7,3},-11/64}, {{6,4},-23/256}, {{5,5},-9/128}, {{4,6},-7/256}, {{3,7},1/64}, {{2,8},23/512}, {{1,9},29/128}, {{0,10},-5/16}, {{8,0},-15403/512}, {{7,1},23/1}, {{6,2},123/32}, {{5,3},161/64}, {{4,4},65/128}, {{3,5},9/64}, {{2,6},-3/32}, {{1,7},-81/32}, {{0,8},1383/512}, {{6,0},1597/12}, {{5,1},-987/8}, {{4,2},-219/32}, {{3,3},-11/4}, {{2,4},-3/32}, {{1,5},21/8}, {{0,6},-127/48}, {{4,0},-1013/16}, {{3,1},251/4}, {{2,2},1/2}, {{0,4},1/16}, {{2,0},1/2}, {{1,1},-1/2}},
  12 -> {{{12,0},1/12}, {{11,1},-21/1024}, {{10,2},-21/2048}, {{9,3},-35/3072}, {{8,4},-35/4096}, {{7,5},-5/512}, {{6,6},-25/3072}, {{5,7},-5/512}, {{4,8},-35/4096}, {{3,9},-35/3072}, {{2,10},-21/2048}, {{1,11},-21/1024}, {{0,12},281/6144}, {{10,0},-35/4}, {{9,1},1365/256}, {{8,2},725/512}, {{7,3},149/128}, {{6,4},7/16}, {{5,5},73/256}, {{4,6},3/32}, {{3,7},17/128}, {{2,8},87/512}, {{1,9},171/128}, {{0,10},-207/128}, {{8,0},9245/64}, {{7,1},-15403/128}, {{6,2},-3627/256}, {{5,3},-1075/128}, {{4,4},-587/512}, {{3,5},-63/128}, {{2,6},-63/256}, {{1,7},-1059/128}, {{0,8},4407/512}, {{6,0},-5033/12}, {{5,1},1597/4}, {{4,2},233/16}, {{3,3},263/48}, {{2,4},1/8}, {{1,5},85/16}, {{0,6},-85/16}, {{4,0},509/4}, {{3,1},-1013/8}, {{2,2},-9/16}, {{1,3},-1/8}, {{0,4},1/16}, {{2,0},-1/2}, {{1,1},1/2}},
  13 -> {{{14,0},33/14336}, {{0,14},-33/14336}, {{12,0},-2387/2048}, {{11,1},1/2}, {{10,2},193/1024}, {{9,3},93/512}, {{8,4},209/2048}, {{7,5},23/256}, {{6,6},25/512}, {{5,7},7/256}, {{4,8},1/2048}, {{3,9},-23/512}, {{2,10},-67/1024}, {{1,11},-65/256}, {{0,12},801/2048}, {{10,0},39391/640}, {{9,1},-175/4}, {{8,2},-4375/512}, {{7,3},-405/64}, {{6,4},-441/256}, {{5,5},-15/16}, {{4,6},-45/256}, {{3,7},17/64}, {{2,8},273/512}, {{1,9},27/4}, {{0,10},-9777/1280}, {{8,0},-20991/32}, {{7,1},9245/16}, {{6,2},3087/64}, {{5,3},1709/64}, {{4,4},157/64}, {{3,5},33/64}, {{2,6},-27/64}, {{1,7},-837/32}, {{0,8},1717/64}, {{6,0},62377/48}, {{5,1},-5033/4}, {{4,2},-121/4}, {{3,3},-87/8}, {{2,4},-1/8}, {{1,5},85/8}, {{0,6},-511/48}, {{4,0},-4083/16}, {{3,1},509/2}, {{2,2},5/8}, {{0,4},1/16}, {{2,0},1/2}, {{1,1},-1/2}},
  14 -> {{{14,0},-1/14}, {{13,1},33/2048}, {{12,2},33/4096}, {{11,3},9/1024}, {{10,4},27/4096}, {{9,5},15/2048}, {{8,6},25/4096}, {{7,7},25/3584}, {{6,8},25/4096}, {{5,9},15/2048}, {{4,10},27/4096}, {{3,11},9/1024}, {{2,12},33/4096}, {{1,13},33/2048}, {{0,14},-85/2048}, {{12,0},99/8}, {{11,1},-7161/1024}, {{10,2},-4089/2048}, {{9,3},-883/512}, {{8,4},-3007/4096}, {{7,5},-35/64}, {{6,6},-225/1024}, {{5,7},-9/64}, {{4,8},-303/4096}, {{3,9},-129/512}, {{2,10},-709/2048}, {{1,11},-2143/1024}, {{0,12},5623/2048}, {{10,0},-12555/32}, {{9,1},39391/128}, {{8,2},11391/256}, {{7,3},3913/128}, {{6,4},1547/256}, {{5,5},761/256}, {{4,6},117/256}, {{3,7},143/128}, {{2,8},27/16}, {{1,9},8049/256}, {{0,10},-8751/256}, {{8,0},91823/32}, {{7,1},-20991/8}, {{6,2},-2501/16}, {{5,3},-2647/32}, {{4,4},-653/128}, {{3,5},-15/8}, {{2,6},-59/64}, {{1,7},-2597/32}, {{0,8},10559/128}, {{6,0},-15929/4}, {{5,1},62377/16}, {{4,2},1981/32}, {{3,3},173/8}, {{2,4},5/32}, {{1,5},341/16}, {{0,6},-341/16}, {{4,0},4089/8}, {{3,1},-4083/8}, {{2,2},-11/16}, {{1,3},-1/8}, {{0,4},1/16}, {{2,0},-1/2}, {{1,1},1/2}},
  15 -> {{{16,0},-429/262144}, {{0,16},429/262144}, {{14,0},1287/1024}, {{13,1},-1/2}, {{12,2},-793/4096}, {{11,3},-193/1024}, {{10,4},-453/4096}, {{9,5},-209/2048}, {{8,6},-255/4096}, {{7,7},-25/512}, {{6,8},-95/4096}, {{5,9},-1/2048}, {{4,10},75/4096}, {{3,11},67/1024}, {{2,12},331/4096}, {{1,13},281/1024}, {{0,14},-955/2048}, {{12,0},-342199/3072}, {{11,1},297/4}, {{10,2},16533/1024}, {{9,3},39467/3072}, {{8,4},4259/1024}, {{7,5},1397/512}, {{6,6},1225/1536}, {{5,7},157/512}, {{4,8},-59/1024}, {{3,9},-1111/1024}, {{2,10},-1651/1024}, {{1,11},-7363/512}, {{0,12},17709/1024}, {{10,0},1498279/640}, {{9,1},-62775/32}, {{8,2},-54145/256}, {{7,3},-4415/32}, {{6,4},-2537/128}, {{5,5},-5719/640}, {{4,6},-51/64}, {{3,7},185/64}, {{2,8},1323/256}, {{1,9},17853/128}, {{0,10},-18925/128}, {{8,0},-6264305/512}, {{7,1},91823/8}, {{6,2},7859/16}, {{5,3},8085/32}, {{4,4},667/64}, {{3,5},61/32}, {{2,6},-7/4}, {{1,7},-3981/16}, {{0,8},128761/512}, {{6,0},24231/2}, {{5,1},-47787/4}, {{4,2},-4017/32}, {{3,3},-1033/24}, {{2,4},-5/32}, {{1,5},341/8}, {{0,6},-2047/48}, {{4,0},-16369/16}, {{3,1},4089/4}, {{2,2},3/4}, {{0,4},1/16}, {{2,0},1/2}, {{1,1},-1/2}}
|>;
gpPsi[p_, xs_, xp_] := Total[(#[[2]] xs^#[[1, 1]] xp^#[[1, 2]]) & /@ GPF[p]];   (* g_p^psi(xs,xp) *)

PolyDetIpPsi[M_, radii_, pmax_] := Table[
   Sum[gpPsi[p, M[[s]] radii[[s]], M[[s + 1]] radii[[s]]], {s, Length[radii]}], {p, 1, pmax}];

PolyDetXMaxPsi[M_, radii_] := Max@Table[
   Max[Abs[M[[s]] radii[[s]]], Abs[M[[s + 1]] radii[[s]]]], {s, Length[radii]}];

(* vectorised per-multipole ln R_mu^psi for l = 0..Ns-1, mu = l + dim/2 - 1 (two Bessel
   calls; the D-general orders are HALF-INTEGER in odd dim -- the chain itself is
   D-blind).  The per-segment T^s = Minv(M_{s+1}) . Mmat(M_s) at R_s, z=|M|R,
   e=sign(M), Mmat = {{I_mu, K_mu}, {-e I_{mu+1}, e K_{mu+1}}},
   Minv = (M R){{e K_{mu+1},-K_mu},{e I_{mu+1},I_mu}}. *)
lnRallPsi[M_, radii_, Ns_, prec_, dim_: 4] := Module[
  {MM, rr, Nm, sgn, av, bv, nus, orders, allargs, om, am, Itab, Ktab, No},
  MM = SetPrecision[M, prec];  rr = SetPrecision[radii, prec];  Nm = Length[radii];  sgn = Sign[MM];
  av = Table[Abs[MM[[s]]] rr[[s]], {s, Nm}];       (* |M_s|   R_s *)
  bv = Table[Abs[MM[[s + 1]]] rr[[s]], {s, Nm}];   (* |M_{s+1}| R_s *)
  nus = Table[l + dim/2 - 1, {l, 0, Ns - 1}];      (* mu_l; consecutive, spaced by 1 *)
  orders = Append[nus, nus[[-1]] + 1];  No = Ns + 1;   (* row i = mu_i, row i+1 = mu_i + 1 *)
  allargs = Join[av, bv];
  om = Transpose[ConstantArray[orders, Length@allargs]];  am = ConstantArray[allargs, No];
  Itab = BesselI[om, am];  Ktab = BesselK[om, am];
  Quiet[Table[Module[{nu = nus[[i]], Tlist},        (* row i = nu, row i+1 = nu+1 *)
    Tlist = Table[Module[
       {es = sgn[[s]], es1 = sgn[[s + 1]], MR = MM[[s + 1]] rr[[s]],
        I0a = Itab[[i, s]], I1a = Itab[[i + 1, s]], K0a = Ktab[[i, s]], K1a = Ktab[[i + 1, s]],
        I0b = Itab[[i, Nm + s]], I1b = Itab[[i + 1, Nm + s]], K0b = Ktab[[i, Nm + s]], K1b = Ktab[[i + 1, Nm + s]]},
      MR {{es1 K1b I0a + es K0b I1a, es1 K1b K0a - es K0b K1a},
          {es1 I1b I0a - es I0b I1a, es1 I1b K0a + es I0b K1a}}], {s, Nm}];
    2 Re[Log[(Abs[MM[[-1]]]/Abs[MM[[1]]])^nu (Dot @@ Reverse[Tlist])[[1, 1]]]]], {i, Length@nus}],
   {General::munfl}]];

(* ---- FrozenWall fermion band ----
   Case-b thin wall: exact Bessels ONLY for the plateau spinor at R_0 and the FV
   decomposition (4 per mode, N-independent); across the wall shells the
   first-order Dirac system u' = A u, A = {{l/rho, -M}, {-M, -(l+D-1)/rho}}, is
   frozen per sub-interval at the LOG-MEAN c = Log[b/a]/h of 1/rho (which makes
   the trace part EXACT: the WEIGHTLESS-component measure factor telescopes to (R_0/R_N)^(1/2),
   pulled out analytically) -> P = cosh(q h) + sinh(q h)/q B, with
   B = {{(mu+1/2) c, -M}, {-M, -(mu+1/2) c}} and q^2 = (mu+1/2)^2 c^2 + M^2.
   q^2 > 0 ALWAYS (M real): the fermion frozen band is never oscillatory.  The
   spinor is continuous at matchings (first order), so shells just chain.
   Signed-mass bookkeeping identical to lnRallPsi (z = |M| rho, e = sign(M):
   regular pair (I_mu, -e I_{mu+1}), growing amplitude
   A_N = z (c1 K_{mu+1} - e c2 K_mu), the Minv row of the exact chain). *)
lnRallPsiFrozen[a___] := Quiet[lnRallPsiFrozenC[a], {General::munfl}];
lnRallPsiFrozenC[M_, radii_, dim_, Ns_, lminF_, nsub_, precB_: 40] := Module[
  {MM, rr, m1a, mNa, e1, eN, R0, RN, nus, x0, xN, I0v, I1v, K0v, K1v,
   w2v, lnu0, lnKN, lnMeas, Mw, subgrid},
  MM = N[M];  rr = N[radii];
  m1a = Abs[SetPrecision[M[[1]], precB]];  mNa = Abs[SetPrecision[M[[-1]], precB]];
  e1 = Sign[N[M[[1]]]];  eN = Sign[N[M[[-1]]]];
  R0 = SetPrecision[radii[[1]], precB];    RN = SetPrecision[radii[[-1]], precB];
  nus = Table[l + dim/2 - 1, {l, lminF, Ns - 1}];
  x0 = m1a R0;  xN = mNa RN;
  I0v = BesselI[nus, x0];  I1v = BesselI[nus + 1, x0];     (* plateau: orders mu, mu+1 *)
  K0v = BesselK[nus, xN];  K1v = BesselK[nus + 1, xN];     (* FV: same pair *)
  w2v = N[-e1 I1v/I0v];                                    (* c2/c1 at R_0, O(1) *)
  lnu0 = N[Log[I0v]];
  lnKN = N[Log[K0v]];                                      (* K_mu magnitude, additive *)
  (* WEIGHTLESS components (the package chain convention, pair (I_mu, -e I_{mu+1})
     with no rho^-(D/2-1) weight): their system has rows {mu/rho, -M}, {-M, -(mu+1)/rho}
     -> trace -1/rho (the weight shifts only the trace; B is the trace-free part),
     so the exact telescoped measure factor is (R_0/R_N)^(1/2) *)
  lnMeas = N[-(1/2) Log[RN/R0]];
  Mw = MM[[2 ;; -2]];                                      (* wall-shell SIGNED masses *)
  subgrid = Table[Module[{a, b, h},
      a = rr[[s]] + (rr[[s + 1]] - rr[[s]]) (j - 1)/nsub;
      b = rr[[s]] + (rr[[s + 1]] - rr[[s]]) j/nsub;  h = b - a;
      {h, Log[b/a]/h}],
    {s, 1, Length[rr] - 1}, {j, nsub}];                    (* shell s: [rr[s], rr[s+1]], mass Mw[[s]] *)
  Table[Module[{mu = N[nus[[i]]], w1 = 1., w2, h, c, q, ch, sh, bd, lnscale = 0., mag, Kr},
     w2 = w2v[[i]];
     Do[Do[
        {h, c} = subgrid[[s, j]];
        bd = (mu + 0.5) c;
        q = Sqrt[bd^2 + Mw[[s]]^2];
        ch = Cosh[q h];  sh = Sinh[q h]/q;
        {w1, w2} = {(ch + sh bd) w1 - sh Mw[[s]] w2, -sh Mw[[s]] w1 + (ch - sh bd) w2},
        {j, nsub}], {s, Length[Mw]}];
     mag = Max[Abs[w1], Abs[w2]];
     If[mag > 10.^15, w1 /= mag; w2 /= mag; lnscale = Log[mag]];
     Kr = N[K1v[[i]]/K0v[[i]]];                            (* K_{mu+1}/K_mu, O(1) *)
     (* ln|r_mu| = mu ln(|M_N|/|M_1|) + ln|A_N|;  ln R^psi = 2 ln|r_mu| *)
     2 (nus[[i]] N[Log[mNa/m1a]] + N[Log[xN]] + lnKN[[i]] + lnu0[[i]] + lnMeas + lnscale +
        Log[Abs[w1 Kr - eN w2]])],
   {i, Length[nus]}]];

PolyDetLnRPsi[nu_, M_, radii_, prec_: 60] := First@Quiet@Module[{},
  (* single nu via the same chain (small wrapper around the vectorised core at one order) *)
  Module[{MM = SetPrecision[M, prec], sgn, Tlist},
    sgn = Sign[MM];
    Tlist = Table[Module[
       {es = sgn[[s]], es1 = sgn[[s + 1]], MR = MM[[s + 1]] radii[[s]],
        a = Abs[MM[[s]]] SetPrecision[radii[[s]], prec], b = Abs[MM[[s + 1]]] SetPrecision[radii[[s]], prec]},
      MR {{es1 BesselK[nu + 1, b] BesselI[nu, a] + es BesselK[nu, b] BesselI[nu + 1, a],
           es1 BesselK[nu + 1, b] BesselK[nu, a] - es BesselK[nu, b] BesselK[nu + 1, a]},
          {es1 BesselI[nu + 1, b] BesselI[nu, a] - es BesselI[nu, b] BesselI[nu + 1, a],
           es1 BesselI[nu + 1, b] BesselK[nu, a] + es BesselI[nu, b] BesselK[nu + 1, a]}}], {s, Length[radii]}];
    {2 Re[Log[(Abs[MM[[-1]]]/Abs[MM[[1]]])^nu (Dot @@ Reverse[Tlist])[[1, 1]]]]}]];

(* D-general assembly.  The chain and the
   I_p^psi tables are D-INDEPENDENT; "Dimension" enters via the orders mu = l+D/2-1,
   the degeneracy weight deg_psi/2 = 2^(floor(D/2)-1) Binomial[l+D-2, l]
   (D=4: (l+1)(l+2) = nu(nu+1); D=3: l+1 = mu+1/2), and the add-back:
     D=4: Z^psi(p) = zeta(p-2)+zeta(p-1), convergent p>=4 (p<=3 = the renormalised running);
     D=3: Hurwitz at the half-integer offset,
          I_2 (1/2) zeta(2,1/2) + sum_{p>=3} I_p [zeta(p-1,1/2) + zeta(p,1/2)/2];
          the dropped pieces are the scaleless power divergence (I_1 constant) and the
          O(1/N) jump-artifact log  I_2 + I_1/2 = -(1/4) sum R_s^2 (dM_s)^2  (scheme:
          the continuum D=3 determinant is finite and scale-free). *)
PolyDetSigmaPsi::dim = "Fermion add-back implemented for Dimension 3 and 4 only (got `1`).";
PolyDetSigmaPsi::frozen = "\"Method\" -> \"FrozenWall\" needs a case-b-like chain (nonzero end masses, ordered radii); falling back to \"Exact\".";
Options[PolyDetSigmaPsi] = {WorkingPrecision -> Automatic, "Nsum" -> Automatic, "Pmax" -> 7, "Target" -> Automatic, "Dimension" -> 4, "Method" -> "Exact", "NSub" -> 8, "SubtractFrom" -> Automatic};
PolyDetSigmaPsi[M_, radii_, OptionsPattern[]] := Module[
  {dim, prec, pmax, tgt, xs, Ns, meth, lstar, precLow, lam, Ip, sub, wdeg, lnRv, direct, addback},
  dim = OptionValue["Dimension"];
  If[!MemberQ[{3, 4}, dim], Message[PolyDetSigmaPsi::dim, dim]; Return[$Failed]];
  pmax = OptionValue["Pmax"];  tgt = OptionValue["Target"];
  xs = PolyDetXMaxPsi[M, radii];
  (* "Method"/"NSub"/"SubtractFrom": the FrozenWall fast band and the offset (Hurwitz)
     assembly, exactly as for the scalar PolyDetSigma (see Scalar.wl); the fermion band
     is never oscillatory (q^2 = (mu+1/2)^2 c^2 + M^2 > 0), and there are no zero modes,
     so all l >= 2 run frozen with l = 0, 1 on the exact chain for uniformity. *)
  meth = Replace[OptionValue["Method"], Automatic :>
     If[NumericQ[tgt] && tgt >= 1/1000 && xs >= 30 && Length[radii] >= 3 &&
        N[Abs[M[[1]]]] > 0 && N[Abs[M[[-1]]]] > 0 && N[radii[[1]]/radii[[-1]]] >= 1/2,
        "FrozenWall", "Exact"]];
  If[meth === "FrozenWall" && ! (N[Abs[M[[1]]]] > 0 && N[Abs[M[[-1]]]] > 0),
     Message[PolyDetSigmaPsi::frozen]; meth = "Exact"];
  Ns = Replace[OptionValue["Nsum"], Automatic :>     (* deg^psi ~ nu^(dim-2): tail law as the scalar *)
       If[NumericQ[tgt], Max[Ceiling[xs] + 2, Ceiling[xs tgt^(-1/(pmax + 2 - dim))]], Max[100, Ceiling[3 xs]]]];
  prec = Replace[OptionValue[WorkingPrecision], Automatic :>
         With[{cancel = (pmax - 1) Log10[Max[xs, 2]]},      (* I_p^psi ~ x^(p+1): cancellation depth *)
           If[Ns <= 7/2 xs && xs <= 35 && cancel <= 10, MachinePrecision, Ceiling[Max[40, 11/5 xs, cancel + 25]]]]];
  lstar = Replace[OptionValue["SubtractFrom"], {
     Automatic :> If[xs >= 30, Max[0, Ceiling[Ceiling[xs] - (dim/2 - 1)]], 0],
     n_?NumericQ /; n > 0 :> Max[0, Ceiling[n - (dim/2 - 1)]], _ :> 0}];
  lstar = Min[lstar, Ns];
  Ip = PolyDetIpPsi[M, radii, pmax];
  sub[v_] := Sum[Ip[[p]]/v^p, {p, 1, pmax}];
  wdeg[l_] := 2^(Floor[dim/2] - 1) Binomial[l + dim - 2, l];          (* deg_psi/2 *)
  lnRv = If[meth === "FrozenWall",
     precLow = Max[50, Ceiling[11/5 xs]];
     SetPrecision[
       Join[lnRallPsi[M, radii, Min[2, Ns], precLow, dim],
            If[Ns > 2, lnRallPsiFrozen[M, radii, dim, Ns, 2, OptionValue["NSub"]], {}]],
       If[lstar > 0, 40, Max[40, (pmax - 1) Log10[Max[xs, 2]] + 25]]],
     lnRallPsi[M, radii, Ns, prec, dim]];
  direct = Sum[wdeg[l] (lnRv[[l + 1]] - If[l >= lstar, sub[l + dim/2 - 1], 0]), {l, 0, Ns - 1}];
  (* add-back: full spectral zetas at lstar = 0; Hurwitz tails + explicit low-sum
     compensation of the not-added-back (divergent/dropped) orders at lstar > 0
     (value-identical: the zeta-split identity, as in the scalar offset assembly) *)
  lam[p_] := Sum[wdeg[l]/(l + dim/2 - 1)^p, {l, 0, lstar - 1}];        (* finite low sums *)
  addback = If[dim == 4,
     If[lstar == 0,
        Sum[Ip[[p]] (Zeta[p - 2] + Zeta[p - 1]), {p, 4, pmax}],        (* Z^psi(p), convergent p>=4 *)
        Sum[Ip[[p]] (HurwitzZeta[p - 2, lstar + 1] + HurwitzZeta[p - 1, lstar + 1]), {p, 4, pmax}] -
          Sum[Ip[[p]] lam[p], {p, 1, 3}]],
     If[lstar == 0,
        Ip[[2]] HurwitzZeta[2, 1/2]/2 +
          Sum[Ip[[p]] (HurwitzZeta[p - 1, 1/2] + HurwitzZeta[p, 1/2]/2), {p, 3, pmax}],
        Ip[[2]] (HurwitzZeta[2, lstar + 1/2]/2 - Sum[1/(l + 1/2), {l, 0, lstar - 1}]) +
          Sum[Ip[[p]] (HurwitzZeta[p - 1, lstar + 1/2] + HurwitzZeta[p, lstar + 1/2]/2), {p, 3, pmax}] -
          Ip[[1]] lam[1]]];
  direct + addback];

Options[PolyDetRenPsi] = Join[{"Scale" -> 1}, Options[PolyDetSigmaPsi]];
PolyDetRenPsi[M_, radii_, opts : OptionsPattern[]] := Module[{dim, sig, beta},
  dim = OptionValue["Dimension"];
  sig = PolyDetSigmaPsi[M, radii, FilterRules[{opts}, Options[PolyDetSigmaPsi]]];
  If[sig === $Failed, Return[$Failed]];
  (* D=4: beta^psi = I_2^psi + I_3^psi running; D=3 is scale-free (odd-D dim reg:
     no log divergence -- the polygon's jump artifact is dropped by the Sigma scheme) *)
  If[dim == 4,
     beta = With[{Ip = PolyDetIpPsi[M, radii, 3]}, Ip[[2]] + Ip[[3]]];
     sig + beta Log[OptionValue["Scale"]],
     sig]];

(* ---- fermion gradient moment I_g^psi (the wave-function running the polygon discards) ----
   The first-order polygonal fermion chain is delta-free at matchings, so beta^psi =
   I_2^psi + I_3^psi is PURE POTENTIAL (the gradient (dm_psi)^2 = (M'[phi] phibar')^2 never
   enters the transfer chain).  The missing piece of the fermion determinant's running is the
   continuum GRADIENT moment -- the fermion analogue of the gauge I_g (Rnu_gauge.tex) --
     I_g^psi = -Int rho^{D-1} (d/drho M(phibar(rho)))^2 drho ,
   fixed by the same wave-function bridge gamma_phi Int (d phibar)^2 = -1/4 I_g^psi with
   gamma_phi = 2 y^2/16 pi^2.  The full fermion determinant then runs with BOTH moments,
   dSigma^psi/dln mu = beta^psi + I_g^psi, exactly as the gauge runs with I_3^eff + I_g.
   Adding I_g^psi restores the Baratella cancellation the polygon alone cannot (see
   Examples/RGImprovedRateYukawa.wl).  Mfun is the SIGNED fermion mass as a function of the
   field, M(phi) = m + y phi (a pure function, matching PolyDetDecayRate's "Mass"); bf a
   FindBounce BounceFunction.  Needs the smooth profile bf["Bounce"], sampled + interpolated
   (NIntegrate of M'[prof[rho]] prof'[rho] on the raw piecewise profile returns 0). *)
Options[PolyDetIgPsi] = {"SamplePoints" -> 4000, "RadiusFactor" -> 5/4};
PolyDetIgPsi[bf_, Mfun_, OptionsPattern[]] := Module[
   {dim, prof, Rfar, rs, ip, mdot},
   dim = bf["Dimension"];
   prof = bf["Bounce"][[1]];
   Rfar = OptionValue["RadiusFactor"] Max[bf["Radii"]];
   rs = Range[0., Rfar, Rfar/OptionValue["SamplePoints"]];
   ip = Interpolation[Transpose[{rs, prof /@ rs}], InterpolationOrder -> 3];
   mdot[rho_] := Mfun'[ip[rho]] ip'[rho];                 (* d/drho M(phibar(rho)) *)
   -Quiet@NIntegrate[rho^(dim - 1) mdot[rho]^2, {rho, 0, Rfar}, MaxRecursion -> 20]];
