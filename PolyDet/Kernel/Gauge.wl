(* ::Package:: *)

(* PolyDet` -- Abelian gauge boson (U(1) Higgs) section.  Loaded by Kernel/PolyDet.wl from inside
   PolyDet`Private` (Get'd by the loader); NOT a standalone package -- it has no
   BeginPackage/Begin of its own and relies on the loader's context.
   Transverse + coupled diagonal-Goldstone transfer chain, ghost collapse,
   PolyDet*Gauge*. Depends on the scalar section (PolyDetLnR/Ip/XMax). *)


(* ============================================================================
   Abelian gauge boson (U(1) Higgs) polygonal determinant, D=4.
   Transverse T modes = scalar GY of the gauge mass m_A^2; coupled diagonal-Goldstone
   D^pm,a block = 3-channel (2 at j=0) transfer chain of consecutive orders {j,j+2,j+1},
   gradient jump g (phi_{s+1}-phi_s) at each matching.  Ghost collapse +2j(j+2) - 2(j+1)^2
   = -2 leaves Sigma_gauge = -2 Sum lnR^(T) + Sum (j+1)^2 lndet R^(Da).
   ============================================================================ *)

gaugePrec[j_, nseg_] := Ceiling[Max[70, 6.5 (j + 2) + 0.9 nseg + 25]];

PolyDetXMaxGauge[mA2_, radii_] := PolyDetXMax[mA2, radii];

(* transverse: scalar GY with the gauge mass, orbital l=j -> order nu=j+1 *)
PolyDetLnRGaugeT[j_, mA2_, radii_, prec_: 60] := PolyDetLnR[j + 1, mA2, radii, prec];

(* per-channel boundary matrix (A,B) -> (psi, psi') at radius R, mass M, order nu.
   D-general weight psi = (A I_nu + B K_nu)/rho^w, w = dim/2 - 1 (default D=4, w=1):
   psi' = M I_{nu-1}/R^w - (nu+w)/R^{w+1} I_nu. *)
gaugeEbnd[nu_, R_, M_, w_: 1] := {{BesselI[nu, M R]/R^w, BesselK[nu, M R]/R^w},
   {M BesselI[nu - 1, M R]/R^w - (nu + w)/R^(w + 1) BesselI[nu, M R],
    -M BesselK[nu - 1, M R]/R^w - (nu + w)/R^(w + 1) BesselK[nu, M R]}};
gaugeEbndMat[nus_, R_, Ms_, w_: 1] := Module[{e = MapThread[gaugeEbnd[#1, R, #2, w] &, {nus, Ms}]},
   ArrayFlatten[{{DiagonalMatrix[e[[All, 1, 1]]], DiagonalMatrix[e[[All, 1, 2]]]},
                 {DiagonalMatrix[e[[All, 2, 1]]], DiagonalMatrix[e[[All, 2, 2]]]}}]];
(* ANALYTIC per-component 2x2 inverse of gaugeEbndMat.  Each gaugeEbnd has det = -1/R^{2w+1}
   exactly (= -1/R^{D-1}; the Bessel Wronskian I_{nu-1}K_nu + I_nu K_{nu-1} = 1/(M R)), so
   [[a,b],[c,d]]^-1 = -R^{2w+1} [[d,-b],[-c,a]].
   Replaces the numerical Inverse[], whose conditioning blows up at high nu (entries span e^{+-2 nu eta}
   as I_nu underflows / K_nu overflows) and corrupts the summand past nu~40 even at high precision;
   the analytic form keeps the diagonal block accurate to nu ~ 200+ for the optimal-truncation demo. *)
gaugeEbndInv[nus_, R_, Ms_, w_: 1] := Module[{e = MapThread[gaugeEbnd[#1, R, #2, w] &, {nus, Ms}], iv},
   iv = (-R^(2 w + 1)) {{#[[2, 2]], -#[[1, 2]]}, {-#[[2, 1]], #[[1, 1]]}} & /@ e;
   ArrayFlatten[{{DiagonalMatrix[iv[[All, 1, 1]]], DiagonalMatrix[iv[[All, 1, 2]]]},
                 {DiagonalMatrix[iv[[All, 2, 1]]], DiagonalMatrix[iv[[All, 2, 2]]]}}]];
(* MASSLESS-FV boundary (symmetric false vacuum, m_A(FV) = 0): the power-law limit of
   gaugeEbnd, basis psi = (A rho^nu + B rho^-nu)/rho^w, det = -2 nu R^(-2w-1), with its
   analytic inverse.  Used only at the LAST matching, for the massless gauge channels. *)
gaugeEbndInv0[nu_, R_, w_: 1] := {{(nu + w) R^(w - nu)/(2 nu),  R^(w - nu + 1)/(2 nu)},
                                  {(nu - w) R^(nu + w)/(2 nu), -R^(nu + w + 1)/(2 nu)}};
(* mixed-channel inverse boundary for a massless-FV final matching: power-law inverse for
   the gauge ("A") channels, the analytic Bessel inverse for the massive Goldstone ("a") *)
gaugeEbndInvML[nus_, R_, Ms1_, labs_, w_: 1] := Module[{iv},
   iv = Table[If[labs[[a]] === "A", gaugeEbndInv0[nus[[a]], R, w],
        Module[{e = gaugeEbnd[nus[[a]], R, Ms1[[a]], w]},
          (-R^(2 w + 1)) {{e[[2, 2]], -e[[1, 2]]}, {-e[[2, 1]], e[[1, 1]]}}]], {a, Length[nus]}];
   ArrayFlatten[{{DiagonalMatrix[iv[[All, 1, 1]]], DiagonalMatrix[iv[[All, 1, 2]]]},
                 {DiagonalMatrix[iv[[All, 2, 1]]], DiagonalMatrix[iv[[All, 2, 2]]]}}]];

(* matching jump matrix [[1,0],[d C,1]] (psi continuous, psi' jumps by the mdot_A delta) *)
gaugeJmat[d_, Cm_] := Module[{n = Length[Cm]},
   ArrayFlatten[{{IdentityMatrix[n], 0 IdentityMatrix[n]}, {d Cm, IdentityMatrix[n]}}]];
(* channel data: orders, mass labels (A=gauge, a=Goldstone), coupling matrix C (coeff of g*dphi).
   D-general: label l = j, middle order nu = j + dim/2 - 1, channels
   {nu-1, nu+1, nu}; Clebsches c1 = 2 Sqrt[j/(2j+dim-2)], c2 = 2 Sqrt[(j+dim-2)/(2j+dim-2)]
   (c1^2 + c2^2 = 4 in any D; dim=4 -> 2j/(j+1), 2(j+2)/(j+1)).  j=0: the 2-channel {D^+ at nu+1, a at nu}. *)
gaugeChans[j_, dim_: 4] := If[j == 0,
   {{dim/2, dim/2 - 1}, {"A", "a"}, {{0, 2}, {2, 0}}},
   {{j + dim/2 - 2, j + dim/2, j + dim/2 - 1}, {"A", "A", "a"},
    Module[{c1 = 2 Sqrt[j/(2 j + dim - 2)], c2 = 2 Sqrt[(j + dim - 2)/(2 j + dim - 2)]},
     {{0, 0, -c1}, {0, 0, c2}, {-c1, c2, 0}}]}];

(* coupled diagonal block ln det R^(Da)_j.  dmA = per-matching signed jump g (phi_{s+1}-phi_s)
   (length nsh-1); coupOn flag drops it (decoupled validation limit -> sum of 3 scalar ratios). *)
gaugeLnRDaCore[j_, mA2_, ma2_, dmA_, radii_, prec_, coupOn_, dim_: 4] := Module[
   {nus, lab, Cm, n, w = dim/2 - 1, nsh = Length[radii] + 1, mlfv, mA2p, ma2p, rp, dp, Ms, Tl, Tt, A11, norm},
   {nus, lab, Cm} = gaugeChans[j, dim]; n = Length[nus];
   mlfv = TrueQ[mA2[[-1]] == 0];        (* massless gauge boson in the FV (symmetric FV) *)
   mA2p = SetPrecision[mA2, prec]; ma2p = SetPrecision[ma2, prec];
   rp = SetPrecision[radii, prec]; dp = SetPrecision[dmA, prec];
   Ms[s_] := Sqrt[(If[# == "A", mA2p[[s]], ma2p[[s]]] & /@ lab)];
   Tl = Table[If[mlfv && s == nsh - 1,
              gaugeEbndInvML[nus, rp[[s]], Ms[s + 1], lab, w],
              gaugeEbndInv[nus, rp[[s]], Ms[s + 1], w]] .
              gaugeJmat[If[coupOn, dp[[s]], 0], Cm] . gaugeEbndMat[nus, rp[[s]], Ms[s], w], {s, nsh - 1}];
   Tt = Dot @@ Reverse[Tl]; A11 = Tt[[1 ;; n, 1 ;; n]];
   (* massless "A" channels are normalized against the massless free operator:
      LogGamma[nu+1] + nu Log[2/m_1] per channel, as the scalar massless-FV chain *)
   norm = Sum[If[mlfv && lab[[i]] === "A",
      LogGamma[nus[[i]] + 1] + nus[[i]] Log[2/Ms[1][[i]]],
      nus[[i]] Log[Ms[nsh][[i]]/Ms[1][[i]]]], {i, n}];
   Re[norm + Log[Det[A11]]]];
(* ---- general-R_xi chi-reduction ----
   The correct general-xi diagonal block is NOT the naive mass matrix (which puts the longitudinal
   at orders j,j+2 where the ghost can't cancel it) but the chi-REDUCTION: the longitudinal is a
   SCALAR chi at order nu=j+1 -- the SAME order as the would-be Goldstone a and the ghost -- giving
   a clean 2-channel same-nu {chi, a} chain (masses xi mA^2, xi mA^2 + V_aa) with the gradient jump
   [f'_chi]=2 dmA f_a, [f'_a]=2 dmA f_chi, plus the decoupled transverse-D vector P (mass mA^2).
   At xi=1 it reproduces the default 3-channel block; for general xi the ghost (order j+1, mass
   xi mA^2) cancels chi+a-2gh per multipole down to the V_aa residual.  gaugeEbndMat is the scalar
   (Gmat) boundary, so the chi-a chain reuses the same machinery with orders {j+1,j+1} and the
   rotated coupling {{0,2},{2,0}}. *)
gaugeLnRChiCore[j_, mchi2_, ma2_, dmA_, radii_, prec_, gradOn_] := Module[
   {nus = {j + 1, j + 1}, Cm = {{0, 2}, {2, 0}}, nsh = Length[radii] + 1, mlfv,
    mchi2p, ma2p, rp, dp, Ms, Tl, Tt, A11, norm},
   mlfv = TrueQ[mchi2[[-1]] == 0];      (* massless FV: the chi channel (mass xi m_A^2)
      turns long range for any xi; the a channel stays massive through V_aa(0) *)
   mchi2p = SetPrecision[mchi2, prec]; ma2p = SetPrecision[ma2, prec];
   rp = SetPrecision[radii, prec]; dp = SetPrecision[dmA, prec];
   Ms[s_] := Sqrt[{mchi2p[[s]], ma2p[[s]]}];
   Tl = Table[If[mlfv && s == nsh - 1,
              gaugeEbndInvML[nus, rp[[s]], Ms[s + 1], {"A", "a"}, 1],
              gaugeEbndInv[nus, rp[[s]], Ms[s + 1]]] .
              gaugeJmat[If[gradOn, dp[[s]], 0], Cm] . gaugeEbndMat[nus, rp[[s]], Ms[s]], {s, nsh - 1}];
   Tt = Dot @@ Reverse[Tl]; A11 = Tt[[1 ;; 2, 1 ;; 2]];
   norm = If[mlfv,
      LogGamma[nus[[1]] + 1] + nus[[1]] Log[2/Ms[1][[1]]] + nus[[2]] Log[Ms[nsh][[2]]/Ms[1][[2]]],
      Sum[nus[[i]] Log[Ms[nsh][[i]]/Ms[1][[i]]], {i, 2}]];
   Re[norm + Log[Det[A11]]]];
(* ---- FrozenWall gauge diagonal block ----
   Case-b thin wall: the n channels (orders {nu-1, nu+1, nu}, masses {mA, mA, ma}) DECOUPLE
   within shells -- the dmA coupling lives ONLY in the matching jumps -- so in the
   u_i = rho^((D-1)/2) psi_i variables the wall band is n independent frozen scalar boxes
   (kappa_i^2 = m_i^2 + (nu_i^2 - 1/4) c^2, c = Log[b/a]/h the log-mean, exact eikonal)
   chained with the EXACT jump U' += d_s Cm.U at each matching (u is continuous and the
   (D-1)/(2 rho) weight shift is continuous, so the jump takes the same form in u as in psi).
   Exact Bessels only at the plateau edge (n regular log-derivatives) and the FV
   decomposition (n Wronskian rows, W[uI,uK] = -1): 4n per multipole, N-independent,
   machine precision in the band.  Tachyonic wall ma^2 < 0 just makes kappa imaginary
   (cosh -> cos, exact per box).  Returns ln det R^(Da)_j, = gaugeLnRDaCore up to the
   frozen band error ~1/x_max^3. *)
gaugeLnRDaFrozen[a___] := Quiet[gaugeLnRDaFrozenC[a], {General::munfl}];
gaugeLnRDaFrozenC[j_, mA2_, ma2_, dmA_, radii_, nsub_, dim_: 4, precB_: 40] := Module[
   {nus, lab, Cm, n, R0, RN, m1, mN, x0, xN, sig0, lnN0, sigK, lnKN, mw, rr, dd,
    U, Up, ch, sh, kap, a, b, h, c, mag, lnscale = 0., nsh},
   {nus, lab, Cm} = gaugeChans[j, dim];  n = Length[nus];
   rr = N[radii];  dd = N[dmA];  nsh = Length[radii];
   R0 = SetPrecision[radii[[1]], precB];  RN = SetPrecision[radii[[-1]], precB];
   m1 = Sqrt[SetPrecision[(If[# == "A", mA2[[1]], ma2[[1]]] & /@ lab), precB]];
   mN = Sqrt[SetPrecision[(If[# == "A", mA2[[-1]], ma2[[-1]]] & /@ lab), precB]];
   x0 = m1 R0;  xN = mN RN;
   sig0 = N[1/(2 R0) + m1 BesselI[nus - 1, x0]/BesselI[nus, x0] - nus/R0];
   lnN0 = N[Log[Sqrt[R0]] + Log[BesselI[nus, x0]]];
   sigK = N[1/(2 RN) - mN BesselK[nus - 1, xN]/BesselK[nus, xN] - nus/RN];
   lnKN = N[Log[Sqrt[RN]] + Log[BesselK[nus, xN]]];
   mw = Table[N[(If[# == "A", mA2[[s + 1]], ma2[[s + 1]]] & /@ lab)], {s, nsh - 1}];  (* wall-shell masses^2 *)
   U = N[IdentityMatrix[n]];  Up = DiagonalMatrix[sig0];
   Up = Up + dd[[1]] Cm . U;                                (* jump at the plateau edge R_0 *)
   Do[
     Do[
       a = rr[[s]] + (rr[[s + 1]] - rr[[s]]) (q - 1)/nsub;
       b = rr[[s]] + (rr[[s + 1]] - rr[[s]]) q/nsub;  h = b - a;  c = Log[b/a]/h;
       kap = Sqrt[mw[[s]] + (nus^2 - 1/4) c^2 + 0. I];      (* n-vector; complex if tachyonic *)
       ch = Cosh[kap h];  sh = Sinh[kap h]/kap;
       {U, Up} = {ch U + sh Up, (kap^2 sh) U + ch Up},      (* per-channel = per-row scaling *)
       {q, nsub}];
     Up = Up + dd[[s + 1]] Cm . U;                          (* exact matching jump at rr[[s+1]] *)
     mag = Max[Abs[U], Abs[Up]];
     If[mag > 10.^15, U /= mag; Up /= mag; lnscale += Log[mag]],
     {s, nsh - 1}];
   Re[Sum[nus[[i]] Log[mN[[i]]/m1[[i]]], {i, n}] + Total[lnN0] + Total[lnKN] +
      n lnscale + Log[Det[Up - sigK U]]]];

(* ---- DIRECT general-xi diagonal block (D=4, massive FV) -------------------------
   The coupled (D^-, D^+, a) block built directly, NOT via the chi-decomposition (whose
   same-order value coupling is O(1)-wrong per multipole -- the a channel couples to the
   longitudinal DERIVATIVE, not its value).  Within a shell the 6-dim (4-dim at j=0)
   solution space is closed form: transverse pair T(Z_{j+1}(m_A rho)), longitudinal pair
   G(Z_{j+1}(Sqrt[xi] m_A rho)) [the Sqrt[xi] mass shift], a-channel Z_{j+1}(m_a rho);
   the xi-deformation makes matching ASYMMETRIC, [f'_pair] = xi 2 d_s Lhat f_a,
   [f'_a] = 2 d_s Lhat.f_pair, Lhat = (-c_-, c_+).  Validated against verify_gauge.lndetR_Da
   (xi=1, 1e-40), independent ODE integration (xi = 1/2, 2), and the algebraic gauge-mode
   identity (python/verify_gauge_xi_direct.py, 7/7).  cm, cp = Sqrt[j/(2j+2)],
   Sqrt[(j+2)/(2j+2)] (the "sum = 1" Clebsches; = gaugeChans' c1, c2 halved). *)
gaugeEbndXiJ[j_, xi_, m2A_, m2a_, R0_, prec_] := Module[
   {cm, cp, m, mx, ma, R, zm, zx, za, E},
   cm = Sqrt[j/(2 j + 2)]; cp = Sqrt[(j + 2)/(2 j + 2)];
   m = Sqrt[SetPrecision[m2A, prec]]; mx = Sqrt[SetPrecision[xi, prec]] m;
   ma = Sqrt[SetPrecision[m2a, prec]]; R = SetPrecision[R0, prec];
   zm = m R; zx = mx R; za = ma R;
   E = ConstantArray[0 R, {6, 6}];   (* rows fm,fp,fa,dm,dp,da; cols T+ T- L+ L- a+ a- *)
   E[[1, 1]] = cp m BesselI[j, zm];  E[[1, 2]] = -cp m BesselK[j, zm];
   E[[1, 3]] = -cm mx BesselI[j, zx]; E[[1, 4]] = cm mx BesselK[j, zx];
   E[[2, 1]] = cm m BesselI[j + 2, zm]; E[[2, 2]] = -cm m BesselK[j + 2, zm];
   E[[2, 3]] = cp mx BesselI[j + 2, zx]; E[[2, 4]] = -cp mx BesselK[j + 2, zx];
   E[[3, 5]] = BesselI[j + 1, za];   E[[3, 6]] = BesselK[j + 1, za];
   E[[4, 1]] = cp m^2 BesselI[j + 1, zm] + cp m j BesselI[j, zm]/R;
   E[[4, 2]] = cp m^2 BesselK[j + 1, zm] - cp m j BesselK[j, zm]/R;
   E[[4, 3]] = -cm mx^2 BesselI[j + 1, zx] - cm mx j BesselI[j, zx]/R;
   E[[4, 4]] = -cm mx^2 BesselK[j + 1, zx] + cm mx j BesselK[j, zx]/R;
   E[[5, 1]] = cm m^2 BesselI[j + 1, zm] - cm m (j + 2) BesselI[j + 2, zm]/R;
   E[[5, 2]] = cm m^2 BesselK[j + 1, zm] + cm m (j + 2) BesselK[j + 2, zm]/R;
   E[[5, 3]] = cp mx^2 BesselI[j + 1, zx] - cp mx (j + 2) BesselI[j + 2, zx]/R;
   E[[5, 4]] = cp mx^2 BesselK[j + 1, zx] + cp mx (j + 2) BesselK[j + 2, zx]/R;
   E[[6, 5]] = ma (BesselI[j, za] + BesselI[j + 2, za])/2;
   E[[6, 6]] = -ma (BesselK[j, za] + BesselK[j + 2, za])/2;
   E];
gaugeEbndXi0[xi_, m2A_, m2a_, R0_, prec_] := Module[{mx, ma, R, zx, za, E},
   mx = Sqrt[SetPrecision[xi, prec]] Sqrt[SetPrecision[m2A, prec]];
   ma = Sqrt[SetPrecision[m2a, prec]]; R = SetPrecision[R0, prec];
   zx = mx R; za = ma R;
   E = ConstantArray[0 R, {4, 4}];   (* rows fp,fa,dp,da; cols L+ L- a+ a-  (cp=1, cm=0) *)
   E[[1, 1]] = mx BesselI[2, zx];    E[[1, 2]] = -mx BesselK[2, zx];
   E[[2, 3]] = BesselI[1, za];       E[[2, 4]] = BesselK[1, za];
   E[[3, 1]] = mx^2 BesselI[1, zx] - mx 2 BesselI[2, zx]/R;
   E[[3, 2]] = mx^2 BesselK[1, zx] + mx 2 BesselK[2, zx]/R;
   E[[4, 3]] = ma (BesselI[0, za] + BesselI[2, za])/2;
   E[[4, 4]] = -ma (BesselK[0, za] + BesselK[2, za])/2;
   E];
gaugeDetNxi[j_, m_, mx_] := Module[{cm = Sqrt[j/(2 j + 2)], cp = Sqrt[(j + 2)/(2 j + 2)]},
   (m mx/(j! (j + 2)! 2^(2 j + 2))) (cp^2 m^j mx^(j + 2) + cm^2 mx^j m^(j + 2))];

(* ---- massless FALSE-VACUUM boundary of the direct general-xi block (symmetric FV,
   Coulomb phase) ----  At the FV both gauge channels are massless; the pair operator is
   still xi-deformed, but its indicial powers are the integers {j, j+2, -j, -(j+2)}
   (xi-INDEPENDENT) with xi-dependent eigenvectors.  gaugeIndN is the 2x2 indicial matrix
   (M f'' = B for f = v rho^p at rho=1); gaugePairVecML the null vector normalised to
   DOMINANT component 1 (f_- for p = +-j, f_+ for p = +-(j+2)) -- the m -> 0 limit of the
   massive T/L columns, for which the prefactor is UNCHANGED (no extra massless factor,
   prefML = 0; python/verify_gauge_xi_massless.py 4/4). *)
gaugeIndN[p_, xi_, j_] := Module[
   {cm = Sqrt[j/(2 j + 2)], cp = Sqrt[(j + 2)/(2 j + 2)], kap = 1/xi - 1, Mi, Mm, col, vm, vp, dm, dp, s, u, tp, Bm, Bp, mv},
   Mi = {{1 - (1 - xi) cm^2, (1 - xi) cm cp}, {(1 - xi) cm cp, 1 - (1 - xi) cp^2}};
   Mm = Inverse[Mi];
   col[{vm_, vp_}] := Module[{dmm = p vm, dpp = p vp, ss, uu, tpp, bm, bp, mvv},
      ss = -cm (dmm - j vm) + cp (dpp + (j + 2) vp);
      uu = cm j vm + cp (j + 2) vp;  tpp = (cm j dmm + cp (j + 2) dpp) - uu;
      bm = -dmm + j^2 vm + kap cm (tpp + (j + 1) ss);
      bp = -dpp + (j + 2)^2 vp - kap cp (tpp - (j + 1) ss);
      mvv = Mm . {vm, vp};
      {p (p - 1) mvv[[1]] - bm, p (p - 1) mvv[[2]] - bp}];
   Transpose[{col[{1, 0}], col[{0, 1}]}]];
gaugePairVecML[p_, xi_, j_, dom_] := Module[{N2 = gaugeIndN[p, xi, j], v},
   v = If[Abs[N2[[1, 1]]] + Abs[N2[[1, 2]]] >= Abs[N2[[2, 1]]] + Abs[N2[[2, 2]]],
          {-N2[[1, 2]], N2[[1, 1]]}, {-N2[[2, 2]], N2[[2, 1]]}];
   v/v[[dom]]];
(* massless FV boundary matrix (rows fm,fp,fa,dm,dp,da): cols 1,2 = growing gauge p=+j,+(j+2);
   3,4 = decaying p=-j,-(j+2); 5,6 = a (I,K).  j=0: only the L (p=+-2) gauge + a. *)
gaugeEbndXiML[j_, xi_, m2a_, R0_, prec_] := Module[
   {ma = Sqrt[SetPrecision[m2a, prec]], R = SetPrecision[R0, prec], za, pv, E, pw, ii},
   za = ma R;
   pv = {{j, 1}, {j + 2, 2}, {-j, 1}, {-(j + 2), 2}};   (* {power, dominant-component} *)
   E = ConstantArray[0 R, {6, 6}];
   Do[pw = pv[[k, 1]]; With[{v = gaugePairVecML[pw, xi, j, pv[[k, 2]]]},
       E[[1, k]] = v[[1]] R^pw;            E[[2, k]] = v[[2]] R^pw;
       E[[4, k]] = v[[1]] pw R^(pw - 1);   E[[5, k]] = v[[2]] pw R^(pw - 1)], {k, 4}];
   E[[3, 5]] = BesselI[j + 1, za];  E[[3, 6]] = BesselK[j + 1, za];
   E[[6, 5]] = ma (BesselI[j, za] + BesselI[j + 2, za])/2;
   E[[6, 6]] = -ma (BesselK[j, za] + BesselK[j + 2, za])/2;
   E];

gaugeLnRDaXiDirect[j_, mA2_, ma2_, dmA_, radii_, xi_, prec_] := Module[
   {mlfv, nn, regcols, growcols, ma2xi, mAl, nsh = Length[radii] + 1, cm, cp, state, Es, En, vals, d, s,
    A, m1, mx1, mh, mxh, ma1, mahat, pref, last},
   mlfv = TrueQ[mA2[[-1]] == 0];              (* symmetric FV: gauge massless at infinity *)
   cm = Sqrt[j/(2 j + 2)]; cp = Sqrt[(j + 2)/(2 j + 2)];
   ma2xi = Table[xi mA2[[t]] + (ma2[[t]] - mA2[[t]]), {t, nsh}];   (* general-xi a-mass *)
   mAl = Sqrt[SetPrecision[mA2, prec]];
   {nn, regcols} = If[j == 0, {2, {1, 3}}, {3, {1, 3, 5}}];
   state = IdentityMatrix[2 nn][[All, regcols]];
   Do[
     last = (s == nsh - 1);
     Es = If[j == 0, gaugeEbndXi0[xi, mA2[[s]], ma2xi[[s]], radii[[s]], prec],
                     gaugeEbndXiJ[j, xi, mA2[[s]], ma2xi[[s]], radii[[s]], prec]];
     En = Which[
        last && mlfv && j > 0, gaugeEbndXiML[j, xi, ma2xi[[s + 1]], radii[[s]], prec],
        last && mlfv && j == 0, gaugeEbndXiML[0, xi, ma2xi[[s + 1]], radii[[s]], prec][[{2, 3, 5, 6}, {2, 4, 5, 6}]],
        j == 0, gaugeEbndXi0[xi, mA2[[s + 1]], ma2xi[[s + 1]], radii[[s]], prec],
        True, gaugeEbndXiJ[j, xi, mA2[[s + 1]], ma2xi[[s + 1]], radii[[s]], prec]];
     vals = Es . state;
     d = SetPrecision[dmA[[s]], prec];
     If[j == 0,
        vals[[3]] += xi 2 d vals[[2]]; vals[[4]] += 2 d vals[[1]],
        vals[[4]] += xi 2 d (-cm) vals[[3]];
        vals[[5]] += xi 2 d cp vals[[3]];
        vals[[6]] += 2 d (-cm vals[[1]] + cp vals[[2]])];
     state = LinearSolve[En, vals],
     {s, nsh - 1}];
   m1 = mAl[[1]]; mh = mAl[[-1]]; mx1 = Sqrt[xi] m1; mxh = Sqrt[xi] mh;
   ma1 = Sqrt[ma2xi[[1]]]; mahat = Sqrt[ma2xi[[-1]]];
   If[mlfv,
      (* growing FV modes: gauge power p=+j,+(j+2) (or +2 at j=0) + a growing (I);
         prefactor drops the FV gauge norm (prefML=0), keeps origin detN + a-channel *)
      growcols = If[j == 0, {1, 3}, {1, 2, 5}];
      A = state[[growcols, All]];
      pref = If[j == 0, Log[mahat/ma1] - Log[mx1^3/8],
         (j + 1) Log[mahat/ma1] - Log[gaugeDetNxi[j, m1, mx1]]],
      A = state[[regcols, All]];
      pref = If[j == 0, 3 Log[mxh/mx1] + Log[mahat/ma1],
         (j + 1) Log[mahat/ma1] + Log[gaugeDetNxi[j, mh, mxh]] - Log[gaugeDetNxi[j, m1, mx1]]]];
   Re[pref + Log[Det[A]]]];

(* legacy chi-reduced general-xi block (kept for reference / the xi=1 P-difference; NOT
   used for xi != 1 physics -- superseded by gaugeLnRDaXiDirect, see PolyDet_v1_check.md
   2026-07-05b/c).  lndetR_Da(xi) = lnR[P] + lndetR_chi(xi). *)
gaugeLnRDaXi[j_, mA2_, ma2_, dmA_, radii_, xi_, prec_] := Module[{Vaa = ma2 - mA2},
   (gaugeLnRDaCore[j, mA2, ma2, dmA, radii, prec, True] -
      gaugeLnRChiCore[j, mA2, ma2, dmA, radii, prec, True]) +
   gaugeLnRChiCore[j, xi mA2, xi mA2 + Vaa, dmA, radii, prec, True]];

PolyDetLnRGaugeDa::dimxi = "General Xi is implemented for Dimension 4 only (got Dimension `1`, Xi `2`).";
PolyDetLnRGaugeDa::mlxi = "A massless false vacuum (mA2[[-1]] == 0) is implemented for Xi == 1 only (got Xi `1`).";
Options[PolyDetLnRGaugeDa] = {"Xi" -> 1, WorkingPrecision -> Automatic, "Dimension" -> 4};
PolyDetLnRGaugeDa[j_, mA2_, ma2_, dmA_, radii_, OptionsPattern[]] := Module[
   {xi = OptionValue["Xi"], dim = OptionValue["Dimension"], prec},
   prec = Replace[OptionValue[WorkingPrecision], Automatic :> gaugePrec[j, Length[radii] + 1]];
   Which[
     xi == 1, gaugeLnRDaCore[j, mA2, ma2, dmA, radii, prec, True, dim],
     (* gaugeLnRDaXiDirect handles both a massive FV (detN prefactor) and a massless FV
        (mA2[[-1]] == 0 -> the power-law pair boundary gaugeEbndXiML, prefML = 0) *)
     dim == 4, gaugeLnRDaXiDirect[j, mA2, ma2, dmA, radii, xi, prec],
     True, Message[PolyDetLnRGaugeDa::dimxi, dim, xi]; $Failed]];

(* effective UV (Laurent) trace: 2 I_p^{mA2} + I_p^{ma2} (2 transverse + Goldstone scalars) *)
PolyDetIpGauge[mA2_, ma2_, radii_, pmax_] :=
   2 PolyDetIp[mA2, radii, pmax] + PolyDetIp[ma2, radii, pmax];

(* per-multipole physical summand (ghost collapse) and its scalar (potential) Laurent subtraction *)
gaugeSummand[j_, mA2_, ma2_, dmA_, radii_, prec_] :=
   -2 PolyDetLnR[j + 1, mA2, radii, prec] +
   (j + 1)^2 gaugeLnRDaCore[j, mA2, ma2, dmA, radii, prec, True];
gaugeSubScalar[IpA_, Ipa_, j_, pmax_] := Module[{f},
   f[Ip_, nu_] := Sum[Ip[[p]]/nu^p, {p, pmax}];
   -2 f[IpA, j + 1] + (j + 1)^2 (If[j >= 1, f[IpA, j], 0] + f[IpA, j + 2] + f[Ipa, j + 1])];

(* CLOSED-FORM leading gradient summand (j+1)^2 C_j  (the gauge gradient under
   analytic control -- NO fit, NO Bessels).  The D^pm-a mixing is a SECOND-ORDER (matching-pair)
   effect: C_j = ln det R^(Da)_on - off = -sum_{s,t} dmA_s dmA_t [c1^2 G_j G_{j+1} + c2^2 G_{j+2} G_{j+1}]
   with the massless GY Green's fns G_mu(r,r') = (1/2mu) r<^(mu+1/2) r>^(1/2-mu).  The exact Clebsch
   (c1^2 = 2j/(j+1), c2^2 = 2(j+2)/(j+1)) collapses c1^2/(4j(j+1)) = c2^2/(4(j+2)(j+1)) = 1/(2(j+1)^2),
   leaving the per-pair kernel (r< r>)(1+w^2) w^(2j+1), w = r</r> = min/max radius.  Limits:
     diagonal s=t (w=1) -> c0 = -sum_s dmA_s^2 R_s^2, the O(1/N) finite-N self-energy (-> 0 continuum);
     off-diagonal coherent sum's 1/nu coefficient -> the continuum moment -g^2 Int rho^3 phibar'^2 = I_g.
   Validated vs the exact gaugeLnRDaCore on/off to ~1e-3 at j~100 (the O(1/nu) residual is the mass
   correction, the trap-free fixed-argument Bessel series I_mu(MR)K_mu(MR)=(1/2mu)[1-(MR)^2/2/(mu^2-1)+..],
   NOT the turning-point-singular Debye form).  dmA = per-matching jump g(phi_{s+1}-phi_s). *)
gaugeGradLead[j_, dmA_, radii_] := Module[{nsh = Length[radii] + 1, R = radii, d = dmA},
   -1/2 Sum[With[{rl = Min[R[[s]], R[[t]]], rh = Max[R[[s]], R[[t]]]},
       d[[s]] d[[t]] (rl rh) (1 + (rl/rh)^2) (rl/rh)^(2 j + 1)], {s, nsh - 1}, {t, nsh - 1}]];

(* D-GENERAL closed-form gradient pair sum.  Second order in the jumps, with the massless
   Green's functions G_mu = r_<^{l_mu} r_>^{-(l_mu+D-2)}/(2 mu) and the D-general
   Clebsches c_-^2 = 2l/nu, c_+^2 = 2(l+D-2)/nu:
     C_l = -(1/2 nu^2) sum_{s,t} d_s d_t (R_< R_>)
           [ (nu-D/2+1)/(nu-1) w^(2 nu - 1) + (nu+D/2-1)/(nu+1) w^(2 nu + 1) ],
   w = R_</R_>, nu = l+D/2-1.  Returns C_l (UNweighted).  dim=4: both weights -> 1 and
   nu^2 C_l = gaugeGradLead (exact reduction, unit test "gauge-gradpair-D4");
   w=1 diagonal: the -(1/2nu^2)[2+(4-D)/(nu^2-1)] sum d^2 R^2 self-energy.
   part: "Full" | "Diag" (s=t, the O(1/N) jump artifact) | "Coh" (s!=t, the convergent
   physical gradient content). *)
gaugeGradPair[l_, dmA_, radii_, dim_: 4, part_: "Full"] := Module[
   {nu = l + dim/2 - 1, wm, wp, n = Length[radii], tot},
   wm = (nu - dim/2 + 1)/(nu - 1);  wp = (nu + dim/2 - 1)/(nu + 1);
   tot = Sum[If[(part === "Diag" && s != t) || (part === "Coh" && s == t), 0,
      With[{a = Min[radii[[s]], radii[[t]]], b = Max[radii[[s]], radii[[t]]]},
        dmA[[s]] dmA[[t]] (a b) (wm (a/b)^(2 nu - 1) + wp (a/b)^(2 nu + 1))]],
     {s, n}, {t, n}];
   -tot/(2 nu^2)];

(* ---- D=3 renormalised gauge assembly ----
   xi=1 only.  deg_T(l;3) = deg_H(l;3) = 2l+1, so transverse + ghost collapse
   to MINUS ONE scalar tower at m_A^2:
     Sigma = [l=0: lndet R^Da_0 - 2 lnR_{1/2}(m_A)]  (raw, single finite term)
           + sum_{l>=1} (2l+1)[lndet R^Da_l - lnR_{l+1/2}(m_A)]  (subtracted).
   Per-l subtraction: the exact per-channel scalar Laurents at orders {nu-1, nu+1, nu} minus
   the collapsed tower at nu (nu = l+1/2), PLUS the FULL closed-form gradient pair sum
   (2l+1) gaugeGradPair[l, ..., 3, "Full"].  The pair sum splits: the DIAGONAL (s=t) part,
   (2l+1) C_l^diag = -(2 + 1/(nu^2-1)) (sum_s dmA_s^2 R_s^2)/nu, is the O(1/N) matching-jump
   self-energy (the D=3 analogue of the D=4 c0) -- log-divergent in the l-sum, dropped with no
   add-back (the apparent D=3 log divergence is pure polygon artifact, the continuum is
   scale-free); the COHERENT (s!=t) part is CONVERGENT physical gradient content and is
   restored exactly (a pure-power l-sum, no Bessels).  Convergent Hurwitz add-backs of the
   potential Laurents (l>=1 sums, nu = l+1/2):
     A1(p) = sum (2l+1)/(nu-1)^p = 2 HurwitzZeta(p-1,1/2) + 2 HurwitzZeta(p,1/2)
     A2(p) = sum (2l+1)/(nu+1)^p = 2 HurwitzZeta(p-1,5/2) - 2 HurwitzZeta(p,5/2)
     A3(p) = sum (2l+1)/nu^p     = 2 HurwitzZeta(p-1,3/2)          (p >= 3),
   plus the p=1 convergent remainder sum_{l>=1} 4/(nu^2-1) = 16/3 on I_1^{m_A} (the p=1
   constants are the scaleless power divergence, dropped; I_2 = 0 identically).
   Net trace I_p^eff(D=3) = I_p^{m_A} + I_p^{m_a} (2 gauge + Goldstone - collapsed tower). *)
gaugeSigmaD3[mA2_, ma2_, dmA_, radii_, pmax_, nsumOpt_, wpOpt_] := Module[
   {xs, ns, prec, IpA, Ipa, dmAp, rp, fL, summ0, direct, addback, cohterm, ll, nu},
   xs = PolyDetXMaxGauge[mA2, radii];
   ns = Replace[nsumOpt, Automatic :> Ceiling[2.4 xs + 6]];
   prec = Replace[wpOpt, Automatic :> gaugePrec[ns, Length[radii] + 1]];
   IpA = PolyDetIp[mA2, radii, pmax]; Ipa = PolyDetIp[ma2, radii, pmax];
   dmAp = SetPrecision[dmA, prec];  rp = SetPrecision[radii, prec];
   fL[Ip_, v_] := Sum[Ip[[p]]/v^p, {p, pmax}];
   summ0 = gaugeLnRDaCore[0, mA2, ma2, dmA, radii, prec, True, 3] -
           2 PolyDetLnR[1/2, mA2, radii, prec];
   direct = summ0 + Sum[nu = j + 1/2;
      (2 j + 1) (gaugeLnRDaCore[j, mA2, ma2, dmA, radii, prec, True, 3] -
                 PolyDetLnR[nu, mA2, radii, prec]) -
      (2 j + 1) (fL[IpA, nu - 1] + fL[IpA, nu + 1] + fL[Ipa, nu] - fL[IpA, nu]) -
      (2 j + 1) gaugeGradPair[j, dmAp, rp, 3, "Full"], {j, 1, ns}];
   addback = 16/3 IpA[[1]] + Sum[
      IpA[[p]] (2 HurwitzZeta[p - 1, 1/2] + 2 HurwitzZeta[p, 1/2] +
                2 HurwitzZeta[p - 1, 5/2] - 2 HurwitzZeta[p, 5/2] -
                2 HurwitzZeta[p - 1, 3/2]) +
      Ipa[[p]] 2 HurwitzZeta[p - 1, 3/2], {p, 3, pmax}];
   (* restore the convergent coherent gradient over ALL l >= 1 (pure powers, cheap) *)
   ll = 1;
   While[ll <= 100000,
      cohterm = (2 ll + 1) gaugeGradPair[ll, dmAp, rp, 3, "Coh"];
      addback += cohterm;
      If[Abs[cohterm] < 10^-30, Break[]];
      ll++];
   <|"Sigma" -> direct + addback, "I1eff" -> IpA[[1]] + Ipa[[1]], "I3eff" -> IpA[[3]] + Ipa[[3]],
     "Ig" -> Missing["ScaleFreeD3"], "Nsum" -> ns, "GradientMethod" -> "ClosedPairSum"|>];

(* gradient coefficient I_g: intercept of (j+1)^3 (lndet R^(Da)_on - off) over a high-j window.
   The general-xi value is I_g(xi) = -(3-xi)/2 int rho^3 mdot_A^2 = (3-xi)/2 I_g(xi=1)  -- a simple
   prefactor (the fit gives the xi=1 moment, I_g(1) = -int rho^3 mdot_A^2); this is the gauge
   wave-function gamma_phi(xi) = (3-xi) g^2/(4pi)^2. *)
Options[PolyDetGaugeIg] = {"Window" -> {22, 28}, WorkingPrecision -> Automatic, "Xi" -> 1};
PolyDetGaugeIg[mA2_, ma2_, dmA_, radii_, OptionsPattern[]] := Module[
   {win = OptionValue["Window"], xi = OptionValue["Xi"], prec, cj},
   prec = Replace[OptionValue[WorkingPrecision], Automatic :> gaugePrec[win[[2]], Length[radii] + 1]];
   cj = Table[{1./(j + 1), (j + 1)^3 (gaugeLnRDaCore[j, mA2, ma2, dmA, radii, prec, True] -
                                      gaugeLnRDaCore[j, mA2, ma2, dmA, radii, prec, False])},
              {j, win[[1]], win[[2]]}];
   (3 - xi)/2 LinearModelFit[N[cj], x, x]["BestFitParameters"][[1]]];   (* (3-xi)/2 x I_g(1) *)

(* ---- non-Abelian factorisation ----
   On the bounce the scalar sits along a fixed direction, so the gauge mass matrix M_A^2 = C phibar^2
   with C constant: its eigenbasis is rho-INDEPENDENT and the determinant FACTORISES into a product of
   the Abelian PolyDets of this section,
     ln det_NA = sum_a mult_a ln det_U(1)(m_a^2 = c_a phibar^2),
   one per broken channel (effective coupling g_a = sqrt(c_a)); unbroken c_a=0 = massless spectators
   (the photon) contribute 0.  The group enters the RATE only via the orbit volume V_G = Vol(G/H)
   (= 2 pi^2 for the SM doublet), supplied ONCE, NOT in Sigma.  E.g. the SM electroweak gauge
   determinant = 2 det_U(1)(m_W) + det_U(1)(m_Z) + det_U(1)(0):
     PolyDetSigmaGaugeNonAbelian[bf, {{g/2, 2}, {Sqrt[g^2+g'^2]/2, 1}, {0, 1}}, "Vaa" -> vaa]. *)
(* NB Options[PolyDetSigmaGaugeNonAbelian] is set BELOW, after Options[PolyDetSigmaGauge] exists
   (it Joins them) -- defining it here would capture an empty option list (load order). *)
PolyDetSigmaGaugeNonAbelian[bf_, couplings_List, opts : OptionsPattern[]] := Module[
   {vaa = OptionValue["Vaa"], chans},
   chans = Table[With[{g = c[[1]], n = c[[2]]},
      {g, n, If[g == 0, 0,
         PolyDetSigmaGauge[Sequence @@ PolyDetInputsGauge[bf, g, vaa],
            Sequence @@ FilterRules[{opts}, Options[PolyDetSigmaGauge]]]["Sigma"]]}],
      {c, couplings}];
   <|"Sigma" -> Total[(#[[2]] #[[3]]) & /@ chans], "Channels" -> chans|>];

(* assembled, renormalised gauge determinant by optimal truncation.  Subtract the potential
   Laurent through pmax + the leading gradient I_g/(j+1); restore the convergent potential tail
   via zeta(-1) I_1 + zeta(p-2) (p>=4).  p=3 (I_3 log) and I_g (gradient log) are the renormalised
   counterterms (running of S_R, dropped here at mu=mu_ref).  Returns <|Sigma,I1eff,I3eff,Ig,Nsum|>. *)
PolyDetSigmaGauge::dim = "Gauge assembly implemented for Dimension 3 and 4 only (got `1`).";
PolyDetSigmaGauge::dimxi = "General Xi is implemented for Dimension 4 only (got Dimension `1`, Xi `2`).";
PolyDetSigmaGauge::frozen = "\"Method\" -> \"FrozenWall\" needs Xi == 1, Dimension == 4 and a case-b-like chain (positive plateau and FV masses); falling back to \"Exact\".";
PolyDetSigmaGauge::mlfv = "A massless false vacuum (mA2[[-1]] == 0) is implemented for Xi == 1, Dimension == 4 only (got Xi `1`, Dimension `2`).";
PolyDetSigmaGauge::mlzr = "\"ZeroRemoval\" (the massless-FV orientation zero mode) is implemented at Xi == 1 only; proceeding without removal.";
Options[PolyDetSigmaGauge] = {"Pmax" -> 6, "Nsum" -> Automatic, WorkingPrecision -> Automatic,
   "Ig" -> Automatic, "Xi" -> 1, "Dimension" -> 4, "Method" -> "Exact", "NSub" -> 8,
   "ZeroRemoval" -> False};
(* deferred from above: now that Options[PolyDetSigmaGauge] exists, the non-Abelian wrapper
   inherits all of them (+ "Vaa") so Ig/Nsum/Pmax/Xi/WP pass through without OptionValue::nodef *)
Options[PolyDetSigmaGaugeNonAbelian] = Join[{"Vaa" -> (0 &)}, Options[PolyDetSigmaGauge]];
PolyDetSigmaGauge[mA2_, ma2_, dmA_, radii_, OptionsPattern[]] := Module[
   {pmax = OptionValue["Pmax"], xi = OptionValue["Xi"], igOpt = OptionValue["Ig"],
    dim = OptionValue["Dimension"], useClosed, meth, nsubF, padP, precLow, lnRT, lnDaF,
    xs, ns, prec, IpA, Ipa, Ieff, Ig, igOut, method, gradSub,
    Vaa, Ipchi, Ipac, dIeff, fL, sub, summ, conv, add},
   If[!MemberQ[{3, 4}, dim], Message[PolyDetSigmaGauge::dim, dim]; Return[$Failed]];
   If[dim == 3 && xi != 1, Message[PolyDetSigmaGauge::dimxi, dim, xi]; Return[$Failed]];
   (* massless gauge boson in the FV (symmetric FV, Coulomb phase): supported at D = 4 for
      ANY xi -- the diagonal block uses the power-law FV boundary (gaugeEbndXiML), the
      transverse/ghost towers the scalar massless PolyDetLnR, the potential trace the
      x^0 -> 1 guard, and the gradient the closed-form pair sum.  The multipole cutoff must
      see the LARGEST argument, which sits in the massive Goldstone tower, so xs takes both
      mass lists. *)
   If[TrueQ[mA2[[-1]] == 0] && dim != 4,
      Message[PolyDetSigmaGauge::mlfv, xi, dim]; Return[$Failed]];
   If[dim == 3, Return[gaugeSigmaD3[mA2, ma2, dmA, radii, pmax,
      OptionValue["Nsum"], OptionValue[WorkingPrecision]]]];
   xs = If[TrueQ[mA2[[-1]] == 0],
      Max[PolyDetXMaxGauge[mA2, radii], PolyDetXMax[ma2, radii]],
      PolyDetXMaxGauge[mA2, radii]];
   (* "Method": "Exact" (default) | "FrozenWall" | Automatic -- the thin-wall fast band
      (gaugeLnRDaFrozen for the diagonal block + the scalar lnRallFrozen for the
      transverse tower), xi = 1 / D = 4 only; l = 0, 1 stay on the exact chain. *)
   meth = Replace[OptionValue["Method"], Automatic :>
      If[xi == 1 && dim == 4 && xs >= 30 && Length[radii] >= 3 &&
         N[mA2[[1]]] > 0 && N[ma2[[1]]] > 0 && N[mA2[[-1]]] > 0 && N[ma2[[-1]]] > 0 &&
         N[radii[[1]]/radii[[-1]]] >= 1/2, "FrozenWall", "Exact"]];
   If[meth === "FrozenWall" && ! (xi == 1 && dim == 4 &&
         N[mA2[[1]]] > 0 && N[ma2[[1]]] > 0 && N[mA2[[-1]]] > 0 && N[ma2[[-1]]] > 0),
      Message[PolyDetSigmaGauge::frozen]; meth = "Exact"];
   (* unphysical turning point x_max(xi) = sqrt(xi) x_max(1) -> the cost scales as sqrt(xi) *)
   ns = Replace[OptionValue["Nsum"], Automatic :> Ceiling[2.4 Max[xs, Sqrt[xi] xs] + 6]];
   prec = Replace[OptionValue[WorkingPrecision], Automatic :> gaugePrec[ns, Length[radii] + 1]];
   IpA = PolyDetIp[mA2, radii, pmax]; Ipa = PolyDetIp[ma2, radii, pmax]; Ieff = 2 IpA + Ipa;
   (* gradient subtraction: CLOSED-FORM polygon pair sum (j+1)^2 C_j = gaugeGradLead (DEFAULT, fit-free,
      ns-CONVERGENT) for the physical xi=1 gauge; else the supplied / continuum-moment single term
      I_g/(j+1) (e.g. the rate, which passes the heat-kernel moment).  Closed: the gradient log + the
      finite-N self-energy c0 is the counterterm, dropped at mu_ref (the OLD I_g/(j+1) left a c0(ns)
      drift -> non-convergent); the running I_3^eff + I_g is presented separately, I_g the moment
      -g^2 Int rho^3 phibar'^2 (needs the bounce profile, hence Missing here). *)
   (* general-xi closed-form gradient: (1+xi)/2 x the xi=1 pair sum gaugeGradLead.  The
      DIAGONAL-block gradient's self-energy constant scales EXACTLY as c_0(xi)=(1+xi)/2 c_0(1)
      (verified 1e-16, python/verify_gauge_xi_assembly.py), and the coherent part matches
      gaugeGradLead to the same order as at xi=1 -> Sigma is ns-CONVERGENT for every xi
      (the old I_g/(j+1) moment left the c_0(xi) drift).  At xi=1 this is exactly gaugeGradLead.
      NB this closed form is the DIAGONAL gradient (for ns-convergence); the PHYSICAL running
      I_g = gamma_phi = (3-xi)/2 (-g^2 Int rho^3 phibar'^2) also has a transverse piece and is
      the profile-based moment reported separately (Missing here). *)
   useClosed = (igOpt === Automatic);
   If[useClosed,
      gradSub[j_] := (1 + xi)/2 gaugeGradLead[j, dmA, radii]; igOut = Missing["BounceProfile"]; method = "Closed",
      Ig = Replace[igOpt, Automatic :> PolyDetGaugeIg[mA2, ma2, dmA, radii, "Xi" -> xi]];
      gradSub[j_] := Ig/(j + 1); igOut = Ig; method = "Moment"];
   (* massless-FV ORIENTATION zero removal ("ZeroRemoval" -> True): with a symmetric
      false vacuum the coupled j=0 block has an EXACT zero mode (the rigid rotation
      dressed by the ghost compensator) for every xi; the raw det R^(Da)_0 is lifted
      only by the polygon's tail truncation.  Replace it by the reduced determinant
      R_0' = mhat_a^2 d/d(mu^2) det R^(Da)_0 (interior masses of both channels offset,
      gradient jumps and FV masses fixed; McKane-Tarlie: R_0' carries the DRESSED mode
      norm, and the rate then gains the orbit measure (Q~/2pi)^(1/2) Vol per broken
      generator).  Implemented at xi = 1. *)
   zrem = TrueQ[OptionValue["ZeroRemoval"]] && TrueQ[mA2[[-1]] == 0];
   If[zrem && xi != 1, Message[PolyDetSigmaGauge::mlzr]; zrem = False];
   If[zrem,
      r0v = Module[{del = SetPrecision[10^-5, prec], sh, lp, lm},
         sh[m2_, d_] := Join[m2[[;; -2]] + d, {m2[[-1]]}];
         lp = gaugeLnRDaCore[0, sh[mA2, del], sh[ma2, del], dmA, radii, prec, True, dim];
         lm = gaugeLnRDaCore[0, sh[mA2, -del], sh[ma2, -del], dmA, radii, prec, True, dim];
         ma2[[-1]] (E^lp - E^lm)/(2 del)]];
   If[xi == 1,                                              (* xi=1 : Feynman, the physical/optimal gauge *)
      sub[j_] := gaugeSubScalar[IpA, Ipa, j, pmax];
      If[meth === "FrozenWall",
         (* frozen band: transverse tower = the scalar frozen band at mass mA^2 (order
            j+1 = scalar l-index j), diagonal block = gaugeLnRDaFrozen; j = 0, 1 exact.
            Band values are machine -> SetPrecision-pad for the subtraction assembly
            (pmax = 6: the ~x_max^(pmax-3) cancellation vs the zeta add-back). *)
         nsubF = OptionValue["NSub"];
         padP = Max[40, (pmax - 3) Log10[Max[xs, 2]] + 25];
         precLow = gaugePrec[2, Length[radii] + 1];
         lnRT = SetPrecision[Join[
            Table[PolyDetLnR[l + 1, mA2, radii, precLow], {l, 0, Min[1, ns]}],
            If[ns >= 2, lnRallFrozen[mA2, radii, 4, ns + 1, 2, nsubF], {}]], padP];
         lnDaF[j_] := SetPrecision[If[j <= 1,
            gaugeLnRDaCore[j, mA2, ma2, dmA, radii, precLow, True],
            gaugeLnRDaFrozen[j, mA2, ma2, dmA, radii, nsubF]], padP];
         summ[j_] := -2 lnRT[[j + 1]] + (j + 1)^2 lnDaF[j],
         summ[j_] := If[zrem && j == 0,
            -2 PolyDetLnR[1, mA2, radii, prec] + Log[Abs[r0v/ma2[[-1]]]],
            gaugeSummand[j, mA2, ma2, dmA, radii, prec]]],
      (* general xi with the DIRECT coupled block (ghost at xi mA^2, diagonal via
         gaugeLnRDaXiDirect).  The summand is written as the xi=1 assembly + (xi-dependent
         Delta): the ghost shift -2(j+1)^2[lnR(xi mA^2) - lnR(mA^2)] plus the diagonal shift
         (j+1)^2[direct(xi) - Core(1)] combine to the exact direct assembly
           2 j(j+2) lnR^T(mA^2) - 2(j+1)^2 lnR^gh(xi mA^2) + (j+1)^2 lndetR^Da,direct(xi),
         so xi=1 reduces exactly.  The Laurent subtraction is the physical effective trace
           Ieff(xi) = 3 I_p^{mA2} - I_p^{xi mA2} + I_p^{xi mA2 + V_aa}
         (the ghost -2 I_p(xi mA2) cancels the unphysical longitudinal+Goldstone xi^2 growth;
         I_1^eff is exactly xi-INDEPENDENT -- the gauge-invariant quadratic divergence -- while
         I_3^eff carries the xi-dependent Goldstone-mass cross-term + the wave function gamma_phi
         via the gradient; see python/verify_gauge_xi_assembly.py).  The gradient gradSub for
         xi != 1 uses the (3-xi)/2 I_g(1) moment (single-term; the closed-form pair sum at general
         xi is not yet ported, so present the running coefficients, not the bare Sigma). *)
      Vaa = ma2 - mA2;
      Ipchi = PolyDetIp[xi mA2, radii, pmax]; Ipac = PolyDetIp[xi mA2 + Vaa, radii, pmax];
      dIeff = IpA - Ipchi + Ipac - Ipa; Ieff = Ieff + dIeff;
      fL[Ip_, nu_] := Sum[Ip[[p]]/nu^p, {p, pmax}];
      summ[j_] := gaugeSummand[j, mA2, ma2, dmA, radii, prec] +
         -2 (j + 1)^2 (PolyDetLnR[j + 1, xi mA2, radii, prec] - PolyDetLnR[j + 1, mA2, radii, prec]) +
         (j + 1)^2 (gaugeLnRDaXiDirect[j, mA2, ma2, dmA, radii, xi, prec] -
                    gaugeLnRDaCore[j, mA2, ma2, dmA, radii, prec, True]);
      sub[j_] := gaugeSubScalar[IpA, Ipa, j, pmax] + (j + 1)^2 fL[dIeff, j + 1]];
   conv = Sum[summ[j] - sub[j] - gradSub[j], {j, 0, ns}];
   add = Zeta[-1] Ieff[[1]] + Sum[Zeta[p - 2] Ieff[[p]], {p, 4, pmax}];
   <|"Sigma" -> conv + add, "I1eff" -> Ieff[[1]], "I3eff" -> Ieff[[3]], "Ig" -> igOut, "Nsum" -> ns, "Method" -> meth,
     "GradientMethod" -> method|>];

(* extract {mA2, ma2, dmA, radii} from a BounceFunction (mirrors PolyDetInputs' shell structure) *)
PolyDetInputsGauge[bf_, g_, Vaa_] := Module[{nodes, rd, caseB, seg, phiShell, rA},
   nodes = Flatten[bf["Path"]]; rd = bf["Radii"]; caseB = Abs[N@rd[[1]]] >= 10^-10;
   seg = (Most[nodes] + Rest[nodes])/2;
   phiShell = If[caseB, Join[{nodes[[1]]}, seg, {nodes[[-1]]}], Join[seg, {nodes[[-1]]}]];
   rA = If[caseB, rd, Rest[rd]];
   {g^2 phiShell^2, g^2 phiShell^2 + (Vaa /@ phiShell), g Differences[phiShell], rA}];

