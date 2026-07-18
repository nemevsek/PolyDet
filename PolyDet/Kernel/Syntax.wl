(* ::Package:: *)

(* PolyDet` -- SyntaxInformation section.  Loaded by Kernel/PolyDet.wl from inside
   PolyDet`Private`.  Gives the front end argument-count syntax coloring / error
   highlighting (red underline on wrong arg count) and feeds the Input Assistant.
   "_." = an optional trailing argument (precision); OptionsPattern[] = options. *)

SyntaxInformation[PolyDetLnR]        = {"ArgumentsPattern" -> {_, _, _, _.}};
SyntaxInformation[PolyDetSigma]      = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SyntaxInformation[PolyDetRen]        = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SyntaxInformation[PolyDetIp]         = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[PolyDetInputs]     = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[PolyDetXMax]      = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[PolyDetR1Prime]    = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

SyntaxInformation[PolyDetLnRMat]     = {"ArgumentsPattern" -> {_, _, _, _.}};
SyntaxInformation[PolyDetIpMat]      = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[PolyDetSigmaMat]   = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SyntaxInformation[PolyDetRenMat]     = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SyntaxInformation[PolyDetR1PrimeMat] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SyntaxInformation[PolyDetXMaxMat]   = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[PolyDetRate]       = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SyntaxInformation[PolyDetDecayRate]  = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

SyntaxInformation[PolyDetLnRPsi]     = {"ArgumentsPattern" -> {_, _, _, _.}};
SyntaxInformation[PolyDetIpPsi]      = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[PolyDetSigmaPsi]   = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SyntaxInformation[PolyDetRenPsi]     = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
SyntaxInformation[PolyDetXMaxPsi]   = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[PolyDetIgPsi]      = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

SyntaxInformation[PolyDetLnRGaugeT]  = {"ArgumentsPattern" -> {_, _, _, _.}};
SyntaxInformation[PolyDetLnRGaugeDa] = {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}};
SyntaxInformation[PolyDetIpGauge]    = {"ArgumentsPattern" -> {_, _, _, _}};
SyntaxInformation[PolyDetGaugeIg]    = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};
SyntaxInformation[PolyDetSigmaGauge] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};
SyntaxInformation[PolyDetXMaxGauge] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[PolyDetInputsGauge]= {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[PolyDetSigmaGaugeNonAbelian] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
