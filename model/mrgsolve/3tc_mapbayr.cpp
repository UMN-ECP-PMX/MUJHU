
$PROB
 Lamivudine (3TC) triphosphate (3TCtp) popPK model 
  
$SET end = 48, ss_cmt = "cent", outvar = "depot,cent,centetp,C2,C4,C7,C5,C8,C6,C9,PAR,MET,DV"    

$CMT @annotated 
depot   : 3tc depot [ADM]
cent    : 3tc plasma central compartment [OBS]
periph  : 3tc plasma peripheral compartment
centv   : 3tc vaginal tissue compartment
cente   : 3tc cervical tissue compartment
centr   : 3tc rectal tissue compartment
centvtp : 3tctp vaginal tissue compartment
centetp : 3tctp cervical tissue compartment [OBS] 
centrtp : 3tctp rectal tissue compartment
trans1  : 3tc gut transit compartment 1 
trans2  : 3tc gut transit compartment 2 
trans3  : 3tc gut transit compartment 3 
trans4  : 3tc gut transit compartment 4 
trans5  : 3tc gut transit compartment 5 
trans6  : 3tc gut transit compartment 6 
trans7  : 3tc gut transit compartment 7 

$PARAM @annotated
  ETA1: 0 : Ka (/hour)
  ETA2: 0 : Vc (L)
  ETA3: 0 : CLttvvtp (L/hour)
  ETA4: 0 : CLvvtp (L/hour)
  ETA5: 0 : CLttve (L/hour)
  ETA6: 0 : CLvetp (L/hour)
  ETA7: 0 : CLttvr (L/hour)
  ETA8: 0 : Kg (/hour)
  
$THETA @annotated
  0.649 	  : Ka (/hour)
  72.3 	    : Vc (L)
  122 	    : Vp (L)
  6.06	    : Q (L/hour)
  18.9	    : CLtt (L/hour)
  0.00131 	: Fv
  0.000069	: Fe
  0.00401		: Fr
  0.325		  : Fvt
  0.0101		: CLttvvtp (L/hour)
  0.09 	    : Vv (L)
  0.09    	: Vvtp (L)
  0.0399 	  : CLVvtp (L/hour)
  1.00  		: Fet
  0.000947  : CLttve (L/hour)
  0.01	    : Ve (L)
  0.01	    : Vetp (L)
  0.0082		: CLVetp (L/hour)
  0.0107	  : Frt
  0.0223		: CLttvr (L/hour)
  0.17	    : Vr (L)
  0.17 	    : Vrtp (L)
  0.647	    : CLvrtp (L/hour)
  0.0724		: Kg (/hour)
  1     		: Kga (/hour)
  1		      : Kgr (/hour)

$OMEGA 
  0.236879973          // IIV Ka
  0.041165935          // IIV Vc
  0.136878659          // IIV CLttvvtp
  0.147730975          // IIV CLvvtp
  0.111841866          // IIV CLttve
  1.31350886           // IIV CLvetp
  0.563550892          // IIV CLttvr
  0.399715947          // IIV Kg

//example: no correlation proportional error between plasma 3tc and cervical tissue 3tctp
$SIGMA 
  0.099856 // proportional error on plasma 3tc
  0.000    // additive error on plasma 3tc
  0.338724 // proportional error on cervical tissue 3tctp
  0.000    // additive error on cervical tissue 3tctp
// reminder: sigma values can be recorded in multiple $SIGMA blocks

$MAIN
  double Fv = THETA(6);					
  double Fe = THETA(7); 					
  double Fr = THETA(8); 				
  
  // Plasma
  double TVKa = THETA(1); 	
  double Ka   = TVKa*exp(ETA1+ETA(1));
  double K1T2 = Ka; 
  
  double TVVc = THETA(2);	
  double Vc   = TVVc*exp(ETA2+ETA(2));
  double TVVp = THETA(3);
  double Vp   = TVVp; 
  double TVQ  = THETA(4);
  double Q    = TVQ;
  double K2T3 = Q/Vc;
  double K3T2 = Q/Vp;
  
  double TVCLtt = THETA(5);	
  double CLtt   = TVCLtt;
  double CLt    = (1-Fv-Fe-Fr)*CLtt;
  double K2T0   = CLt/Vc;
  
  double S2 = Vc;
  
  //vaginal tissue
  double CLv  = Fv*CLtt;
  double K2T4 = CLv/Vc;
  
  double Fvt        = THETA(9);
  double TVCLttvvtp = THETA(10);
  double CLttvvtp   = TVCLttvvtp*exp(ETA3+ETA(3));
  
  double TVVv = THETA(11);
  double Vv   = TVVv;
  double CLvv = (1-Fvt)*CLttvvtp;
  double CLtv = Fvt*CLttvvtp;
  
  double TVVvtp   = THETA(12);
  double Vvtp     = TVVvtp;
  double TVCLvvtp = THETA(13);
  double CLvvtp   = TVCLvvtp*exp(ETA4+ETA(4));
  double K4T0     = CLvv/Vv;
  double K4T7     = CLtv/Vv;
  double K7T0     = CLvvtp/Vvtp;
  
  double S4 = Vv;
  double S7 = Vvtp;
  
  //Cervical tissue
  double CLe  = Fe*CLtt;
  double K2T5 = CLe/Vc;
  
  double Fet      = THETA(14);
  double TVCLttve = THETA(15);
  double CLttve   = TVCLttve*exp(ETA5+ETA(5));
  
  double TVVe = THETA(16);
  double Ve   = TVVe;
  double CLve = (1-Fet)*CLttve;
  double CLte = Fet*CLttve;
  
  double TVVetp   = THETA(17);
  double Vetp     = TVVetp;
  double TVCLvetp = THETA(18);
  double CLvetp   = TVCLvetp*exp(ETA6+ETA(6));
  double K5T0     = CLve/Ve;
  double K5T8     = CLte/Ve;
  double K8T0     = CLvetp/Vetp;
  
  double S5 = Ve;
  double S8 = Vetp;
  
  //Rectal tissue
  double CLr  = Fr*CLtt;
  double K2T6 = CLr/Vc;
  
  double Frt      = THETA(19);
  double TVCLttvr = THETA(20);
  double CLttvr   = TVCLttvr*exp(ETA7+ETA(7));
  
  double TVVr = THETA(21);
  double Vr   = TVVr;
  double CLvr = (1-Frt)*CLttvr;
  double CLtr = Frt*CLttvr;
  
  double TVVrtp   = THETA(22);
  double Vrtp     = TVVrtp;
  double TVCLvrtp = THETA(23);
  double CLvrtp   = TVCLvrtp;
  double K6T0     = CLvr/Vr;
  double K6T9     = CLtr/Vr;
  double K9T0     = CLvrtp/Vrtp;
  
  double S6 = Vr;
  double S9 = Vrtp;
  
  // Transit compartments
  double TVKg  = THETA(24);
  double Kg    = TVKg*exp(ETA8+ETA(8));
  double TVKga = THETA(25);
  double Kga   = TVKga;
  double TVKgr = THETA(26);
  double Kgr   = TVKgr;
  
  double K1T10  = Kg;
  double K10T11 = Kg;
  double K11T12 = Kg;
  double K12T13 = Kg;
  double K13T14 = Kg;
  double K14T15 = Kg;
  double K15T16 = Kg;
  double K16T6  = Kga;
  double K16T0  = Kgr;

$ODE
  dxdt_depot   = -K1T2*depot - K1T10*depot;
  dxdt_cent    =  K1T2*depot - K2T3*cent + K3T2*periph - K2T0*cent - K2T4*cent - K2T5*cent - K2T6*cent;
  dxdt_periph  =  K2T3*cent - K3T2*periph; 
  dxdt_centv   =  K2T4*cent - K4T0*centv - K4T7*centv; 
  dxdt_cente   =  K2T5*cent - K5T0*cente - K5T8*cente;
  dxdt_centr   =  K2T6*cent - K6T0*centr - K6T9*centr + K16T6*trans7;
  dxdt_centvtp =  K4T7*centv - K7T0*centvtp; 
  dxdt_centetp =  K5T8*cente - K8T0*centetp; 
  dxdt_centrtp =  K6T9*centr - K9T0*centrtp; 
  dxdt_trans1  =  K1T10*depot - K10T11*trans1; 
  dxdt_trans2  =  K10T11*trans1 - K11T12*trans2; 
  dxdt_trans3  =  K11T12*trans2 - K12T13*trans3; 
  dxdt_trans4  =  K12T13*trans3 - K13T14*trans4; 
  dxdt_trans5  =  K13T14*trans4 - K14T15*trans5; 
  dxdt_trans6  =  K14T15*trans5 - K15T16*trans6; 
  dxdt_trans7  =  K15T16*trans6 - K16T6*trans7 - K16T0*trans7; 

$TABLE
  double C2 = cent/S2;    // IPRED 3TC in plasma            (fmol/L)
  double C4 = centv/S4;   // IPRED 3TC in vaginal tissue    (fmol/L)
  double C7 = centvtp/S7; // IPRED 3TCtp in vaginal tissue  (fmol/L)
  double C5 = cente/S5;   // IPRED 3TC in cervical tissue   (fmol/L)
  double C8 = centetp/S8; // IPRED 3TCtp in cervical tissue (fmol/L)
  double C6 = centr/S6;   // IPRED 3TC in rectal tissue     (fmol/L)
  double C9 = centrtp/S9; // IPRED 3TCtp in rectal tissue   (fmol/L)
  
  double PAR = C2*(1+EPS(1))+EPS(2); // OBS 3TC in plasma
  double MET = C8*(1+EPS(3))+EPS(4); // OBS 3TCtp in cervical tissue
  
  double DV = PAR ;
  if(self.cmt == 8) DV = MET ; 
  // reminder: use "self.cmt" to internaly refer to a compartment in a mrgsolve model code. 
  
$CAPTURE C2 C4 C7 C5 C8 C6 C9 PAR MET DV
  
  
