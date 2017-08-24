#include "FIT.h"

using namespace std;

bool debug_Fit  = false;
bool debug_Fit_Glo  = false;
bool debug_FitRes  = false;

void FIT::FIT_simple(int nrPnts, double *ax, double *ay, float &m, float &q, float &chi2){
  
  // 1 or 2 SUPERLAYERS SIMPLE FIT
  // liner fit: 	x = a + m*y --> SL	(nrVar, nrPnts)
  if(debug_Fit) printf("...Starting FIT_simple...\n");
  
  //weights
  Double_t er = 1.; // 1./TMath::Power(450.,2);
  Double_t ae[nrPnts];
  for(int i=0;i<nrPnts;i++)  ae[i]=double(er);
  
  const Int_t nrVar  = 2;
  if(debug_Fit) printf("N.Hit to fit %d\n",nrPnts);
  
  
  // Make the vectors 'Use" the data : they are not copied, the vector data
  // pointer is just set appropriately
  TVectorD x; x.Use(nrPnts,ax);
  TVectorD y; y.Use(nrPnts,ay);
  TVectorD e; e.Use(nrPnts,ae);
  
  TMatrixD A(nrPnts,nrVar);
  TMatrixDColumn(A,0) = 1.0;
  TMatrixDColumn(A,1) = x;
  if(debug_Fit) {
    printf("X: ");
    for(int i=0;i<nrPnts;i++)
      printf(" %.1f, ",x(i));
    printf("\nY: ");
    for(int i=0;i<nrPnts;i++)
      printf(" %.1f, ",y(i));
    printf("\n");
  }
  // first bring the weights in place
  TMatrixD Aw = A;
  TVectorD yw = y;
  if(debug_Fit) printf("A.GetNrows() %d\n",A.GetNrows());
  for (Int_t irow = 0; irow < A.GetNrows(); irow++) {
    TMatrixDRow(Aw,irow) *= 1/e(irow);
    yw(irow) /= e(irow);
  }
  
  TDecompSVD svd(Aw);
  Bool_t ok;
  TVectorD *c_svd = NULL;
  c_svd= new TVectorD(svd.Solve(yw,ok));
  
  m=0.;
  q=0.;
  if(ok){
    m=(*c_svd)(1);
    q=(*c_svd)(0); 
  }
  if(debug_Fit) printf("Parametri: m=%.2f,  q=%.1f \n",(*c_svd)(1),(*c_svd)(0));
  
  chi2=0;
  if(ok) 
    for(int i=0;i<nrPnts;i++){
      chi2 += TMath::Power((m*x(i)+q)-y(i),2);
      if(debug_Fit) printf("...chi2 %.3f\n",chi2);
    }
  if(nrPnts>2) 
    chi2=chi2/(nrPnts-2);
  if(debug_Fit) printf("Sum(res2)/ndf: %.3f\n",chi2);
  
  delete c_svd;
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FIT::FIT_t0(double sigmaPhi, double sigmaTimes, int nrPnts, int nrMinPnts, int nrVarv, double *ax, double *ay, int *ac, pair<double,double> &m, pair<double,double> &q, pair<double,double> &t0, double &chi2, int &NPT, bool &okFit){
  
  // 2 SUPERLAYERS t0 FIT
  // xi = fi + ci*v(ti-t0) = fi + ci*v*ti - ci*z0 = x_hit -ci*z0
  // xi = x coordinates along chamber FE, fi = wire coordinate, ci = left/right code
  // v = drift velocity (fixed), ti = drift time, 
  // t0 = time of muon passage respect to the clock 
  // xi = a + m*y --> SL   (nrVarv, nrPnts)
  // redo the fit in case one residual is out of sigmaTimes sigma

  if(debug_Fit) printf("...Starting FIT_t0...\n");
  
  Bool_t reDoFit = true;  
  Int_t numFit = 0;
  okFit=false;
  static const Int_t numMaxFit = 5;
  TVectorD *c_svd = NULL; 
  Int_t totParMax = 0;
  Int_t totPar = 0;
  Int_t nrVar=0;
  Int_t nPF=0;
  bool fitT0=true;
  //reset errors
  Double_t X_err = 0.;
  Double_t X_Phi_err = 0.;
  Double_t Phi_err = 0.;
  Double_t Phi_T0_err = 0.;
  Double_t X_T0_err = 0.;
  Double_t T0_err = 0.;
  Double_t res0[nrPnts];
  Double_t resSum=0.;
  Double_t resSum2=0.;
  //  Double_t meanSum2=0.;
  Double_t meanRes2 = 0.;
  Double_t totMeanRes2 = 0.;
  double vdrift = 0.00547; //ALTEA ??
  //  Double_t sigmaPhi = 400.; //TODO
  //   static const Float_t sigmaTimes = 3.;         //reject hits if |residual| > sigmaTimes * sigma
  
  // weights
  // Double_t er = 1.; // 1./TMath::Power(450.,2);
  // Double_t er = 1/TMath::Power(0.03,2);
  Double_t er = 1/TMath::Power(0.1,2); //ALTEA
  Double_t ae[nrPnts];
  for(int i=0;i<nrPnts;i++)  ae[i]=double(er);
  
  if(debug_Fit) printf("N.Hit to fit %d\n",nrPnts);
  
  while(reDoFit && numFit < numMaxFit)
    {
      // reset number of points to fit
      nPF = 0;
      
      // Count this fit N. of points
      for(int i=0;i<nrPnts;i++) 
	if(ax[i]) 
	  nPF +=1;
      
      totPar = nrVarv;
      
      if(debug_Fit)
	cout << "Total points to fit are " << nPF << ", total parameters are " << totPar << endl;	
      
      totParMax = totPar;
      
      //check nPF > parameter number
      if(nPF <= totPar)
	{
	  if(debug_Fit)
	    cout << "Don't fit : total num.points <= num.parameters ! " << endl;
	  //go to the end of while loop and exit - if there was a previous loop with okFit keep the latest results
      if(nPF<3)
          okFit=false;
	  break;
	}
      
      //check: if nPF>=minimum number of points for each requested SL
      Bool_t breakFlag = false;
      
      if(debug_Fit)
	cout << "Check segment, if points >= nrMinPoints " << endl; 
      if(debug_Fit)
	cout << "nrMinPoints = " << nrMinPnts << endl; 
      
      if(nrVarv)
	if( nPF < TMath::Max(nrVarv,nrMinPnts))
	  {
	    if(debug_Fit)
	      cout << "Segment Phi has: " << nPF << " points < " << TMath::Max(nrVarv,nrMinPnts)<< " nrMinPoints ... reject segment" << endl; 
	    breakFlag = true;
	  }
      
      //if one requested segment doesn't have the minimum points then end fit loop
      if(breakFlag)
	{
	  //okFit=false => reject the event, okFit=true => keep the latest result
	  okFit = false; 
	  break;
	}
      
      //count fit number
      numFit += 1; 
      
      //reset chi2 and residuals
      for(Int_t k=0; k<nrPnts; k++)
	{
	  res0[k] = 0.;
	}
      NPT=0;
      chi2 = 0.;
      resSum = 0.;
      resSum2 = 0.;
      meanRes2 = 0.;
      totMeanRes2 = 0.;
      
      if(debug_Fit)
	cout << "Fit number " << numFit << endl;
      
      nrVar = nrVarv;
      Bool_t acceptSeg = false;
      
      if(debug_Fit)
	cout << "Examing segment for fit... " << endl; 	
      
      // ** perform checks on segment for proper fit
      // 1. check that parameter number is 2 at least...
      if(nrVar<2)
	{
	  if(debug_Fit)
	    cout << "Don't fit segment: min parameter number should be 2!" << endl;
	  continue;
	}
      
      // 2. check: if mu crosses all cells at the same side of the wire do not fit t0!
      Double_t cSum = 0;
      Double_t pSum = 0;
      
      for (int ji=0; ji<nrPnts; ji++)
	{
	  if(ax[ji]!=0){
	    cSum += ac[ji];
	    pSum += 1;
	  }
	}
      if(TMath::Abs(cSum) == pSum)
	{
	  if(debug_Fit)
	    cout << "Don't fit t0 in segment: all cell at same side of wires!" << endl;
	  nrVar = 2;
	  fitT0=false;
	}
      
      //ok, can do the fit
      acceptSeg = true;
      
      if(debug_Fit){
	for(Int_t i=0; i<nrPnts; i++)
	  cout << " ax[" <<  i << "]=" << ax[i] << "  "; 
	cout << endl;
	for(Int_t i=0; i<nrPnts; i++)
	  cout << " ay[" <<  i << "]=" << ay[i] << "  "; 
	cout << endl;
      }
      
      // define matrix
      Int_t matrixRows = nrVar;
      
      TMatrixD A(matrixRows,matrixRows);
      TVectorD y(matrixRows);
      Int_t lr = matrixRows-1; //last row/column index
      
      
      // compute sums
      // this sums go over hits of all segments
      Double_t Scc = 0.;
      Double_t Sx = 0.;
      Double_t Sxy = 0.;
      Double_t Sy = 0.;
      Double_t Syy = 0.;
      Double_t Sc = 0.;
      Double_t Scy = 0.;
      Double_t Scx = 0.;
      Double_t Sw = 0.;
      
      // first row to insert matrix elements for this segment
      Int_t fr = 0;
      
      for (int j=0; j<nrPnts; j++)
	if(ax[j])
	  {
	    Sx  	= Sx + ae[j]*ax[j];
	    Sxy 	= Sxy + ae[j]*ax[j]*ay[j];
	    Sy  	= Sy + ae[j]*ay[j];
	    Syy 	= Syy + ae[j]*ay[j]*ay[j];
	    Sw	        = Sw + ae[j];
	    if(nrVar>=3)
	      {
		Scc 	= Scc + ae[j]*ac[j]*ac[j];
		Scy 	= Scy + ae[j]*ac[j]*ay[j];
		Sc 	= Sc + ae[j]*ac[j];
		Scx 	= Scx + ae[j]*ac[j]*ax[j];
	      }
	  } //end loop
      
      if(nrVar>=2)
	{
	  A(fr,fr)=Syy;	        A(fr,fr+1)=Sy;			
	  A(fr+1,fr)=Sy;	A(fr+1,fr+1)=Sw;
	  
	  if(nrVar==2)
	    {
	      y(fr)=Sxy;
	      y(fr+1)=Sx;
	    }
	}
      
      if(nrVar==3)
	{
	  Int_t lrnew = lr;
	  
	                                                A(fr,lrnew)=Scy;
	                                                A(fr+1,lrnew)=Sc;
	  A(lrnew,fr)=Scy;	A(lrnew,fr+1)=Sc;	A(lrnew,lrnew)=Scc;
	  
	  y(fr)=Sxy;
	  y(fr+1)=Sx;
	  y(lrnew)=Scx;
	}
      
      if(debug_Fit)
	{	
	  cout << "matrice da invertire "<<matrixRows<<" x "<<matrixRows<<":"<< endl; 
	  for(Int_t r=0; r<A.GetNrows(); r++)
	    { 
	      for(Int_t ci=0; ci<A.GetNcols(); ci++)
		cout << " " << A(r,ci);
	      cout << endl;  
	    }
	  cout << "vector: " << endl; 
	  for(Int_t ci=0; ci<A.GetNrows(); ci++)
	    cout << " " << y(ci);
	  cout << endl;
	}	
      
      // solve matrix only if segment was accepted for the fit
      if(acceptSeg)
	{
	  // numerically preferred method
	  if(debug_Fit)
	    cout << "** Fit is accepted: solve through SVD" << endl;
	  
	  // first bring the weights in place
	  TMatrixD Aw = A;
	  TVectorD yw = y;
	  for (Int_t irow = 0; irow < A.GetNrows(); irow++) 
	    {
	      TMatrixDRow(Aw,irow) *= 1;///e(irow);
		yw(irow) /= 1;//e(irow);
	    }
	  TDecompSVD svd(Aw);
	  Bool_t oksolve;
	  c_svd = new TVectorD(svd.Solve(yw,oksolve));
	  
	  if(debug_Fit)
	    {
	      cout << "Solution array: " << endl;
	      for(Int_t i=0; i<matrixRows; i++)
		cout << (*c_svd)(i) << "  ";
	      cout << endl;
	    }
	  
	  //invert matrix to get errors
	  Double_t det1;
	  TMatrixD Aerr(matrixRows,matrixRows);
	  Aerr = A;
	  Aerr.InvertFast(&det1);
	  if(debug_Fit)
	    {	
	      cout << "Parameter errors : sqrt(M^-1) " << endl; 
	      for(Int_t r=0; r<Aerr.GetNrows(); r++)
		{ 
		  for(Int_t ci=0; ci<Aerr.GetNcols(); ci++)
		    cout << " " << TMath::Sqrt(Aerr(r,ci));
		  cout << endl;  
		}
	    }
	  
	  Phi_err = TMath::Sqrt(Aerr(0,0));
	  X_Phi_err = TMath::Sqrt(Aerr(0,1));
	  X_err = TMath::Sqrt(Aerr(1,1));
	  if(nrVar>2)
	    {
	      Phi_T0_err = TMath::Sqrt(Aerr(0,2));
	      X_T0_err = TMath::Sqrt(Aerr(1,2));
	      T0_err = TMath::Sqrt(Aerr(2,2));
	    }	
	  
	  //compute t0
	  
	  if(nrVar>2)
	    t0.first =  (*c_svd)(2) / vdrift;
	  else t0.first=0.;
	  
	  if(debug_Fit)
	    cout << "t0 " << t0.first << "  vdrift " << vdrift << endl;
	  
	  reDoFit = false;
	  
	  //chi2 for each segment
	  //res = single hit residual
	  //resSum = residual sum in micron
	  //resSum2 = sqrt(square res sum/DOF) --> this is sigma in mm
	  
	  //compute residuals in segment if segment was in the fit....
	  if(nrVar>=2){
	    for(Int_t pi=0; pi<nrPnts; pi++)
	      if(ax[pi])
		{
		  Double_t res = ax[pi] - (ac[pi]*vdrift*t0.first+(*c_svd)(0)*ay[pi]+(*c_svd)(1));
		  chi2 += TMath::Power((res/sigmaPhi),2);
		  resSum += res * 10000;			//micron
		  resSum2 += TMath::Power(res * 10,2);	//mm
		  res0[pi] = res * 10000;
		}
	    if( (nPF-nrVarv)!=0 )
	      meanRes2 = TMath::Sqrt(resSum2 / (nPF-nrVarv) );
	    
	    if( (nPF-nrVar)!=0 )
	      chi2 /= (nPF-nrVar);
	  }	  
	  
	  NPT=nPF;
	  for(int i=0;i<nrPnts;i++)
	    if(ax[i]==0) 
	      ac[i]=0;
	  
	  Double_t totresSum2 = resSum2;
	  totMeanRes2 = TMath::Sqrt( totresSum2 / (nPF - nrVar) );
	  
	  if(debug_Fit && debug_FitRes)
	    {
	      cout << " resSum = " << resSum << endl; 
	      cout << " resSum2 = " << resSum2 << endl;
	      cout << " meanRes2 = " << meanRes2 << endl;
	      cout << " totMeanRes2=" << totMeanRes2 << endl;
	      cout << " nPF = " << nPF << "  totPar = " << totPar << " ==> NDF = " << (nPF - totPar) << endl; 
	    }	
	  
	  //check if hit residual is less then sigmaTimes Max(sigma,250 micron), otherwise reject and re-fit
	  Double_t evSigma=0.;
	  evSigma= sigmaPhi*10.;
	  
	  //variables for maximum residual search
	  Double_t maxRes = 0.;
	  Int_t maxPI = 0;
	  
	  if(nrVar>=2)
	    {	
	      if(debug_Fit && debug_FitRes)
		cout << "Now check residuals for sigma " << evSigma << endl;
	      for(Int_t pi=0; pi<nrPnts; pi++)
		if(ax[pi])
		  {
		    Double_t res = ax[pi] - (ac[pi]*(vdrift)*(t0.first)+(*c_svd)(0)*ay[pi]+(*c_svd)(1));
		    if(debug_Fit && debug_FitRes)
		      {
			cout << "Hit position " << ax[pi] - ac[pi]*(vdrift)*(t0.first);
			cout << ",  Fitted position " << (*c_svd)(0)*ay[pi]+(*c_svd)(1);
			cout << ",  Residual (mm) point " << pi << " -> " << res*10 << endl;
			cout << "  t0 " << t0.first << endl;
		      }
		    if( (TMath::Abs(res*10)) > (sigmaTimes * evSigma) )
		      {
			if(debug_Fit && debug_FitRes)
			  cout << "Point: " << (res*10) << 
			    " is " << sigmaTimes << " sigma out " << endl;
			if( TMath::Abs(res) > maxRes)
			  {
			    maxRes = TMath::Abs(res);
			    maxPI = pi;
			  }	
			
			if(numFit < numMaxFit)
			  nPF = nPF - 1;	
			reDoFit = true;
			if(debug_Fit && debug_FitRes) 
			  printf("===> reDoFit\n");
		      } // close if(res>sigmaTimes*evSigma)
		  }//end points loop
	    } // close if(nrVarv>=2)
	  
	  
	  //now reject max residual
	  if(maxRes != 0.)
	    {
	      if(debug_Fit)
		cout << "Rejected max res. and re-do fit: " << (maxRes*10) << " is " << sigmaTimes << " sigma out " << endl;
	      
	      ax[maxPI] = 0;	//set the wire at 0
	      if(numFit < numMaxFit)
		nPF = nPF - 1;	
	      reDoFit = true;
	      if(debug_Fit) 
		printf("   ===> reDoFit\n");
	    }
	  
	  if(reDoFit){
	    delete c_svd;
	    c_svd=NULL;
	  }
	 
	  okFit = true;
	  
	  Aw.Delete();
	  yw.Delete();
	  Aerr.Delete();
	  
	}//end if acceptSeg
      else
	{
	  if(debug_Fit)
	    cout << "** Fit not computed ! " << endl;
	  okFit = false;
	  reDoFit = false;
	  
	  delete c_svd;
	  c_svd=NULL;
	  
	}
      A.Delete();
      y.Delete();
      
    }//end while reDoFit
  
  m.first=-999.; q.first=-999.;
  if(okFit)
    {
      //return solution    
      if(nrVar>=3 )
        t0 =  make_pair( (*c_svd)(2) / vdrift , T0_err / vdrift);

        m = make_pair( (*c_svd)(0) , Phi_err);
        q = make_pair( (*c_svd)(1), X_err);

      delete c_svd;
      c_svd=NULL;
      
      if(debug_Fit)
        {
          cout << "Solution after fit loop:" << endl;
          cout << "x = " << m.first << " * y + " << q.first << endl;
          cout << "t0 = " << t0.first << endl;
          for(int i=0;i<nrPnts;i++) cout<<ac[i];
          cout<<endl;
          cout << "NPT = "<< NPT <<", chi2 = "<< chi2<< endl;
          cout << "ok fit flag: " << okFit << endl;
        }
    } // close if(okFit)
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FIT::FIT_global(double *sigmaPhi,int *nrPnts, int *nrMinPnts, int *nrVarv, vector<vector<double> > axf, vector<vector<double> > at, vector<vector<double> > ay, vector<vector<int> > ac, double *m, double *sm, double *q, double *sq, double &t0, double &st0, double &vdrift, double *chi2, int *NPT, bool &okFit, int &fit_CH){
  
  // 2 SUPERLAYERS t0 FIT
  // xi = fi + ci*v(ti-t0) = fi + ci*v*ti - ci*z0
  // xi = x coordinates along chamber FE, fi = wire coordinate, ci = left/right code
  // v = drift velocity (fixed), ti = drift time, 
  // t0 = time of muon passage respect to the clock 
  // xi = a + m*y --> phi chamber 1     (nrVar[0], nrPnts[0])
  // xi = b + n*y --> theta chamber 1	(nrVar[1], nrPnts[1])
  // xi = c + p*y --> phi chamber 2	(nrVar[2], nrPnts[2])
  // xi = d + q*y --> theta chamber 2	(nrVar[3], nrPnts[3])
  // parameters: slope and intercept of the 2 chambers, t0, vdrift --(10 parameters)
  // redo the fit in case one residual is out of sigmaTimes sigma
  
  if(debug_Fit_Glo) printf("...Starting FIT_global...\n");
  
  Bool_t reDoFit = true;  
  Int_t numFit = 0;
  okFit=false;
  static const Int_t numMaxFit = 5;
  TVectorD * c_svd = 0; 
  
  //max number of parameters
  Int_t totParMax = 0;
  Int_t totPar = 0;
  Int_t nrVar=0;
  Int_t nPF[4]={0,0,0,0};
  Int_t nrSeg=4;
  
  //reset errors
  Double_t X_err = 0.;
  Double_t X_Phi_err = 0.;
  Double_t X_Z_err = 0.;
  Double_t X_Theta_err = 0.;
  Double_t X_T0_err = 0.;
  
  Double_t Phi_err = 0.;
  Double_t Phi_Z_err = 0.;
  Double_t Phi_Theta_err = 0.;
  Double_t Phi_T0_err = 0.;
  
  Double_t Z_err = 0.;
  Double_t Z_Theta_err = 0.;
  Double_t Z_T0_err = 0.;
  
  Double_t Theta_err = 0.;
  Double_t Theta_T0_err = 0.;
  
  Double_t T0_err = 0.;
  
  Double_t resSum[4]={0.,0.,0.,0.};
  Double_t resSum2[4]={0.,0.,0.,0.};
  Double_t meanRes2[4] = {0.,0.,0.,0.};
  Double_t totMeanRes2 = 0.;
  
  vdrift = 0.00547; //<------- DRIFT VELOCITY
  static const Float_t sigmaTimes = 3.;         //reject hits if |residual| > sigmaTimes * sigma
  
  //weights
  // Double_t er = 1.; // 1./TMath::Power(450.,2);
  Double_t er = 1/TMath::Power(0.03,2);
  Double_t ae[8];
  for(int i=0;i<8;i++)  ae[i]=double(er);
  
  if(debug_Fit_Glo) 
    for(int i=0;i<nrSeg;i++) 
      printf("N.Hit to fit %d\n",nrPnts[i]);
  
  while(reDoFit && numFit < numMaxFit)
    {
      if(debug_Fit_Glo)
	{
	  cout << "\n Go with the fit attemp " << numFit << " ! number of parameters: " << endl;
	  cout << "ch1 phi " << nrVarv[0] << " theta " << nrVarv[1] << endl; 	
	  cout << "ch2 phi " << nrVarv[2] << " theta " << nrVarv[3] << endl; 	
	}
      
      //vdrift fit flag
      Bool_t vdriftFit = false;
      
      //reset number of points to fit
      for(Int_t sr=0; sr<nrSeg; sr++)
	nPF[sr] = 0;
      
      //count this fit number of points
      for(Int_t sr=0; sr<nrSeg; sr++)
  	{
	  if(nrSeg >= (sr+1))
	    {	
	      if(debug_Fit_Glo)
		cout << "Now counting fitting points " << endl;
	      for(Int_t pi=0; pi<nrPnts[sr]; pi++)
		{
		  if(axf[sr][pi] && nrVarv[sr]>=2)
		    nPF[sr] += 1;
		}
	    }	
	  if(debug_Fit_Glo)
	    cout << "Points to fit in segment " << sr << " are " << nPF[sr] << endl;	
  	}	
      
      Int_t totPF = 0;
      totPar = 0;
      
      //count total number of points and total parameters to fit
      for(Int_t sr=0; sr<nrSeg; sr++)
	{
	  totPF += nPF[sr];
	  totPar += TMath::Min(2,nrVarv[sr]);
	}
      if(nrVarv[0]>=3 || nrVarv[1]>=3 || nrVarv[2]>=3 || nrVarv[3]>=3 )
	totPar += 1;
      if(nrVarv[0]==4 || nrVarv[1]==4 || nrVarv[2]==4 || nrVarv[3]==4 ){
	totPar += 1;
	vdriftFit = true;
      }
      
      if(debug_Fit_Glo)
	cout << "total points to fit are " << totPF << ", total parameters are " << totPar << endl;	
      
      totParMax = totPar;
      
      //check totPF > parameter number
      if(totPF <= totPar)
        {
	  if(debug_Fit_Glo)
	    cout << "Don't fit : total num.points <= num.parameters ! " << endl;
	  //go to the end of while loop and exit - if there was a previous loop with okFit keep the latest results
	  break;
	}
      
      //check: if nPF[sr]>=minimum number of points for each requested SL
      Bool_t breakFlag = false;
      for(Int_t s=0; s<nrSeg; s++)
	{
	  //if this segment was requested in the fit
	  if(nrVarv[s])
	    {
	      if(debug_Fit_Glo)
		cout << "Check segment " << s << ": if points >= nrMinPoints" << endl; 
	      if(debug_Fit_Glo)
		cout << "nrMinPoints = " << nrMinPnts[s] << endl; 
	      
	      //check SL
	      if( nPF[s] < TMath::Max(nrVarv[s],nrMinPnts[s]) )
		{
		  if(debug_Fit_Glo)
		    cout << "Segment " << s << " has: " << nPF[s] << " points < " << TMath::Max(nrVarv[s],nrMinPnts[s])<< " nrMinPoints ... reject segment" << endl; 
		  breakFlag = true;
		}
	    }
	}
      
      //if one requested segment doesn't have the minimum points then end fit loop
      if(breakFlag)
	{
	  //okFit=false => reject the event, okFit=true => keep the latest result
	  okFit = false; 
	  break;
	}
      
      //count fit number
      numFit += 1; 
      
      //reset chi2 and residuals
      Double_t res0[nrPnts[0]]; Double_t res1[nrPnts[1]]; Double_t res2[nrPnts[2]]; Double_t res3[nrPnts[3]];
      for(Int_t j=0; j<nrSeg; j++)
	for(Int_t k=0; k<nrPnts[j]; k++)
	  {
	    res0[k] = 0.;
	    res2[k] = 0.;
	    res1[k] = 0.;
	    res3[k] = 0.;
	  }
      
      //reset chi2 and residuals
      for(Int_t k=0; k<nrSeg; k++)
  	{
	  chi2[k] = 0.;
	  resSum[k] = 0.;
	  resSum2[k] = 0.;
	  meanRes2[k] = 0.;
	  NPT[k]=0;
	}
      totMeanRes2 = 0.;
      
      if(debug_FitRes)
	cout << "Fit number " << numFit << endl;
      
      // define matrix
      Int_t matrixRows = totPar;
      //Int_t matrixRows = nrVar;
      
      TMatrixD A(matrixRows,matrixRows);
      TVectorD y(matrixRows);
      Int_t lr = matrixRows-1; //last row/column index
      
      //this sums go over hits of all segments
      Double_t Scc = 0.;
      Double_t Scct = 0.;
      Double_t Scctt = 0.;  
      Double_t Sfc = 0.;
      Double_t Sfct = 0.;
      
      Double_t ax[4][8];	//hit positions
      for(int i=0;i<4;i++)
	for(int j=0;j<8;j++)
	  ax[i][j]=0.;
     
      Bool_t acceptSeg = false;
      for(Int_t s=0; s<nrSeg; s++)
  	{
	  if(debug_Fit_Glo)
	    cout << "Examing segment " << s << " for fit... " << endl; 	
	  
	  // ** perform checks on segment for proper fit
	  // 1. check that parameter number is 2 at least...
	  if(nrVarv[s]<2)
	    {
	      if(debug_Fit_Glo)
		cout << "Don't fit segment " << s << ": min parameter number should be 2!" << endl;
	      continue;
	    }
	  
	  // 2. check: if mu crosses all cells at the same side of the wire do not fit t0 and vdrift!
	  Double_t cSum = 0;
	  Double_t pSum = 0;
	  
	  for (int ji=0; ji<nrPnts[s]; ji++)
	    {
	      if(axf[s][ji]!=0)
		{
		  cSum += ac[s][ji];
		  pSum += 1;
		}
	    } 
	  if(TMath::Abs(cSum) == pSum)
	    {
	      if(debug_Fit_Glo)
		cout << "Don't fit t0 in segment " << s << ": all cell at same side of wires!" << endl;
	      nrVarv[s] = 2;
	    }
	  
	  Int_t nrVar = nrVarv[s];
	  
	  //ok, can do the fit
	  acceptSeg = true;
	  
	  // hit position x1 .... x8
	  for(Int_t i=0; i<nrPnts[s]; i++)
	    ax[s][i] = axf[s][i]+ac[s][i]*at[s][i]*vdrift;  
	  
	  //test hit position
	  if(debug_Fit_Glo)
	    {		
	      for(Int_t i=0; i<nrPnts[s]; i++)
		cout << " ax[" << s << "][" << i << "]=" << ax[s][i] << "  "; 
	      cout << endl;
	      for(Int_t i=0; i<nrPnts[s]; i++)
		cout << " ay[" << s << "][" << i << "]=" << ay[s][i] << "  "; 
	      cout << endl;
	    }
	  
	  // compute sums
	  // this sums go over hits of all segments
	  Double_t Sx = 0.;
	  Double_t Sxy = 0.;
	  Double_t Sy = 0.;
	  Double_t Syy = 0.;
	  Double_t Sc = 0.;
	  Double_t Scy = 0.;
	  Double_t Sct = 0.;
	  Double_t Scty = 0.;
	  Double_t Sf = 0.;
	  Double_t Sfy = 0.;
	  Double_t Sw = 0.;
	  
	  // first row to insert matrix elements for this segment
	  Int_t fr = 0;
	  if(fit_CH==1 || fit_CH==12) fr=2*s; 
	  if(fit_CH==2) fr=2*(s-2); 
	  
	  for (int j=0; j<nrPnts[s]; j++)
	    {
	      //add wire to sum only if it's not empty
	      if(axf[s][j])
		{
		  Sx  	= Sx + ae[j]*ax[s][j];
		  Sxy 	= Sxy + ae[j]*ax[s][j]*ay[s][j];
		  Sy  	= Sy + ae[j]*ay[s][j];
		  Syy 	= Syy + ae[j]*ay[s][j]*ay[s][j];
		  Sw	= Sw + ae[j];
		  if(nrVar>=3)
		    {
		      Scc 	= Scc + ae[j]*ac[s][j]*ac[s][j];
		      Scy 	= Scy + ae[j]*ac[s][j]*ay[s][j];
		      Sc 	= Sc + ae[j]*ac[s][j];
		      
		      Sf	= Sf + ae[j]*axf[s][j];
		      Sfc	= Sfc + ae[j]*axf[s][j]*ac[s][j];
		      Sfy	= Sfy + ae[j]*axf[s][j]*ay[s][j];
		      
		      Scty 	= Scty + ae[j]*ac[s][j]*at[s][j]*ay[s][j];
		      Sct 	= Sct + ae[j]*ac[s][j]*at[s][j];
		      Scct 	= Scct + ae[j]*ac[s][j]*ac[s][j]*at[s][j];
		    }
		  if(nrVar>=4)
		    {
		      //TODO: sistemare matrice!!!
		      Sfct	= Sfct + ae[j]*axf[s][j]*ac[s][j]*at[s][j];
		      Scctt 	= Scctt + ae[j]*ac[s][j]*ac[s][j]*at[s][j]*at[s][j];
		    }
		}//end if
	    } //end loop points
	  
	  if(nrVar>=2)
	    {
	      A(fr,fr)=Syy;	A(fr,fr+1)=Sy;			
	      A(fr+1,fr)=Sy;	A(fr+1,fr+1)=Sw;
	      
	      if(nrVar==2)
		{
		  y(fr)=Sxy;
		  y(fr+1)=Sx;
		}
	    }
	  
	  if(nrVar==3)
	    {
	      Int_t lrnew = lr;
	      if(vdriftFit == true)
		lrnew = lr - 1;
	      
	                                                A(fr,lrnew)=Scy;
	                                                A(fr+1,lrnew)=Sc;
	      A(lrnew,fr)=Scy;	   A(lrnew,fr+1)=Sc;	A(lrnew,lrnew)=Scc;

	      y(fr)=Sfy + vdrift * Scty;	
	      y(fr+1)=Sf  + vdrift * Sct;
	      y(lrnew)=Sfc + vdrift * Scct;
	    }
	  
	  if(nrVar==4)
	    // TODO!!!
	    {
	                                                  A(fr,lr-1)=Scy;
	                                                  A(fr+1,lr-1)=Sc;
	      A(lr-1,fr)=Scy;       A(lr-1,fr+1)=Sc;      A(lr-1,lr-1)=Scc;
	      
	      
	      
	                                                     A(fr,lr)=-Scty;
	                                                     A(fr+1,lr)=-Sct;	
	                                                     A(lr-1,lr)=-Scct;	
	      A(lr,fr)=Scty;	A(lr,fr+1)=Sct;	A(lr,lr-1)=Scct;	A(lr,lr)=-Scctt;	
	      
	      y(fr)=Sfy;	
	      y(fr+1)=Sf;
	      y(lr-1)=Sfc;
	      y(lr)=Sfct;
	    }
	} //end loop segment
      
      
      if(debug_Fit_Glo)
	{	
	  //cout << "matrice da invertire "<<matrixRows<<" x "<<matrixRows<<":"<< endl; 
	  for(Int_t r=0; r<A.GetNrows(); r++)
	    { 
	      for(Int_t ci=0; ci<A.GetNcols(); ci++)
		cout << " " << A(r,ci);
	      cout << endl;
	    }
	  cout << "vector: " << endl; 
	  for(Int_t ci=0; ci<A.GetNrows(); ci++)
	    cout << " " << y(ci);
	  cout << endl;
	}	
      
      // solve matrix only if minimum 1 segment was accepted for the fit
      if(acceptSeg)
	{
	  // numerically  preferred method
	  if(debug_Fit_Glo)
	    cout << " Fit is accepted: solve through SVD" << endl;
	  
	  // first bring the weights in place
	  TMatrixD Aw = A;
	  TVectorD yw = y;
	  for (Int_t irow = 0; irow < A.GetNrows(); irow++) 
	    {
	      TMatrixDRow(Aw,irow) *= 1;///e(irow);
		yw(irow) /= 1;//e(irow);
	    }
	  TDecompSVD svd(Aw);
	  Bool_t oksolve;
	  c_svd = new TVectorD(svd.Solve(yw,oksolve));
	  
	  if(debug_Fit_Glo)
	    {
	      cout << "Solution array: " << endl;
	      for(Int_t i=0; i<matrixRows; i++)
		cout << (*c_svd)(i) << "  ";
	      cout << endl;
	    }
	  
	  // *************************************************************
	  // TODO: aggiustarla x i 4 seg!!!! 
	  
	  //invert matrix to get errors
	  Double_t det1;
	  TMatrixD Aerr(matrixRows,matrixRows);
	  Aerr = A;
	  Aerr.InvertFast(&det1);
	  if(debug_Fit_Glo)
	    {	
	      cout << "Parameter errors : sqrt(M^-1) " << endl; 
	      for(Int_t r=0; r<Aerr.GetNrows(); r++)
		{ 
		  for(Int_t ci=0; ci<Aerr.GetNcols(); ci++)
		    cout << " " << TMath::Sqrt(Aerr(r,ci));
		  cout << endl;  
		}
	    }
	  
	  Phi_err = TMath::Sqrt(Aerr(0,0));
	  X_Phi_err = TMath::Sqrt(Aerr(0,1));
	  Phi_Theta_err = TMath::Sqrt(Aerr(0,2));
	  Phi_Z_err = TMath::Sqrt(Aerr(0,3));
	  X_err = TMath::Sqrt(Aerr(1,1));
	  X_Theta_err = TMath::Sqrt(Aerr(1,2));
	  X_Z_err = TMath::Sqrt(Aerr(1,3));
	  Theta_err = TMath::Sqrt(Aerr(2,2));
	  Z_Theta_err = TMath::Sqrt(Aerr(2,3));
	  Z_err = TMath::Sqrt(Aerr(3,3));
	  if(nrVar>2)
	    {
	      Phi_T0_err = TMath::Sqrt(Aerr(0,4));
	      X_T0_err = TMath::Sqrt(Aerr(1,4));
	      Theta_T0_err = TMath::Sqrt(Aerr(2,4));
	      Z_T0_err = TMath::Sqrt(Aerr(3,4));
	      T0_err = TMath::Sqrt(Aerr(4,4));
	    }	

	  
	  // TODO...controllare nrVarv nel caso di hit tutte con lo stesso segno!!!
	  if((nrVarv[0]==4 || nrVarv[1]==4 || nrVarv[2]==4 || nrVarv[3]==4) && (*c_svd)(totPar-1) )
	    vdrift =  (*c_svd)(totPar-1); 
	  
	  if( nrVarv[0]>=3 || nrVarv[1]>=3 || nrVarv[2]>=3 || nrVarv[3]>=3)
	    if( nrVarv[0]==4 || nrVarv[1]==4 || nrVarv[2]==4 || nrVarv[3]==4 )
	      t0 =  (*c_svd)(totPar-2) / vdrift;
	    else 
	      t0 =  (*c_svd)(totPar-1) / vdrift;
	 
	  	  
	  if(debug_Fit_Glo)
	    cout << "t0 " << t0 << ",  vdrift " << vdrift << endl;
	  
	  reDoFit = false;
	  
	  //check vdrift and redo the fit in case -> equal tdrifts in cell
	  if(vdrift < TMath::Power(10,-3))
	    {
	      if(debug_Fit_Glo)
		cout << "vdrift " << vdrift << " not accetable, redo the fit without vrift!" << endl;
	      reDoFit = true;
	      numFit -= 1;
	      for(Int_t rf=0; rf<nrSeg; rf++)
		if(nrVarv[rf]==4)
		  nrVarv[rf]=3;
	    }	 
	  
	  //chi2 for each segment
	  //res = single hit residual
	  //resSum = residual sum in micron
	  //resSum2 = sqrt(square res sum/DOF) --> this is sigma in mm
	  for(Int_t sr=0; sr<nrSeg; sr++)
	    {
	      //compute residuals in segment if segment was in the fit....
	      if(nrSeg>=(sr+1) &&  nrVarv[sr]>=2)
		{	
		  //0: ch 1 phi;     1: ch1 theta;    2: ch2 phi;    3: ch2 theta
		  for(Int_t pi=0; pi<nrPnts[sr]; pi++)
		    {
		      if(axf[sr][pi])
			{
			  int ss = 0;
			  if(fit_CH==1 || fit_CH==12) ss = sr;
			  if(fit_CH==2) ss = sr-2;
			  Double_t res = axf[sr][pi]+ac[sr][pi]*vdrift*at[sr][pi] - 
			    (ac[sr][pi]*vdrift*t0+(*c_svd)(ss*2)*ay[sr][pi]+(*c_svd)(ss*2+1));
			  chi2[sr] += TMath::Power((res/sigmaPhi[sr]),2);
			  
			  resSum[sr] += res * 10000;			//micron
			  resSum2[sr] += TMath::Power(res * 10,2);	//mm
			  if(sr==0)
			    {
			      res0[pi] = res * 10000;
			    }
			  if(sr==1)
			    {
			      res1[pi] = res * 10000;
			    }
			  if(sr==2)
			    {
			      res2[pi] = res * 10000;
			    }
			  if(sr==3)
			    {
			      res3[pi] = res * 10000;    
			    }
			}
		    }
		  //if( (nrPnts[sr]-nrVarv[sr]) )
		  if( (nPF[sr]-nrVarv[sr]) )
		    {
// 		      meanRes2[sr] = TMath::Sqrt(resSum2[sr] / (nrPnts[sr]-nrVarv[sr]) );
// 		      chi2[sr] = (chi2[sr] / (nrPnts[sr]-nrVarv[sr]) );
		      meanRes2[sr] = TMath::Sqrt(resSum2[sr] / (nPF[sr]-nrVarv[sr]) );
		      chi2[sr] = (chi2[sr] / (nPF[sr]-nrVarv[sr]) );
		    }
		}
	      
	      NPT[sr]=nPF[sr];
	      for(int i=0;i<nrPnts[sr];i++)
		if(axf[sr][i]==0) 
		  ac[sr][i]=0;
	    }
	  
	  Double_t totresSum2 = resSum2[0] + resSum2[1] + resSum2[2] + resSum2[3];
	  //Double_t totnrPnt = nrPntsv[0] + nrPntsv[1] + nrPntsv[2] + nrPntsv[3];  
	  totMeanRes2 = TMath::Sqrt( totresSum2 / (totPF - totPar) );
	  
	  if(debug_Fit_Glo && debug_FitRes)
	    {
	      for(int j=0; j<nrSeg; j++)
		{
		  cout << " resSum[" << j << "]=" << resSum[j] << endl; 
		  cout << " resSum2[" << j << "]=" << resSum2[j] << endl;
		  cout << " meanRes2[" << j << "]=" << meanRes2[j] << endl;
		}
	      cout << " totMeanRes2=" << totMeanRes2 << endl;
	      cout << " totPF = " << totPF << "  totPar = " << totPar << " ==> NDF = " << (totPF - totPar) << endl; 
	    }	
	  
	  //check if hit residual is less then sigmaTimes Max(sigma,250 micron), otherwise reject and re-fit
	  Double_t evSigma[4] = {0.,0.,0.,0.};
	  for(int j=0; j<nrSeg; j++)
	    evSigma[j] = sigmaPhi[j]*10.;
	  
	  //variables for maximum residual search
	  Double_t maxRes = 0.;
	  Int_t maxSR = 0;
	  Int_t maxPI = 0;
	  
	  for(Int_t sr=0; sr<nrSeg; sr++)
	    {
	      if(nrSeg >= (sr+1) &&  nrVarv[sr]>=2)
		{	
		  if(debug_FitRes)
		    cout << "Now check residuals for sigma " << evSigma[sr] << endl;
		  for(Int_t pi=0; pi<nrPnts[sr]; pi++)
		    {
		      if(axf[sr][pi])
			{
			  int ss=0;
			  if(fit_CH==1 || fit_CH==12) ss=sr;
			  if(fit_CH==2) ss=sr-2;
			  Double_t res = axf[sr][pi]+ac[sr][pi]*vdrift*at[sr][pi] - 
			    (ac[sr][pi]*vdrift*t0+(*c_svd)(ss*2)*ay[sr][pi]+(*c_svd)(ss*2+1));
			  if(debug_FitRes)
			    {
			      cout << "Hit position " << axf[sr][pi]+ac[sr][pi]*vdrift*at[sr][pi] - ac[sr][pi]*vdrift*t0;
			      cout << ",  Fitted position " << (*c_svd)(ss*2)*ay[sr][pi]+(*c_svd)(ss*2+1);
			      cout << ",  Residual (mm) point " << pi << " -> " << res*10 << endl;
			      
			      cout << "  t0 " << t0 << endl;
			    }
			  if( TMath::Abs(res*10) > sigmaTimes * evSigma[sr] )
			    {
			      if(debug_FitRes || debug_Fit_Glo)
				cout << "Point: " << (res*10) << 
				  " is " << sigmaTimes << " sigma out " << endl;
			      //aw[sr][pi] = 0;
			      if( TMath::Abs(res) > maxRes)
				{
				  maxRes = TMath::Abs(res);
				  maxSR = sr;
				  maxPI = pi;
				}	
			      
			      if(numFit < numMaxFit)
				nPF[sr] = nPF[sr] - 1;	
			      reDoFit = true;
			    }
			}
		    }//end points loop
		}
	    }//end segment loop
	  
	  //now reject max residual
	  if(maxRes != 0.)
	    {
	      if(debug_FitRes || debug_Fit_Glo)
		cout << "Rejected max res. and re-do fit: " << (maxRes*10) << " is " << sigmaTimes << " sigma out " << endl;
	      
	      axf[maxSR][maxPI] = 0;	//set the wire at 0
	      if(numFit < numMaxFit)
		nPF[maxSR] = nPF[maxSR] - 1;	
	      reDoFit = true;
	      if(debug_Fit_Glo) 
		printf("   ===> reDoFit\n");
	    }
	  
	  if(reDoFit){
	    delete c_svd;
	    c_svd=NULL;
	  }
	  
	  okFit = true;
	  
	  Aw.Delete();
	  yw.Delete();
	  Aerr.Delete();
	  
	}//end if acceptSeg
      else
	{
	  if(debug_Fit_Glo )
	    cout << "** Fit not computed ! " << endl;
	  okFit = false;
	  reDoFit = false;

	  delete c_svd;
	  c_svd=NULL;
	}
      A.Delete();
      y.Delete();
      
    }//end while reDoFit
  
  for(int i=0;i<4;i++){
    m[i]=-999.; q[i]=-999.;
    sm[i]=-999.; sq[i]=-999.;
  }
  if(okFit)
    {
      //return solution
      if( (nrVarv[0]==4 || nrVarv[1]==4 || nrVarv[2]==4 || nrVarv[3]==4) && (*c_svd)(totPar-1))
	vdrift =  (*c_svd)(totPar-1); 
      
      if(nrVarv[0]>=3 || nrVarv[1]>=3 || nrVarv[2]>=3 || nrVarv[3]>=3 )
	if( nrVarv[0]==4 || nrVarv[1]==4 || nrVarv[2]==4 || nrVarv[3]==4 )
	  t0 =  (*c_svd)(totPar-2) / vdrift;
	else
	  t0 =  (*c_svd)(totPar-1) / vdrift;
      st0 = T0_err / vdrift;
      
      if(fit_CH==1 || fit_CH==12) {
	m[0] = (*c_svd)(0);
	q[0] = (*c_svd)(1);
	sm[0] = Phi_err;
	sq[0] = X_err;
      }
      
      if(nrSeg>=2 && totParMax>=3)
	{	
	  if(fit_CH==1 || fit_CH==12){
	    m[1] = (*c_svd)(2);
	    q[1] = (*c_svd)(3);
	    sm[1] = Theta_err;
	    sq[1] = Z_err;

	    
	  }
	}
      
      if(nrSeg>=3 && totParMax>=4)
	{
	  if(fit_CH==2) {
	    m[2] = (*c_svd)(0);
	    q[2] = (*c_svd)(1);
	    sm[2] = Phi_err;
	    sq[2] = X_err;
	  }
	  if(fit_CH==12) {
	    m[2] = (*c_svd)(4);
	    q[2] = (*c_svd)(5);
	    sm[2] = Phi_err;
	    sq[2] = X_err;
	  }
	}
      
      if(nrSeg>=4 && totParMax>=5)
	{
	  if(fit_CH==2) {
	    m[3] = (*c_svd)(2);
	    q[3] = (*c_svd)(3);
	    sm[3] = Theta_err;
	    sq[3] = Z_err;
	  }
	  if(fit_CH==12) {
	    m[3] = (*c_svd)(6);
	    q[3] = (*c_svd)(7);
	    sm[3] = Theta_err;
	    sq[3] = Z_err;
	  }
	  
	}
      
      delete c_svd;
      c_svd=NULL;
      
      if(debug_Fit_Glo)
	{
	  cout << "Solution after global fit loop:" << endl;
	  cout << "ch1_phi    -->  m = " << m[0] << ",  a = "<< q[0] << endl;
	  cout << "           -->  sm= " << sm[0] << ",  sa= "<< sq[0] << endl;
	  cout << "ch1_theta  -->  n = " << m[1] << ",  b = "<< q[1] << endl;
	  cout << "           -->  sn= " << sm[1] << ",  sb= "<< sq[1] << endl;
	  cout << "ch2_phi    -->  p = " << m[2] << ",  c = "<< q[2] << endl;
	  cout << "           -->  sp= " << sm[2] << ",  sc= "<< sq[2] << endl;
	  cout << "ch2_theta  -->  q = " << m[3] << ",  d = "<< q[3] << endl;
	  cout << "           -->  sq= " << sm[3] << ",  sd= "<< sq[3] << endl;
	  cout << "global t0  -->  " << t0 << endl;
	  cout << "      st0  -->  " << st0 << endl;
	  cout << "global drift velocity --> " << vdrift << endl;
	  cout << "NPT[0] = "<< NPT[0] <<", chi2[0] = "<< chi2[0]<< endl;
	  cout << "NPT[1] = "<< NPT[1] <<", chi2[1] = "<< chi2[1]<< endl;
	  cout << "NPT[2] = "<< NPT[2] <<", chi2[2] = "<< chi2[2]<< endl;
	  cout << "NPT[3] = "<< NPT[3] <<", chi2[3] = "<< chi2[3]<< endl;
	  cout << "ok fit flag " << okFit << endl;
	}		
    }
  
  
  return;
}

