#ifndef topkinsolver_h
#define topkinsolver_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "TSystem.h"
#include "TObject.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#endif


/**
   @short this structure olds the kinematics solution found for a TTbar system
*/
struct TTbarSolution_t
{
  TLorentzVector pl1,pb1,pn1,pt1;
  TLorentzVector pl2,pb2,pn2,pt2;            
};
typedef std::vector<TTbarSolution_t> TTbarSolutionCollection_t;

//some utility functions
bool mTTbarOrder( TTbarSolution_t a, TTbarSolution_t b);
double getMT(TLorentzVector &pbl, TLorentzVector &pn);
std::vector<double> getMT2( TTbarSolution_t &sol);

/**
   @short the algorithm actual work is done here...
 */
class TopKinSolver
{
public:
      
  /**
     @short CTOR
  */
  TopKinSolver(double ex=0.01, double ey=0.01, double ez=0.01, 
	       double maxDeltaMw=3, double maxDeltaMt=2,
	       int maxIter=30)
    : status_(true), status2_(true), error_("None"), 
      eVec_(ex,ey,ez), 
      maxDeltaMw_(maxDeltaMw), maxDeltaMt_(maxDeltaMt),
      maxIter_(maxIter)
  {
	
  }
      
  /**
     @short getters
  */
  bool status()                   { return status_; }
  std::string what()              { return error_; }
  TVector3 getConstraint()        { return constraint_; }

  /**
     @short setters
  */
  void setConstraint(TVector3 constraint) { constraint_=constraint; }

      
  /**
     @short wraps up the procedure of guessing and solving the equations
  */
  std::vector<TTbarSolution_t> findSolutions(TLorentzVector pl1, TLorentzVector pb1, 
					     TLorentzVector pl2, TLorentzVector pb2, 
					     TVector3 constraint,
					     double mw)
  {
    std::vector<TTbarSolution_t> sols;
	  
    //base guess for neutrino direction (collinear to lepton)	 
    TVector3 baseDir = pl1.Vect().Unit();
    TVector3 orthDir = baseDir.Orthogonal();
    TVector3 baseGuess = baseDir*(mw/(2*pl1.Mag()));
	  
    int nphix(6), nphiy(6);
    for(double phix=0; phix< 2*TMath::Pi(); phix += 2*TMath::Pi()/nphix)
      {
	TVector3 preGuess = baseGuess;
	preGuess.Rotate( phix, orthDir );
	TVector3 newOrthDir = preGuess.Cross( orthDir ).Unit();
	      
	for(double phiy=0; phiy< 2*TMath::Pi(); phiy += 2*TMath::Pi()/nphiy)
	  {
	    TVector3 guess = preGuess;
	    guess.Rotate( phiy, newOrthDir );
		  
	    //solve the equations
	    TLorentzVector theSol = nSolve(pl1,pb1,pl2,pb2,constraint,guess,mw);
	    if(!status()) continue;
		  
	    //check if result is different from the one obtained
	    bool isRedundant(false);
	    for(std::vector<TTbarSolution_t>::iterator it=sols.begin();
		it != sols.end();
		it++)
	      {
		double dR = it->pn1.DeltaR(theSol);
		TLorentzVector diff = it->pn1 - theSol;
		if(diff.Vect().Mag() <  eVec_.Mag() || dR<0.1) { isRedundant=true; break;}
	      }
	    if(isRedundant) continue;

	    //save result
	    TTbarSolution_t newSol;
	    newSol.pl1= pl1; newSol.pb1=pb1; newSol.pn1= theSol;                  newSol.pt1 = newSol.pl1+newSol.pb1+newSol.pn1;
	    newSol.pl2 =pl2; newSol.pb2=pb2; newSol.pn2 = constrainedMom(theSol); newSol.pt2 = newSol.pl2+newSol.pb2+newSol.pn2;
	    sols.push_back( newSol );

	    //cross check: if new solution find the displaced one for consistency should get X=(0,0,0)
	    // 		    TLorentzVector theSol2 = constrainedMom( theSol );
	    // 		    std::vector<TVector3> solsX = nSolve2(pl1,pb1,theSol,pl2,pb2,theSol2, mw);
	    // 		    for(std::vector<TVector3>::iterator it = solsX.begin();
	    // 		    it != solsX.end();
	    // 		    it++)
	    // 		    {
	    // 		    TVector3 nuVec = theSol.Vect() + *it;
	    // 		    TLorentzVector nu(nuVec,nuVec.Mag());
	    // 		    sols.push_back( nu );
	  }
      }
	  
    sort(sols.begin(),sols.end(), mTTbarOrder );
    return sols;
  }
      
  /**
   */
  TLorentzVector nSolve(TLorentzVector pl1, TLorentzVector pb1, 
			TLorentzVector pl2, TLorentzVector pb2, 
			TVector3 constraint,
			TVector3 guess,
			double mw)
  {
    status_= true;
    error_ = "";
    int nIter(0);
    bool hasConverged=false;
	
    TLorentzVector curMom = TLorentzVector(guess,guess.Mag());
    constraint_=constraint;

    do{

      nIter++;
      TVector3 fVec = f(pl1,pb1,pl2,pb2,curMom,mw);
      TMatrixD jacobian = dfdp( pl1, pb1, pl2, pb2, curMom);

      if(status_)
	{
	  TMatrixD mtr = inverse( jacobian, status_ );

	  if(status_)
	    {
	      TVector3 deltaMom(0,0,0);
	      deltaMom.SetX( -mtr(0,0)*fVec.X()-mtr(0,1)*fVec.Y()-mtr(0,2)*fVec.Z() );
	      deltaMom.SetY( -mtr(1,0)*fVec.X()-mtr(1,1)*fVec.Y()-mtr(1,2)*fVec.Z() );
	      deltaMom.SetZ( -mtr(2,0)*fVec.X()-mtr(2,1)*fVec.Y()-mtr(2,2)*fVec.Z() );
		  
	      TVector3 nextMom = curMom.Vect() + deltaMom;
	      curMom = TLorentzVector(nextMom,nextMom.Mag());
		  
	      if(deltaMom.Mag() < eVec_.Mag() ) hasConverged=true;
	    }
	}
    }while(!hasConverged && nIter<maxIter_ && status_);

    if(!hasConverged)
      {

	TLorentzVector otherNu = constrainedMom( curMom );
	TLorentzVector w1 = pl1+curMom;
	TLorentzVector t1 = w1+pb1;
	TLorentzVector w2 = pl2+otherNu;
	TLorentzVector t2 = w2+pb2;
	    
	double deltaMw1 = TMath::Abs( w1.M() - mw );
	double deltaMw2 = TMath::Abs( w2.M() - mw );
	double deltaMt = TMath::Abs( t1.M() - t2.M() );
	    
	// if(deltaMt > maxDeltaMt_)
	// 	      {		
	// 		std::cout << " l1=(" << pl1.Px() << " ; " << pl1.Py() << " ; " << pl1.Pz() << ") "
	// 			  << " n1=(" << curMom.Px() << " ; " << curMom.Py() << " ; " << curMom.Pz() << ") "
	// 			  << " b1=(" << pb1.Px() << " ; " << pb1.Py() << " ; " << pb1.Pz() << ") " << endl
	// 			  << " l2=(" << pl2.Px() << " ; " << pl2.Py() << " ; " << pl2.Pz() << ") "
	// 			  << " n2=(" << otherNu.Px() << " ; " << otherNu.Py() << " ; " << otherNu.Pz() << ") "
	// 			  << " b2=(" << pb2.Px() << " ; " << pb2.Py() << " ; " << pb2.Pz() << ") " << endl << endl
	// 			  << "t1=" << t1.E() << ";" << t1.P() << " " << t1.M()<< endl
	// 			  << "t2=" << t2.E() << ";" << t2.P() << " " << t2.M() << endl;
	// 	      }
	    
	if(deltaMw1 > maxDeltaMw_ || deltaMw2 > maxDeltaMw_ || deltaMt > maxDeltaMt_ )
	  {
	    status_ = false;
	    error_ = "Failed to converge";
	  }
      }

    return curMom;
	
  }
      
      
  /**
   */
  std::vector<TVector3> nSolve2(TLorentzVector pl1, TLorentzVector pb1, TLorentzVector pn1,
				TLorentzVector pl2, TLorentzVector pb2, TLorentzVector pn2,
				double mw)
  {
    std::vector<TVector3> theSols;
	  
    //1st decay
    TVector3 nu1 = pn1.Vect();
    TVector3 w1  = nu1+pl1.Vect();
    TVector3 t1  = w1+pb1.Vect();
    double el1   = pl1.E();
    double eb1   = pb1.E();
	  
    //2nd decay
    TVector3 nu2 = pn2.Vect();
    TVector3 w2  = nu2+pl2.Vect();
    TVector3 t2  = w2+pb2.Vect();
    double el2   = pl2.E();
    double eb2   = pb2.E();
	  
    //scan in a rectangular grid
    double delta=10;
    for(int ix=-10; ix<=10; ix++)
      {
	for(int iy=-10; iy<=10; iy++)
	  {
	    for(int iz=-10; iz<=10; iz++)
	      {
		      
		//initial guess
		TVector3 px(ix,iy,iz);
		px *= delta;
		      
		int nIter=0;
		bool hasConverged=false;
		status2_=true;
		do{
			
		  nIter++;
		  TVector3 gVec = g(t1,w1,nu1,el1,eb1,t2,w2,nu2,el2,eb2,px,mw);
		  TMatrixD jacobian = dgdx(t1,w1,nu1,eb1,t2,w2,nu2,eb2,px,mw);
			
		  if(status2_)
		    {
		      TMatrixD mtr = inverse( jacobian, status2_ );
			    
		      if(status2_)
			{
			  TVector3 deltaX(0,0,0);
			  deltaX.SetX( -mtr(0,0)*gVec.X()-mtr(0,1)*gVec.Y()-mtr(0,2)*gVec.Z() );
			  deltaX.SetY( -mtr(1,0)*gVec.X()-mtr(1,1)*gVec.Y()-mtr(1,2)*gVec.Z() );
			  deltaX.SetZ( -mtr(2,0)*gVec.X()-mtr(2,1)*gVec.Y()-mtr(2,2)*gVec.Z() );
				
			  TVector3 nextX = px + deltaX;
				
			  if(deltaX.Mag() < eVec_.Mag() ) hasConverged=true;
			}
		    }
		}while(!hasConverged && nIter <maxIter_ && status2_);
		      
		if(nIter>maxIter_) hasConverged=false;
		if(!hasConverged || !status2_) continue;
		      
		//don't care about redundant solutions
		if( px.Mag() < eVec_.Mag() ) continue;
		bool isRedundant(false);
		for(std::vector<TVector3>::iterator it = theSols.begin();
		    it != theSols.end();
		    it++)
		  {
		    TVector3 diff = *it - px;
		    if(diff.Mag() < eVec_.Mag()) { isRedundant=true; break; }
		  }
		if(isRedundant) continue;
		      
		theSols.push_back( px );			     
	      }
	  }
      }
	  
    return theSols;
  }

      
  /**
     @short the kinematics based on the invariant masses
  */
  TVector3 f(TLorentzVector &pl1, TLorentzVector &pb1, 
	     TLorentzVector &pl2, TLorentzVector &pb2, 
	     TLorentzVector &pnu, double &mw)
  {
    TVector3 fVec(0,0,0);
	
    TLorentzVector wRec1 = pl1+pnu;
    fVec.SetX( wRec1.M2()-pow(mw,2) );

    TLorentzVector wRec2 = pl2+constrainedMom(pnu);
    fVec.SetY( wRec2.M2()-pow(mw,2) );

    TLorentzVector topRec = wRec1+pb1;
    TLorentzVector antiTopRec = wRec2+pb2;
    fVec.SetZ( topRec.M2() - antiTopRec.M2() );

    return fVec;
  }

  /**
     @short the kinematics based on the energies
  */
  TVector3 g(TVector3 &pt1, TVector3 &pw1, TVector3 &pn1, double el1, double eb1,
	     TVector3 &pt2, TVector3 &pw2, TVector3 &pn2, double el2, double eb2,
	     TVector3 &px,  double mw)
  {
    TVector3 gVec(0,0,0);
	
    //1st decay
    TVector3 nu1  = pn1+px;
    TVector3 w1   = pw1+px;
    TVector3 top1 = pt1+px;

    //2nd decay
    TVector3 nu2  = pn2-px;
    TVector3 w2   = pw2-px;
    TVector3 top2 = pt2-px;

    gVec.SetX( sqrt(pow(mw,2)+w1.Mag2()) -el1 -nu1.Mag() );
    gVec.SetY( sqrt(pow(mw,2)+w2.Mag2()) -el2 -nu2.Mag() );
    gVec.SetZ( sqrt(pow( sqrt(pow(mw,2)+w1.Mag2()) +eb1,2) - top1.Mag2())
	       -sqrt(pow( sqrt(pow(mw,2)+w2.Mag2()) +eb2,2) - top2.Mag2()) );
	
    return gVec;
  }
      
      
  /**
     @short the jacobian matrix for f(pnu)
  */
  TMatrixD dfdp(TLorentzVector &pl1, TLorentzVector &pb1, 
		TLorentzVector &pl2, TLorentzVector &pb2, 
		TLorentzVector &pnu)
  {
    TMatrixD dfdpMtr(3,3); dfdpMtr.Zero();
	
    //check
    TLorentzVector pnu2 = constrainedMom(pnu);
    if(pnu.E()==0 || pnu2.E()==0)
      {
	status_ =false;
	error_  = "Null energy found";
	return dfdpMtr;
      }

    //1st decay
    double alpha1 = 2*(pl1.E()+pnu.E())/pnu.E();
    double beta1  = 2*(pl1.E()+pnu.E()+pb1.E())/pnu.E();
    TVector3 w1ThreeMom = pl1.Vect()+pnu.Vect();
    TVector3 top1ThreeMom = w1ThreeMom + pb1.Vect();

    //2nd decay	
    double alpha2 = 2*(pl2.E()+pnu2.E())/pnu2.E();
    double beta2  = 2*(pl2.E()+pnu2.E()+pb2.E())/pnu2.E();
    TVector3 w2ThreeMom = pl2.Vect()+pnu2.Vect();
    TVector3 top2ThreeMom = w2ThreeMom + pb2.Vect();

    //fill matrix elements
    dfdpMtr(0,0) =  alpha1*pnu.X() - 2*w1ThreeMom.X();
    dfdpMtr(0,1) =  alpha1*pnu.Y() - 2*w1ThreeMom.Y();
    dfdpMtr(0,2) =  alpha1*pnu.Z() - 2*w1ThreeMom.Z();

    dfdpMtr(1,0) =  -alpha2*pnu2.X() + 2*w2ThreeMom.X();
    dfdpMtr(1,1) =  -alpha2*pnu2.Y() + 2*w2ThreeMom.Y();
    dfdpMtr(1,2) =  -alpha2*pnu2.Z() + 2*w2ThreeMom.Z();

    dfdpMtr(2,0) =  beta1*pnu.X() + beta2*pnu2.X() - 2*top1ThreeMom.X() - 2*top2ThreeMom.X();
    dfdpMtr(2,1) =  beta1*pnu.Y() + beta2*pnu2.Y() - 2*top1ThreeMom.Y() - 2*top2ThreeMom.Y();
    dfdpMtr(2,2) =  beta1*pnu.Z() + beta2*pnu2.Z() - 2*top1ThreeMom.Z() - 2*top2ThreeMom.Z();
	
    return dfdpMtr;
  }

  /**
     @short the jacobian matrix for g(X)
  */
  TMatrixD dgdx(TVector3 &pt1, TVector3 &pw1, TVector3 &pn1,double eb1,
		TVector3 &pt2, TVector3 &pw2, TVector3 &pn2,double eb2,
		TVector3 &px,  double mw)
  {
    TMatrixD dgdxMtr(3,3); dgdxMtr.Zero();

	
    //1st decay
    TVector3 nu1  = pn1+px;
    TVector3 w1   = pw1+px;
    TVector3 top1 = pt1+px;
    double alpha1 = sqrt( pow(mw,2)+w1.Mag2() );
    double beta1  = alpha1+eb1;
    double gamma1 = sqrt( pow(beta1,2)-top1.Mag2());
	
    //2nd decay
    TVector3 nu2  = pn2-px;
    TVector3 w2   = pw2-px;
    TVector3 top2 = pt2-px;
    double alpha2 = sqrt( pow(mw,2)+w2.Mag2() );
    double beta2  = alpha2+eb2;
    double gamma2 = sqrt( pow(beta2,2)-top2.Mag2());
	
    //check
    if(alpha1==0 || nu1.Mag()==0 || gamma1==0
       || alpha2==0 || nu2.Mag()==0 || gamma2==0)
      {
	status2_ = false;
	return dgdxMtr;
      }

	
    //fill matrix elements
    dgdxMtr(0,0) =  w1.X()/alpha1 - nu1.X()/nu1.Mag();
    dgdxMtr(0,1) =  w1.Y()/alpha1 - nu1.Y()/nu1.Mag();
    dgdxMtr(0,2) =  w1.Z()/alpha1 - nu1.Z()/nu1.Mag();

    dgdxMtr(1,0) =  -w2.X()/alpha2 + nu2.X()/nu2.Mag();
    dgdxMtr(1,1) =  -w2.Y()/alpha2 + nu2.Y()/nu2.Mag();
    dgdxMtr(1,2) =  -w2.Z()/alpha2 + nu2.Z()/nu2.Mag();

    dgdxMtr(2,0) = ((beta1/alpha1)*w1.X()-top1.X())/beta1-(-(beta2/alpha2)*w2.X()+top2.X())/beta2;
    dgdxMtr(2,1) = ((beta1/alpha1)*w1.Y()-top1.Y())/beta1-(-(beta2/alpha2)*w2.Y()+top2.Y())/beta2;
    dgdxMtr(2,2) = ((beta1/alpha1)*w1.Z()-top1.Z())/beta1-(-(beta2/alpha2)*w2.Z()+top2.Z())/beta2;
	
    return dgdxMtr;
  }

      
  /**
     @short 3x3 matrix inversion
  */
  TMatrixD inverse(TMatrixD &a, bool &theStatus)
  {
    TMatrixD aInv(3,3); aInv.Zero();
	
    //compute the determinant
    double det=a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))
      -a(0,1)*(a(1,0)*a(2,2)-a(1,2)*a(2,0))
      +a(0,2)*(a(1,0)*a(2,1)-a(1,1)*a(2,0));
    if(det==0) 
      {
	theStatus =false;
	return aInv;
      }
	
	
    //invert the matrix
    aInv(0,0) = a(1,1)*a(2,2)-a(1,2)*a(2,1);
    aInv(0,1) = a(0,2)*a(2,1)-a(0,1)*a(2,2);
    aInv(0,2) = a(0,1)*a(1,2)-a(0,2)*a(1,1);
	
    aInv(1,0) = a(1,2)*a(2,0)-a(1,0)*a(2,2);
    aInv(1,1) = a(0,0)*a(2,2)-a(0,2)*a(2,0);
    aInv(1,2) = a(0,2)*a(1,0)-a(0,0)*a(1,2);
	
    aInv(2,0) = a(1,0)*a(2,1)-a(1,1)*a(2,0);
    aInv(2,1) = a(0,1)*a(2,0)-a(0,0)*a(2,1);
    aInv(2,2) = a(0,0)*a(1,1)-a(0,1)*a(1,0);

    //return the result
    aInv *= 1./det;
    return aInv;
	
  }
      
  /**
     @short returns the constrained of the second neutrino
  */
  TLorentzVector constrainedMom(TLorentzVector &mom)
  {
    TVector3 threeMomPrime = constraint_ - mom.Vect();
    TLorentzVector momPrime( threeMomPrime, threeMomPrime.Mag() );
    return momPrime;
  }

  /**
     @short DTOR
  */
  ~TopKinSolver()
  {
  }

private:

  //the status
  bool status_,status2_;
      
  //the error
  std::string error_;
      
  //convergence parameters
  TVector3 eVec_;
  double maxDeltaMw_,maxDeltaMt_;
      
  //the number of iterations
  int maxIter_;
      
  //the constraint
  TVector3 constraint_;
      
};

#endif

