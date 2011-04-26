#include "LIP/Top/interface/TopKinSolver.h"

/**
   @short can be used to sort a collection of solutions by increasing M(TTbar)
*/
bool mTTbarOrder( TTbarSolution_t a, TTbarSolution_t b)
{
  TLorentzVector ttbar_a = a.pt1+a.pt2;
  TLorentzVector ttbar_b = b.pt1+b.pt2;
  return ttbar_a.M() > ttbar_b.M();
}

/**
   @short computes the transverse mass of a top decay
   @param pbl the sum of the b-jet + lepton system
   @param pn the neutrino
   @return the transverse mass
*/
double getMT(TLorentzVector &pbl, TLorentzVector &pn)
{
  double pTbl_pTn=pbl.Px()*pn.Px()+pbl.Py()*pn.Py();
  double mt = pow(pbl.M(),2)+2*(pbl.Et()*pn.Et()-pTbl_pTn);
  if( mt < 0) return -TMath::Sqrt(-mt);
  return TMath::Sqrt(mt);
}

/**
   @short returns the minimum transverse mass of a ttbar pair
   @return the maximum mT       
*/
std::vector<double> getMT2( TTbarSolution_t &sol)
{
  TLorentzVector pbl1 = sol.pl1+sol.pb1;
  TLorentzVector pn1 = sol.pn1;
  double mt21 = getMT(pbl1,pn1);
  
  TLorentzVector pbl2 = sol.pl2+sol.pb2;
  TLorentzVector pn2 =sol.pn2;
  double mt22 = getMT(pbl2,pn2);
  
  std::vector<double> mt2sort(2,0);
  mt2sort[0]=(fabs(mt21) < fabs(mt22) ? mt22 : mt21);
  mt2sort[1]=(fabs(mt21) < fabs(mt22) ? mt21 : mt22);
  return mt2sort;
}
