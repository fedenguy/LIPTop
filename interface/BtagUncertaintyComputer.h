#ifndef btagunccomp_h
#define btagunccomp_h
//from A. Rizzi
namespace btag
{
  class Weight 
  {
  public:
    struct JetInfo {
      JetInfo(float mceff,float datasf) : eff(mceff), sf(datasf) {}
      float eff;
      float sf;
    };
    
    Weight(int jmin, int jmax) : maxTags(jmax), minTags(jmin) {}
      
      bool filter(int t) { return (t >= minTags && t <= maxTags); }
      float weight(vector<JetInfo> jets, int tags);
  
  private:
      int maxTags;
      int minTags;
  };
  
  float BTagWeight::weight(vector<JetInfo> jets, int tags)
  {
    if(!filter(tags)) return 0;
    int njets=jets.size();
    int comb= 1 << njets;
    float pMC=0;
    float pData=0;
    for(int i=0;i < comb; i++)
      {
	float mc=1.;
	float data=1.;
	int ntagged=0;
	for(int j=0;j<njets;j++)
	  {
	    bool tagged = ((i >> j) & 0x1) == 1;
	    if(tagged) 
	      {
		ntagged++;
		mc*=jets[j].eff;
		data*=jets[j].eff*jets[j].sf;
	      }
	    else
	     {
	       mc*=(1.-jets[j].eff);
	       data*=(1.-jets[j].eff*jets[j].sf);
	     }
	  }       
	
	if(filter(ntagged))
	  {
	    std::cout << mc << " " << data << endl;
	    pMC+=mc;
	    pData+=data;
	  }
      }
    
    if(pMC==0) return 0; 
    return pData/pMC;
  }
}

#endif
