#ifndef __USERhss_H
#define __USERhss_H


#ifndef __CONMODEL_H
#include "conmodel.h"
#endif

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

class UserhssModel : public ConstitutiveModel {
  public:
    // static
    enum ModelNum { mnUserhssModel=324 };
    // Creators
    EXPORT UserhssModel(bool bRegister=true);
    // Accessors
    virtual const char *Keyword(void) const { return("userhss"); }
    virtual const char *Name(void) const { return("User-hs-small"); }
    virtual const char **Properties(void) const;
    virtual const char **States(void) const;
    virtual double GetProperty(unsigned ul) const;
    virtual ConstitutiveModel *Clone(void) const { return(new UserhssModel()); }
    virtual double ConfinedModulus(void) const { return(dBulk + d4d3*dShear); }
    virtual double ShearModulus(void) const { return(dShear); }
    virtual double BulkModulus(void) const { return(dBulk); }
    virtual double SafetyFactor(void) const { return(10.0); }
    virtual unsigned Version(void) const { return(2); }
	virtual bool SupportsHystDamp() const {return true;}

    // Manipulators
    virtual void SetProperty(unsigned ul,const double &dVal);
    virtual const char *Copy(const ConstitutiveModel *m);
    virtual const char *Initialize(unsigned uDim,State *ps);
    virtual const char *Run(unsigned uDim,State *ps);
	virtual const char *SaveRestore(ModelSaveObject *mso);
	virtual void HDampInit(const double dHMult);

  private:
	  double dCohesion, dFriction, dDilatation, dM_index, dP_ref, dE50_ref, dEur_ref, dEoed_ref, dRatiofail, dPoisson, dK0, dG0_ref, dR07;//input
	  double dBulk, dShear, dvSH, dvVH, dPc, sh11, sh12, sh13, sh22, sh23, sh33, gama_hist;//invoke
	  double dFcos, dFsin, dCFcot, dFsinCV, dDsin, dDsinMob, dEur, dE50, dQf, dQa, dQ, dQ_, dP, H, alpha, dMM2, dShearaaa, dFsinMob, dVolumeaaa, dG0,dShearMin;//globle
	  Eigen::Matrix3d SHist;
};



#endif
// EOF
