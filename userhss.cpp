#include "userhss.h"
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Variables used by all model objects. Hence only one copy is maintained for all objects
static const double d2d3 = 2.0 / 3.0;
static const double dPi  = 3.141592653589793238462643383279502884197169399;
static const double dDegRad = dPi / 180.0;
// static const double drr = 0.385;
static const double sqrt3 = 1.732050808;

// Plasticity Indicators
static const unsigned long mShearNow    = 0x01;   /* state logic */
static const unsigned long mTensionNow  = 0x02;
static const unsigned long mShearPast   = 0x04;
static const unsigned long mTensionPast = 0x08;
static const unsigned long mVolumeNow   = 0x100;
static const unsigned long mVolumePast  = 0x200;

// One Static instance isneccessary as a part of internal registration process of the model with FLAC/FLAC3D
static UserhssModel userhssModel(true);

UserhssModel::UserhssModel(bool bRegister)
          :ConstitutiveModel(mnUserhssModel,bRegister), dCohesion(0.0), dFriction(0.0), dDilatation(0.0), dM_index(0.0), dP_ref(0.0), dE50_ref(0.0),dEur_ref(0.0), dEoed_ref(0.0), dRatiofail(0.0), dPoisson(0.0), dK0(0.0), dG0_ref(0.0), dR07(0.0),
	dBulk(0.0), dShear(0.0), dvSH(0.0), dvVH(0.0), dPc(0.0), sh11(0.0), sh12(0.0), sh13(0.0), sh22(0.0), sh23(0.0), sh33(0.0), gama_hist(0.0),
	dFcos(0.0), dFsin(0.0), dCFcot(0.0), dFsinCV(0.0), dDsin(0.0), dDsinMob(0.0), dEur(0.0), dE50(0.0), dQf(0.0), dQa(0.0), dQ(0.0), dQ_(0.0), dP(0.0), alpha(0.0), dMM2(0.0), dShearaaa(0.0), dFsinMob(0.0), dVolumeaaa(0.0), dG0(0.0), dShearMin(0.0){  }

	  
const char **UserhssModel::Properties(void) const {
  static const char *strKey[] = {
    "cohesion","friction","dilatation","m_index","p_ref","e50_ref","eur_ref","eoed_ref","ratiofail","poisson","k0","g0_ref","r07",
	"bulk_cur","shear_cur","vsh","vvh","pc","sh11","sh12","sh13","sh22","sh23","sh33","r_hist"
	"qf","ashear","sinmob","avolume",0
  };
  return(strKey);
}


const char **UserhssModel::States(void) const {
  static const char *strKey[] = {
    "shear-n","tension-n","shear-p""tension-p","volume-n","volume-p",0
  };
  return(strKey);
}


/*   *Note:Maintain order of property input/output
*/
double UserhssModel::GetProperty(unsigned ul) const {
  switch (ul) {
    case 1:  return(dCohesion);
    case 2:  return(dFriction);
    case 3:  return(dDilatation);
    case 4:  return(dM_index);
    case 5:  return(dP_ref);
    case 6:  return(dE50_ref);
	case 7:  return(dEur_ref);
	case 8:  return(dEoed_ref);
	case 9:  return(dRatiofail);
	case 10: return(dPoisson);
	case 11: return(dK0);
	case 12:return(dG0_ref);
	case 13:return(dR07);

	case 14: return(dBulk);
	case 15: return(dShear);
	case 16: return(dvSH);
	case 17: return(dvVH);
	case 18: return(dPc);
	case 19: return(sh11);
	case 20: return(sh12);
	case 21: return(sh13);
	case 22: return(sh22);
	case 23: return(sh23);
	case 24: return(sh33);
	case 25: return(gama_hist);

    case 26: return(dQf);
	case 27: return(dShearaaa);
	case 28: return(dFsinMob);
	case 29: return(dVolumeaaa);
  }
  return(0.0);
}


void UserhssModel::SetProperty(unsigned ul,const double &dVal) {
  switch (ul) {
    case 1:  dCohesion   = dVal;  break;
	case 2:  dFriction   = dVal;  break;
	case 3:  dDilatation = dVal;  break;
	case 4:  dM_index    = dVal;  break;
	case 5:  dP_ref      = dVal;  break;
	case 6:  dE50_ref    = dVal;  break;
	case 7:  dEur_ref    = dVal;  break;
	case 8:  dEoed_ref   = dVal;  break;
	case 9:  dRatiofail  = dVal;  break;
	case 10: dPoisson    = dVal;  break;
	case 11: dK0         = dVal;  break;
	case 12: dG0_ref     = dVal;  break;
	case 13: dR07        = dVal;  break;

    case 14: dBulk       = dVal;  break;
	case 15: dShear      = dVal;  break;
	case 16: dvSH        = dVal;  break;
	case 17: dvVH        = dVal;  break;
	case 18: dPc         = dVal;  break;
	case 19: sh11        = dVal;  break;
	case 20: sh12        = dVal;  break;
	case 21: sh13        = dVal;  break;
	case 22: sh22        = dVal;  break;
	case 23: sh23        = dVal;  break;
	case 24: sh33        = dVal;  break;
	case 25: gama_hist   = dVal;  break;
	case 26: dQf = dVal; break;
	case 27: dShearaaa = dVal; break;
	case 28: dFsinMob = dVal; break;
	case 29: dVolumeaaa = dVal; break;
  }
}


const char *UserhssModel::Copy(const ConstitutiveModel *cm) {
  //Detects type mismatch error and returns error string. otherwise returns 0
	const char *str = ConstitutiveModel::Copy(cm);
	if (str) return(str);
	UserhssModel *mm = (UserhssModel *)cm;
	dCohesion   = mm->dCohesion;
	dFriction   = mm->dFriction;
	dDilatation = mm->dDilatation;
	dM_index    = mm->dM_index;
	dP_ref      = mm->dP_ref;
	dE50_ref    = mm->dE50_ref;
	dEur_ref    = mm->dEur_ref;
	dEoed_ref   = mm->dEoed_ref;
	dRatiofail  = mm->dRatiofail;
	dPoisson    = mm->dPoisson;
	dK0         = mm->dK0;
	dG0_ref     = mm->dG0_ref;
	dR07        = mm->dR07;

	dBulk       = mm->dBulk;
	dShear      = mm->dShear;
	dvSH        = mm->dvSH;
	dvVH        = mm->dvVH;
	dPc         = mm->dPc; 
	sh11        = mm->sh11;
	sh12        = mm->sh12;
	sh13        = mm->sh13;
	sh22        = mm->sh22;
	sh23        = mm->sh23;
	sh33        = mm->sh33;
	gama_hist   = mm->gama_hist;

	dQf = mm->dQf;
	dShearaaa = mm->dShearaaa;
	dFsinMob = mm->dFsinMob;
	dVolumeaaa = mm->dVolumeaaa;
	return(0);
}


const char *UserhssModel::Initialize(unsigned uDim,State *ps) {
	if ((uDim != 2) && (uDim != 3)) return("Illegal dimension in userhs constitutive model");
	dFcos = cos(dFriction*dDegRad);
	dFsin = sin(dFriction*dDegRad);
	dDsin = sin(dDilatation*dDegRad);
	double dFcot = (dFsin < 0.01) ? 100 : dFcos / dFsin;
	dCFcot = dCohesion*dFcot;
	dFsinCV = (dFsin - dDsin) / (1 - dFsin*dDsin);
	/* --- relating to stress --- */
	Axes aDir;
	double dPrinMin, dPrinMid, dPrinMax, sdif = 0.0, psdif = 0.0;
	int icase = 0;
	bool bFast = ps->stnS.Resoltopris(&dPrinMin, &dPrinMid, &dPrinMax, &aDir, uDim, &icase, &sdif, &psdif);
	double dS3val = -dPrinMax;
	double dS3toPref = (dCFcot + dS3val) / (dCFcot + dP_ref);
	if (dS3toPref < 0.01)dS3toPref = 0.01;// to make sure the mudulus is not too small
	dQ = dPrinMax - dPrinMin;
	alpha = (3.0 + dFsin) / ((3.0 + dFsin));
	dQ_ = -(dPrinMin + (alpha - 1.0)*dPrinMid - alpha*dPrinMax);
	dP = -(dPrinMin + dPrinMid + dPrinMax) / 3.0;
	double dKs = dEur_ref / (3.0 - 6.0*dPoisson);
	double dKc = dEoed_ref*(1.0 + 2.0*dK0) / 3.0;
	H = (dKs / dKc>1.1) ? dKc*dKs / (dKs - dKc) : 11.0*dKc;
	dEur = dEur_ref*pow(dS3toPref, dM_index);
	//dShear = dEur / (2.0 + 2.0*dPoisson);
	//dBulk = dEur / (3.0 - 6.0*dPoisson);
	dE50 = dE50_ref*pow(dS3toPref, dM_index);
	dQf = (dCFcot + dS3val)*2.0*dFsin / (1.0 - dFsin);
	dQa = dQf / dRatiofail;
	dFsinMob = dQ / (-dPrinMax - dPrinMin + 2.0*dCFcot);
	dDsinMob = (dFsinMob - dFsinCV) / (1.0 - dFsinMob*dFsinCV);
	if (dDsinMob < 0.0)dDsinMob = 0.0;
	double MM = 6.0*dFsin / (3.0 - dFsin);
	dMM2 = MM*MM;

	dG0 = dG0_ref*pow(dS3toPref, dM_index);
	SHist(0, 0) = sh11;
	SHist(0, 1) = sh12;
	SHist(0, 2) = sh13;
	SHist(1, 0) = sh12;
	SHist(1, 1) = sh22;
	SHist(1, 2) = sh23;
	SHist(2, 0) = sh13;
	SHist(2, 1) = sh23;
	SHist(2, 2) = sh33;
	dShear = dG0;
	dBulk = dShear*2.0*(1.0 + dPoisson) / (3.0 - 6.0*dPoisson);
	return(0);
}


const char *UserhssModel::Run(unsigned uDim, State *ps) {
	if ((uDim != 2) && (uDim != 3)) return("Illegal dimension in userhs constitutive model");
	if (ps->dHystDampMult > 0.0) HDampInit(ps->dHystDampMult);
	static double dGama, dEpsilon;
	/* ---plasticity indicator:
	store 'now' info. as 'past' and turn 'now' info off ---*/
	if (ps->mState & mShearNow)
		ps->mState = (unsigned long)(ps->mState | mShearPast);
	ps->mState = (unsigned long)(ps->mState & ~mShearNow);
	if (ps->mState & mTensionNow)
		ps->mState = (unsigned long)(ps->mState | mTensionPast);
	ps->mState = (unsigned long)(ps->mState & ~mTensionNow);
	if (ps->mState & mVolumeNow)
		ps->mState = (unsigned long)(ps->mState | mVolumePast);
	ps->mState = (unsigned long)(ps->mState & ~mVolumeNow);
	int iPlas_T = 0;
	int iPlas_S = 0;
	int iPlas_V = 0;
	int iPlas_MC = 0;
	/* --- hardening:
	initialize stacks to calculate hardening parameters for zone --- */
	if (!ps->bySubZone) {
		dGama = 0.0;
		dEpsilon = 0.0;
	}


	/* --- calculate the mudulus --- */
	Eigen::Matrix3d ss;// actual deviatoric strain tensor rate
	ss(0, 0) = ps->stnE.d11 - (ps->stnE.d11 + ps->stnE.d22 + ps->stnE.d33) / 3.0;
	ss(1, 1) = ps->stnE.d22 - (ps->stnE.d11 + ps->stnE.d22 + ps->stnE.d33) / 3.0;
	ss(2, 2) = ps->stnE.d33 - (ps->stnE.d11 + ps->stnE.d22 + ps->stnE.d33) / 3.0;
	ss(0, 1) = ps->stnE.d12;
	ss(0, 2) = ps->stnE.d13;
	ss(1, 2) = ps->stnE.d23;
	ss(1, 0) = ps->stnE.d12;
	ss(2, 0) = ps->stnE.d13;
	ss(2, 1) = ps->stnE.d23;
	/* --- calculate the orthogonal eigenvetors and eigenvalues of deviatoric strain tensor --- */
	double gama_hist1;
	if (ss.norm()>0.0) /* --- when ss != 0, update gama_hist --- */
	{
		Eigen::EigenSolver<Eigen::Matrix3d> es(ss);
		Eigen::Matrix3d D = es.pseudoEigenvalueMatrix();
		Eigen::Matrix3d V = es.pseudoEigenvectors();
		/* --- transform the strain history tensor to current strain space --- */
		SHist = V.transpose()*SHist*V;
		Eigen::Matrix3d eH = ss*SHist;
		gama_hist1 = sqrt3*eH.norm() / ss.norm();
		/* --- calculate the diagonal transformtion matrix T --- */
		Eigen::Matrix3d T;
		T = Eigen::Matrix3d::Identity();
		if (SHist(0, 0) + 1.0 > 0.0) {
			double T11 = sqrt(SHist(0, 0) + 1.0);
			double u11 = (SHist(0, 0)*D(0, 0) < 0) ? 0.0 : 1.0;
			T(0, 0) = (1 + u11*(T11 - 1.0)) / T11;
		}
		if (SHist(1, 1) + 1.0 > 0.0) {
			double T22 = sqrt(SHist(1, 1) + 1.0);
			double u22 = (SHist(1, 1)*D(1, 1)<0) ? 0.0 : 1.0;
			T(1, 1) = (1 + u22*(T22 - 1.0)) / T22;
		}
		if (SHist(2, 2) + 1.0 > 0.0) {
			double T33 = sqrt(SHist(2, 2) + 1.0);
			double u33 = (SHist(2, 2)*D(2, 2) < 0) ? 0.0 : 1.0;
			T(2, 2) = (1 + u33*(T33 - 1.0)) / T33;
		}
		/* --- update the strain history H --- */
		SHist = T*(SHist + Eigen::Matrix3d::Identity())*T + ss - Eigen::Matrix3d::Identity();
		/* --- update gama_hist --- */
		eH = ss*SHist;
		gama_hist = sqrt3*eH.norm() / ss.norm();
		dShear = dG0 / (gama_hist - gama_hist1)*(gama_hist / (1.0 + sqrt3*gama_hist / 2.0 / dR07) - gama_hist1 / (1.0 + sqrt3*gama_hist1 / 2.0 / dR07));
		dBulk = dShear*2.0*(1.0 + dPoisson) / (3.0 - 6.0*dPoisson);
	}
	


	double dE1 = dBulk + d4d3*dShear;
	double dE2 = dBulk - d2d3*dShear;
	double dG2 = 2.0*dShear;
	/* --- Trial elastic stresses --- */
	double dE11 = ps->stnE.d11;
	double dE22 = ps->stnE.d22;
	double dE33 = ps->stnE.d33;
	ps->stnS.d11 += dE11 * dE1 + (dE22 + dE33) * dE2;
	ps->stnS.d22 += dE22 * dE1 + (dE11 + dE33) * dE2;
	ps->stnS.d33 += dE33 * dE1 + (dE11 + dE22) * dE2;
	ps->stnS.d12 += ps->stnE.d12 * dG2;
	ps->stnS.d13 += ps->stnE.d13 * dG2;
	ps->stnS.d23 += ps->stnE.d23 * dG2;
    Axes aDir;
	double dPrinMin, dPrinMid, dPrinMax, sdif = 0.0, psdif = 0.0;
	int icase = 0;
	bool bFast = ps->stnS.Resoltopris(&dPrinMin, &dPrinMid, &dPrinMax, &aDir, uDim, &icase, &sdif, &psdif);
	double dPrinMinCopy = dPrinMin;
	double dPrinMidCopy = dPrinMid;
	double dPrinMaxCopy = dPrinMax;

	/* --- stress state on elastic hypothesis --- */
	double dQI = dPrinMax - dPrinMin;
	double dQI_ = -(dPrinMin + (alpha - 1.0)*dPrinMid - alpha*dPrinMax);
	double dPI = -(dPrinMin + dPrinMid + dPrinMax) / 3.0;
	/* --- tesion function --- */
	double dFTension = dPrinMax - dCFcot;
	/* --- MC function --- */
	double dFMC = dQf - dQI;
	/* --- shear function --- */
	double dFShear = dQa*dQI / (dE50*(dQa - dQI)) - 2.0*dQI / dEur - dvSH;
	/* --- volume function --- */
	double dFVolume = dQI_*dQI_ / dMM2 + dPI*dPI - dPc*dPc;
    double dMClamda = 0.0;
	double dSlamda = 0.0;
	double dVlamda = 0.0;

	/* --- stress corection --- */
	if (dFTension > 0.0) {/* --- tension corection --- */
		iPlas_T = 1;
		ps->mState = (unsigned long)(ps->mState | mTensionNow);
		double dTco = dE2*dFTension / dE1;
		dPrinMin -= dTco;
		dPrinMid -= dTco;
		dPrinMax -= dFTension;
	}
	else {
		if (dFMC <= 0.0 && dFVolume > 0.0) {/* --- MC yield && volume yield --- */
			iPlas_MC = 1;
			iPlas_V = 1;
			ps->mState |= (unsigned long)(ps->mState | mShearNow | mVolumeNow);
			double dNDil = (1.0 + dDsinMob) / (1.0 - dDsinMob);
			double dA1 = dE1 + dE2*(-dNDil);
			double dA2 = dE2*(1 - dNDil);
			double dA3 = dE1*(-dNDil) + dE2;
			double dC1 = -dE1*(2.0*dQ_ / dMM2 + d2d3*dP) - dE2*(-2.0*dQ_ / dMM2 + 2.0* d2d3*dP);
			double dC2 = -dE1*(2.0*dQ_*(alpha - 1.0) / dMM2 + d2d3*dP) - dE2*(2.0*dQ_*(1.0 - alpha) / dMM2 + 2.0* d2d3*dP);
			double dC3 = -dE1*(-2.0*dQ_*alpha / dMM2 + d2d3*dP) - dE2*(2.0*dQ_*alpha / dMM2 + 2.0* d2d3*dP);
			double dAAA = dA1 + (alpha - 1.0)*dA2 - alpha*dA3;
			double daaa = (dA1 + dA2 + dA3) / 3.0;
			double dCCC = dC1 + (alpha - 1.0)*dC2 - alpha*dC3;
			double dccc = (dC1 + dC2 + dC3) / 3.0;
			double num_mc = dFMC*(2.0*dQI_*dCCC / dMM2 + 2.0*dPI*dccc) + dFVolume*(dC1 - dC3);
			double num_v = -(dA1 - dA3)*dFVolume - dFMC*(2.0*dQI_*dAAA / dMM2 + 2.0*dPI*daaa);
			double den = (dA1 - dA3)*(2.0*dQI_*dCCC / dMM2 + 2.0*dPI*dccc) - (2.0*dQI_*dAAA / dMM2 + 2.0*dPI*daaa)*(dC1 - dC3);
			dMClamda = num_mc / den;
			dVlamda = num_v / den;
			dPrinMin -= (dMClamda*dA1 + dVlamda*dC1);
			dPrinMid -= (dMClamda*dA2 + dVlamda*dC2);
			dPrinMax -= (dMClamda*dA3 + dVlamda*dC3);
		}
		if (dFMC <= 0.0 && dFVolume <= 0.0) {/* --- MC yield --- */
			iPlas_MC = 1;
			ps->mState |= (unsigned long)(ps->mState | mShearNow);
			double dNDil = (1.0 + dDsinMob) / (1.0 - dDsinMob);
			double dA1 = dE1 + dE2*(-dNDil);
			double dA2 = dE2*(1 - dNDil);
			double dA3 = dE1*(-dNDil) + dE2;
			double num_mc = dFMC;
			double den = dA1 - dA3;
			dMClamda = num_mc / den;
			dPrinMin -= dMClamda*dA1;
			dPrinMid -= dMClamda*dA2;
			dPrinMax -= dMClamda*dA3;
		}
		if (dFMC > 0.0 && dFShear > 0.0 && dFVolume > 0.0) {/* --- shear yield && volume yield --- */
			iPlas_S = 1;
			iPlas_V = 1;
			ps->mState |= (unsigned long)(ps->mState | mShearNow | mVolumeNow);
			double dB1 = dE1*(-0.5 + 0.5*dDsinMob) + dE2*(0.5 + 0.5*dDsinMob);
			double dB2 = dE2*dDsinMob;
			double dB3 = dE1*(0.5 + 0.5*dDsinMob) + dE2*(-0.5 + 0.5*dDsinMob);
			double dC1 = -dE1*(2.0*dQ_ / dMM2 + d2d3*dP) - dE2*(-2.0*dQ_ / dMM2 + 2.0* d2d3*dP);
			double dC2 = -dE1*(2.0*dQ_*(alpha - 1.0) / dMM2 + d2d3*dP) - dE2*(2.0*dQ_*(1.0 - alpha) / dMM2 + 2.0* d2d3*dP);
			double dC3 = -dE1*(-2.0*dQ_*alpha / dMM2 + d2d3*dP) - dE2*(2.0*dQ_*alpha / dMM2 + 2.0* d2d3*dP);
			double dBBB = dB1 + (alpha - 1.0)*dB2 - alpha*dB3;
			double dbbb = (dB1 + dB2 + dB3) / 3.0;
			double dCCC = dC1 + (alpha - 1.0)*dC2 - alpha*dC3;
			double dccc = (dC1 + dC2 + dC3) / 3.0;
			double num_s = -dFShear*(2.0*dQI_*dCCC / dMM2 + 2.0*dPI*dccc) / (dQa*dQa / (dQa - dQI) / (dQa - dQI) / dE50 - 2.0 / dEur) + dFVolume*(dC1 - dC3);
			double num_v = -(dB1 - dB3)*dFVolume + dFShear*(2.0*dQI_*dBBB / dMM2 + 2.0*dPI*dbbb) / ((dQa*dQa / (dQa - dQI) / (dQa - dQI) / dE50 - 2.0 / dEur));
			double den = (dB1 - dB3)*(2.0*dQI_*dCCC / dMM2 + 2.0*dPI*dccc) - (2.0*dQI_*dBBB / dMM2 + 2.0*dPI*dbbb)*(dC1 - dC3);
			dSlamda = num_s / den;
			dVlamda = num_v / den;
			dPrinMin -= (dSlamda*dB1 + dVlamda*dC1);
			dPrinMid -= (dSlamda*dB2 + dVlamda*dC2);
			dPrinMax -= (dSlamda*dB3 + dVlamda*dC3);
		}
		if (dFMC > 0.0 && dFShear > 0 && dFVolume <= 0.0) {/* --- shear yield --- */
			iPlas_S = 1;
			ps->mState |= (unsigned long)(ps->mState | mShearNow);
			double dB1 = dE1*(-0.5 + 0.5*dDsinMob) + dE2*(0.5 + 0.5*dDsinMob);
			double dB2 = dE2*dDsinMob;
			double dB3 = dE1*(0.5 + 0.5*dDsinMob) + dE2*(-0.5 + 0.5*dDsinMob);
			double num_s = -dFShear;
			double den = (dB1 - dB3)*(dQa*dQa / (dQa - dQI) / (dQa - dQI) / dE50 - 2.0 / dEur) - 1.0 + dDsinMob;
			dSlamda = num_s / den;
			dPrinMin -= dSlamda*dB1;
			dPrinMid -= dSlamda*dB2;
			dPrinMax -= dSlamda*dB3;
			double aaaaa = 0.0;
		}
		if (dFMC > 0.0 && dFShear <= 0.0 && dFVolume > 0.0) {/* --- volume yield --- */
			iPlas_V = 1;
			ps->mState |= (unsigned long)(ps->mState | mVolumeNow);
			double dC1 = -dE1*(2.0*dQ_ / dMM2 + d2d3*dP) - dE2*(-2.0*dQ_ / dMM2 + 2.0* d2d3*dP);
			double dC2 = -dE1*(2.0*dQ_*(alpha - 1.0) / dMM2 + d2d3*dP) - dE2*(2.0*dQ_*(1.0 - alpha) / dMM2 + 2.0* d2d3*dP);
			double dC3 = -dE1*(-2.0*dQ_*alpha / dMM2 + d2d3*dP) - dE2*(2.0*dQ_*alpha / dMM2 + 2.0* d2d3*dP);
			double dCCC = dC1 + (alpha - 1.0)*dC2 - alpha*dC3;
			double dccc = (dC1 + dC2 + dC3) / 3.0;
			double num_v = -dFVolume;
			double den = 2.0*dQI_*dCCC / dMM2 + 2.0*dPI*dccc;
			dVlamda = num_v / den;
			dPrinMin -= dVlamda*dC1;
			dPrinMid -= dVlamda*dC2;
			dPrinMax -= dVlamda*dC3;
		}
	}
	dQI = dPrinMax - dPrinMin;
	dQf = (dCFcot - dPrinMax)*2.0*dFsin / (1.0 - dFsin);
	dFMC = dQf - dQI;
	if (dFMC < 0.0 ) {// --- rejudge MC ---
		iPlas_MC = 1;
		ps->mState |= (unsigned long)(ps->mState | mShearNow);
		double dNDil = (1.0 + dDsinMob) / (1.0 - dDsinMob);
		double dA1 = dE1 + dE2*(-dNDil);
		double dA2 = dE2*(1 - dNDil);
		double dA3 = dE1*(-dNDil) + dE2;
		double num_mc = dFMC;
		double den = dA1 - dA3;
		dMClamda = num_mc / den;
		dPrinMin -= dMClamda*dA1;
		dPrinMid -= dMClamda*dA2;
		dPrinMax -= dMClamda*dA3;
	}
	
	/* --- update the variables related to hardening rule:
	                                                          --- */
	if (iPlas_T+ iPlas_S+ iPlas_V+ iPlas_MC) {
		ps->stnS.Resoltoglob(dPrinMin, dPrinMid, dPrinMax, aDir, dPrinMinCopy, dPrinMidCopy, dPrinMaxCopy, uDim, icase, sdif, psdif, bFast);
		/* --- hadRdening padRameter accumulation --- */
		if (iPlas_S == 1) {
			/* ---                      shear padRameter --- */
			double dGamaval = dSlamda*(1.0 - dDsinMob);
			dGama += dGamaval*ps->dSubZoneVolume;
		}
		if (iPlas_V == 1) {
			/* ---                      volume padRameter --- */
			double dEpsilonval = -2.0*dP*dVlamda;
			dEpsilon += dEpsilonval*ps->dSubZoneVolume;
		}
	}
	/* --- plastic padRameter incrementation and properties update --- */
	if (ps->bySubZone == ps->byTotSubZones - 1) {
		double dAux = 1.0 / ps->dZoneVolume;
		if (ps->byOverlay == 2) dAux *= 0.5;
		dvSH += dGama*dAux;// update shear plastic strain
		double ddvVH = dEpsilon*dAux;
		dvVH += dEpsilon*dAux;// update volume plastic strain
		double dPc_pref = (dCFcot - dPc) / (dCFcot + dP_ref);
		dPc += H*pow(dPc_pref, dM_index)*ddvVH;
		double aaaaaaa = 0.0;
		}
	/* --- update the parameters related to stresses in the initialize() function 
	                                                                                     ---*/
	double dS3val = -dPrinMax;
	double dS3toPref = (dCFcot + dS3val) / (dCFcot + dP_ref);
	if (dS3toPref < 0.2)dS3toPref = 0.2;// to make sure the mudulus is not too small
	dQ = dPrinMax - dPrinMin;

	//1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111//
	dShearaaa = fabs(dQa*dQ / (dE50*(dQa - dQ)) - 2.0*dQ / dEur - dvSH);
	if (iPlas_MC == 1) dShearaaa = dQf - dQ;
	double bbbb = 0.0;
	//1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111//

	dQ_ = -(dPrinMin + (alpha - 1.0)*dPrinMid - alpha*dPrinMax);
	dP = -(dPrinMin + dPrinMid + dPrinMax) / 3.0;

	//22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222//
	dVolumeaaa = dQ_*dQ_ / dMM2 + dP*dP - dPc*dPc;
	//22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222//

	dEur = dEur_ref*pow(dS3toPref, dM_index);
	dE50 = dE50_ref*pow(dS3toPref, dM_index);
	dQf = (dCFcot + dS3val)*2.0*dFsin / (1.0 - dFsin);
	dQa = dQf / dRatiofail;
	dFsinMob = dQ / (-dPrinMax - dPrinMin + 2.0*dCFcot);
	dDsinMob = (dFsinMob - dFsinCV) / (1.0 - dFsinMob*dFsinCV);
	if (dDsinMob < 0.0)dDsinMob = 0.0;
	return(0);
}
  
 const char *UserhssModel::SaveRestore(ModelSaveObject *mso) {
	 const char *str = ConstitutiveModel::SaveRestore(mso);
	 if (str) return(str);
	 mso->Initialize(28, 0);
	 mso->Save(0, dCohesion);
	 mso->Save(1, dFriction);
	 mso->Save(2, dDilatation);
	 mso->Save(3, dM_index);
	 mso->Save(4, dP_ref);
	 mso->Save(5, dE50_ref);
	 mso->Save(6, dEur_ref);
	 mso->Save(7, dEoed_ref);
	 mso->Save(8, dRatiofail);
	 mso->Save(9, dPoisson);
	 mso->Save(10, dK0);
	 mso->Save(11, dG0_ref);
	 mso->Save(12, dR07);

	 mso->Save(13, dBulk);
	 mso->Save(14, dShear);
	 mso->Save(15, dvSH);
	 mso->Save(16, dvVH);
	 mso->Save(17, dPc);
	 mso->Save(18, sh11);
	 mso->Save(19, sh12);
	 mso->Save(20, sh13);
	 mso->Save(21, sh22);
	 mso->Save(22, sh23);
	 mso->Save(23, sh33);
	 mso->Save(24, gama_hist);

	 mso->Save(25, dQf);
	 mso->Save(26, dShearaaa);
	 mso->Save(27, dFsinMob);
	 mso->Save(28, dVolumeaaa);
	 return(0);
}


void UserhssModel::HDampInit(const double dHMult)
{
	/* --- it means that the hystertic damping is out of consideration --- */
}