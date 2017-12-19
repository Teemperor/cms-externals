#ifndef __PARTICLE_DEF__
#define __PARTICLE_DEF__
/**
  Single particle class with basic methods such as boost, rotation,
  and angle calculation.
*/

#include <cstdio>
#include <iostream>
#include <math.h>

namespace TauSpinner {

const double ELECTRON_MASS = 0.0005111;
const double MUON_MASS     = 0.105659;
const double TAU_MASS      = 1.777;

class Particle
{
public:
/*
	Constructors
*/
  Particle():_px(0),_py(0),_pz(0),_e(0),_pdgid(0)                                                      {}
	Particle(double px, double py, double pz, double e, int id):_px(px),_py(py),_pz(pz),_e(e),_pdgid(id) {}
  
public:
/*
	Accessors
*/
	double px() const { return _px; }
	double py() const { return _py; }
	double pz() const { return _pz; }
	double e () const { return _e; }
	int pdgid() const { return _pdgid; }

	// Invariant mass. If mass is negative then -sqrt(-mass) is returned
	double recalculated_mass()  const
	{
		double p2 = _px*_px + _py*_py + _pz*_pz;
		double e2 = _e*_e;
		double m2 = e2 - p2;

		if ( m2 < 0.0 ) return -sqrt( -m2 );

		return sqrt( m2 );
	}

	void   setPx(double px) { _px = px; }
	void   setPy(double py) { _py = py; }
	void   setPz(double pz) { _pz = pz; }
	void   setE (double e ) { _e  = e;  }
	void   print()
	{
		if(_pdgid) printf("%4d: %15.8e %15.8e %15.8e %15.8e  | %15.8e\n",
		                  pdgid(),px(),  py(),  pz(),  e(),    recalculated_mass());
		else       printf(" SUM: %15.8e %15.8e %15.8e %15.8e  | %15.8e\n",
		                          px(),  py(),  pz(),  e(),    recalculated_mass());
	}
public:
/*
	These methods work on above accessors - they don't have to be changed
	even if above methods change.
*/
	double getAnglePhi();
	double getAngleTheta();
	void   rotateXZ(double theta);
	void   rotateXY(double theta);
	void   boostAlongZ(double pz, double e);
	void   boostToRestFrame(Particle &p);
	void   boostFromRestFrame(Particle &p);
private:
  double _px,_py,_pz,_e;
  int    _pdgid;
};

inline double Particle::getAnglePhi()
{
	// conventions as in ANGFI(X,Y) of tauola.f and PHOAN1 of photos.f
	// but now 'phi' in name define that it is rotation in px py
  
	double buf = 0.;
  
	if(fabs(py())<fabs(px()))
	{
		buf = atan( fabs(py()/px()) );
		if(px()<0.) buf = M_PI-buf;
	}
	else buf = acos( px()/sqrt(px()*px()+py()*py()) );
  
	if(py()<0.)	buf = 2.*M_PI-buf;

	return buf;
}

inline double Particle::getAngleTheta()
{
	// conventions as in ANGXY(X,Y) of tauola.f or PHOAN2 of photos.f
	// but now 'theta' in name define that it is rotation in pz px 
	// note that first argument of PHOAN2 was usually z
  
	double buf = 0.;
  
	if(fabs(px())<fabs(pz()))
	{
		buf = atan( fabs(px()/pz()) );
		if(pz()<0.) buf = M_PI-buf;
	}
	else buf = acos( pz()/sqrt(pz()*pz()+px()*px()) );

	return buf;
}

inline void Particle::rotateXZ(double theta)
{
	// as PHORO2 of photos.f
	double pX=px();
	double pZ=pz();
	setPx( cos(theta)*pX + sin(theta)*pZ);
	setPz(-sin(theta)*pX + cos(theta)*pZ);
}

inline void Particle::rotateXY(double theta)
{
	// as PHORO3 of photos.f
	double pX=px();
	double pY=py();
	setPx( cos(theta)*pX - sin(theta)*pY);
	setPy( sin(theta)*pX + cos(theta)*pY);
}

inline void Particle::boostAlongZ(double p_pz, double p_e)
{
	// as PHOBO3 of photos.f
	double m=sqrt(p_e*p_e-p_pz*p_pz);
	double buf_pz=pz();
	double buf_e=e();

	setPz((p_e *buf_pz + p_pz*buf_e)/m);
	setE ((p_pz*buf_pz + p_e *buf_e)/m);
}

inline void Particle::boostToRestFrame(Particle &p)
{
	double p_len = sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz());
	double phi   = p.getAnglePhi();
	p.rotateXY(-phi );
	double theta = p.getAngleTheta();
	p.rotateXY( phi );

	//Now rotate coordinates to get boost in Z direction.
	rotateXY(-phi  );
	rotateXZ(-theta);
	boostAlongZ(-1*p_len,p.e());
	rotateXZ( theta);
	rotateXY( phi  );
}

inline void Particle::boostFromRestFrame(Particle &p)
{
	double p_len = sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz());
	double phi   = p.getAnglePhi();
	p.rotateXY(-phi );
	double theta = p.getAngleTheta();
	p.rotateXY( phi );

	//Now rotate coordinates to get boost in Z direction.
	rotateXY(-phi  );
	rotateXZ(-theta);
	boostAlongZ(p_len,p.e());
	rotateXZ( theta);
	rotateXY( phi  );
}

} // namespace TauSpinner
#endif
