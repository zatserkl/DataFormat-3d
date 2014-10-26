// Andriy Zatserklyaniy, April 17, 2014

#ifndef Track_h
#define Track_h

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TClonesArray.h>
#include <TPaletteAxis.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <iomanip>
#include <cstdio>
#include <algorithm>

using std::cout;     using std::endl;

class BeamSpot {
public:
  Double_t v_;
  Double_t t_;
  Double_t u_;
  Double_t r_;      // beam spot sigma
  BeamSpot(Double_t v, Double_t t, Double_t u, Double_t r=15): v_(v), t_(t), u_(u), r_(r) {}   // r = 15 mm
};

class Vertex2D {
public:
  Double_t x_;
  Double_t u_;
  Double_t r_;      // vertex sigma
  Vertex2D(Double_t x, Double_t u, Double_t r=15): x_(x), u_(u), r_(r) {}   // r = 15 mm
};

//-------------------------------- CRay.h begin ----------------------------------

class Hit2D: public TObject {
public:
   Double_t u_;
   Double_t pos_;
   Hit2D(): TObject(), u_(0), pos_(0) {}
   Hit2D(Double_t u, Double_t pos): TObject(), u_(u), pos_(pos) {}
   Hit2D(const Hit2D& hit): TObject(hit), u_(hit.u_), pos_(hit.pos_) {
      //cout<< "copy ctor Hit2D" <<endl;
   }
   ~Hit2D() {
      //cout<< "dtor Hit2D" <<endl;
   }

   ClassDef(Hit2D, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class Hit2D;
#endif

class Hit3D: public TObject {
public:
   const Hit2D* hitv_;
   const Hit2D* hitt_;
   Hit3D(): TObject(), hitv_(0), hitt_(0) {}
   Hit3D(const Hit2D* hitv, const Hit2D* hitt): TObject(), hitv_(hitv), hitt_(hitt) {}
   Hit3D(const Hit3D& hit3D): TObject(hit3D), hitv_(hit3D.hitv_), hitt_(hit3D.hitt_) {
      //cout<< "copy ctor Hit3D" <<endl;
   }
   ~Hit3D() {
      //cout<< "dtor Hit3D" <<endl;
   }

   ClassDef(Hit3D, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class Hit3D;
#endif

class SensorHit: public Hit2D {
   friend std::ostream& operator<<(std::ostream&, const SensorHit&);
public:
   Int_t sensorId_;
   Int_t nfirst_;                         // the first strip in the cluster (0..383)
   Int_t nstrips_;                        // the number of strips in the cluster
public:
   SensorHit(): Hit2D()
      , sensorId_(-1)
      , nfirst_(-1)
      , nstrips_(0)
   {
      //cout<< "ctor default SensorHit" <<endl;
   }
   SensorHit(Int_t sensorId, Int_t nfirst, Int_t nstrips, Double_t u, Double_t pos): Hit2D(u,pos)
      , sensorId_(sensorId)
      , nfirst_(nfirst)
      , nstrips_(nstrips)
   {}
   SensorHit(const SensorHit& hit): Hit2D(hit)
                                    , sensorId_(hit.sensorId_)
                                    , nfirst_(hit.nfirst_)
                                    , nstrips_(hit.nstrips_)
   {
      //cout<< "copy ctor SensorHit" <<endl;
   }
   ~SensorHit() {
      //cout<< "dtor SensorHit" <<endl;
   }
   ClassDef(SensorHit, 5);
};

std::ostream& operator << (std::ostream& os, const SensorHit& hit) {
   return os
     << &hit
     << " sensorId " << std::setw(4) << hit.sensorId_
     << " nfirst " << std::setw(3) << hit.nfirst_
     << " nstrips " << hit.nstrips_
     << " u " << hit.u_
     << " pos " << hit.pos_;
}

#ifdef __MAKECINT__
#pragma link C++ class SensorHit;
#endif

class SensorHit3D: public TObject {
   friend std::ostream& operator<<(std::ostream&, const SensorHit3D&);
public:
   const SensorHit *hitv_;
   const SensorHit *hitt_;
public:
   SensorHit3D(): TObject(), hitv_(0), hitt_(0) {}
   SensorHit3D(const SensorHit* hitv, const SensorHit* hitt): TObject(), hitv_(hitv), hitt_(hitt) {}
   SensorHit3D(const SensorHit3D& hit3d): TObject(hit3d)
   {
      //cout<< "copy ctor SensorHit3D" <<endl;
      hitv_ = hit3d.hitv_;
      hitt_ = hit3d.hitt_;
   }
   ~SensorHit3D() {
      //cout<< "dtor SensorHit3D" <<endl;
   }
   ClassDef(SensorHit3D, 1);
};

std::ostream& operator << (std::ostream& os, const SensorHit3D& hit) {
   return os
     << "hitv_"
     << hit.hitv_
     << "\n"
     << "hitt_"
     << hit.hitt_;
}

#ifdef __MAKECINT__
#pragma link C++ class SensorHit3D;
#endif

class CRay2D: public TObject {
friend std::ostream& operator<<(std::ostream&, const CRay2D&);
public:
   Double_t x_;               // 2D radiant: u_ = 0
   Double_t cx_, cu_;	        // direction cosines
public:
   void clear() {
      x_ = 0;
      cx_ = cu_ = 0;
   }
   CRay2D(): TObject() {
      clear();
   }
   CRay2D(const CRay2D& ray2d): TObject(ray2d) {
      //cout<< "copy ctor CRay2D" <<endl;
      x_ = ray2d.x_;
      cx_ = ray2d.cx_;
      cu_ = ray2d.cu_;
   }
   CRay2D(Double_t x1, Double_t u1, Double_t x2, Double_t u2): TObject() {
      Double_t dx = x2 - x1;
      Double_t du = u2 - u1;
      Double_t alpha = TMath::ATan2(dx, du);
      cx_ = TMath::Sin(alpha);                         // cos(pi/2 - alpha) = sin(alpha)
      cu_ = TMath::Cos(alpha);
      // propagate to the plane u = 0
      Double_t p = -u1 / cu_;               // (0 - hit1.u_)/cu_
      x_ = x1 + p*cx_;
   }
   CRay2D(const Hit2D* hit1, const Hit2D* hit2): TObject() {
      Double_t dx = hit2->pos_ - hit1->pos_;
      Double_t du = hit2->u_ - hit1->u_;
      Double_t alpha = TMath::ATan2(dx, du);
      cx_ = TMath::Sin(alpha);                         // cos(pi/2 - alpha) = sin(alpha)
      cu_ = TMath::Cos(alpha);
      // propagate to the plane u = 0
      Double_t p = -hit1->u_ / cu_;          // (0 - hit1->u_)/cu_
      x_ = hit1->pos_ + p*cx_;
   }
   ~CRay2D() {
      //cout<< "dtor CRay2D" <<endl;
   }
   Double_t Theta() const {return TMath::ACos(cu_);}
   static Double_t Angle(const CRay2D* track1, const CRay2D* track2) {     // standalone function
      Double_t scalar = track1->cx_*track2->cx_ + track1->cu_*track2->cu_;
      Double_t theta = TMath::Pi()/2.;
      if (TMath::Abs(scalar) < 1.) theta = TMath::ACos(scalar);
      return theta;
   }
   Double_t at(Double_t uplane) const {
      const Double_t eps = 1e-7;
      Double_t x_proj = 1./eps;
      if (TMath::Abs(cu_) > eps) {
      	 Double_t p = uplane/cu_;        // (uplane - 0)/cu_
      	 x_proj = x_ + cx_*p;
      }
      // cout<< "CRay2D::at: x_proj = " << x_proj <<endl;
      return x_proj;
   }

   ClassDef(CRay2D, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class CRay2D;
#endif

std::ostream& operator << (std::ostream& os, const CRay2D& ray) {
   os << "x_ = " << ray.x_ << " cx_ = " << ray.cx_ << " cu_ = " << ray.cu_;
   return os;
}

class CRay: public TObject {
public:
   Double_t v_, t_;           // radiant: u_ = 0
   Double_t cv_, ct_, cu_;	   // direction cosines
public:
   void clear() {
      v_ = t_ = 0;
      cv_ = ct_ = cu_ = 0;
   }
   CRay(): TObject() {
      clear();
   }
   CRay(const CRay& ray): TObject(ray) {
      //cout<< "copy ctor CRay" <<endl;
      v_ = ray.v_;
      t_ = ray.t_;
      cv_ = ray.cv_;
      ct_ = ray.ct_;
      cu_ = ray.cu_;
   }
   CRay(Double_t v1, Double_t t1, Double_t u1, Double_t v2, Double_t t2, Double_t u2): TObject()
   {
      Double_t dv = v2 - v1;
      Double_t dt = t2 - t1;
      Double_t du = u2 - u1;
      cv_ = TMath::Sin(TMath::ATan2(dv, du));                         // cos(pi/2 - alpha) = sin(alpha)
      ct_ = TMath::Sin(TMath::ATan2(dt, du));
      cu_ = TMath::Cos(TMath::ATan2(TMath::Sqrt(dv*dv + dt*dt), du));
      // propagate to the plane u = 0
      Double_t p = -u1 / cu_;               // (0 - hit1.u_)/cu_
      v_ = v1 + p*cv_;
      t_ = t1 + p*ct_;
   }
   CRay(const CRay2D* v_ray, const CRay2D* t_ray): TObject() {
      // ratio of the 2D directional cosines (will be used later)
      Double_t v_xu = v_ray->cx_/v_ray->cu_;
      Double_t t_xu = t_ray->cx_/t_ray->cu_;
      // 3D directional cosine for u
      cu_ = 1./TMath::Sqrt(1. + (v_xu*v_xu + t_xu*t_xu));
      // 3D directional cosines for v and t
      cv_ = v_xu * cu_;
      ct_ = t_xu * cu_;
      // v- and t-components of the radiant
      v_ = v_ray->x_;
      t_ = t_ray->x_;
   }
   CRay(const SensorHit3D* hit1, const SensorHit3D* hit2): TObject() {
      CRay2D v_ray(hit1->hitv_, hit2->hitv_);
      CRay2D t_ray(hit1->hitt_, hit2->hitt_);
      // ratio of the 2D directional cosines (will be used later)
      Double_t v_xu = v_ray.cx_/v_ray.cu_;
      Double_t t_xu = t_ray.cx_/t_ray.cu_;
      // 3D directional cosine for u
      cu_ = 1./TMath::Sqrt(1. + (v_xu*v_xu + t_xu*t_xu));
      // 3D directional cosines for v and t
      cv_ = v_xu * cu_;
      ct_ = t_xu * cu_;
      // v- and t-components of the radiant
      v_ = v_ray.x_;
      t_ = t_ray.x_;
   }
   ~CRay() {
      //cout<< "dtor CRay" <<endl;
   }
   void at(Double_t uplane, Double_t& v_proj, Double_t& t_proj) const {
      const Double_t eps = 1e-7;
      v_proj = t_proj = 1/eps;
      if (TMath::Abs(cu_) < eps) return;
      Double_t p = uplane/cu_;        // (uplane - 0)/cu_
      v_proj = v_ + cv_*p;
      t_proj = t_ + ct_*p;
   }
   Double_t V(Double_t uplane) const {
      const Double_t eps = 1e-7;
      Double_t v_proj = 1./eps;
      if (TMath::Abs(cu_) > eps) {
      	 Double_t p = uplane/cu_;        // (uplane - 0)/cu_
      	 v_proj = v_ + cv_*p;
      }
      return v_proj;
   }
   Double_t T(Double_t uplane) const {
      const Double_t eps = 1e-7;
      Double_t t_proj = 1./eps;
      if (TMath::Abs(cu_) > eps) {
      	 Double_t p = uplane/cu_;        // (uplane - 0)/cu_
      	 t_proj = t_ + ct_*p;
      }
      return t_proj;
   }
   Double_t static Angle(const CRay& ray1, const CRay& ray2) {
      Double_t scalar = ray1.cv_*ray2.cv_ + ray1.ct_*ray2.ct_ + ray1.cu_*ray2.cu_;
      return scalar < 1.? TMath::ACos(ray1.cv_*ray2.cv_ + ray1.ct_*ray2.ct_ + ray1.cu_*ray2.cu_): 0;
   }
   Double_t Angle(const CRay& ray) const {
      Double_t scalar = cv_*ray.cv_ + ct_*ray.ct_ + cu_*ray.cu_;
      return scalar < 1.? TMath::ACos(cv_*ray.cv_ + ct_*ray.ct_ + cu_*ray.cu_): 0;
   }
   Double_t Theta() const {return TMath::ACos(cu_);}

   ClassDef(CRay, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class CRay;
#endif

//-------------------------------- CRay.h end ----------------------------------

class Track2D: public CRay2D {
friend std::ostream& operator<<(std::ostream&, const Track2D&);
  // CRay2D with 2 hits
public:
   const SensorHit* hit1_;
   const SensorHit* hit2_;
   Bool_t owner;                    // owner of hits (created them with copy constructor)
public:
   Track2D(): CRay2D(), hit1_(0), hit2_(0), owner(kFALSE) {
      //cout<< "ctor default Track2D" <<endl;
   }
   Track2D(const Track2D& track2D): CRay2D(track2D), owner(kTRUE) {
      //cout<< "copy ctor Track2D" <<endl;
      hit1_ = new SensorHit(*track2D.hit1_);
      hit2_ = new SensorHit(*track2D.hit2_);
      // assign hit pointers from the track2D provided
      //-- hit1_ = track2D.hit1_;
      //-- hit2_ = track2D.hit2_;
   }
   Track2D(const SensorHit* hit1, const SensorHit* hit2): CRay2D(hit1,hit2), owner(kFALSE) {
      hit1_ = hit1;
      hit2_ = hit2;
   }
   Track2D(Double_t x1, Double_t u1, Double_t x2, Double_t u2): CRay2D(x1,u1, x2,u2), owner(kFALSE) {}
   ~Track2D() {
      //cout<< "dtor Track2D" <<endl;
      // do NOT delete the hits

      if (owner) {
         delete hit1_;
         delete hit2_;
      }

      //cout<< "dtor Track2D: return" <<endl;
   }
   bool hitSensor(Int_t sensorId) const {
      return hit1_->sensorId_ == sensorId || hit2_->sensorId_ == sensorId;
   }

   ClassDef(Track2D, 3)
};

#ifdef __MAKECINT__
#pragma link C++ class Track2D;
#endif

class Track3D: public CRay {
friend std::ostream& operator<<(std::ostream&, const Track3D&);
  // CRay3D with 2 hits
public:
   const SensorHit3D* hit1_;
   const SensorHit3D* hit2_;
public:
   Track3D(): CRay(), hit1_(0), hit2_(0) {
      //cout<< "ctor default Track3D" <<endl;
   }
   Track3D(const Track3D& track3D): CRay(track3D) {
      //cout<< "copy ctor Track3D" <<endl;
      hit1_ = new SensorHit3D(*track3D.hit1_);
      hit2_ = new SensorHit3D(*track3D.hit2_);
      // assign hit pointers from the track3D provided   TODO: check if that applicable
      //-- hit1_ = track3D.hit1_;
      //-- hit2_ = track3D.hit2_;
   }
   Track3D(const SensorHit3D* hit1, const SensorHit3D* hit2): CRay(hit1,hit2) {
      hit1_ = hit1;
      hit2_ = hit2;
   }
   //--n/a-- Track3D(Double_t x1, Double_t u1, Double_t x2, Double_t u2): CRay(x1,u1, x2,u2) {}
   ~Track3D() {
      //cout<< "dtor Track3D" <<endl;
      // do NOT delete the hits
      delete hit1_;
      delete hit2_;
      //cout<< "dtor Track3D: return" <<endl;
   }
   bool hitSensor(Int_t sensorId) const {
      return false
         || hit1_->hitv_->sensorId_ == sensorId
         || hit1_->hitt_->sensorId_ == sensorId
         || hit2_->hitv_->sensorId_ == sensorId
         || hit2_->hitt_->sensorId_ == sensorId
      ;
   }
   bool hitVSensor(Int_t sensorId) const {
      return false
         || hit1_->hitv_->sensorId_ == sensorId
         || hit2_->hitv_->sensorId_ == sensorId
      ;
   }
   bool hitTSensor(Int_t sensorId) const {
      return false
         || hit1_->hitt_->sensorId_ == sensorId
         || hit2_->hitt_->sensorId_ == sensorId
      ;
   }

   ClassDef(Track3D, 3)
};

#ifdef __MAKECINT__
#pragma link C++ class Track3D;
#endif

std::ostream& operator << (std::ostream& os, const Track2D& track2D) {
   os << "x_ = " << track2D.x_ << " cx_ = " << track2D.cx_ << " cu_ = " << track2D.cu_ << " hit1 " << *track2D.hit1_ << " hit2 " << *track2D.hit2_;
   return os;
}

class Track: public CRay {
  // CRay with v- and t-tracks
public:
   const Track2D* vTrack_;
   const Track2D* tTrack_;
public:
   Track(): CRay(), vTrack_(0), tTrack_(0) {
      //cout<< "ctor default Track" <<endl;
   }
   Track(const Track& track): CRay(track) {
      //cout<< "copy ctor Track" <<endl;
      // vTrack_ = track.vTrack_;
      // tTrack_ = track.tTrack_;
      vTrack_ = new Track2D(*track.vTrack_);
      tTrack_ = new Track2D(*track.tTrack_);
   }
   Track(const Track2D* vTrack, const Track2D* tTrack): CRay(vTrack,tTrack), vTrack_(vTrack), tTrack_(tTrack) {}
   ~Track() {
      //cout<< "dtor Track" <<endl;
      delete vTrack_;
      delete tTrack_;
   }
   bool hitSensor(Int_t sensorId) const {
      return vTrack_->hitSensor(sensorId) || tTrack_->hitSensor(sensorId);
   }

   ClassDef(Track, 4);
};

std::ostream& operator << (std::ostream& os, const Track& track) {
   os << "vTrack: " << *track.vTrack_ << "\ntTrack: " << *track.tTrack_ <<endl;
   return os;
}

#ifdef __MAKECINT__
#pragma link C++ class Track;
#endif

class SuperTrack: public TObject {
public:
   const Track* itrack_;
   const Track* otrack_;
   Double_t angle;
   Double_t d;                            // distance between hits in the plane u = 0
   // TODO: add vcal and tcal
public:
   SuperTrack(): TObject(), itrack_(0), otrack_(0), angle(0), d(0) {
      //cout<< "ctor default SuperTrack" <<endl;
   }
   SuperTrack(const SuperTrack& superTrack): TObject(superTrack) {
      //cout<< "copy ctor SuperTrack" <<endl;
      //-- itrack_ = superTrack.itrack_;
      //-- otrack_ = superTrack.otrack_;
      itrack_ = new Track(*superTrack.itrack_);
      otrack_ = new Track(*superTrack.otrack_);
      angle = superTrack.angle;
      d = superTrack.d;
   }
   SuperTrack(const Track* itrack, const Track* otrack): TObject(), itrack_(itrack), otrack_(otrack) {
      angle = Angle();
      d = Distance();
   }
   ~SuperTrack() {
      //cout<< "dtor SuperTrack: itrack_ = " << itrack_ << " otrack_ = " << otrack_ <<endl;
      delete itrack_;
      delete otrack_;
      //cout<< "dtor SuperTrack: return" <<endl;
   }
   void Clear(Option_t*) {
      //cout<< "SuperTrack::Clear" <<endl;
      delete itrack_;   itrack_ = 0;
      delete otrack_;   otrack_ = 0;
      angle = 0;
      d = 0;
   }
   Double_t Angle() const {
      Double_t scalar = itrack_->cv_*otrack_->cv_ + itrack_->ct_*otrack_->ct_ + itrack_->cu_*otrack_->cu_;
      //cout<< "scalar = " << scalar <<endl;
      return scalar < 1.? TMath::ACos(scalar): 0;
   }
   Double_t Distance() const {
      // Distance between the hits in the plane u = 0. Just distance between the radiants.
      Double_t dv = otrack_->v_ - itrack_->v_;
      Double_t dt = otrack_->t_ - itrack_->t_;
      return TMath::Sqrt(dv*dv + dt*dt);
   }
   Double_t V(Double_t u) const {
      if (u > 0) return otrack_->V(u);
      else return itrack_->V(u);
   }
   Double_t T(Double_t u) const {
      if (u > 0) return otrack_->T(u);
      else return itrack_->T(u);
   }
   void at(Double_t u, Double_t& v, Double_t& t) const {
      v = t = 0;
      if (u > 0) {
         v = itrack_->V(u);
         t = itrack_->T(u);
      }
      else {
         v = otrack_->V(u);
         t = otrack_->T(u);
      }
   }
   bool iHitSensor(Int_t sensorId) const {
      return itrack_->hitSensor(sensorId);
   }
   bool oHitSensor(Int_t sensorId) const {
      return otrack_->hitSensor(sensorId);
   }
   // bool HitSensor(Int_t sensorId) const {                      // NB: do I really need this time-consuming method?
   //    return iHitSensor(sensorId) || oHitSensor(sensorId);
   // }
   Bool_t SharedHits(const SuperTrack* superTrack) const {
      if (itrack_->vTrack_->hit1_ == superTrack->itrack_->vTrack_->hit1_) return kTRUE;
      if (itrack_->vTrack_->hit2_ == superTrack->itrack_->vTrack_->hit2_) return kTRUE;
      if (itrack_->tTrack_->hit1_ == superTrack->itrack_->tTrack_->hit1_) return kTRUE;
      if (itrack_->tTrack_->hit2_ == superTrack->itrack_->tTrack_->hit2_) return kTRUE;
      if (otrack_->vTrack_->hit1_ == superTrack->otrack_->vTrack_->hit1_) return kTRUE;
      if (otrack_->vTrack_->hit2_ == superTrack->otrack_->vTrack_->hit2_) return kTRUE;
      if (otrack_->tTrack_->hit1_ == superTrack->otrack_->tTrack_->hit1_) return kTRUE;
      if (otrack_->tTrack_->hit2_ == superTrack->otrack_->tTrack_->hit2_) return kTRUE;
      return kFALSE;
   }
   bool hitSensor(Int_t sensorId) const {
      //return itrack_->hitSensor(sensorId) || itrack_->hitSensor(sensorId);
      return iHitSensor(sensorId) || oHitSensor(sensorId);
   }

   ClassDef(SuperTrack, 6);
};

std::ostream& operator << (std::ostream& os, const SuperTrack& superTrack) {
   os << &superTrack << "\tangle = " << superTrack.angle << " distance = " << superTrack.Distance() << "\nitrack: " << *superTrack.itrack_ << "\notrack: " << *superTrack.otrack_ <<endl;
   return os;
}

#ifdef __MAKECINT__
#pragma link C++ class SuperTrack;
#endif

#endif	// #ifdef Track_h
