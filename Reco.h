// Andriy Zatserklyaniy, April 17, 2014:  production release

#ifndef Reco_h
#define Reco_h

#include "DataFormat.h"
#include "Track.h"
#include "Sensor.h"
#include "Geometry.h"

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <vector>
#include <utility>

using std::cout;    using std::endl;

ClassImp(Hit2D);
ClassImp(SensorHit);
ClassImp(CRay2D);
ClassImp(CRay);
ClassImp(Track2D);
ClassImp(Track);
ClassImp(SuperTrack);

class RecoEvent;

class Reco {
   friend class RecoEvent;
public:
   Bool_t debug_;
   Float_t deltaT_;
   const Geometry* geometry_;
   PCTSensors* pCTSensors_;
   const BeamSpot* beamSpot_;
   const PCTEvent* pCTEvent_;
   Int_t event_;
   Int_t nhitv_[4];
   Int_t nhitt_[4];
   Float_t vreco_[4];
   Float_t treco_[4];
   std::list<const Track2D*> tin_;
   std::list<const Track2D*> tout_;
   std::list<const Track2D*> vin_;
   std::list<const Track2D*> vout_;

   std::list<const Track*> itracks_;
   std::list<const Track*> otracks_;
   std::list<const SuperTrack*> superTracks_;

   std::list<const Track2D*> tTracksIn_;
   std::list<const Track2D*> tTracksOut_;
   std::list<const Track2D*> vTracksIn_;
   std::list<const Track2D*> vTracksOut_;

   static TClonesArray* poolTrack2D_;                 //->
   static TClonesArray* poolTrack_;                   //->
   static TClonesArray* poolSuperTrack_;              //->
public:
   Reco(const Geometry* geometry, const BeamSpot* beamSpot, PCTSensors* pCTSensors, const PCTEvent* pCTEvent, Int_t event, bool debug=kFALSE): debug_(debug)
      , geometry_(geometry), pCTSensors_(pCTSensors), beamSpot_(beamSpot), pCTEvent_(pCTEvent)
      , event_(event)
   {
      //cout<< "Reco::Reco" <<endl;

      for (int ilayer=0; ilayer<4; ++ilayer) {
         vreco_[ilayer] = 0;
         treco_[ilayer] = 0;
      }

      deltaT_ = pCTEvent_->deltaT;
      //cout<< "Reco::Reco: create pool" <<endl;
      //Sensor::CreateHitPool();
      Sensor::ClearPool();

      if (!poolTrack2D_) poolTrack2D_ = new TClonesArray("Track2D", 1024);
      if (!poolTrack_) poolTrack_ = new TClonesArray("Track", 1024);
      if (!poolSuperTrack_) poolSuperTrack_ = new TClonesArray("SuperTrack", 1024);

      //cout<< "Sensor::poolSensorHit_->GetLast()+1 = " << Sensor::poolSensorHit_->GetLast()+1 <<endl;
      //cout<< "poolTrack2D_->GetLast()+1 = " << poolTrack2D_->GetLast()+1 <<endl;
      //cout<< "poolSuperTrack2D_->GetLast()+1 = " << poolSuperTrack2D_->GetLast()+1 <<endl;
      //cout<< "poolSuperTrack_->GetLast()+1 = " << poolSuperTrack_->GetLast()+1 <<endl;

      //----------- pCTSensors_ = new PCTSensors();
      pCTSensors_->clear();

      // get hits for this event
      for (int iFPGA=0; iFPGA<12; ++iFPGA)
      {
         const TrackerFPGA& trackerFPGA = pCTEvent_->trackerFPGA[iFPGA];
         for (unsigned int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
         {
            TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
            for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
               if (debug_) {
                  Int_t strip = 64*trackerChip->address + (63-trackerChip->nfirst[icluster]);
                  cout<< event_ << "\t iFPGA = " << std::setw(2) << iFPGA
                  << " trackerChip->address = " << std::setw(2) << (unsigned) trackerChip->address
                  << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster]
                  << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster]
                  << " strip = " << strip <<endl;
                  // printf("%-8d iFPGA = %2d trackerChip->address = %2d nstrips[%d] = %d nfirst[%d] = %d strip = %d\n",
                  //        event,iFPGA,(unsigned) trackerChip->address,icluster,(unsigned) trackerChip->nstrips[icluster],icluster,(unsigned) trackerChip->nfirst[icluster],strip);
               }
               pCTSensors_->AddCluster(iFPGA, trackerChip->address, trackerChip->nfirst[icluster], trackerChip->nstrips[icluster]);
            }
         }
      }
      // process hits: combine clusters extended over the chip boundary, convert strips to mm
      pCTSensors_->GetHits();

      Int_t n_v_hits = 0;
      for (unsigned ilayer=0; ilayer<4; ++ilayer) {
         n_v_hits += pCTSensors_->v_hits[ilayer].size();
         nhitv_[ilayer] = pCTSensors_->v_hits[ilayer].size();
      }
      Int_t n_t_hits = 0;
      for (unsigned ilayer=0; ilayer<4; ++ilayer) {
         n_t_hits += pCTSensors_->t_hits[ilayer].size();
         //--------------------------------------------------- nhitt_[ilayer] = pCTSensors_->v_hits[ilayer].size();
         nhitt_[ilayer] = pCTSensors_->t_hits[ilayer].size();
      }

      if (debug_) cout<< "n_v_hits = " << n_v_hits << " n_t_hits = " << n_t_hits <<endl;

      GenerateTracksTV();

      if (debug_) cout<< "------------------------ generate Tracks ---" <<endl;
      GenerateTracks();

      if (debug_) cout<< "------------------------ generate super tracks ---" <<endl;
      GenerateSuperTracks();

      if (debug_) cout<< "------------------------ FilterPhaseSpace ---" <<endl;
      FilterPhaseSpace();
   }
   void GenerateTracksTV()
   {
      bool debug0 = debug_;

      std::vector<const Vertex2D*> vertices;

      Vertex2D beamSpotT(beamSpot_->t_, beamSpot_->u_, beamSpot_->r_);
      vertices.push_back(&beamSpotT);

      Int_t layer = 0;

      // debug_ = true;
      GenerateTracksT(vertices, layer, tTracksIn_);
      debug_ = debug0;

      // start from the t-boards to know the position of the t-gaps in the V-boards
      //GenerateTracksT();

      vertices.clear();

      Vertex2D beamSpotV(beamSpot_->v_, beamSpot_->u_, beamSpot_->r_);
      vertices.push_back(&beamSpotV);

      layer = 0;
      GenerateTracksV(vertices, layer, tTracksIn_, vTracksIn_);

      if (debug_) cout<< "vTracksIn_.size() = " << vTracksIn_.size() << " tTracksIn_.size() = " << tTracksIn_.size() <<endl;

      if (vTracksIn_.size() == 0) return;
      if (vTracksIn_.size() > 2) return;
      if (vTracksIn_.size() == 2) {
         bool different_halves = false;
         const Track2D* track1 = vTracksIn_.front();
         const Track2D* track2 = vTracksIn_.back();
         if (true
             && track1->hit1_->sensorId_ > 0
             && track2->hit1_->sensorId_ > 0
             && track1->hit1_->sensorId_ != track2->hit1_->sensorId_
            )
         different_halves = true;
         if (true
             && track1->hit2_->sensorId_ > 0
             && track2->hit2_->sensorId_ > 0
             && track1->hit2_->sensorId_ != track2->hit2_->sensorId_
            )
         different_halves = true;
         if (!different_halves) return;
      }

      layer = 2;
      // Double_t sigma = 10.;   // mm
      Double_t sigma = 1.;   // mm

      if (debug_) cout<< "\nGenerateTracksTV: Generate output tracks" <<endl;

      vertices.clear();

      for (std::list<const Track2D*>::const_iterator it=tTracksIn_.begin(); it!=tTracksIn_.end(); ++it) {
         const Track2D* track2D = *it;
         Vertex2D* vertex = new Vertex2D(track2D->at(0), 0, sigma);
         vertices.push_back(vertex);
      }
      //vertices.push_back(&beamSpotT);

      if (debug_) cout<< "call GenerateTracksT(vertices, " << layer << ", tTracksOut_)" <<endl;
      GenerateTracksT(vertices, layer, tTracksOut_);
      if (debug_) cout<< "tTracksOut_.size() = " << tTracksOut_.size() <<endl;

      vertices.clear();

      for (std::list<const Track2D*>::const_iterator it=vTracksIn_.begin(); it!=vTracksIn_.end(); ++it) {
         const Track2D* track2D = *it;
         Vertex2D* vertex = new Vertex2D(track2D->at(0), 0, sigma);
         vertices.push_back(vertex);
      }
      GenerateTracksV(vertices, layer, tTracksOut_, vTracksOut_);
      if (debug_) cout<< "vTracksOut_.size() = " << vTracksOut_.size() << " tTracksOut_.size() = " << tTracksOut_.size() <<endl;

      vin_.merge(vTracksIn_);
      vout_.merge(vTracksOut_);
      tin_.merge(tTracksIn_);
      tout_.merge(tTracksOut_);
   }
   void GenerateTracksT(std::vector<const Vertex2D*>vertices, Int_t layer_start, std::list<const Track2D*>& tTracks)
   {
      // lists of all hits
      std::list<const SensorHit*> tlist[4];

      for (unsigned ilayer=0; ilayer<4; ++ilayer) {
         for (unsigned ihit=0; ihit<pCTSensors_->t_hits[ilayer].size(); ++ihit) {
            tlist[ilayer].push_back(pCTSensors_->t_hits[ilayer][ihit]);
         }
      }

      for (unsigned ivertex=0; ivertex<vertices.size(); ++ivertex)
      {
         Int_t layer = layer_start;

         const Vertex2D* vertex = vertices[ivertex];

         SensorHit vertexHit_t(0, 0, 0, vertex->u_, vertex->x_);

         if (debug_) cout<< "layer " << layer << " vertexHit_t: " << vertexHit_t <<endl;

         ///////////////////////////////////////////////
         //
         // look at t-hits at the upstream layer
         //
         ///////////////////////////////////////////////

         const Double_t distance_max = 1. + 5.*vertex->r_*(geometry_->ut_[layer+1]-geometry_->ut_[layer])/(geometry_->ut_[layer]-vertex->u_);
         if (debug_) cout<< "GenerateTracksT: distance_max = " << distance_max <<endl;

         if (debug_) cout<< "tlist[" << layer << "].size() = " << tlist[layer].size() <<endl;

         //-- std::list<const SensorHit*>::const_iterator it0 = tlist[layer].begin();   // erase with const_iterator is an error in SL6
         std::list<const SensorHit*>::iterator it0 = tlist[layer].begin();
         while (it0!=tlist[layer].end())
         {
            const SensorHit* hit0 = *it0;
            bool erase_hit = false;                         // hit has been assigned to Track2D and need to be erased from tlist

            Track2D vertex_track(&vertexHit_t, hit0);

            // project the vertex track to the next layer
            Double_t pos = vertex_track.at(geometry_->ut_[layer+1]);

            // find the closest hit in the next layer
            Double_t distance = distance_max;
            const SensorHit* closestHit = 0;

            for (std::list<const SensorHit*>::const_iterator it1=tlist[layer+1].begin(); it1!=tlist[layer+1].end(); ++it1) {
               const SensorHit* hit1 = *it1;
               Double_t d = hit1->pos_ - pos;
               if (TMath::Abs(d) < distance) {
                  closestHit = hit1;
                  distance = d;
               }
            }
            if (closestHit) tlist[layer+1].remove(closestHit);       // remove matched hit from the tlist[layer+1]

            if (!closestHit) {
               // there is no hit for this vertex track in the next layer. One of two possibilities:
               // 1) this is a noise hit
               // 2) the proton hit the gap in the next layer

               // check the gaps in the next layer
               for (int igap=0; igap<3; ++igap) {
                  if (TMath::Abs(pos - pCTSensors_->gap_.tgap_[layer+1][igap]) < pCTSensors_->gap_.width_) {
                     // recover this hit: place a hit in the middle of the gap TODO: fit existing hit and the vertex
                     Double_t sigma_vertex = vertex->r_*(geometry_->ut_[layer+1] - geometry_->ut_[layer])/(geometry_->ut_[layer]-vertex->u_);
                     Double_t sigma_gap = 1./TMath::Sqrt(12.);
                     Double_t w_vertex = 1./(sigma_vertex*sigma_vertex);
                     Double_t w_gap = 1./(sigma_gap*sigma_gap);
                     Double_t pos_mean = (w_vertex*pos + w_gap*pCTSensors_->gap_.tgap_[layer+1][igap]) / (w_vertex + w_gap);
                     Int_t sensorId = -1*(200 + layer*10 + (layer+1));
                     SensorHit* gapHit = new ((*Sensor::poolSensorHit_)[Sensor::poolSensorHit_->GetLast()+1]) SensorHit(sensorId, 0, 0, geometry_->ut_[layer+1], pos_mean);
                     //-- SensorHit* gapHit = new ((*Sensor::poolSensorHit_)[Sensor::poolSensorHit_->GetLast()+1]) SensorHit(sensorId, 0, 0, geometry_->ut_[layer+1], pCTSensors_->gap_.tgap_[layer+1][igap]);
                     closestHit = gapHit;
                     treco_[layer+1] = pos_mean;
                     if (debug_) cout<< "\n--> event_ = " << event_ << " recovered t-hit in downstream layer " << layer+1 << " pos - pCTSensors_->gap_.tgap_[" << layer+1 << "][igap] = " << pos - pCTSensors_->gap_.tgap_[layer+1][igap] << " pos = " << pos << " gapHit: " << *gapHit <<endl;
                     break;
                  }
               }

               if (!closestHit) {
                  // check the dead strips
               }
            }

            if (closestHit) {
               //--debug-- if (debug_) cout<< "found good hit in the downstream layer: erase it from the tlist[" << layer << "]" << endl;
               erase_hit = true;                              // erase the hit at the end of the loop body
               //
               // 1) create Track2D and add it to tTrack2D
               // 2) remove hits from the tlist
               //
               Track2D* track2D = new ((*poolTrack2D_)[poolTrack2D_->GetLast()+1]) Track2D(hit0, closestHit);
               tTracks.push_back(track2D);
            }

            //
            //-- erase the paired hit or increament the iterator
            //
            if (erase_hit) {
               it0 = tlist[layer].erase(it0);
            }
            else ++it0;
         }  // loop over hits in the upstream layer

         ///////////////////////////////////////////////
         //
         // look at t-hits at the downstream layer (if any)
         //
         ///////////////////////////////////////////////

         layer += 1;

         //-- std::list<const SensorHit*>::const_iterator it1 = tlist[layer].begin(); // erase with const_iterator is an error in SL6
         std::list<const SensorHit*>::iterator it1 = tlist[layer].begin();
         while (it1!=tlist[layer].end())
         {
            const SensorHit* hit1 = *it1;
            bool erase_hit = false;                         // hit has been assigned to Track2D and need to be erased from tlist

            Track2D vertex_track(&vertexHit_t, hit1);

            // project the vertex track to the next layer
            Double_t pos = vertex_track.at(geometry_->ut_[layer-1]);

            // find the closest hit in the next layer
            Double_t distance = distance_max;
            const SensorHit* closestHit = 0;

            for (it0=tlist[layer-1].begin(); it0!=tlist[layer-1].end(); ++it0) {       // reuse of std::list<const SensorHit*>::const_iterator it0 
               const SensorHit* hit0 = *it0;
               Double_t d = hit0->pos_ - pos;
               if (TMath::Abs(d) < distance) {
                  closestHit = hit0;
                  distance = d;
               }
            }
            if (closestHit) tlist[layer-1].remove(closestHit);       // remove matched hit from the tlist[layer+1]

            if (!closestHit) {
               // there is no hit for this vertex track in the upstream layer. One of two possibilities:
               // 1) this is a noise hit
               // 2) the proton hit the gap in the next layer

               // check the gaps in the previous layer
               for (int igap=0; igap<3; ++igap) {
                  if (TMath::Abs(pos - pCTSensors_->gap_.tgap_[layer-1][igap]) < pCTSensors_->gap_.width_) {
                     // recover this hit: place a hit in the middle of the gap TODO: fit existing hit and the vertex
                     Double_t sigma_vertex = vertex->r_*(geometry_->ut_[layer] - geometry_->ut_[layer-1])/(geometry_->ut_[layer]-vertex->u_);
                     Double_t sigma_gap = 1./TMath::Sqrt(12.);
                     Double_t w_vertex = 1./(sigma_vertex*sigma_vertex);
                     Double_t w_gap = 1./(sigma_gap*sigma_gap);
                     Double_t pos_mean = (w_vertex*pos + w_gap*pCTSensors_->gap_.tgap_[layer-1][igap]) / (w_vertex + w_gap);
                     Int_t sensorId = -1*(200 + layer*10 + (layer-1));
                     SensorHit* gapHit = new ((*Sensor::poolSensorHit_)[Sensor::poolSensorHit_->GetLast()+1]) SensorHit(sensorId, 0, 0, geometry_->ut_[layer-1], pos_mean);
                     //-- SensorHit* gapHit = new ((*Sensor::poolSensorHit_)[Sensor::poolSensorHit_->GetLast()+1]) SensorHit(sensorId, 0, 0, geometry_->ut_[layer-1], pCTSensors_->gap_.tgap_[layer-1][igap]);
                     closestHit = gapHit;
                     treco_[layer-1] = pos_mean;
                     if (debug_) cout<< "\n--> event_ = " << event_ << " recovered t-hit in upstream layer " << layer-1 << " pos - pCTSensors_->gap_.tgap_[" << layer-1 << "][igap] = " << pos - pCTSensors_->gap_.tgap_[layer-1][igap] << " pos = " << pos << " gapHit: " << *gapHit <<endl;
                     break;
                  }
               }

               if (!closestHit) {
                  // check the dead strips
               }
            }

            if (closestHit) {
               //--debug-- if (debug_) cout<< "found good hit in the upstream layer: erase it from the tlist[" << layer << "]" << endl;
               erase_hit = true;                              // erase the hit at the end of the loop body
               //
               // 1) create Track2D and add it to tTrack2D
               // 2) remove hits from the tlist
               //
               Track2D* track2D = new ((*poolTrack2D_)[poolTrack2D_->GetLast()+1]) Track2D(closestHit, hit1);
               tTracks.push_back(track2D);
            }

            //
            //-- erase the paired hit or increament the iterator
            //
            if (erase_hit) {
               it1 = tlist[layer].erase(it1);
            }
            else ++it1;
         }  // loop over hits in the upstream layer

         if (debug_) {
            //--debug-- cout<< "tTracks.size() = " << tTracks.size() <<endl;
         }
      }
   }
   void GenerateTracksV(std::vector<const Vertex2D*> vertices, Int_t layer_start, const std::list<const Track2D*>& tTracks, std::list<const Track2D*>& vTracks)
   {
      ////////////////////////////////////
      //
      // uses tTracks to find t-gaps
      //
      ////////////////////////////////////

      // lists of all hits
      std::list<const SensorHit*> vlist[4];

      for (unsigned ilayer=0; ilayer<4; ++ilayer) {
         for (unsigned ihit=0; ihit<pCTSensors_->v_hits[ilayer].size(); ++ihit) {
            vlist[ilayer].push_back(pCTSensors_->v_hits[ilayer][ihit]);
         }
      }

      for (unsigned ivertex=0; ivertex<vertices.size(); ++ivertex)
      {
         Int_t layer = layer_start;

         const Vertex2D* vertex = vertices[ivertex];

         SensorHit vertexHit_v(0, 0, 0, vertex->u_, vertex->x_);
         if (debug_) cout<< "layer " << layer << " vertexHit_v: " << vertexHit_v <<endl;

         ///////////////////////////////////////////////
         //
         // look at v-hits at the upstream layer
         //
         ///////////////////////////////////////////////

         const Double_t distance_max = 1. + 5.*vertex->r_*(geometry_->ut_[layer+1]-geometry_->ut_[layer])/(geometry_->ut_[layer]-vertex->u_);
         //--debug-- if (debug_) cout<< "GenerateTracksV: distance_max = " << distance_max <<endl;

         //-- std::list<const SensorHit*>::const_iterator it0 = vlist[layer].begin();   // erase with const_iterator is an error in SL6
         std::list<const SensorHit*>::iterator it0 = vlist[layer].begin();
         while (it0!=vlist[layer].end())
         {
            const SensorHit* hit0 = *it0;
            bool erase_hit = false;                         // hit has been assigned to Track2D and need to be erased from tlist

            Track2D vertex_track(&vertexHit_v, hit0);

            // project the vertex track to the next layer
            Double_t pos = vertex_track.at(geometry_->uv_[layer+1]);        // v-position to look for the hit in the next layer

            // find the closest hit in the next layer
            Double_t distance = distance_max;
            const SensorHit* closestHit = 0;

            for (std::list<const SensorHit*>::const_iterator it1=vlist[layer+1].begin(); it1!=vlist[layer+1].end(); ++it1) {
               const SensorHit* hit1 = *it1;
               Double_t d = hit1->pos_ - pos;
               if (TMath::Abs(d) < distance) {
                  closestHit = hit1;
                  distance = d;
               }
            }
            if (closestHit) vlist[layer+1].remove(closestHit);       // remove matched hit from the vlist[layer+1]

            if (!closestHit) {
               // there is no hit for this vertex track in the next layer. One of two possibilities:
               // 1) this is a noise hit
               // 2) the proton hit the gap in the next layer

               // check the gaps in the NEXT layer

               // loop over tTracks to find is some of them hit a gap
               for (std::list<const Track2D*>::const_iterator it=tTracks.begin(); it!=tTracks.end(); ++it) {
                  const Track2D* t_track = *it;
                  Double_t t_pos = t_track->at(geometry_->uv_[layer+1]);   // t-position to look for the gap
                  bool found_gap = false;
                  for (int igap=0; igap<3; ++igap) {
                     if (TMath::Abs(t_pos - pCTSensors_->gap_.vgap_[layer+1][igap]) < pCTSensors_->gap_.width_) {
                        // recover this hit: place a hit in the middle of the gap
                        Int_t sensorId = -1*(100 + layer*10 + (layer+1));
                        SensorHit* gapHit = new ((*Sensor::poolSensorHit_)[Sensor::poolSensorHit_->GetLast()+1]) SensorHit(sensorId, 0, 0, geometry_->uv_[layer+1], pos);
                        closestHit = gapHit;
                        found_gap = true;
                        vreco_[layer+1] = pos;
                        if (debug_) cout<< "\n--> event_ = " << event_ << " recovered v-hit in downstream layer " << layer+1 << " t_pos - pCTSensors_->gap_.vgap_[layer+1][igap] = " << t_pos - pCTSensors_->gap_.vgap_[layer+1][igap] <<endl;
                        //
                        // TODO: mark this pair of t- and v- tracks
                        //
                        break;
                     }
                  }
                  if (found_gap) break;
               }
            }

            if (closestHit) {
               //--debug-- if (debug_) cout<< "found good hit in the downstream layer: erase it from the vlist[" << layer << "]" << endl;
               erase_hit = true;                                  // erase the hit at the end of the loop body
               //
               // 1) create Track2D and add it to tTrack2D
               // 2) remove hits from the tlist
               //
               Track2D* track2D = new ((*poolTrack2D_)[poolTrack2D_->GetLast()+1]) Track2D(hit0, closestHit);
               vTracks.push_back(track2D);
            }

            //
            //-- erase the paired hit or increament the iterator
            //
            if (erase_hit) {
               it0 = vlist[layer].erase(it0);
            }
            else ++it0;
         }  // loop over hits in the upstream layer

         ///////////////////////////////////////////////
         //
         // look at v-hits at the downstream layer (if any)
         //
         ///////////////////////////////////////////////

         layer += 1;

         //-- std::list<const SensorHit*>::const_iterator it1 = vlist[layer].begin(); // erase with const_iterator is an error in SL6
         std::list<const SensorHit*>::iterator it1 = vlist[layer].begin();
         while (it1!=vlist[layer].end())
         {
            const SensorHit* hit1 = *it1;
            bool erase_hit = false;                         // hit has been assigned to Track2D and need to be erased from tlist

            Track2D vertex_track(&vertexHit_v, hit1);

            // project the vertex track to the next layer
            Double_t pos = vertex_track.at(geometry_->uv_[layer-1]);

            // find the closest hit in the next layer
            Double_t distance = distance_max;
            const SensorHit* closestHit = 0;

            for (it0=vlist[layer-1].begin(); it0!=vlist[layer-1].end(); ++it0) {       // reuse of std::list<const SensorHit*>::const_iterator it0 
               const SensorHit* hit0 = *it0;
               Double_t d = hit0->pos_ - pos;
               if (TMath::Abs(d) < distance) {
                  closestHit = hit0;
                  distance = d;
               }
            }
            if (closestHit) vlist[layer-1].remove(closestHit);       // remove matched hit from the tlist[layer+1]

            if (!closestHit) {
               // there is no hit for this vertex track in the upstream layer. One of two possibilities:
               // 1) this is a noise hit
               // 2) the proton hit the gap in the next layer

               // check the gaps

               // loop over tTracks to find is some of them hit a gap
               for (std::list<const Track2D*>::const_iterator it=tTracks.begin(); it!=tTracks.end(); ++it) {
                  const Track2D* t_track = *it;
                  Double_t t_pos = t_track->at(geometry_->uv_[layer-1]);   // t-position to look for the gap
                  bool found_gap = false;
                  for (int igap=0; igap<3; ++igap) {
                     if (TMath::Abs(t_pos - pCTSensors_->gap_.vgap_[layer-1][igap]) < pCTSensors_->gap_.width_) {
                        // recover this hit: place a hit in the middle of the gap
                        Int_t sensorId = -1*(100 + layer*10 + (layer-1));
                        SensorHit* gapHit = new ((*Sensor::poolSensorHit_)[Sensor::poolSensorHit_->GetLast()+1]) SensorHit(sensorId, 0, 0, geometry_->uv_[layer-1], pos);
                        closestHit = gapHit;
                        found_gap = true;
                        vreco_[layer-1] = pos;
                        if (debug_) cout<< "\n--> event_ = " << event_ << " recovered v-hit in upstream layer " << layer-1 << " t_pos - pCTSensors_->gap_.vgap_[layer-1][igap] = " << t_pos - pCTSensors_->gap_.vgap_[layer-1][igap] <<endl;
                        //
                        // TODO: mark this pair of t- and v- tracks
                        //
                        break;
                     }
                  }
                  if (found_gap) break;
               }
            }

            if (closestHit) {
               //--debug-- if (debug_) cout<< "found good hit in the upstream layer: erase it from the vlist[" << layer << "]" << endl;
               erase_hit = true;                                    // erase the hit at the end of the loop body
               //
               // 1) create Track2D and add it to tTrack2D
               // 2) remove hits from the tlist
               //
               Track2D* track2D = new ((*poolTrack2D_)[poolTrack2D_->GetLast()+1]) Track2D(closestHit, hit1);
               vTracks.push_back(track2D);
            }

            //
            //-- erase the paired hit or increament the iterator
            //
            if (erase_hit) {
               it1 = vlist[layer].erase(it1);
            }
            else ++it1;
         }  // loop over hits in the upstream layer

         if (debug_) {
            //--debug-- cout<< "vTracks.size() = " << vTracks.size() <<endl;
         }
      }
   }
   ~Reco() {
      if (poolTrack2D_) poolTrack2D_->Clear();
      if (poolTrack_) poolTrack_->Clear();
      if (poolSuperTrack_) poolSuperTrack_->Clear();
   }
   void GenerateTracks()         //------------- new algorithm --------------//
   {
      if (debug_) cout<< "GenerateTracks: vin_.size() = " << vin_.size() << " tin_.size() = " << tin_.size() << " vout_.size() = " << vout_.size() << " tout_.size() = " << tout_.size() <<endl;

      // input tracks
      for (std::list<const Track2D*>::const_iterator itt=tin_.begin(); itt!=tin_.end(); ++itt) {
         const Track2D* tTrack = *itt;
         //if (debug_) cout<< "tTrack = " << tTrack <<endl;
         //if (debug_) cout<< "*tTrack->hit1_ = " << *tTrack->hit1_ << " *tTrack->hit2_ = " << *tTrack->hit2_ <<endl;
         for (std::list<const Track2D*>::const_iterator itv=vin_.begin(); itv!=vin_.end(); ++itv) {
            const Track2D* vTrack = *itv;
            //if (debug_) cout<< "vTrack = " << vTrack <<endl;
            //if (debug_) cout<< "*vTrack->hit1_ = " << *vTrack->hit1_ << " *vTrack->hit2_ = " << *vTrack->hit2_ <<endl;
            // TODO: check the proper v-half
            const Track* track = new ((*poolTrack_)[poolTrack_->GetLast()+1]) Track(vTrack, tTrack);
            itracks_.push_back(track);
         }
      }
      if (debug_) cout<< "itracks_.size() = " << itracks_.size() <<endl;

      // output tracks
      for (std::list<const Track2D*>::const_iterator itt=tout_.begin(); itt!=tout_.end(); ++itt) {
         const Track2D* tTrack = *itt;
         for (std::list<const Track2D*>::const_iterator itv=vout_.begin(); itv!=vout_.end(); ++itv) {
            const Track2D* vTrack = *itv;
            // TODO: check the proper v-half
            const Track* track = new ((*poolTrack_)[poolTrack_->GetLast()+1]) Track(vTrack, tTrack);
            otracks_.push_back(track);
         }
      }
      if (debug_) cout<< "otracks_.size() = " << otracks_.size() <<endl;
   }
   void GenerateSuperTracks()
   {
      for (std::list<const Track*>::const_iterator it=itracks_.begin(); it!=itracks_.end(); ++it) {
         const Track* itrack = *it;
         for (std::list<const Track*>::const_iterator ot=otracks_.begin(); ot!=otracks_.end(); ++ot) {
            const Track* otrack = *ot;
            SuperTrack* superTrack = new ((*poolSuperTrack_)[poolSuperTrack_->GetLast()+1]) SuperTrack(itrack, otrack);
            superTracks_.push_back(superTrack);
         }
      }
   }
   void RemoveDuplicates(std::multimap<Double_t, const SuperTrack*>& multimap)
   {
      // remove supertracks down the multimap which share the same hits
      for (std::multimap<Double_t, const SuperTrack*>::iterator it=multimap.begin(); it!=multimap.end(); ++it)
      {
         const SuperTrack* it_superTrack = it->second;
         if (debug_) cout<< "it_superTrack = " << it_superTrack <<endl;
         std::multimap<Double_t, const SuperTrack*>::iterator next = it;   // cannot use const_iterator with old compilers like SL6
         ++next;
         while (next != multimap.end()) {
            const SuperTrack* next_superTrack = next->second;
            if (debug_) cout<< "   *next_superTrack = " << next_superTrack <<endl;
            if (next->second->SharedHits(it->second)) {
               if (debug_) cout<< "   removed" <<endl;
               superTracks_.remove(next->second);
               multimap.erase(next++);
            }
            else ++next;
         }
      }
   }
   void FilterAngle(Double_t angle_max=3.14159265359)
   {
      if (debug_) cout<< "FilterAngle: superTracks_.size() = " << superTracks_.size() <<endl;

      if (superTracks_.size() < 2) return;

      std::multimap<Double_t, const SuperTrack*> mapAngle;

      for (std::list<const SuperTrack*>::const_iterator it=superTracks_.begin(); it!=superTracks_.end(); ++it) {
         const SuperTrack* superTrack = *it;
         if (debug_) if (mapAngle.find(superTrack->angle) != mapAngle.end()) cout<< "found: mapAngle[" << superTrack->angle << "] = " << superTrack << "\t " << *superTrack <<endl;
         //mapAngle[superTrack->angle] = superTrack;
         mapAngle.insert(std::make_pair(superTrack->angle, superTrack));
      }

      if (debug_) {
         cout<< "angle multimap mapAngle.size() = " << mapAngle.size() <<endl;
         for (std::multimap<Double_t, const SuperTrack*>::const_iterator it=mapAngle.begin(); it!=mapAngle.end(); ++it) {
            Double_t angle = it->first;
            cout<< std::distance<std::multimap<Double_t, const SuperTrack*>::const_iterator>(mapAngle.begin(), it) << "\t angle = " << angle <<endl;
         }
      }

      if (debug_) cout<< "loop over the multimap to remove tracks with angle above the angle_max" <<endl;

      // clear list superTracks_ and fill it in accordance with the angle
      superTracks_.clear();
      for (std::multimap<Double_t, const SuperTrack*>::iterator it=mapAngle.begin(); it!=mapAngle.end(); ++it)  // cannot use const_iterator with old compilers
      {
         if (it->first > angle_max) break;
         superTracks_.push_back(it->second);
      }

      if (debug_) cout<< "FilterAngle: after angle filter: superTracks_.size() = " << superTracks_.size() <<endl;

      if (debug_) for (std::list<const SuperTrack*>::const_iterator it=superTracks_.begin(); it!=superTracks_.end(); ++it) {
         const SuperTrack* superTrack = *it;
         cout<< std::distance<std::list<const SuperTrack*>::const_iterator>(superTracks_.begin(), it) << " " << superTrack << " angle = " << superTrack->angle << " " << *superTrack <<endl;
      }

      // make sure that the rest of tracks do not share the same hits
      RemoveDuplicates(mapAngle);

      if (debug_) cout<< "FilterAngle: GenerateSuperTracks after shared hits filter: superTracks_.size() = " << superTracks_.size() <<endl<<endl;

      if (debug_) {
         for (std::list<const SuperTrack*>::const_iterator it=superTracks_.begin(); it!=superTracks_.end(); ++it) {
            const SuperTrack* superTrack = *it;
            int ntrack = std::distance<std::list<const SuperTrack*>::const_iterator>(superTracks_.begin(), it);
            cout<< ntrack << "\t superTrack->angle = " << superTrack->angle << " " << *superTrack <<endl;
         }
      }
   }
   void FilterDistance(Double_t rmax=10.)
   {
      if (debug_) cout<< "FilterDistance: superTracks_.size() = " << superTracks_.size() <<endl;

      if (superTracks_.size() < 2) return;

      // apply filter on the distance between the hits in the plane u = 0

      std::multimap<Double_t, const SuperTrack*> mapCloseTracks;            // for the filter on the distance in the plane u = 0

      for (std::list<const SuperTrack*>::const_iterator it=superTracks_.begin(); it!=superTracks_.end(); ++it) {
         const SuperTrack* superTrack = *it;
         //mapCloseTracks[superTrack->Distance()] = superTrack;
         mapCloseTracks.insert(std::make_pair(superTrack->Distance(), superTrack));
      }

      if (debug_) {
         cout<< "resulting distance multimap from SuperTracks" <<endl;
         for (std::multimap<Double_t, const SuperTrack*>::const_iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end(); ++it) {
            Double_t distance = it->first;
            cout<< std::distance<std::multimap<Double_t, const SuperTrack*>::const_iterator>(mapCloseTracks.begin(), it) << "\t distance = " << distance <<endl;
         }
      }

      if (debug_) cout<< "loop over the multimap to remove tracks with distance above the r" <<endl;

      // clear list superTracks_ and fill it in accordance with the distance
      superTracks_.clear();
      for (std::multimap<Double_t, const SuperTrack*>::iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end(); ++it)  // cannot use const_iterator with old compilers
      {
         if (it->first > rmax) break;
         superTracks_.push_back(it->second);
      }

      if (debug_) cout<< "FilterDistance: after distance filter: superTracks_.size() = " << superTracks_.size() <<endl;

      if (debug_) for (std::list<const SuperTrack*>::const_iterator it=superTracks_.begin(); it!=superTracks_.end(); ++it) {
         const SuperTrack* superTrack = *it;
         cout<< std::distance<std::list<const SuperTrack*>::const_iterator>(superTracks_.begin(), it) << " " << superTrack << " distance = " << superTrack->Distance() << " " << *superTrack <<endl;
      }

      // make sure that the rest of tracks do not share the same hits
      RemoveDuplicates(mapCloseTracks);

      if (debug_) cout<< "FilterDistance: GenerateSuperTracks after shared hits filter: superTracks_.size() = " << superTracks_.size() <<endl<<endl;

      if (debug_) {
         for (std::list<const SuperTrack*>::const_iterator it=superTracks_.begin(); it!=superTracks_.end(); ++it) {
            const SuperTrack* superTrack = *it;
            int ntrack = std::distance<std::list<const SuperTrack*>::const_iterator>(superTracks_.begin(), it);
            cout<< ntrack << "\t superTrack->angle = " << superTrack->angle << " " << *superTrack <<endl;
            //const Track* itrack = superTrack->itrack_;
            //const Track* otrack = superTrack->otrack_;
            //cout<< "itrack->cv_ = " << itrack->cv_ << " itrack->ct_ = " << itrack->ct_ << " itrack->cu_ = " << itrack->cu_ <<endl;
            //cout<< "otrack->cv_ = " << otrack->cv_ << " otrack->ct_ = " << otrack->ct_ << " otrack->cu_ = " << otrack->cu_ <<endl;
         }
      }
   }
   void FilterPhaseSpace(Double_t rmax=10., Double_t angle_max=0.1)
   {
      // cut on the element of the phase space angle*distance together with reasonable cuts on the distance and the angle

      if (debug_) cout<< "FilterPhaseSpace: superTracks_.size() = " << superTracks_.size() <<endl;

      if (superTracks_.size() < 2) return;

      std::multimap<Double_t, const SuperTrack*> mapAngle;

      for (std::list<const SuperTrack*>::const_iterator it=superTracks_.begin(); it!=superTracks_.end(); ++it) {
         const SuperTrack* superTrack = *it;
         if (debug_) if (mapAngle.find(superTrack->angle) != mapAngle.end()) cout<< "found: mapAngle[" << superTrack->angle << "] = " << superTrack << "\t " << *superTrack <<endl;
         //-- if (superTrack->d < rmax) mapAngle.insert(std::make_pair(superTrack->angle, superTrack));
         Double_t dangle = superTrack->angle*superTrack->d;
         if (superTrack->d < rmax && superTrack->angle < angle_max) mapAngle.insert(std::make_pair(dangle, superTrack));
      }

      if (debug_) {
         cout<< "angle multimap mapAngle.size() = " << mapAngle.size() <<endl;
         for (std::multimap<Double_t, const SuperTrack*>::const_iterator it=mapAngle.begin(); it!=mapAngle.end(); ++it) {
            //-- Double_t angle = it->first;
            Double_t dangle = it->first;
            cout<< std::distance<std::multimap<Double_t, const SuperTrack*>::const_iterator>(mapAngle.begin(), it) << "\t dangle = " << dangle <<endl;
         }
      }

      if (debug_) cout<< "loop over the multimap to remove tracks with angle above the angle_max" <<endl;

      // clear list superTracks_ and fill it in accordance with the angle
      superTracks_.clear();
      for (std::multimap<Double_t, const SuperTrack*>::iterator it=mapAngle.begin(); it!=mapAngle.end(); ++it)  // cannot use const_iterator with old compilers
      {
         if (it->first > angle_max) break;
         superTracks_.push_back(it->second);
      }

      if (debug_) cout<< "FilterPhaseSpace: after angle filter: superTracks_.size() = " << superTracks_.size() <<endl;

      if (debug_) for (std::list<const SuperTrack*>::const_iterator it=superTracks_.begin(); it!=superTracks_.end(); ++it) {
         const SuperTrack* superTrack = *it;
         cout<< std::distance<std::list<const SuperTrack*>::const_iterator>(superTracks_.begin(), it) << " angle*d = " << superTrack->angle*superTrack->d << "\t" << *superTrack <<endl;
      }

      // make sure that the rest of tracks do not share the same hits
      RemoveDuplicates(mapAngle);

      if (debug_) cout<< "FilterPhaseSpace: GenerateSuperTracks after shared hits filter: superTracks_.size() = " << superTracks_.size() <<endl<<endl;

      if (debug_) {
         for (std::list<const SuperTrack*>::const_iterator it=superTracks_.begin(); it!=superTracks_.end(); ++it) {
            const SuperTrack* superTrack = *it;
            int ntrack = std::distance<std::list<const SuperTrack*>::const_iterator>(superTracks_.begin(), it);
            cout<< ntrack << "\t superTrack->angle*superTrack->d = " << superTrack->angle*superTrack->d << " " << *superTrack <<endl;
         }
      }
   }
};

class RecoEvent: public TObject {
public:
   Bool_t ok;                       // error flag
   Float_t deltaT;
   TClonesArray* track;             //-> 
   Int_t nt;                        // the number of super tracks
   Float_t a[5];                    // Energy detector channels
   Float_t ped[5];
   Float_t sample[5][16];           // to plot e.g. channel 1:    r->Draw("sample[1][]:Iteration$","Entry$==0") 
   Float_t wepl;
   Int_t nhitv[4];
   Int_t nhitt[4];
   Float_t vreco[4];
   Float_t treco[4];
   RecoEvent(): TObject(), ok(kTRUE), nt(0) {
      //cout<< "ctor RecoEvent: create the track" <<endl;
      track = new TClonesArray("SuperTrack");
      //cout<< "RecoEvent: created" <<endl;
   }
   virtual ~RecoEvent() {
      //cout<< "dtor RecoEvent" <<endl;
      clear();
      delete track;
   }
   void clear() {
      //cout<< "RecoEvent::clear" <<endl;
      deltaT = 0;
      track->Clear("C");            // initiate call of the method SuperTrack::Clear
      nt = 0;
      for (int i=0; i<5; ++i) {
         a[i] = 0;
         ped[i] = 0;
      }
      for (int ilayer=0; ilayer<4; ++ilayer) {
         vreco[ilayer] = 0;
         treco[ilayer] = 0;
      }
      for (int ichan=0; ichan<5; ++ichan) for (int isample=0; isample<16; ++isample) sample[ichan][isample] = 0;
      wepl = 0;
   }
   Float_t SampleSum(Int_t chan, Int_t nfront, Int_t ntail, Double_t pedestal) const {
      if (chan < 0 || chan > 4) return 0;
      // find a position of the maximum
      Int_t imax = 0;
      // assumes that the number of samples is 16
      Float_t sum = sample[chan][imax];
      for (int isample=0; isample<16; ++isample) {
         sum += sample[chan][isample];
         if (sample[chan][isample] > sample[chan][imax]) {
            imax = isample;
         }
      }
      if (sum == 0) return 0;          // there are no samples in this event
      Int_t n1 = imax - nfront;
      if (n1 < 0) n1 = 0;
      Int_t n2 = imax + ntail;
      if (n2 > 15) n2 = 15;

      sum = 0;
      for (int isample=n1; isample<=n2; ++isample) sum += sample[chan][isample];
      sum -= (n2 - n1 + 1)*pedestal;

      return sum;
   }
   Int_t Nvreco() const {
      Int_t nrecovered = 0;
      for (int ilayer=0; ilayer<4; ++ilayer) if (vreco[ilayer] != 0) ++nrecovered;
      return nrecovered;
   }
   Int_t Ntreco() const {
      Int_t nrecovered = 0;
      for (int ilayer=0; ilayer<4; ++ilayer) if (treco[ilayer] != 0) ++nrecovered;
      return nrecovered;
   }
   void SetWEPL(Double_t the_wepl) {wepl = the_wepl;}
   void Extract(const Reco& reco) {
      //cout<< "RecoEvent::Extract" <<endl;
      clear();
      deltaT = reco.deltaT_;
      for (int ilayer=0; ilayer<4; ++ilayer) {
         vreco[ilayer] = reco.vreco_[ilayer];
         treco[ilayer] = reco.treco_[ilayer];
      }
      for (std::list<const SuperTrack*>::const_iterator it=reco.superTracks_.begin(); it!=reco.superTracks_.end(); ++it) {
         const SuperTrack* superTrack = *it;
         //cout<< "Extract: add SuperTrack into the TClonesArray" <<endl;
         new ((*track)[track->GetLast()+1]) SuperTrack(*superTrack);
         ++nt;
      }
      for (int ilayer=0; ilayer<4; ++ilayer) {
         nhitv[ilayer] = reco.nhitv_[ilayer];
         nhitt[ilayer] = reco.nhitt_[ilayer];
      }
      // Energy detector
      //a[0] = reco.pCTEvent_->energyBoard[0].pulse[0];
      //a[1] = reco.pCTEvent_->energyBoard[0].pulse[1];
      //a[2] = reco.pCTEvent_->energyBoard[0].pulse[2];
      //a[3] = reco.pCTEvent_->energyBoard[1].pulse[0];
      //a[4] = reco.pCTEvent_->energyBoard[1].pulse[1];

      /////////////////////////////
      int brd=0;
      ///// //--orig int enrgTag0= thisEvent->Event->Board[0].enrgTag;
      ///// //--orig int enrgTag1= thisEvent->Event->Board[1].enrgTag;
      ///// //		if (enrgTag0 != enrgTag1) cout << "enrg tag mismatch: " << enrgTag0 << " vs " << enrgTag1;
      if (reco.pCTEvent_->energyBoard[brd].numChan>0 || reco.pCTEvent_->energyBoard[brd].numSamples>0) {
        /////   //			if (thisEvent->Event->Board[brd].enrgTag != evtNum % 4) cout << "tag mismatch, energy tag=" << thisEvent->Event->Board[brd].enrgTag << "\n";
        if (reco.pCTEvent_->energyBoard[brd].reduced) {
          a[0] = reco.pCTEvent_->energyBoard[brd].pulse[0];  ped[0] = reco.pCTEvent_->energyBoard[brd].pedestal[0];
          a[1] = reco.pCTEvent_->energyBoard[brd].pulse[1];	ped[1] = reco.pCTEvent_->energyBoard[brd].pedestal[1];
          a[2] =	reco.pCTEvent_->energyBoard[brd].pulse[2];	ped[2] = reco.pCTEvent_->energyBoard[brd].pedestal[2];
        } else {
          a[0] = 0.; a[1] = 0.; a[2] = 0.;
          /// enrgSamp *thisSamp = thisEvent->Event->Board[brd].firstSample;
          //EnergySample* energySample = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(0);
          //if (thisSamp != 0)
          if (reco.pCTEvent_->energyBoard[brd].samples->GetLast()+1 > 0)
          {
            EnergySample* energySample0 = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(0);
            sample[0][0] = energySample0->pulse[0];   // assign the first sample of the RecoEvent::sample for the first energy board
            sample[1][0] = energySample0->pulse[1];
            sample[2][0] = energySample0->pulse[2];
            int ped0 = energySample0->pulse[0]; ped[0] = ped0;
            int ped1 = energySample0->pulse[1]; ped[1] = ped1;
            int ped2 = energySample0->pulse[2]; ped[2] = ped2;
            //while (thisSamp != 0)
            for (int isample=1; isample<reco.pCTEvent_->energyBoard[brd].samples->GetLast()+1; ++isample)
            {
              EnergySample* energySample = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(isample);
              sample[0][isample] = energySample->pulse[0];  // assign the rest of 16 samples of the RecoEvent::sample for the first energy board
              sample[1][isample] = energySample->pulse[1];
              sample[2][isample] = energySample->pulse[2];
              a[0] = a[0] + energySample->pulse[0] - ped0;
              a[1] = a[1] + energySample->pulse[1] - ped1;
              a[2] = a[2] + energySample->pulse[2] - ped2;
            }
          }
        }
      } else {
        a[0] = 0; ped[0] = 0;
        a[1] = 0; ped[1] = 0;
        a[2] = 0; ped[2] = 0;
      }
      brd=1;
      if (reco.pCTEvent_->energyBoard[brd].numChan>0 || reco.pCTEvent_->energyBoard[brd].numSamples>0) {
        /////   //			if (thisEvent->Event->Board[brd].enrgTag != evtNum % 4) cout << "tag mismatch, energy tag=" << thisEvent->Event->Board[brd].enrgTag << "\n";
        if (reco.pCTEvent_->energyBoard[brd].reduced) {
          a[3] = reco.pCTEvent_->energyBoard[brd].pulse[0];  ped[3] = reco.pCTEvent_->energyBoard[brd].pedestal[0];
          a[4] = reco.pCTEvent_->energyBoard[brd].pulse[1];	ped[4] = reco.pCTEvent_->energyBoard[brd].pedestal[1];
          //--no such channel-- PhCh5 =	reco.pCTEvent_->energyBoard[brd].pulse[2];
        } else {
          a[3] = 0.; a[4] = 0.; //--no such channel-- PhCh5 = 0.;
          int ped3 = 0; int ped4 = 0; int ped5 = 0;
          /// enrgSamp *thisSamp = thisEvent->Event->Board[brd].firstSample;
          //EnergySample* energySample = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(0);
          //if (thisSamp != 0)
          if (reco.pCTEvent_->energyBoard[brd].samples->GetLast()+1 > 0)
          {
            EnergySample* energySample0 = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(0);
            sample[3][0] = energySample0->pulse[0];      // assign the first sample of the RecoEvent::sample for the second energy board
            sample[4][0] = energySample0->pulse[1];
            ped3 = energySample0->pulse[0]; ped[3] = ped3;
            ped4 = energySample0->pulse[1]; ped[4] = ped4;
            ped5 = energySample0->pulse[2];
            //while (thisSamp != 0)
            for (int isample=1; isample<reco.pCTEvent_->energyBoard[brd].samples->GetLast()+1; ++isample)
            {
              EnergySample* energySample = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(isample);
              sample[3][isample] = energySample->pulse[0];  // assign the rest of 16 samples of the RecoEvent::sample for the second energy board
              sample[4][isample] = energySample->pulse[1];
              a[3] = a[3] + energySample->pulse[0] - ped3;
              a[4] = a[4] + energySample->pulse[1] - ped4;
              //--no such channel-- PhCh5 = PhCh5 + energySample->pulse[2] - ped5;
            }
          }
        }
      } else {
        a[3] = 0; ped[3] = 0;
        a[4] = 0; ped[4] = 0;
        //--no such channel-- PhCh5 = 0;
      }
      /////////////////////////////
   }

   ClassDef(RecoEvent, 11);
};

#ifdef __MAKECINT__
#pragma link C++ class RecoEvent;
#endif

#endif  // Reco_h
