/*
 * \file EBPedestalOnlineClient.cc
 *
 * $Date: 2010/03/27 20:07:57 $
 * $Revision: 1.155 $
 * \author G. Della Ricca
 * \author F. Cossutti
 *
*/

#include <memory>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DQMServices/Core/interface/DQMStore.h"

#ifdef WITH_ECAL_COND_DB
#include "OnlineDB/EcalCondDB/interface/MonPedestalsOnlineDat.h"
#include "OnlineDB/EcalCondDB/interface/RunCrystalErrorsDat.h"
#include "OnlineDB/EcalCondDB/interface/RunTTErrorsDat.h"
#include "OnlineDB/EcalCondDB/interface/EcalCondDBInterface.h"
#endif

#include "CondTools/Ecal/interface/EcalErrorDictionary.h"

#include "DQM/EcalCommon/interface/EcalErrorMask.h"
#include "DQM/EcalCommon/interface/UtilsClient.h"
#include "DQM/EcalCommon/interface/LogicID.h"
#include "DQM/EcalCommon/interface/Numbers.h"

#include <DQM/EcalBarrelMonitorClient/interface/EBPedestalOnlineClient.h>

EBPedestalOnlineClient::EBPedestalOnlineClient(const edm::ParameterSet& ps) {

  // cloneME switch
  cloneME_ = ps.getUntrackedParameter<bool>("cloneME", true);

  // verbose switch
  verbose_ = ps.getUntrackedParameter<bool>("verbose", true);

  // debug switch
  debug_ = ps.getUntrackedParameter<bool>("debug", false);

  // prefixME path
  prefixME_ = ps.getUntrackedParameter<std::string>("prefixME", "");

  // enableCleanup_ switch
  enableCleanup_ = ps.getUntrackedParameter<bool>("enableCleanup", false);

  // vector of selected Super Modules (Defaults to all 36).
  superModules_.reserve(36);
  for ( unsigned int i = 1; i <= 36; i++ ) superModules_.push_back(i);
  superModules_ = ps.getUntrackedParameter<std::vector<int> >("superModules", superModules_);

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    h03_[ism-1] = 0;

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    meg03_[ism-1] = 0;

    mep03_[ism-1] = 0;

    mer03_[ism-1] = 0;

  }

  expectedMean_ = 200.0;
  discrepancyMean_ = 25.0;
  RMSThreshold_ = 2.0;

}

EBPedestalOnlineClient::~EBPedestalOnlineClient() {

}

void EBPedestalOnlineClient::beginJob(void) {

  dqmStore_ = edm::Service<DQMStore>().operator->();

  if ( debug_ ) std::cout << "EBPedestalOnlineClient: beginJob" << std::endl;

  ievt_ = 0;
  jevt_ = 0;

}

void EBPedestalOnlineClient::beginRun(void) {

  if ( debug_ ) std::cout << "EBPedestalOnlineClient: beginRun" << std::endl;

  jevt_ = 0;

  this->setup();

}

void EBPedestalOnlineClient::endJob(void) {

  if ( debug_ ) std::cout << "EBPedestalOnlineClient: endJob, ievt = " << ievt_ << std::endl;

  this->cleanup();

}

void EBPedestalOnlineClient::endRun(void) {

  if ( debug_ ) std::cout << "EBPedestalOnlineClient: endRun, jevt = " << jevt_ << std::endl;

  this->cleanup();

}

void EBPedestalOnlineClient::setup(void) {

  char histo[200];

  dqmStore_->setCurrentFolder( prefixME_ + "/EBPedestalOnlineClient" );

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( meg03_[ism-1] ) dqmStore_->removeElement( meg03_[ism-1]->getName() );
    sprintf(histo, "EBPOT pedestal quality G12 %s", Numbers::sEB(ism).c_str());
    meg03_[ism-1] = dqmStore_->book2D(histo, histo, 85, 0., 85., 20, 0., 20.);
    meg03_[ism-1]->setAxisTitle("ieta", 1);
    meg03_[ism-1]->setAxisTitle("iphi", 2);

    if ( mep03_[ism-1] ) dqmStore_->removeElement( mep03_[ism-1]->getName() );
    sprintf(histo, "EBPOT pedestal mean G12 %s", Numbers::sEB(ism).c_str());
    mep03_[ism-1] = dqmStore_->book1D(histo, histo, 100, 150., 250.);
    mep03_[ism-1]->setAxisTitle("mean", 1);

    if ( mer03_[ism-1] ) dqmStore_->removeElement( mer03_[ism-1]->getName() );
    sprintf(histo, "EBPOT pedestal rms G12 %s", Numbers::sEB(ism).c_str());
    mer03_[ism-1] = dqmStore_->book1D(histo, histo, 100, 0.,  10.);
    mer03_[ism-1]->setAxisTitle("rms", 1);

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( meg03_[ism-1] ) meg03_[ism-1]->Reset();

    for ( int ie = 1; ie <= 85; ie++ ) {
      for ( int ip = 1; ip <= 20; ip++ ) {

        meg03_[ism-1]->setBinContent( ie, ip, 2. );

      }
    }

    if ( mep03_[ism-1] ) mep03_[ism-1]->Reset();
    if ( mer03_[ism-1] ) mer03_[ism-1]->Reset();

  }

}

void EBPedestalOnlineClient::cleanup(void) {

  if ( ! enableCleanup_ ) return;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( cloneME_ ) {
      if ( h03_[ism-1] ) delete h03_[ism-1];
    }

    h03_[ism-1] = 0;

  }

  dqmStore_->setCurrentFolder( prefixME_ + "/EBPedestalOnlineClient" );

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( meg03_[ism-1] ) dqmStore_->removeElement( meg03_[ism-1]->getName() );
    meg03_[ism-1] = 0;

    if ( mep03_[ism-1] ) dqmStore_->removeElement( mep03_[ism-1]->getName() );
    mep03_[ism-1] = 0;

    if ( mer03_[ism-1] ) dqmStore_->removeElement( mer03_[ism-1]->getName() );
    mer03_[ism-1] = 0;

  }

}

#ifdef WITH_ECAL_COND_DB
bool EBPedestalOnlineClient::writeDb(EcalCondDBInterface* econn, RunIOV* runiov, MonRunIOV* moniov, bool& status) {

  status = true;

  EcalLogicID ecid;

  MonPedestalsOnlineDat p;
  map<EcalLogicID, MonPedestalsOnlineDat> dataset;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( verbose_ ) {
      std::cout << " " << Numbers::sEB(ism) << " (ism=" << ism << ")" << std::endl;
      std::cout << std::endl;
      UtilsClient::printBadChannels(meg03_[ism-1], h03_[ism-1]);
    }

    float num03;
    float mean03;
    float rms03;

    for ( int ie = 1; ie <= 85; ie++ ) {
      for ( int ip = 1; ip <= 20; ip++ ) {

        bool update03;

        update03 = UtilsClient::getBinStatistics(h03_[ism-1], ie, ip, num03, mean03, rms03);

        if ( update03 ) {

          if ( Numbers::icEB(ism, ie, ip) == 1 ) {

            if ( verbose_ ) {
              std::cout << "Preparing dataset for " << Numbers::sEB(ism) << " (ism=" << ism << ")" << std::endl;
              std::cout << "G12 (" << ie << "," << ip << ") " << num03  << " " << mean03 << " " << rms03  << std::endl;
              std::cout << std::endl;
            }

          }

          p.setADCMeanG12(mean03);
          p.setADCRMSG12(rms03);

          if ( UtilsClient::getBinStatus(meg03_[ism-1], ie, ip) ) {
            p.setTaskStatus(true);
          } else {
            p.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQuality(meg03_[ism-1], ie, ip);

          int ic = Numbers::indexEB(ism, ie, ip);

          if ( econn ) {
            ecid = LogicID::getEcalLogicID("EB_crystal_number", Numbers::iSM(ism, EcalBarrel), ic);
            dataset[ecid] = p;
          }

        }

      }
    }

  }

  if ( econn ) {
    try {
      if ( verbose_ ) std::cout << "Inserting MonPedestalsOnlineDat ..." << std::endl;
      if ( dataset.size() != 0 ) econn->insertDataArraySet(&dataset, moniov);
      if ( verbose_ ) std::cout << "done." << std::endl;
    } catch (runtime_error &e) {
      cerr << e.what() << std::endl;
    }
  }

  return true;

}
#endif

void EBPedestalOnlineClient::analyze(void) {

  ievt_++;
  jevt_++;
  if ( ievt_ % 10 == 0 ) {
    if ( debug_ ) std::cout << "EBPedestalOnlineClient: ievt/jevt = " << ievt_ << "/" << jevt_ << std::endl;
  }

  uint64_t bits03 = 0;
  bits03 |= EcalErrorDictionary::getMask("PEDESTAL_ONLINE_HIGH_GAIN_MEAN_WARNING");
  bits03 |= EcalErrorDictionary::getMask("PEDESTAL_ONLINE_HIGH_GAIN_RMS_WARNING");
  bits03 |= EcalErrorDictionary::getMask("PEDESTAL_ONLINE_HIGH_GAIN_MEAN_ERROR");
  bits03 |= EcalErrorDictionary::getMask("PEDESTAL_ONLINE_HIGH_GAIN_RMS_ERROR");

  char histo[200];

  MonitorElement* me;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    sprintf(histo, (prefixME_ + "/EBPedestalOnlineTask/Gain12/EBPOT pedestal %s G12").c_str(), Numbers::sEB(ism).c_str());
    me = dqmStore_->get(histo);
    h03_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h03_[ism-1] );
    if ( meg03_[ism-1] ) meg03_[ism-1]->Reset();
    if ( mep03_[ism-1] ) mep03_[ism-1]->Reset();
    if ( mer03_[ism-1] ) mer03_[ism-1]->Reset();

    for ( int ie = 1; ie <= 85; ie++ ) {
      for ( int ip = 1; ip <= 20; ip++ ) {

        if ( meg03_[ism-1] ) meg03_[ism-1]->setBinContent(ie, ip, 2.);

        bool update03;

        float num03;
        float mean03;
        float rms03;

        update03 = UtilsClient::getBinStatistics(h03_[ism-1], ie, ip, num03, mean03, rms03);

        if ( update03 ) {

          float val;

          val = 1.;
          if ( fabs(mean03 - expectedMean_) > discrepancyMean_ )
            val = 0.;
          if ( rms03 > RMSThreshold_ )
            val = 0.;
          if ( meg03_[ism-1] ) meg03_[ism-1]->setBinContent(ie, ip, val);

          if ( mep03_[ism-1] ) mep03_[ism-1]->Fill(mean03);
          if ( mer03_[ism-1] ) mer03_[ism-1]->Fill(rms03);

        }

      }
    }

  }

#ifdef WITH_ECAL_COND_DB
  if ( EcalErrorMask::mapCrystalErrors_.size() != 0 ) {
    map<EcalLogicID, RunCrystalErrorsDat>::const_iterator m;
    for (m = EcalErrorMask::mapCrystalErrors_.begin(); m != EcalErrorMask::mapCrystalErrors_.end(); m++) {

      if ( (m->second).getErrorBits() & bits03 ) {
        EcalLogicID ecid = m->first;

        if ( strcmp(ecid.getMapsTo().c_str(), "EB_crystal_number") != 0 ) continue;

        int ism = Numbers::iSM(ecid.getID1(), EcalBarrel);
        std::vector<int>::iterator iter = find(superModules_.begin(), superModules_.end(), ism);
        if (iter == superModules_.end()) continue;

        int ic = ecid.getID2();
        int ie = (ic-1)/20 + 1;
        int ip = (ic-1)%20 + 1;

        UtilsClient::maskBinContent( meg03_[ism-1], ie, ip );

      }

    }
  }

  if ( EcalErrorMask::mapTTErrors_.size() != 0 ) {
    map<EcalLogicID, RunTTErrorsDat>::const_iterator m;
    for (m = EcalErrorMask::mapTTErrors_.begin(); m != EcalErrorMask::mapTTErrors_.end(); m++) {

      if ( (m->second).getErrorBits() & bits03 ) {
        EcalLogicID ecid = m->first;

        if ( strcmp(ecid.getMapsTo().c_str(), "EB_trigger_tower") != 0 ) continue;

        int ism = Numbers::iSM(ecid.getID1(), EcalBarrel);
        std::vector<int>::iterator iter = find(superModules_.begin(), superModules_.end(), ism);
        if (iter == superModules_.end()) continue;

        int itt = ecid.getID2();
        int iet = (itt-1)/4 + 1;
        int ipt = (itt-1)%4 + 1;

        for ( int ie = 5*(iet-1)+1; ie <= 5*iet; ie++ ) {
          for ( int ip = 5*(ipt-1)+1; ip <= 5*ipt; ip++ ) {
            UtilsClient::maskBinContent( meg03_[ism-1], ie, ip );
          }
        }

      }

    }
  }
#endif

}

