#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "nexus/NeXusFile.hpp"
#include <vector>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
using namespace std;
const uint32_t NX_MAX_DET =   2;
const uint32_t NX_MAX_PRD =   1;

#include "log.h"
#include "config.h"

const uint32_t MAX_TOF =3626 ;
//const uint32_t MAX_TOF =13 ;
const uint32_t MAX_DET =8000;
//const uint32_t MAX_DET =30;
const uint32_t BIN_DET =80;
double NxsMap[MAX_DET+5][MAX_TOF];
double ErrMap[MAX_DET+5][MAX_TOF];
double TofMap[MAX_TOF];
int DetMap[MAX_DET];
int DetCount[MAX_DET];
int SpectraIdx[MAX_DET];
double DetPositions[MAX_DET][3];

typedef struct __MODULEEVT{
  uint8_t psd;
  uint32_t tof;
}Module_Evt;

typedef struct __HEADER{
  uint32_t subsecond;
  uint8_t  module;
  uint8_t  reserve2;
  uint8_t  type;
  uint8_t  header;
}Pulse_Header;

typedef struct __TIME{
  time_t second;
}Pulse_Time;

typedef struct __EVENT{
  uint8_t evt[7];
  uint8_t psd;
}Event;

typedef struct __EOP{
  uint8_t reserve2;
  uint8_t reserve3;
  uint8_t reserve4;
  uint8_t reserve5;
  uint8_t reserve6;
  uint8_t reserve7;
  uint8_t reserve8;
  uint8_t eop;
}EndOfPulse;


uint32_t MapIdx(uint32_t tofidx, uint32_t detidx){
  return tofidx*MAX_DET+detidx;
}

uint32_t TofIdx(uint32_t mapidx){
  return mapidx/MAX_DET;
}

uint32_t DetIdx(uint32_t mapidx){
  return mapidx%MAX_DET;
}

uint8_t PSDIdx(uint32_t detidx){
  return detidx/BIN_DET;
}

uint8_t PosIdx(uint32_t detidx){
  return detidx%BIN_DET;
}

void Encode_PulseHeader(Pulse_Header* pulseHeader, uint8_t *type, uint8_t *module, uint32_t *subsecond){
  pulseHeader->header    = 0x0;
  pulseHeader->type      = *type;
  pulseHeader->module    = *module;
  pulseHeader->subsecond = *subsecond; 
}

void Decode_PulseHeader(uint64_t *buff, uint32_t *type, uint32_t *module, uint32_t *subsecond ){
  *type      = (uint32_t) (((*buff)>>48)&0xFF);
  *module    = (uint32_t) (((*buff)>>32)&0xFF); 
  *subsecond = (uint32_t) ((*buff)&0xFFFFFFFF);
}

void Encode_PulseTime(Pulse_Time* pulseTime, time_t *second){
  pulseTime->second    = *second; 
}

void Decode_PulseTime(uint64_t *buff, time_t *second ){
  *second = (time_t) (* buff); 
}

void Encode_Event(Event* event, uint8_t *psd, uint32_t *tof, uint32_t *qa, uint32_t *qb){
  event->psd = *psd;
  event->evt[6] = ( (*tof) >>20); 
  event->evt[5] = ( (*tof) >>12); 
  event->evt[4] = ( (*tof) >>4); 
  event->evt[3] = (((*tof) & 0xF  )<<4) + (((*qa) & 0x3C00)>>10); 
  event->evt[2] = (((*qa)  & 0x3FC)>>2); 
  event->evt[1] = (((*qa)  & 0x3  )<<6) + (((*qb) & 0x3F00)>>8); 
  event->evt[0] = ((*qb) & 0xFF); 
}

void Decode_Event(uint64_t *buff, uint32_t *psd, uint32_t *tof, uint32_t *qa, uint32_t *qb ){
  *psd = (uint32_t) (((*buff) >> 56 ) & 0xFF);
  *tof = (uint32_t) (((*buff) >> 28 ) & 0xFFFFFFF);
  *qa  = (uint32_t) (((*buff) >> 14 ) & 0x3FFF);
  *qb  = (uint32_t) (( *buff) & 0x3FFF);
}

void Encode_EOP(EndOfPulse* eop){
  eop->eop = 0xFF;
}

void SaveNexusFile(uint32_t* rebinmap, uint32_t* mmap, uint32_t* midx, std::string nexusfilename ){
  std::cout << "SavenexusFile2 " << __LINE__  << std::endl;
  NeXus::File file("test.nxs",  NXACC_CREATE5);
  std::vector<int> dim;

  //enter group: raw_data_1
  file.makeGroup("raw_data_1","NXentry");
  file.openGroup("raw_data_1","NXentry");

  const std::string name = "CSNS_REFL";
  dim.clear();
  dim.push_back(name.length());
  file.makeData("name",NeXus::CHAR,dim,true);
  file.putData(name.c_str());
  file.closeData();

  const double proton_charge[1]={15.0};
  dim.clear();
  dim.push_back(1);
  file.makeData("proton_charge",NeXus::FLOAT64,dim,true);
  file.putData(proton_charge);
  file.closeData();

  const int run_number[1]={11};
  dim.clear();
  dim.push_back(1);
  file.makeData("run_number",NeXus::INT32,dim,true);
  file.putData(run_number);
  file.closeData();

  const std::string end_time = "2014-04-01T21:19:32";
  dim.clear();
  dim.push_back(end_time.length());
  file.makeData("end_time",NeXus::CHAR,dim,true);
  file.putData(end_time.c_str());
  file.closeData();

  const std::string start_time = "2014-04-02T21:19:32";
  dim.clear();
  dim.push_back(start_time.length());
  file.makeData("start_time",NeXus::CHAR,dim,true);
  file.putData(start_time.c_str());
  file.closeData();

  //enter group: raw_data_1/isis_vms_compat
  file.makeGroup("isis_vms_compat","NXentry");
  file.openGroup("isis_vms_compat","NXentry");

  const int nsp1[1]={5};
  dim.clear();
  dim.push_back(1);
  file.makeData("NSP1",NeXus::INT32,dim,true);
  file.putData(nsp1);
  file.closeData();

  const int udet[6]={1,2,3,4,5,6};
  dim.clear();
  dim.push_back(6);
  file.makeData("UDET",NeXus::INT32,dim,true);
  file.putData(udet);
  file.closeData();

  const int spec[6]={1,2,3,4,5,6};
  dim.clear();
  dim.push_back(6);
  file.makeData("SPEC",NeXus::INT32,dim,true);
  file.putData(spec);
  file.closeData();

  const std::string hdr = "what is HDR";
  dim.clear();
  dim.push_back(hdr.length());
  file.makeData("HDR",NeXus::CHAR,dim,true);
  file.putData(hdr.c_str());
  file.closeData();

  int Irpb[32];
  for(int i = 0; i< 32; i++){
    Irpb[i]=i+10; 
  }
  dim.clear();
  dim.push_back(32);
  file.makeData("IRPB",NeXus::INT32,dim,true);
  file.putData(Irpb);
  file.closeData();

  double Rrpb[32];
  for(int i = 0; i< 32; i++){
    Rrpb[i]=(double)(i+100); 
  }
  dim.clear();
  dim.push_back(32);
  file.makeData("RRPB",NeXus::FLOAT32,dim,true);
  file.putData(Rrpb);
  file.closeData();

  int Spb[64];
  for(int i = 0; i< 64; i++){
    Rrpb[i]=0; 
  }
  dim.clear();
  dim.push_back(64);
  file.makeData("SPB",NeXus::INT32,dim,true);
  file.putData(Spb);
  file.closeData();

  double Rspb[64];
  for(int i = 0; i< 64; i++){
    Rspb[i]=0.0; 
  }
  dim.clear();
  dim.push_back(64);
  file.makeData("RSPB",NeXus::FLOAT32,dim,true);
  file.putData(Rspb);
  file.closeData();

  //close group: raw_data_1/isis_vms_compat
  file.closeGroup();

  //enter group: raw_data_1/sample
  file.makeGroup("sample","NXsample");
  file.openGroup("sample","NXsample");


  const std::string samplename = "";
  dim.clear();
  dim.push_back(samplename.length());
  file.makeData("name",NeXus::CHAR,dim,true);
  file.putData(samplename.c_str());
  file.closeData();

  file.makeData("id",NeXus::CHAR,dim,true);
  file.putData(samplename.c_str());
  file.closeData();

  file.makeData("shape",NeXus::CHAR,dim,true);
  file.putData(samplename.c_str());
  file.closeData();

  file.makeData("type",NeXus::CHAR,dim,true);
  file.putData(samplename.c_str());
  file.closeData();

  const double samplenum[1]={0.0};
  dim.clear();
  dim.push_back(1);
  file.makeData("distance",NeXus::FLOAT32,dim,true);
  file.putData(samplenum);
  file.closeData();

  file.makeData("height",NeXus::FLOAT32,dim,true);
  file.putData(samplenum);
  file.closeData();

  file.makeData("thickness",NeXus::FLOAT32,dim,true);
  file.putData(samplenum);
  file.closeData();

  file.makeData("width",NeXus::FLOAT32,dim,true);
  file.putData(samplenum);
  file.closeData();

  //close group: raw_data_1/sample
  file.closeGroup();

  //enter group: raw_data_1/detector_1
  file.makeGroup("detector_1","NXdata");
  file.openGroup("detector_1","NXdata");

  int CountMap[NX_MAX_PRD][NX_MAX_DET][MAX_TOF];
  float TofIdx[MAX_TOF];
  int PeriodIdx[NX_MAX_PRD];
  int SpectrumIdx[NX_MAX_DET];

  for(int i = 0; i< NX_MAX_PRD; i++){
    for(int j = 0; j< NX_MAX_DET; j++){
      for(int k = 0; k< MAX_TOF; k++){
	CountMap[i][j][k]=rebinmap[k]; 
      }
    }
  }

  for(int i = 0; i< MAX_TOF; i++){
    TofIdx[i]=float(i*8.0+11000); 
    //std::cout << TofIdx[i] << std::endl;
  }

  for(int i = 0; i< NX_MAX_PRD; i++){
    PeriodIdx[i]=i+1; 
  }

  for(int i = 0; i< NX_MAX_DET; i++){
    SpectrumIdx[i]=i+4; 
  }

  dim.clear();
  dim.push_back(NX_MAX_PRD);
  dim.push_back(NX_MAX_DET);
  dim.push_back(MAX_TOF);
  file.makeData("counts",NeXus::INT32,dim,true);
  file.putData(CountMap);
  file.putAttr("units","counts");
  file.putAttr("signal","1");
  file.putAttr("axes","period_index,spectrum_index,time_of_flight");
  file.closeData();

  dim.clear();
  dim.push_back(MAX_TOF);
  file.makeData("time_of_flight",NeXus::FLOAT32,dim,true);
  file.putData(TofIdx);
  file.putAttr("units","microseconds");
  file.putAttr("axis","1");
  file.closeData();

  dim.clear();
  dim.push_back(NX_MAX_DET);
  file.makeData("spectrum_index",NeXus::INT32,dim,true);
  file.putData(SpectrumIdx);
  file.closeData();

  dim.clear();
  dim.push_back(NX_MAX_PRD);
  file.makeData("period_index",NeXus::INT32,dim,true);
  file.putData(PeriodIdx);
  file.closeData();

  //close group: raw_data_1/detector_1
  file.closeGroup();

  //enter group: raw_data_1/monitor_1
  file.makeGroup("monitor_1","NXmonitor");
  file.openGroup("monitor_1","NXmonitor");

  int MonitorMap[1][1][2000];

  for(int i = 0; i< 1; i++){
    for(int j = 0; j< 1; j++){
      for(int k = 0; k< 2000; k++){
	MonitorMap[i][j][k]=mmap[k]; 
      }
    }
  }

  float MIdx[2000];
  for(int i = 0; i< 2000; i++){
    MIdx[i]=float(midx[i]); 
  }

  for(int i = 0; i< 1; i++){
    SpectrumIdx[i]=i+1; 
  }

  dim.clear();
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2000);
  file.makeData("data",NeXus::INT32,dim,true);
  file.putData(MonitorMap);
  file.putAttr("units","counts");
  file.putAttr("signal","1");
  file.putAttr("axes","period_index,spectrum_index,time_of_flight");
  file.closeData();

  dim.clear();
  dim.push_back(2000);
  file.makeData("time_of_flight",NeXus::FLOAT32,dim,true);
  file.putData(MIdx);
  file.putAttr("units","microseconds");
  file.putAttr("axis","1");
  file.closeData();

  dim.clear();
  dim.push_back(1);
  file.makeData("spectrum_index",NeXus::INT32,dim,true);
  file.putData(SpectrumIdx);
  file.closeData();

  dim.clear();
  dim.push_back(1);
  file.makeData("monitor_number",NeXus::INT32,dim,true);
  file.putData(SpectrumIdx);
  file.closeData();

  dim.clear();
  dim.push_back(NX_MAX_PRD);
  file.makeData("period_index",NeXus::INT32,dim,true);
  file.putData(PeriodIdx);
  file.closeData();


  //close group: raw_data_1/monitor_1
  file.closeGroup();

  //enter group: raw_data_1/monitor_2
  file.makeGroup("monitor_2","NXmonitor");
  file.openGroup("monitor_2","NXmonitor");

  for(int i = 0; i< 1; i++){
    for(int j = 0; j< 1; j++){
      for(int k = 0; k< 2000; k++){
	MonitorMap[i][j][k]=mmap[k]; 
      }
    }
  }

  for(int i = 0; i< 1; i++){
    SpectrumIdx[i]=i+2; 
  }

  dim.clear();
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2000);
  file.makeData("data",NeXus::INT32,dim,true);
  file.putData(MonitorMap);
  file.putAttr("units","counts");
  file.putAttr("signal","1");
  file.putAttr("axes","period_index,spectrum_index,time_of_flight");
  file.closeData();

  dim.clear();
  dim.push_back(2000);
  file.makeData("time_of_flight",NeXus::FLOAT32,dim,true);
  file.putData(MIdx);
  file.putAttr("units","microseconds");
  file.putAttr("axis","1");
  file.closeData();

  dim.clear();
  dim.push_back(1);
  file.makeData("spectrum_index",NeXus::INT32,dim,true);
  file.putData(SpectrumIdx);
  file.closeData();

  dim.clear();
  dim.push_back(1);
  file.makeData("monitor_number",NeXus::INT32,dim,true);
  file.putData(SpectrumIdx);
  file.closeData();

  dim.clear();
  dim.push_back(NX_MAX_PRD);
  file.makeData("period_index",NeXus::INT32,dim,true);
  file.putData(PeriodIdx);
  file.closeData();


  //close group: raw_data_1/monitor_2
  file.closeGroup();

  //enter group: raw_data_1/monitor_3
  file.makeGroup("monitor_3","NXmonitor");
  file.openGroup("monitor_3","NXmonitor");

  for(int i = 0; i< 1; i++){
    for(int j = 0; j< 1; j++){
      for(int k = 0; k< 2000; k++){
	MonitorMap[i][j][k]=mmap[k]; 
      }
    }
  }

  for(int i = 0; i< 1; i++){
    SpectrumIdx[i]=i+3; 
  }

  dim.clear();
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(2000);
  file.makeData("data",NeXus::INT32,dim,true);
  file.putData(MonitorMap);
  file.putAttr("units","counts");
  file.putAttr("signal","1");
  file.putAttr("axes","period_index,spectrum_index,time_of_flight");
  file.closeData();

  dim.clear();
  dim.push_back(2000);
  file.makeData("time_of_flight",NeXus::FLOAT32,dim,true);
  file.putData(MIdx);
  file.putAttr("units","microseconds");
  file.putAttr("axis","1");
  file.closeData();

  dim.clear();
  dim.push_back(1);
  file.makeData("spectrum_index",NeXus::INT32,dim,true);
  file.putData(SpectrumIdx);
  file.closeData();

  dim.clear();
  dim.push_back(1);
  file.makeData("monitor_number",NeXus::INT32,dim,true);
  file.putData(SpectrumIdx);
  file.closeData();

  dim.clear();
  dim.push_back(NX_MAX_PRD);
  file.makeData("period_index",NeXus::INT32,dim,true);
  file.putData(PeriodIdx);
  file.closeData();


  //close group: raw_data_1/monitor_3
  file.closeGroup();

  //enter group: raw_data_1/selog
  file.makeGroup("selog","IXselog");
  file.openGroup("selog","IXselog");

  float Values[1] = {0.0};
  std::string read_control = "Detector Mode";
  std::string set_control  = "Detector Select";

  //enter group: raw_data_1/selog/PD1H
  file.makeGroup("PD1H","IXseblock");
  file.openGroup("PD1H","IXseblock");

  Values[0] = 200.;

  read_control = "Height Position";
  dim.clear();
  dim.push_back(read_control.length());
  file.makeData("read_control",NeXus::CHAR,dim,true);
  file.putData(read_control.c_str());
  file.closeData();

  set_control = "Height";
  dim.clear();
  dim.push_back(set_control.length());
  file.makeData("set_control",NeXus::CHAR,dim,true);
  file.putData(set_control.c_str());
  file.closeData();

  dim.clear();
  dim.push_back(1);
  file.makeData("value",NeXus::FLOAT32,dim,true);
  file.putData(Values);
  file.putAttr("units","mm");
  file.closeData();

  dim.clear();
  dim.push_back(1);
  file.makeData("setpoint",NeXus::FLOAT32,dim,true);
  file.putData(Values);
  file.putAttr("units","mm");
  file.closeData();

  //enter group: raw_data_1/selog/PD1H/value_log
  file.makeGroup("value_log","NXlog");
  file.openGroup("value_log","NXlog");

  float Log_Values[1] = {0.0};
  float Log_Time[1]   = {0.0};
  std::string  Log_Name = "PD1H log";

  dim.clear();
  dim.push_back(Log_Name.length());
  file.makeData("name",NeXus::CHAR,dim,true);
  file.putData(Log_Name.c_str());
  file.closeData();

  dim.clear();
  dim.push_back(1);
  file.makeData("value",NeXus::FLOAT32,dim,true);
  file.putData(Values);
  file.putAttr("units","mm");
  file.closeData();

  dim.clear();
  dim.push_back(1);
  file.makeData("time",NeXus::FLOAT32,dim,true);
  file.putData(Log_Time);
  file.putAttr("units","seconds");
  file.putAttr("start","2014-05-10T13:56:10");
  file.closeData();


  //close group: raw_data_1/selog/PD1H/value_log
  file.closeGroup();


  //close group: raw_data_1/selog/PD1H
  file.closeGroup();

  //close group: raw_data_1/selog
  file.closeGroup();

  //close group: raw_data_1
  file.closeGroup();

  // close file
  file.close();
}

void SaveHeaderToBinaryFile(ofstream& fout, uint8_t *type, uint8_t *module, uint32_t *subsecond) {
  /*----------------------------------------------*/
  Pulse_Header* pulseHeader = new Pulse_Header;
  Encode_PulseHeader(pulseHeader, type, module, subsecond);
  fout.write((char*)pulseHeader, sizeof(Pulse_Header));
  delete pulseHeader;
}

void SaveTimeStampToBinaryFile(ofstream& fout, time_t *second){
  /*----------------------------------------------*/
  Pulse_Time* pulseTime = new Pulse_Time;
  Encode_PulseTime(pulseTime, second);
  fout.write((char*)pulseTime, sizeof(Pulse_Time));
  delete pulseTime;
}

void SaveEventToBinaryFile(ofstream& fout, uint8_t *PSD, uint32_t *TOF, uint32_t *QA, uint32_t *QB){
  Event* event = new Event;
  Encode_Event(event, PSD, TOF, QA, QB);
  fout.write((char*)event, sizeof(Event));
  delete event;
}

void SaveEOPToBinaryFile(ofstream& fout){
  /*----------------------------------------------*/
  EndOfPulse* eop = new EndOfPulse;
  Encode_EOP(eop); 
  fout.write((char*)eop, sizeof(EndOfPulse));
  delete eop;
}

uint32_t Get_PositionID(uint32_t qa, uint32_t qb){
  if ((qa+qb)<1){
    //std::cout << "qa " << qa << " qb " << qb << std::endl;
    return 0;
  }
  double R = (double)qa/(qa+qb);
  double P = R*BIN_DET;
  int IntP = (uint32_t)P;
  double D = P - IntP;
  if(D>0.5){
    return (IntP + 1);
  }
  else{
    return IntP;
  }
}

void SaveBinaryFile(uint32_t *cmap, std::string binaryfilename){
  int counts = 0;
  int pulses = 0;
  std::cout << "SaveBinaryFile" << std::endl;
  std::ofstream fout(binaryfilename.c_str(), std::ios::binary); 

  std::time_t UnixTime = std::time(0);  // t is an integer type
  time_t   second    = UnixTime;
  uint32_t subsecond = 0;

  uint8_t type   = 0x0;
  uint8_t module = 0x1;

  SaveHeaderToBinaryFile(fout, &type, &module, &subsecond);
  SaveTimeStampToBinaryFile(fout, &second);
  /*----------------------------------------------*/
  srand((int)time(0));
  std::cout << "SaveEvent" << std::endl;
  for(int i=0; i < MAX_TOF ; i++){
    uint32_t TOF = i;
    for(int j=0; j < MAX_DET ; j++){
      uint8_t PSD = PSDIdx(j);
      uint8_t pos = PosIdx(j);
      double R = ((double)pos)/BIN_DET;
      for(int k=0;k<cmap[MapIdx(i,j)];k++){
	// rand end of one pluse, and begin an new pulse 
	if((rand()%10)<2.0){
	  pulses++;
	  subsecond = 150000000*(pulses%25);
	  if(0==(pulses%25))UnixTime=UnixTime+100000;
	  if(0==(pulses%100000)) std::cout << "SaveEOP " << pulses << " pulses "<< counts  << " events " << std::endl; 
	  SaveEOPToBinaryFile(fout);
	  SaveHeaderToBinaryFile(fout, &type, &module, &subsecond);
	  SaveTimeStampToBinaryFile(fout, &second);

	}
	// one neutron event
	counts +=1;
	uint32_t QB = 0;
	uint32_t QA = 0;
	if(R<0.5){
	  QB = (rand()%1024+1);
	  QA = R*QB/(1-R);
	}
	else{
	  QA = (rand()%1024+1);
	  QB = QA*(1-R)/R;
	}
	SaveEventToBinaryFile(fout, &PSD, &TOF, &QA, &QB);
	// one neutron event
      }
    }
  }
  SaveEOPToBinaryFile(fout);
  std::cout << "SaveEOP " << pulses << " pulses "<< counts  << " events " << std::endl; 
  std::cout << "Save count " << counts << std::endl; 


  /*----------------------------------------------*/
  fout.close();
}


void Map_EventToDetector(uint32_t *dmap, uint32_t *module, uint32_t *psd, uint32_t *tof, uint32_t *qa, uint32_t *qb){
  uint32_t pos_id = Get_PositionID(*qa, *qb);
  uint32_t det_id = (*psd) * BIN_DET + pos_id;
  uint32_t id = (*tof) * MAX_DET + det_id;
  //std::cout << "(" << *psd << "/" << pos_id << "/" << det_id << "/" << *tof << "/"<<id << ")" << std::endl;
  dmap[id] += 1;
  //std::cout << " Map Event: module " << *module << " psd " << *psd 
  //  << " tof " << *tof << " qa "  << *qa << " qb " << *qb  
  //  << " pos id " << pos_id  << " det id " << det_id 
  //  << " Matrix " << id << std::endl; 
}


void LoadSimulationFile(uint32_t* cmap, std::string samplefilename){

    std::cout << "LoadSimulationFile "<< std::endl;
    std::ifstream samplefile(samplefilename.c_str());
    //std::ifstream samplefile("/home/tianhl/workarea/CSNS_SANS_SIM/app/test/test_raw");
    string samplebuff;
    getline(samplefile, samplebuff);
    uint32_t tot=0;

    for (int tofidx=0;tofidx<MAX_TOF ;tofidx++){
      getline(samplefile, samplebuff);
      vector<string> substring;
      //std::cout << "new line " << tofidx << " " << samplebuff << std::endl;
      vector<double> counts;
      //boost::split( substring, samplebuff, boost::is_any_of( "\t " ), boost::token_compress_on );
      boost::split( substring, samplebuff, boost::is_any_of( ";" ), boost::token_compress_on );
      //std::cout <<"Process Line " <<  tofidx << std::endl;
      for(uint32_t detidx = 0; detidx < MAX_DET  ; detidx++){
	cmap[MapIdx(tofidx, detidx)] =atoi(substring[detidx+1].c_str());
	cmap[MapIdx(tofidx, detidx)] =(int)(atoi(substring[detidx+1].c_str()));
	//      std::cout << "tofidx " << tofidx << " detidx " << detidx 
	//	<< " MapIdx " << MapIdx(tofidx,detidx) << std::endl;
	//      std::cout << "rtof   " << TofIdx(MapIdx(tofidx,detidx))  
	//	<< " rdet   " << DetIdx(MapIdx(tofidx,detidx)) << std::endl;
	tot+=(int)(atoi(substring[detidx+1].c_str()));
      }
    }
    std::cout << "total neutron hit count: " << tot << std::endl;
    samplefile.close();
}

void LoadMonitorFile(uint32_t* midx, uint32_t* mmap, std::string samplefilename){
  std::cout << "LoadMonitorFile " << samplefilename << std::endl;
  std::ifstream samplefile(samplefilename.c_str());

  string samplebuff;
  getline(samplefile, samplebuff);

  //std::cout << samplebuff << std::endl; 
  for (int tofidx=0;tofidx<2000 ;tofidx++){
    getline(samplefile, samplebuff);

    vector<string> substring;
    boost::split( substring, samplebuff, boost::is_any_of( ";" ), boost::token_compress_on );
    mmap[tofidx] =atoi(substring[1].c_str());
    midx[tofidx] =atoi(substring[0].c_str());
    std::cout << midx[tofidx] << " " << mmap[tofidx] << std::endl; 
  }
  samplefile.close();

}


uint64_t Decode_RawDataSegment(uint64_t *Buff, uint32_t *dmap, uint32_t size, uint8_t *flag){
  //std::cout << "Enter Decode_RawDataSegment(), buffer size: " << size << std::endl;
  uint64_t count = 0;
  uint64_t *ReadRawData;// = new uint8_t[8]; 
  time_t second;
  uint32_t type, module, subsecond;
  for (uint32_t i = 0; i < size ; i++ ){
    ReadRawData = (uint64_t*)(Buff+i);
    //std::cout << "idx " << i << " " << std::hex << *ReadRawData << std::endl;// << std::endl;
    //std::cout << "zzz: " << std::hex << *ReadRawData << std::endl;
    if ((((*ReadRawData)>>56) == 0x0) && (*flag == 0))  {
      Decode_PulseHeader(ReadRawData, &type, &module, &subsecond);
      //std::cout << " Header: type " << type << " module  " << module << " subsecond " << subsecond << std::endl;
      *flag = 1;
    }
    else if (*flag == 1){
      //int64_t second;
      Decode_PulseTime(ReadRawData, &second);
      //std::cout << " Time: second " << second  <<" "<< ctime((time_t*)&second)    <<  std::endl; 
      *flag = 2; 
      continue;
    }
    if ((((*ReadRawData)>>56) == 0xFF)&&(*flag >= 2)) {
      //std::cout << " EndOfPulse" << std::endl; 
      *flag = 0;
    }
    else if ((*flag == 2)||(*flag == 3)){
      uint32_t psd, tof, qa, qb;
      Decode_Event(ReadRawData, &psd, &tof, &qa, &qb);
      count += 1;
      if ((qa+qb)<1){
	std::cout << "decode qa " << qa << " qb " << qb << std::endl;
      }
      Map_EventToDetector(dmap, &module,&psd, &tof, &qa, &qb );
      *flag = 3;
    }
  }
  return count;
}

void LoadBinaryFile(uint32_t *dmap, std::string binaryfilename){

  /*----------------------------------------------*/
  uint32_t size = 10000;
  uint64_t count = 0;
  uint8_t *flag = new uint8_t;
  *flag = 0;
  size_t buffsize = 0; 
  uint64_t *Buff = new uint64_t[size]; 
  std::ifstream fin(binaryfilename.c_str(), std::ios::binary);
  fin.read((char*)Buff, sizeof(uint64_t)*size);
  buffsize = fin.gcount();  
  count += Decode_RawDataSegment(Buff, dmap, size, flag);
  while (buffsize == (sizeof(uint64_t)*size)){
    //std::cout << "LoadBinaryFile " << fin.gcount() << std::endl;
    fin.read((char*)Buff, sizeof(uint64_t)*size);
    buffsize = fin.gcount();  
    if ((sizeof(uint64_t)*size) == buffsize ){
      count += Decode_RawDataSegment(Buff, dmap, size, flag);
    }
    else{
      //std::cout << "Read file " << buffsize/(sizeof(uint64_t)) << std::endl;
      count += Decode_RawDataSegment(Buff, dmap, buffsize/(sizeof(uint64_t)), flag);
    }
  }
  delete flag;
  //std::cout << " raw data: "<< (uint32_t)ReadRawData[2] << " " << sizeof(time_t) << " " << sizeof(uint32_t)<< std::endl;
  std::cout << " Read Event Count " << count << std::endl;


  fin.close();
}


double Rebin(uint32_t* dmap, uint32_t* rebinmap){
  int sum=0;
  float positionidx=0;

  for(int t=0; t< MAX_TOF ; t++){
    rebinmap[t]=0;
    int integral[100];
    for(int i = 0; i < 100; i++){
      integral[i]=0;
      for(int j = 0; j < 80; j++){
	int det = j+i*100;
	int mapidx = MapIdx(t, det);
	integral[i]+=dmap[mapidx];
      }
    }
    for(int i =0; i < 100; i++){
      rebinmap[t]+=integral[i];
      sum += integral[i];
      positionidx += float(i*integral[i]);
    }

  }
  return  positionidx/sum*2 ;
}

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    std::cout << "Usage: " << argv[0] << "  option.txt" << std::endl;
    return 1;
  }

  uint32_t *cmap = new uint32_t[MAX_TOF*MAX_DET];
  uint32_t *dmap = new uint32_t[MAX_TOF*MAX_DET];
  uint32_t *mmap = new uint32_t[2000];
  uint32_t *midx = new uint32_t[2000];
  uint32_t *rebinmap = new uint32_t[MAX_TOF];

  std::string configfile(argv[1]);
  Config* fConfig = new Config(configfile);

  std::string samplefile  ; 
  std::string binaryfile  ;
  std::string monitorfile ;
  std::string nexusfile   ;

  samplefile  = fConfig->pString("samplefile") ;  
  binaryfile  = fConfig->pString("binaryfile") ; 
  monitorfile = fConfig->pString("monitorfile") ; 
  nexusfile   = fConfig->pString("nexusfile") ;

  std::cout << "read sample file: " << samplefile << std::endl;
  std::cout << "save binary file: " << binaryfile << std::endl;

  //LoadSimulationFile(cmap, samplefile); 
  //for(int i = 0; i< MAX_TOF; i++){
  //  for(int j = 0; j< MAX_DET; j++){
  //    std::cout << cmap[MapIdx(i, j)] << " " ;
  //  }
  //  std::cout << std::endl;
  //}
  //SaveBinaryFile(cmap, binaryfile);
  //LoadMonitorFile(midx, mmap, monitorfile); 
  LoadBinaryFile(dmap, binaryfile);
  //std::cout << std::endl;
  //for(int i = 0; i< MAX_TOF; i++){
  //  for(int j = 0; j< MAX_DET; j++){
  //    std::cout << dmap[MapIdx(i, j)] << " " ;
  //  }
  //  std::cout << std::endl;
  //}

  Rebin(dmap, rebinmap);
  SaveNexusFile(rebinmap ,mmap,midx,nexusfile);

  return 0;
}
