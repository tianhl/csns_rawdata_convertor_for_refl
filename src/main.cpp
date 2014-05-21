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
const uint32_t MAX_TOF = 128;
const uint32_t MAX_DET =   2;
const uint32_t MAX_PRD =   1;

#include "log.h"
#include "config.h"


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

void SaveNexusFile2(){
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

  int CountMap[MAX_PRD][MAX_DET][MAX_TOF];
  float TofIdx[MAX_TOF];
  int PeriodIdx[MAX_PRD];
  int SpectrumIdx[MAX_DET];

  for(int i = 0; i< MAX_PRD; i++){
    for(int j = 0; j< MAX_DET; j++){
      for(int k = 0; k< MAX_TOF; k++){
	CountMap[i][j][k]=j*10+k; 
      }
    }
  }

  for(int i = 0; i< MAX_TOF; i++){
    TofIdx[i]=float(i*5.0+0.5); 
    std::cout << TofIdx[i] << std::endl;
  }

  for(int i = 0; i< MAX_PRD; i++){
    PeriodIdx[i]=i+1; 
  }

  for(int i = 0; i< MAX_DET; i++){
    SpectrumIdx[i]=i+4; 
  }

  dim.clear();
  dim.push_back(MAX_PRD);
  dim.push_back(MAX_DET);
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
  dim.push_back(MAX_DET);
  file.makeData("spectrum_index",NeXus::INT32,dim,true);
  file.putData(SpectrumIdx);
  file.closeData();

  dim.clear();
  dim.push_back(MAX_PRD);
  file.makeData("period_index",NeXus::INT32,dim,true);
  file.putData(PeriodIdx);
  file.closeData();

  //close group: raw_data_1/detector_1
  file.closeGroup();

  //enter group: raw_data_1/monitor_1
  file.makeGroup("monitor_1","NXmonitor");
  file.openGroup("monitor_1","NXmonitor");

  int MonitorMap[1][1][MAX_TOF];

  for(int i = 0; i< 1; i++){
    for(int j = 0; j< 1; j++){
      for(int k = 0; k< MAX_TOF; k++){
	MonitorMap[i][j][k]=1000+k; 
      }
    }
  }

  for(int i = 0; i< 1; i++){
    SpectrumIdx[i]=i+1; 
  }

  dim.clear();
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(MAX_TOF);
  file.makeData("data",NeXus::INT32,dim,true);
  file.putData(MonitorMap);
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
  dim.push_back(MAX_PRD);
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
      for(int k = 0; k< MAX_TOF; k++){
	MonitorMap[i][j][k]=2000+k; 
      }
    }
  }

  for(int i = 0; i< 1; i++){
    SpectrumIdx[i]=i+2; 
  }

  dim.clear();
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(MAX_TOF);
  file.makeData("data",NeXus::INT32,dim,true);
  file.putData(MonitorMap);
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
  dim.push_back(MAX_PRD);
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
      for(int k = 0; k< MAX_TOF; k++){
	MonitorMap[i][j][k]=3000+k; 
      }
    }
  }

  for(int i = 0; i< 1; i++){
    SpectrumIdx[i]=i+3; 
  }

  dim.clear();
  dim.push_back(1);
  dim.push_back(1);
  dim.push_back(MAX_TOF);
  file.makeData("data",NeXus::INT32,dim,true);
  file.putData(MonitorMap);
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
  dim.push_back(MAX_PRD);
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

  ////enter group: raw_data_1/selog/Detector
  //file.makeGroup("Detector","IXseblock");
  //file.openGroup("Detector","IXseblock");


  //dim.clear();
  //dim.push_back(read_control.length());
  //file.makeData("read_control",NeXus::CHAR,dim,true);
  //file.putData(read_control.c_str());
  //file.closeData();

  //dim.clear();
  //dim.push_back(set_control.length());
  //file.makeData("set_control",NeXus::CHAR,dim,true);
  //file.putData(set_control.c_str());
  //file.closeData();

  //dim.clear();
  //dim.push_back(1);
  //file.makeData("value",NeXus::FLOAT32,dim,true);
  //file.putData(Values);
  //file.putAttr("units","mm");
  //file.closeData();

  //dim.clear();
  //dim.push_back(1);
  //file.makeData("setpoint",NeXus::FLOAT32,dim,true);
  //file.putData(Values);
  //file.putAttr("units","mm");
  //file.closeData();

  ////close group: raw_data_1/selog/Detector
  //file.closeGroup();

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

  ////enter group: raw_data_1/selog/PD1A
  //file.makeGroup("PD1A","IXseblock");
  //file.openGroup("PD1A","IXseblock");

  //Values[0] = 0.3;

  //read_control = "PD1 Angle Position";
  //dim.clear();
  //dim.push_back(read_control.length());
  //file.makeData("read_control",NeXus::CHAR,dim,true);
  //file.putData(read_control.c_str());
  //file.closeData();

  //set_control = "PD1 Angle";
  //dim.clear();
  //dim.push_back(set_control.length());
  //file.makeData("set_control",NeXus::CHAR,dim,true);
  //file.putData(set_control.c_str());
  //file.closeData();

  //dim.clear();
  //dim.push_back(1);
  //file.makeData("value",NeXus::FLOAT32,dim,true);
  //file.putData(Values);
  //file.putAttr("units","mm");
  //file.closeData();

  //dim.clear();
  //dim.push_back(1);
  //file.makeData("setpoint",NeXus::FLOAT32,dim,true);
  //file.putData(Values);
  //file.putAttr("units","mm");
  //file.closeData();

  ////close group: raw_data_1/selog/PD1A
  //file.closeGroup();

  //close group: raw_data_1/selog
  file.closeGroup();

  //close group: raw_data_1
  file.closeGroup();

  // close file
  file.close();
}

int main(int argc, char *argv[])
{
  SaveNexusFile2();

  return 0;
}
