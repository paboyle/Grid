/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/Init.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@MacBook-Pro.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

NAMESPACE_BEGIN(Grid);
///////////////////////////////////////////////////////
// Grid Norm logging for repro testing
///////////////////////////////////////////////////////
int FlightRecorder::PrintEntireLog;
int FlightRecorder::ContinueOnFail;
int FlightRecorder::LoggingMode;
int FlightRecorder::ChecksumComms;
int FlightRecorder::ChecksumCommsSend;
int32_t  FlightRecorder::XmitLoggingCounter;
int32_t  FlightRecorder::RecvLoggingCounter;
int32_t  FlightRecorder::CsumLoggingCounter;
int32_t  FlightRecorder::NormLoggingCounter;
int32_t  FlightRecorder::ReductionLoggingCounter;
uint64_t FlightRecorder::ErrorCounter;
std::vector<double> FlightRecorder::NormLogVector;
std::vector<double> FlightRecorder::ReductionLogVector;
std::vector<uint64_t> FlightRecorder::CsumLogVector;
std::vector<uint64_t> FlightRecorder::XmitLogVector;
std::vector<uint64_t> FlightRecorder::RecvLogVector;

void FlightRecorder::ResetCounters(void)
{
  XmitLoggingCounter=0;
  RecvLoggingCounter=0;
  CsumLoggingCounter=0;
  NormLoggingCounter=0;
  ReductionLoggingCounter=0;
}
void FlightRecorder::Truncate(void)
{
  ResetCounters();
  XmitLogVector.resize(0);
  RecvLogVector.resize(0);
  NormLogVector.resize(0);
  CsumLogVector.resize(0);
  ReductionLogVector.resize(0);
}
void FlightRecorder::SetLoggingMode(FlightRecorder::LoggingMode_t mode)
{
  switch ( mode ) {
  case LoggingModePrint:
    SetLoggingModePrint();
    break;
  case LoggingModeRecord:
    SetLoggingModeRecord();
    break;
  case LoggingModeVerify:
    SetLoggingModeVerify();
    break;
  case LoggingModeNone:
    LoggingMode = mode;
    Truncate();
    break;
  default:
    assert(0);
  }
}

void FlightRecorder::SetLoggingModePrint(void)
{
  std::cout << " FlightRecorder: set to print output " <<std::endl;
  Truncate();
  LoggingMode = LoggingModePrint;
}
void FlightRecorder::SetLoggingModeRecord(void)
{
  std::cout << " FlightRecorder: set to RECORD " <<std::endl;
  Truncate();
  LoggingMode = LoggingModeRecord;
}
void FlightRecorder::SetLoggingModeVerify(void)
{
  std::cout << " FlightRecorder: set to VERIFY " << NormLogVector.size()<< " log entries "<<std::endl;
  ResetCounters();
  LoggingMode = LoggingModeVerify;
}
uint64_t FlightRecorder::ErrorCount(void)
{
  return ErrorCounter;
}
void FlightRecorder::NormLog(double value)
{
  uint64_t hex = * ( (uint64_t *)&value );
  if(LoggingMode == LoggingModePrint) {
    std::cerr<<"FlightRecorder::NormLog : "<< NormLoggingCounter <<" "<<std::hex<< hex<<std::dec <<std::endl;
    NormLoggingCounter++;
  }
  if(LoggingMode == LoggingModeRecord) {
    std::cerr<<"FlightRecorder::NormLog RECORDING : "<< NormLoggingCounter <<" "<<std::hex<< hex<<std::dec <<std::endl;
    NormLogVector.push_back(value);
    NormLoggingCounter++;
  }
  if(LoggingMode == LoggingModeVerify) {

    if(NormLoggingCounter < NormLogVector.size()){
      uint64_t hexref  = * ( (uint64_t *)&NormLogVector[NormLoggingCounter] );

      if ( (value != NormLogVector[NormLoggingCounter]) || std::isnan(value) ) {

	std::cerr<<"FlightRecorder::NormLog Oops, I did it again "<< NormLoggingCounter
		 <<std::hex<<" "<<hex<<" "<<hexref<<std::dec<<" "
		 <<std::hexfloat<<value<<" "<< NormLogVector[NormLoggingCounter]<<std::endl;

	std::cerr << " Oops got norm "<< std::hexfloat<<value<<" expect "<<NormLogVector[NormLoggingCounter] <<std::endl;

	fprintf(stderr,"%s:%d Oops, I did it again! Reproduce failure for norm %d/%zu %.16e expect %.16e\n",
		GridHostname(),
		GlobalSharedMemory::WorldShmRank,
		NormLoggingCounter,NormLogVector.size(),
		value, NormLogVector[NormLoggingCounter]); fflush(stderr);

	if(!ContinueOnFail)assert(0); // Force takedown of job
	  
	ErrorCounter++;
      } else {
	if ( PrintEntireLog ) { 
	  std::cerr<<"FlightRecorder::NormLog VALID "<< NormLoggingCounter << std::hex
		   <<" "<<hex<<" "<<hexref
		   <<" "<<std::hexfloat<<value<<" "<< NormLogVector[NormLoggingCounter]<<std::dec<<std::endl;
	}
      }
       
    }
    if ( NormLogVector.size()==NormLoggingCounter ) {
      std::cout << "FlightRecorder:: Verified entire sequence of "<<NormLoggingCounter<<" norms "<<std::endl;
    }
    NormLoggingCounter++;
  }
}
void FlightRecorder::CsumLog(uint64_t hex)
{
  if(LoggingMode == LoggingModePrint) {
    std::cerr<<"FlightRecorder::CsumLog : "<< CsumLoggingCounter <<" "<<std::hex<< hex<<std::dec <<std::endl;
    CsumLoggingCounter++;
  }

  if(LoggingMode == LoggingModeRecord) {
    std::cerr<<"FlightRecorder::CsumLog RECORDING : "<< NormLoggingCounter <<" "<<std::hex<< hex<<std::dec <<std::endl;
    CsumLogVector.push_back(hex);
    CsumLoggingCounter++;
  }

  if(LoggingMode == LoggingModeVerify) {
    
    if(CsumLoggingCounter < CsumLogVector.size()) {

      uint64_t hexref  = CsumLogVector[CsumLoggingCounter] ;

      if ( hex != hexref ) {

        std::cerr<<"FlightRecorder::CsumLog Oops, I did it again "<< CsumLoggingCounter
		 <<std::hex<<" "<<hex<<" "<<hexref<<std::dec<<std::endl;

	fprintf(stderr,"%s:%d Oops, I did it again! Reproduce failure for csum %d %lx expect %lx\n",
		GridHostname(),
		GlobalSharedMemory::WorldShmRank,
		CsumLoggingCounter,hex, hexref);
	fflush(stderr);

	if(!ContinueOnFail) assert(0); // Force takedown of job
	  
	ErrorCounter++;

      } else {

	if ( PrintEntireLog ) { 
	  std::cerr<<"FlightRecorder::CsumLog VALID "<< CsumLoggingCounter << std::hex
		   <<" "<<hex<<" "<<hexref<<std::dec<<std::endl;
	}
      }
    }  
    if ( CsumLogVector.size()==CsumLoggingCounter ) {
      std::cout << "FlightRecorder:: Verified entire sequence of "<<CsumLoggingCounter<<" checksums "<<std::endl;
    }
    CsumLoggingCounter++;
  }
}
void FlightRecorder::ReductionLog(double local,double global)
{
  uint64_t hex_l = * ( (uint64_t *)&local );
  uint64_t hex_g = * ( (uint64_t *)&global );
  if(LoggingMode == LoggingModePrint) {
    std::cerr<<"FlightRecorder::ReductionLog : "<< ReductionLoggingCounter <<" "<< std::hex << hex_l << " -> " <<hex_g<<std::dec <<std::endl;
    ReductionLoggingCounter++;
  }
  if(LoggingMode == LoggingModeRecord) {
    std::cerr<<"FlightRecorder::ReductionLog RECORDING : "<< ReductionLoggingCounter <<" "<< std::hex << hex_l << " -> " <<hex_g<<std::dec <<std::endl;
    ReductionLogVector.push_back(global);
    ReductionLoggingCounter++;
  }
  if(LoggingMode == LoggingModeVerify) {
    if(ReductionLoggingCounter < ReductionLogVector.size()){
      if ( global != ReductionLogVector[ReductionLoggingCounter] ) {
	fprintf(stderr,"%s:%d Oops, MPI_Allreduce did it again! Reproduce failure for norm %d/%zu glb %.16e lcl %.16e expect glb %.16e\n",
		GridHostname(),
		GlobalSharedMemory::WorldShmRank,
		ReductionLoggingCounter,ReductionLogVector.size(),
		global, local, ReductionLogVector[ReductionLoggingCounter]); fflush(stderr);
	
	if ( !ContinueOnFail ) assert(0);

	ErrorCounter++;
      } else {
	if ( PrintEntireLog ) { 
	  std::cerr<<"FlightRecorder::ReductionLog : VALID "<< ReductionLoggingCounter <<" "<< std::hexfloat << local << "-> "<< global <<std::endl;
	}
      }
    }
    if ( ReductionLogVector.size()==ReductionLoggingCounter ) {
      std::cout << "FlightRecorder::ReductionLog : Verified entire sequence of "<<ReductionLoggingCounter<<" norms "<<std::endl;
    }
    ReductionLoggingCounter++;
  }
}
void FlightRecorder::xmitLog(void *buf,uint64_t bytes)
{
  if ( ChecksumCommsSend ){
  uint64_t *ubuf = (uint64_t *)buf;
  if(LoggingMode == LoggingModeNone) return;
#ifdef GRID_SYCL
  uint64_t _xor = svm_xor(ubuf,bytes/sizeof(uint64_t));
  if(LoggingMode == LoggingModePrint) {
    std::cerr<<"FlightRecorder::xmitLog : "<< XmitLoggingCounter <<" "<< std::hex << _xor <<std::dec <<std::endl;
    XmitLoggingCounter++;
  }
  if(LoggingMode == LoggingModeRecord) {
    std::cerr<<"FlightRecorder::xmitLog RECORD : "<< XmitLoggingCounter <<" "<< std::hex << _xor <<std::dec <<std::endl;
    XmitLogVector.push_back(_xor);
    XmitLoggingCounter++;
  }
  if(LoggingMode == LoggingModeVerify) {
    if(XmitLoggingCounter < XmitLogVector.size()){
      if ( _xor != XmitLogVector[XmitLoggingCounter] ) {
	fprintf(stderr,"%s:%d Oops, send buf difference! Reproduce failure for xmit %d/%zu  %lx expect glb %lx\n",
		GridHostname(),
		GlobalSharedMemory::WorldShmRank,
		XmitLoggingCounter,XmitLogVector.size(),
		_xor, XmitLogVector[XmitLoggingCounter]); fflush(stderr);
	
	if ( !ContinueOnFail ) assert(0);

	ErrorCounter++;
      } else {
	if ( PrintEntireLog ) { 
	  std::cerr<<"FlightRecorder::XmitLog : VALID "<< XmitLoggingCounter <<" "<< std::hexfloat << _xor << " "<<  XmitLogVector[XmitLoggingCounter] <<std::endl;
	}
      }
    }
    if ( XmitLogVector.size()==XmitLoggingCounter ) {
      std::cout << "FlightRecorder::ReductionLog : Verified entire sequence of "<<XmitLoggingCounter<<" sends "<<std::endl;
    }
    XmitLoggingCounter++;
  }
#endif
  } else {
    uint64_t word = 1;
    deviceVector<uint64_t> dev(1);
    acceleratorCopyToDevice(&word,&dev[0],sizeof(uint64_t));
    acceleratorCopySynchronise();
#ifndef GRID_COMMS_NONE
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
}
void FlightRecorder::recvLog(void *buf,uint64_t bytes,int rank)
{
  if ( ChecksumComms ){
  uint64_t *ubuf = (uint64_t *)buf;
  if(LoggingMode == LoggingModeNone) return;
#ifdef GRID_SYCL
  uint64_t _xor = svm_xor(ubuf,bytes/sizeof(uint64_t));
  if(LoggingMode == LoggingModePrint) {
    std::cerr<<"FlightRecorder::recvLog : "<< RecvLoggingCounter <<" "<< std::hex << _xor <<std::dec <<std::endl;
    RecvLoggingCounter++;
  }
  if(LoggingMode == LoggingModeRecord) {
    std::cerr<<"FlightRecorder::recvLog RECORD : "<< RecvLoggingCounter <<" "<< std::hex << _xor <<std::dec <<std::endl;
    RecvLogVector.push_back(_xor);
    RecvLoggingCounter++;
  }
  if(LoggingMode == LoggingModeVerify) {
    if(RecvLoggingCounter < RecvLogVector.size()){
      if ( _xor != RecvLogVector[RecvLoggingCounter] ) {
	fprintf(stderr,"%s:%d Oops, recv buf difference! Reproduce failure for recv %d/%zu  %lx expect glb %lx from MPI rank %d\n",
		GridHostname(),
		GlobalSharedMemory::WorldShmRank,
		RecvLoggingCounter,RecvLogVector.size(),
		_xor, RecvLogVector[RecvLoggingCounter],rank); fflush(stderr);
	
	if ( !ContinueOnFail ) assert(0);

	ErrorCounter++;
      } else {
	if ( PrintEntireLog ) { 
	  std::cerr<<"FlightRecorder::RecvLog : VALID "<< RecvLoggingCounter <<" "<< std::hexfloat << _xor << " "<<  RecvLogVector[RecvLoggingCounter] <<std::endl;
	}
      }
    }
    if ( RecvLogVector.size()==RecvLoggingCounter ) {
      std::cout << "FlightRecorder::ReductionLog : Verified entire sequence of "<<RecvLoggingCounter<<" sends "<<std::endl;
    }
    RecvLoggingCounter++;
  }
#endif
  }
}

NAMESPACE_END(Grid);
