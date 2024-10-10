#pragma once

NAMESPACE_BEGIN(Grid);
class FlightRecorder {
 public:
  enum LoggingMode_t {
    LoggingModeNone,
    LoggingModePrint,
    LoggingModeRecord,
    LoggingModeVerify
  };
  
  static int                   LoggingMode;
  static uint64_t              ErrorCounter;
  static const char *                StepName;
  static int32_t               StepLoggingCounter;
  static int32_t               XmitLoggingCounter;
  static int32_t               RecvLoggingCounter;
  static int32_t               CsumLoggingCounter;
  static int32_t               NormLoggingCounter;
  static int32_t               ReductionLoggingCounter;
  static std::vector<uint64_t> XmitLogVector;
  static std::vector<uint64_t> RecvLogVector;
  static std::vector<uint64_t> CsumLogVector;
  static std::vector<double>   NormLogVector;
  static std::vector<double>   ReductionLogVector;
  static int ContinueOnFail;
  static int PrintEntireLog;
  static int ChecksumComms;
  static int ChecksumCommsSend;
  static void SetLoggingModePrint(void);
  static void SetLoggingModeRecord(void);
  static void SetLoggingModeVerify(void);
  static void SetLoggingMode(LoggingMode_t mode);
  static bool StepLog(const char *name);
  static bool NormLog(double value);
  static bool CsumLog(uint64_t csum);
  static void ReductionLog(double lcl, double glbl);
  static void Truncate(void);
  static void ResetCounters(void);
  static uint64_t ErrorCount(void);
  static void xmitLog(void *,uint64_t bytes);
  static void recvLog(void *,uint64_t bytes,int rank);
};
NAMESPACE_END(Grid);

