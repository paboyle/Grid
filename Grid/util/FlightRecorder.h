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

#ifdef GRID_LOG_VIEWS
class ViewLogger {
  struct Entry_t {
    const char* filename;
    int line;
    int index;
    uint64_t head, tail;
  };

public:
  static bool Enabled;
  static std::vector<Entry_t> LogVector;
  static void Begin();
  static void End();
  static void Log(const char* filename, int line, int index, int mode, void* data, uint64_t bytes);
};
#endif
NAMESPACE_END(Grid);

