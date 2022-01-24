#ifndef AOFLAGGER
#define AOFLAGGER

#include <sstream>
#include <iostream>

#include <mutex>

class Logger {
 public:
  enum LoggerLevel {
    NoLevel = 5,
    FatalLevel = 4,
    ErrorLevel = 3,
    WarningLevel = 2,
    InfoLevel = 1,
    DebugLevel = 0
  };

  enum VerbosityLevel { QuietVerbosity, NormalVerbosity, VerboseVerbosity };

  template <enum LoggerLevel Level, bool ToStdErr = false>
  class LogWriter {
   public:
    LogWriter() : _atNewLine(true) {}

    LogWriter &operator<<(const std::string &str) {
      std::lock_guard<std::mutex> lock(_mutex);
      size_t start = 0, end;
      while (std::string::npos != (end = str.find('\n', start))) {
        outputLinePart(str.substr(start, end - start + 1), true);
        start = end + 1;
      }
      outputLinePart(str.substr(start, str.size() - start), false);
      return *this;
    }
    LogWriter &operator<<(const char *str) {
      (*this) << std::string(str);
      return *this;
    }
    LogWriter &operator<<(const char c) {
      std::lock_guard<std::mutex> lock(_mutex);
      outputLinePart(std::string(1, c), c == '\n');
      return *this;
    }
    template <typename S>
    LogWriter &operator<<(const S &str) {
      std::ostringstream stream;
      stream << str;
      (*this) << stream.str();
      return *this;
    }
    void Flush() {
      std::lock_guard<std::mutex> lock(_mutex);
      if (ToStdErr)
        std::cerr.flush();
      else
        std::cout.flush();
    }

   private:
    std::mutex _mutex;
    bool _atNewLine;

    void outputLinePart(const std::string &str, bool endsWithCR) {
      if ((int)_coutLevel <= (int)Level && !str.empty()) {
        if (_atNewLine && _logTime) outputTime(ToStdErr);
        if (ToStdErr)
          std::cerr << str;
        else
          std::cout << str;
        _atNewLine = endsWithCR;
      }
    }
  };

  static void SetVerbosity(VerbosityLevel verbosityLevel);

  static bool IsVerbose() { return _coutLevel == DebugLevel; }

  static void SetLogTime(bool logTime) { _logTime = logTime; }

  static bool LogTime() { return _logTime; }

  static class LogWriter<DebugLevel> Debug;
  static class LogWriter<InfoLevel> Info;
  static class LogWriter<WarningLevel> Warn;
  static class LogWriter<ErrorLevel> Error;
  static class LogWriter<FatalLevel> Fatal;
  static class LogWriter<NoLevel, true> Progress;

 private:
  Logger() {}

  static void outputTime(bool toStdErr);

  static enum LoggerLevel _coutLevel;

  static bool _logTime;
};

class LogReceiver {
 public:
  LogReceiver()
      : Fatal(this), Error(this), Warn(this), Info(this), Debug(this) {}

  template <enum Logger::LoggerLevel Level>
  class LevelReceiver {
   public:
    LevelReceiver(LogReceiver *parent) : _parent(parent) {}
    LevelReceiver &operator<<(const std::string &str) {
      size_t start = 0, end;
      while (std::string::npos != (end = str.find('\n', start))) {
        _parent->output(Level, str.substr(start, end - start + 1));
        start = end + 1;
      }
      _parent->output(Level, str.substr(start, str.size() - start));
      return *this;
    }
    LevelReceiver &operator<<(const char *str) {
      (*this) << std::string(str);
      return *this;
    }
    LevelReceiver &operator<<(const char c) {
      _parent->output(Level, std::string(1, c));
      return *this;
    }
    template <typename S>
    LevelReceiver &operator<<(const S &str) {
      std::ostringstream stream;
      stream << str;
      (*this) << stream.str();
      return *this;
    }

   private:
    LogReceiver *_parent;
  };  // end of class LevelReceiver

  LevelReceiver<Logger::FatalLevel> Fatal;
  LevelReceiver<Logger::ErrorLevel> Error;
  LevelReceiver<Logger::WarningLevel> Warn;
  LevelReceiver<Logger::InfoLevel> Info;
  LevelReceiver<Logger::DebugLevel> Debug;

 protected:
  virtual void output(enum Logger::LoggerLevel level,
                      const std::string &str) = 0;

  void forward(enum Logger::LoggerLevel level, const std::string &str) {
    switch (level) {
      case Logger::FatalLevel:
        Logger::Fatal << str;
        break;
      case Logger::ErrorLevel:
        Logger::Error << str;
        break;
      case Logger::WarningLevel:
        Logger::Warn << str;
        break;
      case Logger::InfoLevel:
        Logger::Info << str;
        break;
      case Logger::DebugLevel:
        Logger::Debug << str;
        break;
      case Logger::NoLevel:
        break;
    }
  }
};

class ForwardingLogReceiver : public LogReceiver {
 protected:
  virtual void output(enum Logger::LoggerLevel level,
                      const std::string &str) final override {
    forward(level, str);
  }
};

#endif
