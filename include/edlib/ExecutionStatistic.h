#ifndef EDLIB_EXECUTIONSTATISTIC_H
#define EDLIB_EXECUTIONSTATISTIC_H

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <chrono>
#include <iostream>
#include <map>
#include <string>
#include <utility>

namespace edlib {

  class ExecutionStatistic {
  public:
    void updateEvent(const std::string& name) {
      double t = time();
      _events[name] = {_events[name].first + t - _events[name].second, t};
    }

    void registerEvent(const std::string& name) {
      _events[name] = {_events[name].first, time()};
    }

    void print() const {
      for (const auto& kv : _events) {
        std::cout << "Event " << kv.first << " take " << kv.second.first << "s." << std::endl;
      }
    }

    std::pair<double, double> event(const std::string& name) const {
      auto it = _events.find(name);
      if (it == _events.end()) return {0.0, 0.0};
      return it->second;
    }

  private:
    std::map<std::string, std::pair<double, double>> _events;

    static double time() {
#ifdef USE_MPI
      return MPI_Wtime();
#else
      return std::chrono::duration_cast<std::chrono::duration<double>>(
                 std::chrono::high_resolution_clock::now().time_since_epoch())
          .count();
#endif
    }
  };

  inline ExecutionStatistic statistics;

}

#endif
