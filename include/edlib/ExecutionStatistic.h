//
// Created by iskakoff on 02/02/17.
//

#ifndef HUBBARD_EXECUTIONSTATISTIC_H
#define HUBBARD_EXECUTIONSTATISTIC_H


#ifdef USE_MPI
#include <mpi.h>
#endif

#include <string>
#include <map>
#include <chrono>

namespace EDLib {
  namespace common {
/**
 * @brief ExecutionStatistic class
 *
 * @author iskakoff
 */
    class ExecutionStatistic {
    public:

      ExecutionStatistic() {}

      /**
       * Update event time
       * @param name - event name
       */
      void updateEvent(const std::string& name) {
        double time1 = time();
        _events[name] = std::make_pair(_events[name].first + time1 - _events[name].second, time1);
      }

      /**
       * register the start point of the event
       *
       * @param name - event name
       */
      void registerEvent(const std::string& name) {
        _events[name] = std::make_pair(_events[name].first, time());
      }

      /**
       * Print all observed events
       */
      void print() {
        for (auto& kv : _events) {
          std::cout <<"Event "<< kv.first << " take " << kv.second.first << "s." << std::endl;
        }
      }

      /**
       * Return event timing pair
       * @param event_name - event name
       * @return event timing
       */
      std::pair<double, double> event(const std::string & event_name) {
        if(_events.find(event_name) != _events.end()) {
          return _events[event_name];
        }
        return std::make_pair(0.0, 0.0);
      };
    private:
      // registered events timing pairs
      // pair.first corresponds to total event time
      // pair.second corresponds to last time when event was happened
      std::map<std::string, std::pair<double, double> > _events;

      double time() const {
#ifdef USE_MPI
        return MPI_Wtime();
#else
        return std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
#endif
      }
    };

    static EDLib::common::ExecutionStatistic statistics;
  }
}



#endif //HUBBARD_EXECUTIONSTATISTIC_H
