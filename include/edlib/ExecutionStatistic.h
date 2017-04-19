//
// Created by iskakoff on 02/02/17.
//

#ifndef HUBBARD_EXECUTIONSTATISTIC_H
#define HUBBARD_EXECUTIONSTATISTIC_H


namespace EDLib {
  namespace common {
/**
 * @brief ExecutionStatistic class
 *
 * @author iskakoff
 */
    class ExecutionStatistic {
    public:


      void updateEvent(const std::string& name, long long newTime) {

      }

      void registerEvent(const std::string& name) {

      }
    private:
      std::map<std::string, int> _events;
    };
  }
}


#endif //HUBBARD_EXECUTIONSTATISTIC_H
