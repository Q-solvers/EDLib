//
// Created by iskakoff on 04/05/17.
//

#include <gtest/gtest.h>
#include <edlib/ExecutionStatistic.h>

TEST(Statistics, UpdateEventTest) {
  EDLib::common::statistics.registerEvent("test");
  std::pair < double, double > timing = EDLib::common::statistics.event("test");
  ASSERT_EQ(0, timing.first);
  // single event
  EDLib::common::statistics.updateEvent("test");
  std::pair < double, double > updated_timing = EDLib::common::statistics.event("test");
  ASSERT_DOUBLE_EQ(updated_timing.second - timing.second, updated_timing.first);
  // repeated event re-registration
  EDLib::common::statistics.registerEvent("test");
  std::pair < double, double > timing2 = EDLib::common::statistics.event("test");
  ASSERT_DOUBLE_EQ(timing2.first, updated_timing.first);
  // repeated event
  EDLib::common::statistics.updateEvent("test");
  std::pair < double, double > updated_timing2 = EDLib::common::statistics.event("test");
  ASSERT_DOUBLE_EQ(updated_timing2.second - timing2.second + updated_timing.first, updated_timing2.first);
}