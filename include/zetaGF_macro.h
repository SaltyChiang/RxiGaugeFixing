#pragma once

#include <cstdlib>

#define zgfMalloc(_dataType, _dataCount) \
  (_dataType *)malloc(sizeof(_dataType) * (_dataCount))

#include <chrono>

#define StartTimeChrono(_num) \
  auto start##_num = std::chrono::high_resolution_clock::now()

#define StopTimeChrono(_num) \
  auto stop##_num = std::chrono::high_resolution_clock::now()

#ifdef _WIN32
#define PrintTimeChrono(_num, _str)                                                                  \
  auto dura##_num = std::chrono::duration_cast<std::chrono::microseconds>(stop##_num - start##_num); \
  printf("Time for calculating %s is %lld microseconds.\n", _str, dura##_num.count())
#else
#define PrintTimeChrono(_num, _str)                                                                  \
  auto dura##_num = std::chrono::duration_cast<std::chrono::microseconds>(stop##_num - start##_num); \
  printf("Time for calculating %s is %ld microseconds.\n", _str, dura##_num.count())
#endif