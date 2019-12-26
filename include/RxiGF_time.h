#pragma once

#include <time.h>
#ifdef _WIN32

#include <windows.h>

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS 11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS 11644473600000000ULL
#endif

struct timezone
{
  long tz_minuteswest; /* minutes W of Greenwich */
  int tz_dsttime;      /* type of dst correction */
};

int gettimeofday(timeval *tv, timezone *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag = 0;

  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);

    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;

    tmpres /= 10; /*convert into microseconds*/
    /*converting file time to unix epoch*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS;
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }

  if (NULL != tz)
  {
    if (!tzflag)
    {
      _tzset();
      tzflag++;
    }
    _get_timezone(&(tz->tz_minuteswest));
    tz->tz_minuteswest /= 60;
    _get_daylight(&(tz->tz_dsttime));
  }

  return 0;
}
#else
#include <sys/time.h>
#endif

#define StartTime(_num)            \
  timeval start##_num, stop##_num; \
  gettimeofday(&start##_num, NULL);

#define StopTime(_num) \
  gettimeofday(&stop##_num, NULL);

#define PrintTime(_num, _str) \
  printf("Time for calculating %s is %ld microseconds.\n", _str, (stop##_num.tv_sec - start##_num.tv_sec) * 1000000 + stop##_num.tv_usec - start##_num.tv_usec);