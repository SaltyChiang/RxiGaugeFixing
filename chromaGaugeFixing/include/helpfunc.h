#pragma once

#include <getopt.h>

template <typename T>
void printMatrix(T &matrix)
{
  for (int i = 0; i < Nc; i++)
  {
    for (int j = 0; j < Nc; j++)
      printf("(%10.6f,%10.6f) ", matrix.elem(i, j).real(), matrix.elem(i, j).imag());
    printf("\n");
  }
}

class OptionParser
{
public:
  OptionParser(int argc, char **argv)
  {
    int c = -1;
    while ((c = getopt_long(argc, argv, shortOpts, longOpts, NULL)) != -1)
    {
      if (c != '?')
      {
        name[c] = true;
        value[c] = optarg;
      }
      else
        printf("No option named %c.\n", c);
    }
  }

  ~OptionParser()
  {
  }

  char *parseOption(const int option)
  {
    if (name[option])
      return value[option];
    else
    {
      printf("No option named %c.\n", option);
      return NULL;
    }
  }

  char *parseOption(const int option, char *defaultVal)
  {
    if (name[option])
      return value[option];
    else
    {
      return defaultVal;
    }
  }

private:
  const char *shortOpts = (char *)"i:o:x:y:z:t:r:f:a:m:d:p:";
  struct option longOpts[8] = {{"input", 1, NULL, 'i'},
                               {"output", 1, NULL, 'o'},
                               {"xi", 1, NULL, 'r'},
                               {"format", 1, NULL, 'f'},
                               {"iter-accuracy", 1, NULL, 'a'},
                               {"iter-max", 1, NULL, 'm'},
                               {"overrelax-do", 1, NULL, 'd'},
                               {"overrelax-param", 1, NULL, 'p'}};
  bool name[128] = {false};
  char *value[128];
};