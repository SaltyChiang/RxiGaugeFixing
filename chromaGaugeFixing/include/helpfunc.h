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
  std::string input, output;
  int x = 24, y = 24, z = 24, t = 64;
  double xi = 0.0;
  double gfAccu = 1e-10;
  int gfMax = 1000;
  bool orDo = false;
  double orPara = 1.7;

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

    parseOption('i', input);
    parseOption('o', output);
    parseOption('x', x);
    parseOption('y', y);
    parseOption('z', z);
    parseOption('t', t);
    parseOption('r', xi);
    parseOption('a', gfAccu);
    parseOption('m', gfMax);
    parseOption('d', orDo);
    parseOption('p', orPara);
  }

  ~OptionParser()
  {
  }

private:
  const char *shortOpts = "i:o:x:y:z:t:r:f:a:m:d:p:";
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

  void parseOption(const int optionId, int &optionValue)
  {
    if (name[optionId])
      optionValue = atoi(value[optionId]);
  }
  void parseOption(const int optionId, bool &optionValue)
  {
    if (name[optionId])
      optionValue = atoi(value[optionId]);
  }
  void parseOption(const int optionId, double &optionValue)
  {
    if (name[optionId])
      optionValue = atof(value[optionId]);
  }
  void parseOption(const int optionId, float &optionValue)
  {
    if (name[optionId])
      optionValue = atof(value[optionId]);
  }
  void parseOption(const int optionId, std::string &optionValue)
  {
    if (name[optionId])
      optionValue = value[optionId];
  }
};